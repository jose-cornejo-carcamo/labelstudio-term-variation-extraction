# =============================================================================
# labelstudio-term-variation
# Extraction of base term / variant pairs from Label Studio JSON exports
# =============================================================================
#
# What this script does
#   Reads the JSON export of ONE annotated text from Label Studio and writes
#   three CSV files listing the annotated base terms, their variants and the
#   type(s) of variation linking them.
#
# Expected annotation scheme (Label Studio labeling interface)
#   - Two span labels: one for base terms ("Terme vedette" or "Terme de base")
#     and one for variants ("Variante").
#   - Relations drawn between a base term and a variant, carrying one or more
#     variation labels (e.g. "VR extension", "VL changement base").
#   - Discontinuous terms: each constituent is annotated separately, then the
#     constituents are attached to one another (parent-child) in the Regions
#     panel, which fills the "parentID" field of the export.
#
# Outputs (written in the same folder as the JSON file)
#   <file>_terminology_analysis_detailed.csv : one row per base term / variant pair
#   <file>_all_terms_with_counts.csv         : every term with its number of occurrences
#   <file>_compound_terms_list.csv           : discontinuous terms rebuilt from parentID
#
# How to use
#   1. Install the packages below if needed.
#   2. Replace the placeholders between [ ] (look for the "[CHANGE]" comments).
#   3. Run the whole script.
# =============================================================================

# --- Required packages -------------------------------------------------------
# If needed, install them once with:
# install.packages(c("jsonlite", "dplyr", "tidyr", "stringr", "purrr"))
library(jsonlite)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# --- Step 1: read the JSON export ---------------------------------------------
# Read a JSON export file from Label Studio
read_json_file <- function(path) {
  data <- fromJSON(path, simplifyVector = FALSE)
  
  if (is.list(data) && !is.null(names(data))) {
    # If it's a named list (single object), wrap in list
    return(list(data))
  }
  return(data)
}

# --- Step 2: extract the annotations and their metadata (task, project) -------
# Extract annotations (results) from a list of LS tasks
extract_annotations <- function(data_list) {
  annotations <- list()
  
  for (item in data_list) {
    anns <- item[["annotations"]]
    if (is.null(anns)) next
    
    for (ann in anns) {
      results <- ann[["result"]]
      if (is.null(results) || length(results) == 0) next
      
      annotations[[length(annotations) + 1]] <- list(
        task = ifelse(is.null(ann[["task"]]), 
                      ifelse(is.null(item[["id"]]), "unknown", item[["id"]]), 
                      ann[["task"]]),
        project = ifelse(is.null(item[["project"]]), 
                         ifelse(is.null(ann[["project"]]), "unknown", ann[["project"]]), 
                         item[["project"]]),
        results = results
      )
    }
  }
  
  return(annotations)
}

# --- Helper used in step 7: removes one trailing punctuation mark that may have
# been included by mistake when selecting a term in the interface
# Clean term text (remove trailing punctuation)
clean_term <- function(text) {
  if (!is.character(text) || is.na(text)) return(text)
  
  text <- str_trim(text)
  if (nchar(text) > 0 && str_sub(text, -1) %in% c("(", ")", "[", "]", ",", ".", ";", ":", "?", "!", "_", "-")) {
    text <- str_trim(str_sub(text, 1, -2))
  }
  return(text)
}

# --- Step 3: separate text segments ("labels") from relations ----------------
# Split one annotation's results into labels and relations
split_labels_relations <- function(results) {
  labels <- Filter(function(r) !is.null(r[["type"]]) && r[["type"]] == "labels", results)
  relations <- Filter(function(r) !is.null(r[["type"]]) && r[["type"]] == "relation", results)
  
  return(list(labels = labels, relations = relations))
}

# --- Step 4: build the table of labels (id, text, label type, parentID,
# start and end positions) ------------------------------------------------------
# Build a dataframe of labels with parent/child info
build_labels_df <- function(labels) {
  if (length(labels) == 0) {
    return(data.frame(
      id = character(),
      text = character(),
      label_type = character(),
      parent_id = character(),
      start_pos = numeric(),
      end_pos = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  rows <- lapply(labels, function(lab) {
    val <- lab[["value"]]
    if (is.null(val[["text"]]) || is.null(val[["labels"]])) return(NULL)
    
    label_type <- if (is.list(val[["labels"]])) val[["labels"]][[1]] else val[["labels"]]
    parent_id <- ifelse(is.null(lab[["parentID"]]), 
                        ifelse(is.null(val[["parentID"]]), NA, val[["parentID"]]), 
                        lab[["parentID"]])
    
    data.frame(
      id = lab[["id"]],
      text = val[["text"]],
      label_type = label_type,
      parent_id = parent_id,
      start_pos = ifelse(is.null(val[["start"]]), NA, val[["start"]]),
      end_pos = ifelse(is.null(val[["end"]]), NA, val[["end"]]),
      stringsAsFactors = FALSE
    )
  })
  
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame(
      id = character(),
      text = character(),
      label_type = character(),
      parent_id = character(),
      start_pos = numeric(),
      end_pos = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  return(bind_rows(rows))
}

# --- Step 5: rebuild discontinuous terms from the parentID field -------------
# Segments attached to the same parent are grouped, ordered by their position
# in the text, and later joined with a space (step 6).
# Find compound terms (labels with children)
find_compound_groups <- function(labels_df) {
  if (nrow(labels_df) == 0) {
    return(list(groups = list(), processed_ids = character()))
  }
  
  groups <- list()
  processed_ids <- character()
  
  for (i in 1:nrow(labels_df)) {
    item <- labels_df[i, ]
    if (item$id %in% processed_ids) next
    
    children <- labels_df[!is.na(labels_df$parent_id) & labels_df$parent_id == item$id, ]
    
    if (nrow(children) > 0) {
      group <- bind_rows(item, children) %>%
        arrange(start_pos)
      groups[[item$id]] <- group
      processed_ids <- c(processed_ids, group$id)
    }
  }
  
  return(list(groups = groups, processed_ids = processed_ids))
}

# --- Step 6: build the terminological registry --------------------------------
# One entry per term (text + label type), with the IDs of all its occurrences.
# Build registry of final terms (compound + simple)
build_term_registry <- function(labels_df, compound_groups, processed_ids, project_id, task_id) {
  id_to_term <- list()
  term_registry <- list()
  
  # Compound terms
  for (parent_id in names(compound_groups)) {
    group <- compound_groups[[parent_id]]
    combined_text <- paste(group$text, collapse = " ")
    combined_label_type <- group$label_type[1]
    term_key <- paste0(combined_text, "|", combined_label_type)
    
    if (is.null(term_registry[[term_key]])) {
      term_registry[[term_key]] <- list(
        projects = character(),
        tasks = character(),
        final_text = combined_text,
        label_type = combined_label_type,
        is_compound = TRUE,
        all_component_ids = character(),
        group_ids = character()
      )
    }
    
    term_registry[[term_key]]$all_component_ids <- c(term_registry[[term_key]]$all_component_ids, group$id)
    term_registry[[term_key]]$group_ids <- c(term_registry[[term_key]]$group_ids, parent_id)
    term_registry[[term_key]]$projects <- unique(c(term_registry[[term_key]]$projects, as.character(project_id)))
    term_registry[[term_key]]$tasks <- unique(c(term_registry[[term_key]]$tasks, as.character(task_id)))
    
    for (gid in group$id) {
      id_to_term[[gid]] <- list(
        final_text = combined_text,
        label_type = combined_label_type,
        is_compound = TRUE,
        parent_id = parent_id,
        component_ids = group$id,
        term_key = term_key
      )
    }
  }
  
  # Simple terms
  simple_items <- labels_df[!(labels_df$id %in% processed_ids), ]
  
  if (nrow(simple_items) > 0) {
    for (i in 1:nrow(simple_items)) {
      item <- simple_items[i, ]
      term_key <- paste0(item$text, "|", item$label_type)
      
      if (is.null(term_registry[[term_key]])) {
        term_registry[[term_key]] <- list(
          projects = character(),
          tasks = character(),
          final_text = item$text,
          label_type = item$label_type,
          is_compound = FALSE,
          all_component_ids = character(),
          group_ids = character()
        )
      }
      
      term_registry[[term_key]]$all_component_ids <- c(term_registry[[term_key]]$all_component_ids, item$id)
      term_registry[[term_key]]$group_ids <- c(term_registry[[term_key]]$group_ids, item$id)
      term_registry[[term_key]]$projects <- unique(c(term_registry[[term_key]]$projects, as.character(project_id)))
      term_registry[[term_key]]$tasks <- unique(c(term_registry[[term_key]]$tasks, as.character(task_id)))
      
      id_to_term[[item$id]] <- list(
        final_text = item$text,
        label_type = item$label_type,
        is_compound = FALSE,
        parent_id = NA,
        component_ids = item$id,
        term_key = term_key
      )
    }
  }
  
  return(list(id_to_term = id_to_term, term_registry = term_registry))
}

# --- Step 7a: convert the relations into a table (from_id, to_id, variation
# labels, direction) --------------------------------------------------------------
# Convert raw relations into a dataframe
build_relations_df <- function(relations) {
  if (length(relations) == 0) {
    return(data.frame(
      from_id = character(),
      to_id = character(),
      relationship = character(),
      direction = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  rows <- lapply(relations, function(rel) {
    if (is.null(rel[["from_id"]]) || is.null(rel[["to_id"]])) return(NULL)
    
    labels <- rel[["labels"]]
    relationship <- if (!is.null(labels) && length(labels) > 0) {
      paste(labels, collapse = "; ")
    } else {
      ""
    }
    
    data.frame(
      from_id = rel[["from_id"]],
      to_id = rel[["to_id"]],
      relationship = relationship,
      direction = ifelse(is.null(rel[["direction"]]), "right", rel[["direction"]]),
      stringsAsFactors = FALSE
    )
  })
  
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame(
      from_id = character(),
      to_id = character(),
      relationship = character(),
      direction = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  return(bind_rows(rows))
}

# --- Step 7b: pair base terms and variants ------------------------------------
# The base term and the variant are identified by their labels, so a relation
# drawn in either direction is paired correctly; the direction recorded in the
# export is kept in the "direction" column. Base terms without any variant are
# also listed (with an empty "Variant" column).
# Build raw result rows (base_term -> variant)
process_relationships <- function(relations_df, id_to_term, term_registry, project_id, task_id) {
  results <- list()
  skipped <- 0
  found <- 0
  seen_bases <- character()
  
  if (nrow(relations_df) > 0) {
    for (i in 1:nrow(relations_df)) {
      rel <- relations_df[i, ]
      
      from_term <- id_to_term[[rel$from_id]]
      to_term <- id_to_term[[rel$to_id]]
      
      if (is.null(from_term) || is.null(to_term)) {
        skipped <- skipped + 1
        next
      }
      
      base <- NULL
      var <- NULL
      
      # [CHANGE IF NEEDED] The label names below ("Terme de base", "Terme vedette",
      # "Variante") must match exactly the span labels defined in your Label
      # Studio labeling interface. If you use other names, replace them here and
      # in the block "Add base terms without variants" further down.
      # Handle both "Terme de base" and "Terme vedette" as base terms
      is_base_from <- from_term$label_type %in% c("Terme de base", "Terme vedette")
      is_base_to <- to_term$label_type %in% c("Terme de base", "Terme vedette")
      
      if (is_base_from && to_term$label_type == "Variante") {
        base <- from_term
        var <- to_term
      } else if (from_term$label_type == "Variante" && is_base_to) {
        base <- to_term
        var <- from_term
      }
      
      if (!is.null(base) && !is.null(var)) {
        found <- found + 1
        base_key <- base$term_key
        var_key <- var$term_key
        seen_bases <- c(seen_bases, base$final_text)
        
        base_component_ids <- unique(term_registry[[base_key]]$all_component_ids)
        var_component_ids <- unique(term_registry[[var_key]]$all_component_ids)
        base_occurrences <- length(unique(term_registry[[base_key]]$group_ids))
        var_occurrences <- length(unique(term_registry[[var_key]]$group_ids))
        
        base_text <- clean_term(base$final_text)
        var_text <- clean_term(var$final_text)
        
        results[[length(results) + 1]] <- data.frame(
          task = as.character(task_id),
          project = as.character(project_id),
          `Terme de base` = base_text,
          `Terme de base ID` = paste(base_component_ids, collapse = ";"),
          `Terme de base Occurrences` = base_occurrences,
          base_is_compound = base$is_compound,
          Variant = var_text,
          `Variant ID` = paste(var_component_ids, collapse = ";"),
          `Variant Occurrences` = var_occurrences,
          variant_is_compound = var$is_compound,
          `Relation Label` = rel$relationship,
          direction = rel$direction,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    }
  }
  
  # [CHANGE IF NEEDED] Same label names as above.
  # Add base terms without variants (handle both "Terme de base" and "Terme vedette")
  for (key in names(term_registry)) {
    info <- term_registry[[key]]
    if (info$label_type %in% c("Terme de base", "Terme vedette") && !(info$final_text %in% seen_bases)) {
      component_ids <- unique(info$group_ids)
      occurrences <- length(unique(info$group_ids))
      
      results[[length(results) + 1]] <- data.frame(
        task = paste(sort(unique(as.character(info$tasks))), collapse = ";"),
        project = paste(sort(unique(as.character(info$projects))), collapse = ";"),
        `Terme de base` = info$final_text,
        `Terme de base ID` = paste(component_ids, collapse = ";"),
        `Terme de base Occurrences` = occurrences,
        base_is_compound = info$is_compound,
        Variant = NA,
        `Variant ID` = NA,
        `Variant Occurrences` = 0,
        variant_is_compound = FALSE,
        `Relation Label` = NA,
        direction = NA,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  
  if (length(results) == 0) {
    return(list(
      df = data.frame(
        task = character(),
        project = character(),
        `Terme de base` = character(),
        `Terme de base ID` = character(),
        `Terme de base Occurrences` = numeric(),
        base_is_compound = logical(),
        Variant = character(),
        `Variant ID` = character(),
        `Variant Occurrences` = numeric(),
        variant_is_compound = logical(),
        `Relation Label` = character(),
        direction = character(),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      found = found,
      skipped = skipped
    ))
  }
  
  return(list(df = bind_rows(results), found = found, skipped = skipped))
}

# --- Summary of all terms with their occurrence counts (second and third CSV) --
# Summarize all final terms with occurrence counts
summarize_terms <- function(term_registry) {
  rows <- lapply(names(term_registry), function(key) {
    info <- term_registry[[key]]
    component_ids <- unique(info$group_ids)
    occurrences <- length(unique(info$group_ids))
    
    data.frame(
      projects = paste(sort(unique(as.character(info$projects))), collapse = ";"),
      tasks = paste(sort(unique(as.character(info$tasks))), collapse = ";"),
      final_text = info$final_text,
      final_id = paste(component_ids, collapse = ";"),
      label_type = info$label_type,
      is_compound = info$is_compound,
      occurrences = occurrences,
      stringsAsFactors = FALSE
    )
  })
  
  if (length(rows) == 0) {
    return(data.frame(
      projects = character(),
      tasks = character(),
      final_text = character(),
      final_id = character(),
      label_type = character(),
      is_compound = logical(),
      occurrences = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  return(bind_rows(rows))
}

# --- Step 8a: remove duplicates ----------------------------------------------
# Rows with the same base term, variant, variant ID and set of variation labels
# are kept once (the order of the labels does not matter).
# Remove duplicates considering relation label order as non-relevant
deduplicate_results <- function(df) {
  initial_count <- nrow(df)
  if (initial_count == 0) return(df)
  
  # Normalize relation labels by sorting
  df <- df %>%
    mutate(
      Relation_Label_Normalized = sapply(`Relation Label`, function(x) {
        if (is.na(x) || x == "") return("")
        labels <- str_trim(unlist(str_split(x, ";")))
        labels <- labels[labels != ""]
        paste(sort(labels), collapse = "; ")
      })
    )
  
  # Remove duplicates
  df_clean <- df %>%
    distinct(`Terme de base`, Variant, `Variant ID`, Relation_Label_Normalized, .keep_all = TRUE) %>%
    select(-Relation_Label_Normalized)
  
  return(df_clean)
}

# --- Step 8b: write the three CSV files --------------------------------------
# Write CSV outputs
write_outputs <- function(final_results, all_terms, out_dir, base_name) {
  # Reorder columns
  col_order <- c(
    "task", "project",
    "Terme de base", "Terme de base ID", "Terme de base Occurrences", "base_is_compound",
    "Variant", "Variant ID", "Variant Occurrences", "variant_is_compound",
    "Relation Label", "direction"
  )
  
  final_results <- final_results %>%
    select(all_of(col_order)) %>%
    arrange(tolower(`Terme de base`))
  
  # [OPTIONAL] Change the suffixes below if you prefer other file names.
  detailed_path <- file.path(out_dir, paste0(base_name, "_terminology_analysis_detailed.csv"))
  all_terms_path <- file.path(out_dir, paste0(base_name, "_all_terms_with_counts.csv"))
  compounds_path <- file.path(out_dir, paste0(base_name, "_compound_terms_list.csv"))
  
  write.csv(final_results, detailed_path, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(all_terms, all_terms_path, row.names = FALSE, fileEncoding = "UTF-8")
  
  compounds <- all_terms %>% filter(is_compound == TRUE)
  if (nrow(compounds) > 0) {
    write.csv(compounds, compounds_path, row.names = FALSE, fileEncoding = "UTF-8")
  }
  
  return(compounds)
}

# --- Main function: runs steps 1 to 8 on one JSON file -------------------------
# Main pipeline: JSON -> CSVs
process_file <- function(json_file) {
  data <- read_json_file(json_file)
  annotations <- extract_annotations(data)
  
  all_results <- list()
  
  for (ann in annotations) {
    split_result <- split_labels_relations(ann$results)
    labels <- split_result$labels
    relations <- split_result$relations
    
    labels_df <- build_labels_df(labels)
    compound_result <- find_compound_groups(labels_df)
    compound_groups <- compound_result$groups
    processed_ids <- compound_result$processed_ids
    
    registry_result <- build_term_registry(labels_df, compound_groups, processed_ids, ann$project, ann$task)
    id_to_term <- registry_result$id_to_term
    term_registry <- registry_result$term_registry
    
    relations_df <- build_relations_df(relations)
    relationship_result <- process_relationships(relations_df, id_to_term, term_registry, ann$project, ann$task)
    
    all_results[[length(all_results) + 1]] <- list(
      df = relationship_result$df,
      term_registry = term_registry
    )
  }
  
  # Combine everything
  combined_results <- bind_rows(lapply(all_results, function(x) x$df))
  combined_results <- deduplicate_results(combined_results)
  
  # Merge all term registries
  all_term_registries <- list()
  for (result in all_results) {
    for (key in names(result$term_registry)) {
      if (is.null(all_term_registries[[key]])) {
        all_term_registries[[key]] <- result$term_registry[[key]]
      } else {
        # Merge projects and tasks
        all_term_registries[[key]]$projects <- unique(c(
          all_term_registries[[key]]$projects,
          result$term_registry[[key]]$projects
        ))
        all_term_registries[[key]]$tasks <- unique(c(
          all_term_registries[[key]]$tasks,
          result$term_registry[[key]]$tasks
        ))
        all_term_registries[[key]]$all_component_ids <- unique(c(
          all_term_registries[[key]]$all_component_ids,
          result$term_registry[[key]]$all_component_ids
        ))
        all_term_registries[[key]]$group_ids <- unique(c(
          all_term_registries[[key]]$group_ids,
          result$term_registry[[key]]$group_ids
        ))
      }
    }
  }
  
  all_terms <- summarize_terms(all_term_registries)
  
  # Write outputs
  # Output folder: the folder of the JSON file.
  # [OPTIONAL] To write elsewhere, replace dirname(json_file) with a folder path,
  # e.g. "[path/to/output/folder]".
  out_dir <- dirname(json_file)
  base_name <- tools::file_path_sans_ext(basename(json_file))
  compounds <- write_outputs(combined_results, all_terms, out_dir, base_name)
  
  cat(sprintf("[%s] Processed %d relations, %d unique terms, %d compounds.\n",
              base_name, nrow(combined_results), nrow(all_terms), nrow(compounds)))
  
  return(list(
    combined_results = combined_results,
    all_terms = all_terms,
    compounds = compounds
  ))
}

# Main execution
# [CHANGE] Replace [your_export_file] with the name of your Label Studio JSON
# export (keep the quotes and the .json extension), e.g. "Text_1.json".
# If the file is not in the working directory, either set the working
# directory first, e.g. setwd("[path/to/your/folder]"), or give the full
# path to the file.
json_file <- "[your_export_file].json"

result <- process_file(json_file)

# [OPTIONAL] To process every JSON file of a folder at once, comment out the
# two lines "json_file <- ..." and "result <- ..." above, then uncomment the
# two lines below and replace [path/to/your/folder]:
# json_files <- list.files("[path/to/your/folder]", pattern = "\\.json$", full.names = TRUE)
# results <- lapply(json_files, process_file)
