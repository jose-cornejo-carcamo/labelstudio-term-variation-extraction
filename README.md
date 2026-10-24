# labelstudio-term-variation-extraction

R pipeline that converts Label Studio JSON exports of annotated terminological variation into CSV tables of base term/variant pairs. It rebuilds discontinuous terms, pairs each variant with its base term whatever the direction of the relation, and removes duplicates while keeping occurrence counts.

The annotation itself is done manually in Label Studio; this script only extracts and structures the resulting data from the JSON export files. It was developed for annotations of intratextual terminological variation in research articles in Earth and Planetary Sciences, and can be used for other domains as long as the annotations follow the same scheme (see below).

## Requirements

- R
- Packages: `jsonlite`, `dplyr`, `tidyr`, `stringr`, `purrr`

```r
install.packages(c("jsonlite", "dplyr", "tidyr", "stringr", "purrr"))
```

## Expected annotation scheme

The script assumes the following setup in the Label Studio labeling interface:

- **Two span labels**: one for base terms (`Terme vedette` or `Terme de base`) and one for variants (`Variante`). If your labels have other names, change them in the script (see the `[CHANGE IF NEEDED]` comments).
- **Relations** drawn between a base term and a variant, carrying one or more variation labels (for example `VR extension` or `VL changement base`). When a variant combines several types of variation, all labels are attached to a single relation.
- **Discontinuous terms**: each constituent is annotated separately. The constituents are then attached to one another (parent–child) by drag and drop in the *Regions* panel, which fills the `parentID` field of the export.

## Usage

1. Export the annotated text from Label Studio in JSON format.
2. Open `Code_for_Extraction_Adapted.R` and replace the placeholders between `[ ]`:
   - `[your_export_file].json`: name of your JSON export (required);
   - `[path/to/your/folder]`: only if the file is not in the working directory or if you process a whole folder (optional).
3. Run the whole script.

The script processes one JSON file at a time. A commented block at the end of the script shows how to process every JSON file in a folder.

## Outputs

The three CSV files (UTF-8) are written in the same folder as the JSON file, with the name of the JSON file as prefix.

| File | Content |
|---|---|
| `<file>_terminology_analysis_detailed.csv` | One row per base term/variant pair: `task`, `project`, `Terme de base`, `Terme de base ID`, `Terme de base Occurrences`, `base_is_compound`, `Variant`, `Variant ID`, `Variant Occurrences`, `variant_is_compound`, `Relation Label`, `direction`. Base terms without any variant are also listed, with an empty `Variant` column. |
| `<file>_all_terms_with_counts.csv` | Every annotated term (base terms and variants): `projects`, `tasks`, `final_text`, `final_id`, `label_type`, `is_compound`, `occurrences`. |
| `<file>_compound_terms_list.csv` | Discontinuous terms rebuilt from `parentID` (same columns as above). Written only if the text contains such terms. |

In these files, `compound` refers to terms rebuilt from several annotated segments (discontinuous terms), not to multiword terms in general.

## How it works

1. **Reading** the JSON export (`read_json_file`).
2. **Extracting** the annotations and their metadata (`extract_annotations`).
3. **Separating** text segments and relations (`split_labels_relations`).
4. **Building** a table of labels with their ID, text, label type, `parentID` and positions (`build_labels_df`).
5. **Rebuilding discontinuous terms**: segments attached to the same parent are grouped, ordered by their position in the text and joined (`find_compound_groups`).
6. **Building** a terminological registry with the occurrences of each term (`build_term_registry`).
7. **Pairing** base terms and variants (`build_relations_df`, `process_relationships`). The base term and the variant are identified by their labels, so a relation drawn in either direction is paired correctly. The direction recorded in the export is kept in the `direction` column.
8. **Removing duplicates** and writing the CSV files (`deduplicate_results`, `write_outputs`). Rows with the same base term, variant, variant ID and set of variation labels are kept once, whatever the order of the labels.

`clean_term` also removes a single trailing punctuation mark that may have been included when selecting a term in the interface.

## License

CC BY 4.0

## Contact

José CORNEJO CÁRCAMO
jose.cornejo.carcamo@gmail.com
