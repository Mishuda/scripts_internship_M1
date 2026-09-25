# M1 Heatmap Development Repository

> **Legacy development snapshot.** The cleaned, recruiter-facing version of these utilities is maintained in [`Mishuda/heatmapper-scripts`](https://github.com/Mishuda/heatmapper-scripts).

This repository contains an earlier development version of the Python scripts I used during my M1 internship for **metabolic reconstruction and comparative genomics of Thermococcales**.

The main task was to transform KEGG module-completeness results into comparative heatmaps grouped by metabolic function and to preprocess significant KofamScan hits into KO lists.

## What is here

- `improved_streamlined.py` — comparative module-completeness heatmaps
- `kofam_output_processor.py` — extraction of significant KO identifiers from KofamScan output
- `pathway_groups.json` — functional grouping of KEGG modules
- `input_data/` — development/example inputs retained from the internship workflow

## Preferred version

For a smaller and more clearly documented version with example usage, see:

**[Metabolic Pathway Completeness Heatmaps →](https://github.com/Mishuda/heatmapper-scripts)**

This repository remains public for provenance of the M1 development work rather than as the canonical version of the tool.
