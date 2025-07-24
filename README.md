# Module Completeness Heatmap Generator

A Python script for generating publication-ready heatmaps of metabolic pathway completeness across archaeal samples, specifically designed for analyzing KEGG module data from the [kegg-pathways-completeness-tool](https://github.com/EBI-Metagenomics/kegg-pathways-completeness-tool).

## Overview

This tool processes TSV files containing module completeness data and creates spatially grouped heatmaps organized by functional categories (Amino acids, Carbohydrates, Nucleotides, Lipids, Energy, Vitamins/Cofactors, and Archaea-specific metabolism).

## Features

- **Automated module grouping** by metabolic function
- **Publication-ready output** in SVG and PDF formats
- **A4-optimized layout** with dynamic font sizing
- **Timestamped output folders** for organized results
- **Transposable heatmaps** (samples vs modules orientation)
- **Customizable completeness thresholds** (default: 50%)
- **Missing value handling** (gray for <50% completeness)

## Sample Data

The repository includes metabolic completeness data for 6 archaeal samples:
- **Oleothermococcus** - hyperthermophilic archaeon
- **Palaeococcus** - hyperthermophilic archaeon  
- **Piezococcus** - piezophilic hyperthermophile
- **Pyrococcus** - hyperthermophilic archaeon
- **Thermococcales** - order-level representative
- **Thermococcus** - hyperthermophilic archaeon

## Key Improvements Made

### 1. Enhanced Module Coverage
Originally included 56 modules, now expanded to **63 modules** by adding missing pathways with ≥50% completeness:

**Added to Amino acids group:**
- `M00034`: Methionine salvage pathway
- `M00621`: Glycine cleavage system  
- `M00978`: Ornithine-ammonia cycle

**Added to Lipids group:**
- `M00098`: Acylglycerol degradation

**Added to Vitamins, Cofactors and polyketides group:**
- `M00122`: Cobalamin biosynthesis (58.33% in Thermococcales)
- `M00127`: Thiamine biosynthesis (prokaryotes, AIR + DXP/tyrosine)
- `M00895`: Thiamine biosynthesis (prokaryotes, AIR + DXP/glycine)

### 2. Organized Output Structure
All outputs are now saved in timestamped folders (`outputs/output_YYYYMMDD_HHMMSS/`) containing:
- **SVG heatmap** (vector format for publications)
- **PDF heatmap** (print-ready format)  
- **CSV data** (raw completeness values)

## Usage

### Basic Usage
```bash
python improved_streamlined.py input_data/*.tsv --groups pathway_groups.json --transpose
```

### Advanced Options
```bash
python improved_streamlined.py input_data/*.tsv \
    --groups pathway_groups.json \
    --transpose \
    --min-completeness 60 \
    --title "Archaeal Metabolic Completeness" \
    --output custom_heatmap.svg
```

### Command Line Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `files` | TSV files to process (one per sample) | Required |
| `--groups` | Pathway groups file (JSON format) | None |
| `--transpose` | Transpose heatmap (samples on Y-axis) | False |
| `--min-completeness` | Minimum completeness threshold (%) | 50 |
| `--output` | Output filename | `module_completeness_heatmap.svg` |
| `--names` | Custom sample names | Use filenames |
| `--figsize` | Figure size in inches (width height) | Auto-calculated |
| `--title` | Custom heatmap title | Auto-generated |

## File Structure

```
heatmapper_script/
├── improved_streamlined.py    # Main script
├── pathway_groups.json        # Module functional groupings  
├── input_data/               # Sample TSV files
│   ├── Oleothermococcus.tsv
│   ├── Palaeococcus.tsv
│   ├── Piezococcus.tsv
│   ├── Pyrococcus.tsv
│   ├── Thermococcales.tsv
│   └── Thermococcus.tsv
├── outputs/                  # Generated heatmaps (gitignored)
└── README.md                # This file
```

## Input Data Format

TSV files should contain three columns:
- `module_accession`: KEGG module ID (e.g., M00001)
- `completeness`: Percentage completeness (0-100)
- `pathway_name`: Module description

Example:
```tsv
module_accession	completeness	pathway_name
M00001	100.0	Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate
M00002	100.0	Glycolysis, core module involving three-carbon compounds
M00122	58.33	Cobalamin biosynthesis, cobyrinate a,c-diamide => cobalamin
```

## Pathway Groups Configuration

The `pathway_groups.json` file organizes modules into functional categories:

- **Amino acids** (31 modules): Amino acid biosynthesis and metabolism
- **Carbohydrates** (21 modules): Central carbon metabolism  
- **Nucleotides** (13 modules): Purine/pyrimidine biosynthesis
- **Lipids** (9 modules): Lipid biosynthesis and metabolism
- **Energy** (2 modules): ATP synthesis and electron transport
- **Vitamins, Cofactors and polyketides** (23 modules): Cofactor biosynthesis
- **Archaea-specific metabolism** (5 modules): Archaea-unique pathways

## Dependencies

```bash
pip install pandas matplotlib seaborn numpy
```

## Output Examples

The script generates:
1. **Heatmap visualization** showing completeness patterns across samples
2. **Spatial grouping** with clear separation between functional categories
3. **Color-coded completeness** (purple gradient for 50-100%, gray for <50%)
4. **Optimized layout** for A4 printing and publication

## Development History

- **Initial version**: Basic heatmap generation
- **Enhanced grouping**: Added spatial separation of functional groups  
- **Missing modules**: Identified and added 7 modules with ≥50% completeness
- **Output organization**: Implemented timestamped folders
- **Publication ready**: A4-optimized layout with proper fonts

## Citation

If you use this tool in your research, please cite the original KEGG pathways completeness tool:
[EBI-Metagenomics/kegg-pathways-completeness-tool](https://github.com/EBI-Metagenomics/kegg-pathways-completeness-tool)

## Author

Developed during M1 internship at IBiTec-S, CEA Saclay
