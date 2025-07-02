#!/usr/bin/env python3
"""
Module Completeness Heatmap Generator

Creates heatmaps for comparing metabolic pathway completeness across samples

This script processes TSV files containing module completeness data where the first 3 columns are: module_accession, completeness, pathway_name;

The tool used to generate these TSV files is https://github.com/EBI-Metagenomics/kegg-pathways-completeness-tool

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import json
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from datetime import datetime

import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = 'none'  # <-- disables path conversion


def create_output_folder():
    """Create a timestamped output folder for this run"""
    # Get current date and time
    now = datetime.now()
    date_str = now.strftime("%Y%m%d")
    time_str = now.strftime("%H%M%S")
    
    # Create base output directory name
    base_name = f"output_{date_str}_{time_str}"
    output_dir = Path("outputs") / base_name
    
    # Create the directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Created output folder: {output_dir}")
    return output_dir

def load_tsv_files(file_paths, sample_names=None):
    """Load multiple TSV files from standard tool output"""
    print(f"Loading {len(file_paths)} TSV files...")
    
    all_data = []
    for i, file_path in enumerate(file_paths):
        # Extract sample name from filename or use provided name
        sample_name = sample_names[i] if sample_names and i < len(sample_names) else Path(file_path).stem
        
        # Load the data - assume standard structure
        df = pd.read_csv(file_path, sep="\t")
        
        # Add sample identifier
        df['sample'] = sample_name
        
        # Truncate pathway_name to text before first comma
        if 'pathway_name' in df.columns:
            df['pathway_name'] = df['pathway_name'].apply(lambda x: str(x).split(',')[0].strip())
        
        all_data.append(df)
        print(f"  Loaded {len(df)} modules for {sample_name}")
    
    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    print(f"Combined data: {len(combined_df)} rows across {combined_df['sample'].nunique()} samples")
    
    return combined_df

def load_pathway_groups(group_file):
    """Load pathway grouping information from JSON file"""
    if not group_file:
        return None
        
    try:
        with open(group_file, 'r') as f:
            data = json.load(f)
            
        # Handle different JSON formats
        if 'groups' in data:
            return data['groups']
        elif isinstance(data, dict) and all(isinstance(v, list) for v in data.values()):
            return data
            
        print(f"Warning: Unrecognized format in {group_file}")
        return None
        
    except Exception as e:
        print(f"Error loading pathway groups: {e}")
        return None

def create_module_info_column(df):
    """Create a combined module_info field for display"""
    df['module_info'] = df['module_accession'] + ': ' + df['pathway_name']
    return df

def apply_module_groups(df, group_file):
    """Assign group to each module based on JSON groupings. Exclude ungrouped modules."""
    groups = load_pathway_groups(group_file)
    if not groups:
        df['module_group'] = 'Ungrouped'
        return df
    module_to_group = {}
    for group, modules in groups.items():
        for m in modules:
            module_to_group[m] = group
    # Only keep rows where module_accession is in the group mapping
    df = df[df['module_accession'].isin(module_to_group.keys())].copy()
    df['module_group'] = df['module_accession'].map(module_to_group)
    return df

def create_pivot_table(df, min_completeness=50, group_order=None):
    """Create a pivot table from the combined data, ordered by group then mean completeness."""
    # Filter to include only modules that are at least min_completeness% complete in at least one sample!
    modules_to_include = []
    for module in df['module_accession'].unique():
        max_completeness = df[df['module_accession'] == module]['completeness'].max()
        if max_completeness >= min_completeness:
            modules_to_include.append(module)
            
    # Filter the dataframe
    filtered_df = df[df['module_accession'].isin(modules_to_include)]
    
    # Create the pivot table with only module_info as index (no group names shown)
    pivot_df = filtered_df.pivot_table(
        index='module_info',
        columns='sample',
        values='completeness'
    )
    
    # Sort by group order but don't show group names
    if group_order is None:
        group_order = filtered_df['module_group'].drop_duplicates().tolist()
    
    # Create a temporary dataframe with both group and module info for sorting
    temp_df = filtered_df[['module_group', 'module_info']].drop_duplicates()
    temp_df['group_order'] = temp_df['module_group'].apply(
        lambda x: group_order.index(x) if x in group_order else len(group_order)
    )
    temp_df = temp_df.sort_values(['group_order', 'module_info'])
    
    # Reorder the pivot table based on the sorted module order
    ordered_modules = temp_df['module_info'].tolist()
    pivot_df = pivot_df.reindex([m for m in ordered_modules if m in pivot_df.index])
    
    print(f"Created pivot table with {len(pivot_df)} modules × {len(pivot_df.columns)} samples")
    return pivot_df

def create_heatmap(pivot_df, output_file=None, figsize=None, title=None, transpose=False):
    """Create and display a heatmap from the pivot table with custom colors and gray for <50%."""
    # Transpose if requested (useful when modules >> samples)
    if transpose:
        pivot_df = pivot_df.T
        print(f"Transposed heatmap: {pivot_df.shape[0]} samples × {pivot_df.shape[1]} modules")
    # Mask values below 50
    masked_data = pivot_df.mask(pivot_df < 50)
    # Calculate appropriate figure size if not specified
    if figsize is None:
        n_modules = len(pivot_df)
        n_samples = len(pivot_df.columns)
        
        # More balanced calculation that prevents extremely tall heatmaps
        module_factor = min(0.25, 10/n_modules)  # Limit how much height per module
        sample_factor = 0.8  # Width per sample
        
        fig_width = max(10, n_samples * sample_factor)
        fig_height = max(8, min(30, n_modules * module_factor))  # Cap height at 30 inches
        
        # Ensure aspect ratio isn't too extreme
        aspect_ratio = fig_height / fig_width
        if aspect_ratio > 2.5:  # If height is more than 2.5x width
            fig_height = fig_width * 2.5
            
        figsize = (fig_width, fig_height)
        print(f"Auto-sized figure: {fig_width:.1f}\" × {fig_height:.1f}\"")
    
    # Create figure
    plt.figure(figsize=figsize)
    
    # Custom color scheme
    color_scheme = [ '#B39DDB', "#8D68CE", "#614D9B", '#311B92']
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", color_scheme)
    custom_cmap.set_bad(color='#BDBDBD')  # Gray for masked values
    ax = sns.heatmap(
        masked_data,
        annot=True,
        fmt=".1f",
        cmap=custom_cmap,
        linewidths=0.5,
        vmin=50,  # Only color values 50–100
        vmax=100,
        cbar_kws={'label': 'Completeness (%)'}
    )
    # Add a horizontal line at 50% on the colorbar
    cbar = ax.collections[0].colorbar
    cbar.ax.hlines(50, *cbar.ax.get_xlim(), colors='black', linewidth=2, linestyles='--')
    cbar.ax.text(0.5, 50, '50% cutoff', color='black', ha='center', va='bottom', fontsize=10, weight='bold', rotation=0, backgroundcolor='white')
    
    # Add title and labels
    plt.title(title or "Module Completeness Across Samples", fontsize=14, pad=20)
    # Set axis labels depending on transpose
    if transpose:
        plt.xlabel("Modules", fontsize=12)
        plt.ylabel("Samples", fontsize=12)
    else:
        plt.ylabel("Modules", fontsize=12)
        plt.xlabel("Samples", fontsize=12)
    
    # Rotate x-axis labels
    plt.xticks(rotation=45, ha='right')
    
    # Adjust layout
    if len(pivot_df) > 15 or len(pivot_df.columns) > 8:
        plt.subplots_adjust(left=0.2)  # More space for y-labels
    
    # Save if output file specified
    if output_file:
        plt.savefig(output_file, dpi=300, format='svg', bbox_inches='tight')
        print(f"Heatmap saved to {output_file} (SVG format)")
        # Also save as PDF
        pdf_file = str(Path(output_file).with_suffix('.pdf'))
        plt.savefig(pdf_file, dpi=300, format='pdf', bbox_inches='tight')
        print(f"Heatmap also saved to {pdf_file} (PDF format)")
    plt.show()
    return plt.gcf()

def main():
    parser = argparse.ArgumentParser(description="Generate module completeness heatmap from TSV files")
    
    parser.add_argument('files', nargs='+', help='TSV files to process (one per sample)')
    parser.add_argument('--names', nargs='+', help='Custom names for samples (optional)')
    parser.add_argument('--output', '-o', default='module_completeness_heatmap.svg', 
                        help='Output filename (default: module_completeness_heatmap.svg)')
    parser.add_argument('--min-completeness', type=float, default=50, 
                        help='Minimum completeness threshold for inclusion (default: 50%%)')
    parser.add_argument('--groups', help='Pathway groups file (JSON format)')
    parser.add_argument('--figsize', nargs=2, type=int, help='Figure size in inches (width height)')
    parser.add_argument('--title', help='Custom title for the heatmap')
    parser.add_argument('--transpose', action='store_true',
                        help='Transpose the heatmap (put samples on Y-axis, modules on X-axis)')
    
    args = parser.parse_args()
    
    # Process the data
    try:
        # 0. Create timestamped output folder
        output_dir = create_output_folder()
        
        # 1. Load the data from TSV files
        combined_df = load_tsv_files(args.files, args.names)
        
        # 2. Create module info column
        combined_df = create_module_info_column(combined_df)
        
        # 2b. Assign groups if provided
        if args.groups:
            combined_df = apply_module_groups(combined_df, args.groups)
            group_order = list(load_pathway_groups(args.groups).keys())
        else:
            group_order = None
        
        # 3. Create pivot table
        pivot_df = create_pivot_table(combined_df, min_completeness=args.min_completeness, group_order=group_order)
        
        # 4. Create and show heatmap
        figsize = tuple(args.figsize) if args.figsize else None
        # Save output files in the timestamped folder
        output_file = output_dir / args.output
        create_heatmap(pivot_df, output_file, figsize, args.title, transpose=args.transpose)
        
        # 5. Save the data
        data_file = output_dir / Path(args.output).with_suffix('.csv')
        pivot_df.to_csv(data_file)
        print(f"Heatmap data saved to {data_file}")
        
        print(f"Heatmap generation completed successfully!")
        print(f"All files saved in: {output_dir}")
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()