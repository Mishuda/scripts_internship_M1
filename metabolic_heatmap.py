#!/usr/bin/env python3
"""
Multi-Genera Metabolic Pathway Heatmap Generator
Creates publication-ready heatmaps for comparing metabolic pathway completeness across genera
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
from typing import List, Dict, Tuple, Optional
import warnings
from matplotlib.colors import LinearSegmentedColormap
warnings.filterwarnings('ignore')

class MetabolicHeatmapGenerator:
    def __init__(self):
        self.data = {}  # genus_name: dataframe
        self.completeness_threshold = 75.0
        self.genus_order = None
          # Create custom blue-pink-violet colormap
        self.custom_colormap = self._create_custom_colormap()
        
    def _create_custom_colormap(self):
        """Create a custom blue-pink-violet colormap with stark transitions"""
        colors = [
            '#0D1B2A',  # Very dark navy blue (low completeness)
            '#1B4F93',  # Deep blue
            '#8E44AD',  # Rich purple/violet  
            '#E91E63',  # Bright magenta/pink
            '#6A1B9A'   # Deep violet (high completeness)
        ]
        
        # Create the colormap with fewer intermediate steps for starker transitions
        custom_cmap = LinearSegmentedColormap.from_list(
            'blue_pink_violet', colors, N=128
        )
        return custom_cmap
        
    def load_files(self, file_paths: List[str], genus_names: Optional[List[str]] = None):
        """Load multiple TSV files, one per genus"""
        self.data = {}
        
        for i, file_path in enumerate(file_paths):
            path = Path(file_path)
            if not path.exists():
                print(f"Warning: File {file_path} not found, skipping...")
                continue
                
            # Extract genus name
            if genus_names and i < len(genus_names):
                genus_name = genus_names[i]
            else:
                genus_name = path.stem  # filename without extension
                
            try:
                df = pd.read_csv(file_path, sep='\t')
                self.data[genus_name] = df
                print(f"Loaded {len(df)} pathways for {genus_name}")
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
                
        print(f"Successfully loaded {len(self.data)} genera")
        
    def set_completeness_threshold(self, threshold: float):
        """Set the completeness threshold (0-100)"""
        self.completeness_threshold = max(0, min(100, threshold))
        print(f"Completeness threshold set to {self.completeness_threshold}%")
        
    def set_genus_order(self, order: List[str]):
        """Set custom order for genera on X-axis"""
        available_genera = set(self.data.keys())
        valid_order = [g for g in order if g in available_genera]
        missing_genera = available_genera - set(valid_order)
        
        if missing_genera:
            print(f"Warning: Adding missing genera to end of order: {missing_genera}")
            valid_order.extend(sorted(missing_genera))
            
        self.genus_order = valid_order
        print(f"Genus order set: {' -> '.join(self.genus_order)}")
        
    def prepare_data(self, group_by_class: bool = True, sort_by: str = 'pathway') -> Tuple[pd.DataFrame, List[str], List[str]]:
        """Prepare data matrix for heatmap"""
        if not self.data:
            raise ValueError("No data loaded. Please load TSV files first.")
            
        # Get all modules that meet threshold in at least one genus
        all_modules = set()
        module_info = {}
        
        for genus_name, df in self.data.items():
            for _, row in df.iterrows():
                completeness = float(row.get('completeness', 0))
                module_id = row.get('module_accession', '')
                
                if completeness >= self.completeness_threshold:
                    all_modules.add(module_id)
                    
                if module_id not in module_info:
                    module_info[module_id] = {
                        'name': row.get('pathway_name', module_id),
                        'class': row.get('pathway_class', 'Unknown'),
                        'short_name': self._truncate_name(row.get('pathway_name', module_id))
                    }
        
        # Sort modules
        modules_list = list(all_modules)
        if group_by_class:
            modules_list.sort(key=lambda x: (
                module_info[x]['class'],
                module_info[x]['short_name']
            ))
        else:
            modules_list.sort(key=lambda x: module_info[x]['short_name'])
            
        # Set genus order
        if self.genus_order is None:
            genera_list = sorted(self.data.keys())
        else:
            genera_list = self.genus_order
            
        # Create data matrix
        matrix = []
        for module_id in modules_list:
            row = []
            for genus in genera_list:
                if genus in self.data:
                    genus_df = self.data[genus]
                    module_row = genus_df[genus_df['module_accession'] == module_id]
                    if not module_row.empty:
                        completeness = float(module_row.iloc[0]['completeness'])
                        row.append(completeness if completeness >= self.completeness_threshold else np.nan)
                    else:
                        row.append(np.nan)
                else:
                    row.append(np.nan)
            matrix.append(row)
            
        # Create DataFrame
        pathway_labels = [module_info[m]['short_name'] for m in modules_list]
        heatmap_df = pd.DataFrame(matrix, index=pathway_labels, columns=genera_list)
        
        # Create class labels for grouping
        class_labels = []
        current_class = None
        for module_id in modules_list:
            module_class = module_info[module_id]['class']
            if module_class != current_class:
                class_labels.append(module_class.split(';')[-1].strip() if ';' in module_class else module_class)
                current_class = module_class
            else:
                class_labels.append('')
                
        return heatmap_df, pathway_labels, class_labels
        
    def _truncate_name(self, name: str, max_length: int = 50) -> str:
        """Truncate pathway names for display"""
        if len(name) <= max_length:
            return name
        # Try to truncate at comma or semicolon
        for sep in [',', ';', ' =>']:
            if sep in name:
                truncated = name.split(sep)[0].strip()
                if len(truncated) <= max_length:
                    return truncated        # Fallback to character limit
        return name[:max_length-3] + '...'
        
    def create_heatmap(self, 
                      figsize: Tuple[int, int] = (12, 16),
                      cmap: str = 'blue_pink_violet',
                      group_by_class: bool = True,
                      show_class_labels: bool = True,
                      save_path: Optional[str] = None,
                      dpi: int = 300,
                      format: str = 'png') -> plt.Figure:
        """Create the metabolic pathway heatmap"""
        
        # Prepare data
        heatmap_df, pathway_labels, class_labels = self.prepare_data(group_by_class=group_by_class)
        
        if heatmap_df.empty:
            raise ValueError("No pathways meet the completeness threshold")
            
        print(f"Creating heatmap with {len(heatmap_df)} pathways and {len(heatmap_df.columns)} genera")
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Use custom colormap if specified, otherwise use the provided one
        if cmap == 'blue_pink_violet':
            colormap = self.custom_colormap
        else:
            colormap = cmap
        
        # Create heatmap
        mask = heatmap_df.isna()
        sns.heatmap(heatmap_df, 
                   mask=mask,
                   cmap=colormap,
                   center=50,
                   vmin=self.completeness_threshold,
                   vmax=100,
                   cbar_kws={'label': 'Pathway Completeness (%)', 'shrink': 0.8},
                   xticklabels=True,
                   yticklabels=True,
                   square=False,
                   linewidths=0.5,
                   linecolor='white',
                   ax=ax)
        
        # Customize appearance
        ax.set_xlabel('Genera', fontsize=14, fontweight='bold')
        ax.set_ylabel('Metabolic Pathways', fontsize=14, fontweight='bold')
        ax.set_title(f'Metabolic Pathway Completeness Comparison\n(Threshold: ≥{self.completeness_threshold}%)', 
                    fontsize=16, fontweight='bold', pad=20)
        
        # Rotate x-axis labels
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
        plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
        
        # Add class labels if requested
        if show_class_labels and group_by_class:
            current_y = 0
            for i, label in enumerate(class_labels):
                if label:  # Non-empty class label
                    ax.text(-0.5, current_y + 0.5, label, 
                           ha='right', va='center', fontsize=9, 
                           fontweight='bold', color='blue',
                           transform=ax.get_yaxis_transform())
                    # Add separator line
                    if current_y > 0:
                        ax.axhline(y=current_y, color='gray', linestyle='-', alpha=0.3, linewidth=1)
                current_y = i + 1
        
        # Add threshold line to colorbar
        cbar = ax.collections[0].colorbar
        cbar.ax.axhline(y=self.completeness_threshold, color='red', linestyle='--', linewidth=2)
        cbar.ax.text(1.05, self.completeness_threshold, f'{self.completeness_threshold}%', 
                    transform=cbar.ax.get_yaxis_transform(), 
                    va='center', ha='left', color='red', fontweight='bold')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save if requested
        if save_path:
            plt.savefig(save_path, dpi=dpi, format=format, bbox_inches='tight')
            print(f"Heatmap saved to {save_path}")
            
        return fig
        
    def print_summary(self):
        """Print summary statistics"""
        if not self.data:
            print("No data loaded")
            return
            
        print(f"\n{'='*50}")
        print("METABOLIC PATHWAY ANALYSIS SUMMARY")
        print(f"{'='*50}")
        print(f"Number of genera: {len(self.data)}")
        print(f"Completeness threshold: {self.completeness_threshold}%")
        
        all_pathways = set()
        for genus_name, df in self.data.items():
            pathways_above_threshold = sum(1 for _, row in df.iterrows() 
                                         if float(row.get('completeness', 0)) >= self.completeness_threshold)
            print(f"  {genus_name}: {len(df)} total pathways, {pathways_above_threshold} above threshold")
            all_pathways.update(df['module_accession'].tolist())
            
        print(f"Total unique pathways across all genera: {len(all_pathways)}")
        print(f"Genus order: {' -> '.join(self.genus_order or sorted(self.data.keys()))}")

def main():
    parser = argparse.ArgumentParser(description='Generate metabolic pathway comparison heatmaps')
    parser.add_argument('files', nargs='+', help='TSV files to process (one per genus)')
    parser.add_argument('--names', nargs='+', help='Custom names for genera (optional)')
    parser.add_argument('--threshold', type=float, default=75.0, help='Completeness threshold (default: 75)')
    parser.add_argument('--order', nargs='+', help='Custom order for genera on X-axis')
    parser.add_argument('--output', '-o', default='metabolic_heatmap.png', help='Output filename')
    parser.add_argument('--figsize', nargs=2, type=int, default=[12, 16], help='Figure size (width height)')
    parser.add_argument('--cmap', default='blue_pink_violet', help='Colormap (default: blue_pink_violet)')
    parser.add_argument('--no-grouping', action='store_true', help='Disable grouping by pathway class')
    parser.add_argument('--no-class-labels', action='store_true', help='Hide pathway class labels')
    parser.add_argument('--dpi', type=int, default=300, help='Output DPI (default: 300)')
    parser.add_argument('--format', default='png', choices=['png', 'pdf', 'svg'], help='Output format')
    
    args = parser.parse_args()
    
    # Create generator
    generator = MetabolicHeatmapGenerator()
    
    # Load data
    generator.load_files(args.files, args.names)
    
    # Set parameters
    generator.set_completeness_threshold(args.threshold)
    if args.order:
        generator.set_genus_order(args.order)
    
    # Print summary
    generator.print_summary()
    
    # Create heatmap
    try:
        fig = generator.create_heatmap(
            figsize=tuple(args.figsize),
            cmap=args.cmap,
            group_by_class=not args.no_grouping,
            show_class_labels=not args.no_class_labels,
            save_path=args.output,
            dpi=args.dpi,
            format=args.format
        )
        plt.show()
    except Exception as e:
        print(f"Error creating heatmap: {e}")

if __name__ == "__main__":
    main()

# Example usage as a script:
"""
# Basic usage
python metabolic_heatmap.py genus1.tsv genus2.tsv genus3.tsv

# With custom names and threshold
python metabolic_heatmap.py file1.tsv file2.tsv file3.tsv --names Escherichia Bacillus Pseudomonas --threshold 80

# With custom order and output
python metabolic_heatmap.py *.tsv --order Bacillus Escherichia Pseudomonas --output comparison.pdf --format pdf

# Full customization
python metabolic_heatmap.py *.tsv --threshold 75 --figsize 15 20 --cmap viridis --dpi 600 --output metabolic_analysis.png
"""

# Example usage as a module:
"""
from metabolic_heatmap import MetabolicHeatmapGenerator

# Create generator
gen = MetabolicHeatmapGenerator()

# Load files
gen.load_files(['genus1.tsv', 'genus2.tsv', 'genus3.tsv'], ['Genus1', 'Genus2', 'Genus3'])

# Set parameters
gen.set_completeness_threshold(75)
gen.set_genus_order(['Genus2', 'Genus1', 'Genus3'])  # Custom X-axis order

# Create heatmap
fig = gen.create_heatmap(figsize=(14, 20), cmap='blue_pink_violet', save_path='my_heatmap.png')
"""
