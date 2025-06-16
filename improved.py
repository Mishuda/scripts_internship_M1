#!/usr/bin/env python3
"""
Optimized Multi-Genera Metabolic Pathway Heatmap Generator
Creates publication-ready heatmaps for comparing metabolic pathway completeness across genera
with manual Y-axis customization and SVG output optimization
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from typing import List, Dict, Tuple, Optional
import warnings
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from collections import defaultdict
import json

warnings.filterwarnings('ignore')

class MetabolicHeatmapGenerator:
    def __init__(self):        
        self.data = {}  # genus_name: dataframe
        self.completeness_threshold = 75.0
        self.genus_order = None
        self.pathway_order = None
        self.pathway_groups = None  # Manual pathway groupings
        self.group_labels = None    # Custom group labels
        # Updated color scheme: light pink -> light blue -> deep violet
        self.color_scheme = ['#FFEBEE', '#E3F2FD', '#B39DDB', '#673AB7', '#4527A0', '#311B92'] 
    
    def _create_colormap(self) -> LinearSegmentedColormap:
        """Create custom colormap from the single pink-blue-violet scheme"""
        return LinearSegmentedColormap.from_list('pink_blue_violet', self.color_scheme, N=256)
        
    def load_files(self, file_paths: List[str], genus_names: Optional[List[str]] = None):
        """Load multiple TSV files, one per genus with enhanced error handling"""
        self.data = {}
        
        for i, file_path in enumerate(file_paths):
            path = Path(file_path)
            if not path.exists():
                print(f"  Warning: File {file_path} not found, skipping...")
                continue
                
            # Extract genus name
            if genus_names and i < len(genus_names):
                genus_name = genus_names[i]
            else:
                genus_name = path.stem.replace('_', ' ').title()
                
            try:
                # Try different separators
                for sep in ['\t', ',', ';']:
                    try:
                        df = pd.read_csv(file_path, sep=sep)
                        if len(df.columns) > 1:
                            break
                    except:
                        continue
                
                # Standardize column names
                df.columns = df.columns.str.lower().str.replace(' ', '_')
                
                # Map common column name variations
                column_mapping = {
                    'module_id': 'module_accession',
                    'pathway_id': 'module_accession', 
                    'accession': 'module_accession',
                    'pathway': 'pathway_name',
                    'name': 'pathway_name',
                    'description': 'pathway_name',
                    'class': 'pathway_class',
                    'category': 'pathway_class',
                    'group': 'pathway_class',
                    'completeness_percent': 'completeness',
                    'completion': 'completeness',
                    'complete': 'completeness'
                }
                
                df.rename(columns=column_mapping, inplace=True)
                
                # Ensure required columns exist
                required_cols = ['module_accession', 'pathway_name', 'completeness']
                missing_cols = [col for col in required_cols if col not in df.columns]
                if missing_cols:
                    print(f" Error: Missing columns {missing_cols} in {file_path}")
                    continue
                
                # Clean and validate data
                df['completeness'] = pd.to_numeric(df['completeness'], errors='coerce')
                df = df.dropna(subset=['completeness'])
                df['completeness'] = df['completeness'].clip(0, 100)
                
                self.data[genus_name] = df
                print(f" Loaded {len(df)} pathways for {genus_name}")
                
            except Exception as e:
                print(f" Error loading {file_path}: {e}")
                
        print(f" Successfully loaded {len(self.data)} genera")
        
    def load_pathway_groups_from_file(self, groups_file: str):
        """Load pathway groups from JSON or text file"""
        try:
            path = Path(groups_file)
            
            if path.suffix.lower() == '.json':
                # Load from JSON file
                with open(groups_file, 'r') as f:
                    data = json.load(f)
                    if 'groups' in data:
                        self.pathway_groups = data['groups']
                    # else: # This was part of the original logic structure, but could be added if direct dict support is needed
                        # if isinstance(data, dict) and all(isinstance(v, list) for v in data.values()):
                        #    self.pathway_groups = data
                    if 'labels' in data:
                        self.group_labels = data['labels']
                        
            else:
                # Load from text file format:
                # GROUP_NAME:
                #   pathway1
                #   pathway2
                # ANOTHER_GROUP:
                #   pathway3
                groups = {}
                current_group = None
                
                with open(groups_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#'):
                            continue
                        
                        if line.endswith(':'):
                            current_group = line[:-1].strip()
                            groups[current_group] = []
                        elif current_group and line.startswith((' ', '\t')):
                            pathway = line.strip()
                            groups[current_group].append(pathway)
                
                self.pathway_groups = groups
                
            print(f" Loaded pathway groups: {list(self.pathway_groups.keys())}")
            
        except Exception as e:
            print(f" Error loading pathway groups from {groups_file}: {e}")

    # anywhere it says "pathway" it means "KEGG module" actually
    def set_pathway_groups(self, groups: Dict[str, List[str]], group_labels: Optional[Dict[str, str]] = None):
        """Manually set pathway groups for Y-axis organization
        
        Args:
            groups: Dict mapping group names to lists of pathway IDs/names
            group_labels: Optional dict mapping group names to display labels
        """
        self.pathway_groups = groups
        self.group_labels = group_labels or {}
        print(f" Set pathway groups: {list(groups.keys())}")
        
    def set_completeness_threshold(self, threshold: float):
        """Set the completeness threshold (0-100)"""
        self.completeness_threshold = max(0, min(100, threshold))
        print(f" Completeness threshold set to {self.completeness_threshold}%")
        
    def set_genus_order(self, order: List[str]):
        """Set custom order for genera on X-axis"""
        available_genera = set(self.data.keys())
        valid_order = [g for g in order if g in available_genera]
        missing_genera = available_genera - set(valid_order)
        
        if missing_genera:
            print(f"  Warning: Adding missing genera to end of order: {missing_genera}")
            valid_order.extend(sorted(missing_genera))
            
        self.genus_order = valid_order
        print(f" Genus order set: {' → '.join(self.genus_order)}")
    def set_pathway_order(self, order: List[str]):
        """Set custom order for pathways on Y-axis"""
        self.pathway_order = order
        print(f" Custom pathway order set ({len(order)} pathways)")
        
    def _get_pathway_group(self, pathway_id: str, pathway_name: str) -> str:
        """Determine pathway group based on manual groupings"""
        if self.pathway_groups:
            # First try with just the module ID (for JSON format groups)
            for group_name, pathway_list in self.pathway_groups.items():
                # Extract just the module ID if it's in format "MXXXXX: Description"
                clean_id = pathway_id.split(':')[0].strip() if ':' in pathway_id else pathway_id
                if clean_id in pathway_list or pathway_id in pathway_list or pathway_name in pathway_list:
                    return group_name
          # Fallback to class from data
        for genus_data in self.data.values():
            pathway_row = genus_data[
                (genus_data['module_accession'] == pathway_id) | 
                (genus_data['pathway_name'] == pathway_name)
            ]
            if not pathway_row.empty and 'pathway_class' in pathway_row.columns:
                class_name = pathway_row.iloc[0]['pathway_class']
                # Extract last part of hierarchical class name
                if ';' in class_name:
                    return class_name.split(';')[-1].strip()
                return class_name
        
        return 'Other'
        
    def _truncate_name(self, name: str, module_id: str = '', max_length: int = 50) -> str:
        """Intelligently truncate pathway names for display with module ID"""
        # Format with module ID if provided
        if module_id:
            prefix = f"{module_id}: "
            available_length = max_length - len(prefix)
        else:
            prefix = ""
            available_length = max_length
        
        # If the full name with prefix fits, return it
        full_display = prefix + name
        if len(full_display) <= max_length:
            return full_display
            
        # Try to truncate at meaningful separators
        separators = [' => ', ', ', '; ', ' / ', ' - ']
        for sep in separators:
            if sep in name:
                parts = name.split(sep)
                if len(prefix + parts[0]) <= max_length:
                    return prefix + parts[0].strip()
        
        # Truncate at word boundary
        words = name.split()
        truncated = ''
        for word in words:
            test_length = len(prefix + truncated + (' ' + word if truncated else word) + '...')
            if test_length <= max_length:
                truncated += (' ' + word) if truncated else word
            else:
                break
        
        if truncated:
            return prefix + truncated + '...'
        else:
            # Last resort: truncate the name to fit
            return prefix + name[:available_length-3] + '...' if available_length > 3 else module_id
        
    def prepare_data(self, group_by_class: bool = True) -> Tuple[pd.DataFrame, List[str], Dict[str, List[int]]]:
        """Prepare data matrix for heatmap with manual grouping"""
        if not self.data:
            raise ValueError("No data loaded. Please load TSV files first.")        # Collect all pathways that meet threshold in at least one genus
        # Custom: Only include modules that are at least 50% complete in at least one genus
        # and respect the grouping from the JSON
        # Collect all pathways and their max completeness
        pathway_info = {}
        pathway_groups = defaultdict(list)
        all_pathways = {}  # pathway_id -> {genus: completeness, ...}
        for genus_name, df in self.data.items():
            for _, row in df.iterrows():
                completeness = float(row.get('completeness', 0))
                pathway_id = row.get('module_accession', '')
                pathway_name = row.get('pathway_name', pathway_id)
                if pathway_id not in all_pathways:
                    all_pathways[pathway_id] = {}
                all_pathways[pathway_id][genus_name] = completeness
                if pathway_id not in pathway_info:
                    pathway_info[pathway_id] = {
                        'name': pathway_name,
                        'short_name': self._truncate_name(pathway_name, pathway_id),
                        'group': self._get_pathway_group(pathway_id, pathway_name),
                        'max_completeness': completeness
                    }
                else:
                    pathway_info[pathway_id]['max_completeness'] = max(
                        pathway_info[pathway_id]['max_completeness'], completeness)
        # Only include modules that are >= 50% in at least one genus
        modules_to_include = set()
        for pathway_id, genus_completeness in all_pathways.items():
            if max(genus_completeness.values()) >= 50:
                modules_to_include.add(pathway_id)        # Respect the grouping from the JSON
        grouped_modules = []
        if self.pathway_groups:
            for group in self.pathway_groups:
                for pid in self.pathway_groups[group]:
                    # Look for matching modules in modules_to_include
                    # Check both exact match and just the module ID part
                    matching_modules = [m for m in modules_to_include if 
                                      m == pid or 
                                      m.split(':')[0].strip() == pid]
                    for module in matching_modules:
                        if module not in grouped_modules:
                            grouped_modules.append(module)
        # Add any remaining modules that meet the threshold but are not in a group
        for pid in modules_to_include:
            if pid not in grouped_modules:
                grouped_modules.append(pid)
        # Ask user for confirmation before generating SVG
        print("The following modules will be included in the heatmap (at least 50% complete in one genus):")
        for pid in grouped_modules:
            print(f"{pid}: {pathway_info[pid]['name']}")
        input("Press Enter to continue and generate the SVG, or Ctrl+C to abort...")
        # Now build the matrix as before, but only for grouped_modules
        genera_list = self.genus_order or sorted(self.data.keys())
        matrix_data = []
        pathway_labels = []
        for pathway_id in grouped_modules:
            row_data = []
            for genus in genera_list:
                if genus in self.data:
                    genus_df = self.data[genus]
                    pathway_row = genus_df[genus_df['module_accession'] == pathway_id]
                    if not pathway_row.empty:
                        completeness = float(pathway_row.iloc[0]['completeness'])
                        row_data.append(completeness)                    
                        else:
                        
                        row_data.append(np.nan)
                else:
                    row_data.append(np.nan)
            matrix_data.append(row_data)
            # Use the formatted short_name as the label but keep track of the original module ID for grouping
            pathway_labels.append(pathway_info[pathway_id]['short_name'])
        heatmap_df = pd.DataFrame(matrix_data, index=pathway_labels, columns=genera_list)
        
        # Create group positions for separators
        group_positions = {}
        if group_by_class:
            current_group = None
            for i, pathway_label in enumerate(pathway_labels):
                # Get the original module ID from the formatted label
                module_id = pathway_label.split(':', 1)[0].strip()
                pathway_group = pathway_info[module_id]['group']                if pathway_group != current_group:
                    if current_group is not None:
                        if current_group not in group_positions:
                            group_positions[current_group] = []
                        if pathway_group not in group_positions:
                            group_positions[pathway_group] = []
                        group_positions[pathway_group].append(i)
                    else:
                        if pathway_group not in group_positions:
                            group_positions[pathway_group] = []
                        group_positions[pathway_group].append(i)
                    current_group = pathway_group
        
        return heatmap_df, pathway_labels, group_positions
        
    def create_heatmap(self, 
                      figsize: Tuple[int, int] = (12, 16),
                      group_by_class: bool = True,
                      show_group_separators: bool = True,
                      show_group_labels: bool = True,
                      save_path: Optional[str] = None,
                      dpi: int = 300,
                      font_size: Dict[str, int] = None) -> plt.Figure:
        """Create optimized metabolic pathway heatmap with SVG output"""
        
        # Default font sizes
        if font_size is None:
            font_size = {
                'title': 16,
                'xlabel': 14, 
                'ylabel': 14,
                'xticks': 10,
                'yticks': 8,
                'colorbar': 10,
                'group_labels': 9
            }
          # Prepare data
        heatmap_df, pathway_labels, group_positions = self.prepare_data(group_by_class=group_by_class)
        
        if heatmap_df.empty:
            raise ValueError("No pathways meet the completeness threshold")
            
        print(f" Creating heatmap with {len(heatmap_df)} pathways and {len(heatmap_df.columns)} genera")
        
        # Calculate square cell dimensions based on number of pathways and genera
        n_pathways = len(heatmap_df)
        n_genera = len(heatmap_df.columns)
        
        # Make cells approximately square by adjusting figure size
        cell_size = 0.4  # Size of each cell in inches
        fig_width = max(8, n_genera * cell_size + 8)  # Extra space for labels and colorbar
        fig_height = max(6, n_pathways * cell_size + 4)  # Extra space for title and labels
        
        if show_group_labels and group_by_class:
            fig_width += 3  # Extra space for group labels
            
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor='white')
          # Create heatmap
        cmap = self._create_colormap()
        
        # Create mask for missing values
        mask = heatmap_df.isna()        # Plot heatmap with square aspect ratio
        im = ax.imshow(heatmap_df.values, 
                      cmap=cmap, 
                      aspect='equal',  # Make cells square
                      vmin=0,  # Show all values from 0 to 100
                      vmax=100,
                      interpolation='nearest')
          # Handle missing values with light gray
        if mask.any().any():
            ax.imshow(np.where(mask.values, 1, np.nan), 
                     cmap=ListedColormap(['#f0f0f0']), 
                     aspect='equal',  # Keep square cells
                     vmin=0, vmax=1,
                     interpolation='nearest')
        
        # Add visible white grid lines between all cells
        # Vertical lines between columns
        for i in range(len(heatmap_df.columns) - 1):
            ax.axvline(x=i + 0.5, color='white', linewidth=1.5, alpha=1.0)
        
        # Horizontal lines between rows
        for i in range(len(heatmap_df.index) - 1):
            ax.axhline(y=i + 0.5, color='white', linewidth=1.5, alpha=1.0)
          # Set ticks and labels with improved formatting
        ax.set_xticks(range(len(heatmap_df.columns)))
        ax.set_xticklabels(heatmap_df.columns, rotation=60, ha='right', fontsize=font_size['xticks'], 
                          fontweight='normal')
        
        ax.set_yticks(range(len(heatmap_df.index)))
        ax.set_yticklabels(heatmap_df.index, fontsize=font_size['yticks'], fontweight='normal')
          # Add group separators and labels
        if show_group_separators and group_by_class and group_positions:
            previous_end = 0
            
            for group_name, positions in group_positions.items():
                if positions:
                    start_pos = positions[0]
                    # Add thicker separator line for groups (more prominent than cell lines)
                    if start_pos > 0:
                        ax.axhline(y=start_pos-0.5, color='white', linewidth=4.0, alpha=1.0)
                    
                    # Add group label
                    if show_group_labels:
                        # Find the middle position of this group
                        next_group_start = min([p[0] for p in group_positions.values() if p[0] > start_pos], default=len(pathway_labels))
                        group_middle = (start_pos + next_group_start) / 2 - 0.5
                        
                        display_label = self.group_labels.get(group_name, group_name) if self.group_labels else group_name
                        
                        ax.text(len(heatmap_df.columns) + 0.2, group_middle, display_label, 
                               fontsize=font_size['group_labels'], 
                               fontweight='bold',
                               va='center', ha='left',
                               color='#2E86C1',
                               rotation=0)
        
        # Colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8, aspect=30, pad=0.02)
        cbar.set_label('Pathway Completeness (%)', fontsize=font_size['colorbar'], fontweight='bold')
        cbar.ax.tick_params(labelsize=font_size['colorbar'])
        
        # Labels and title
        ax.set_xlabel('Genera', fontsize=font_size['xlabel'], fontweight='bold')
        ax.set_ylabel('Metabolic Pathways', fontsize=font_size['ylabel'], fontweight='bold')
        
        title = 'Metabolic Pathway Completeness Comparison'
        ax.set_title(title, fontsize=font_size['title'], fontweight='bold', pad=20)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save as SVG with optimization for Inkscape
        if save_path:
            if not save_path.lower().endswith('.svg'):
                save_path = save_path.rsplit('.', 1)[0] + '.svg'
            
            # SVG-specific optimizations
            plt.rcParams['svg.fonttype'] = 'none'  # Keep text editable
            plt.rcParams['font.family'] = 'sans-serif'
            
            fig.savefig(save_path, 
                       format='svg', 
                       dpi=dpi, 
                       bbox_inches='tight',
                       facecolor='white',
                       edgecolor='none')
            
            print(f" Heatmap saved to {save_path} (optimized for Inkscape)")
        
        return fig
    
    def export_pathway_list(self, output_file: str, above_threshold_only: bool = True):
        """Export list of pathways for manual grouping"""
        if not self.data:
            print(" No data loaded")
            return
            
        pathways = set()
        for genus_name, df in self.data.items():
            for _, row in df.iterrows():
                completeness = float(row.get('completeness', 0))
                if not above_threshold_only or completeness >= self.completeness_threshold:
                    pathway_id = row.get('module_accession', '')
                    pathway_name = row.get('pathway_name', pathway_id)
                    pathways.add(f"{pathway_id}\t{pathway_name}")
        
        with open(output_file, 'w') as f:
            f.write("# Pathway list for manual grouping\n")
            f.write("# Format: MODULE_ID\tPATHWAY_NAME\n\n")
            for pathway in sorted(pathways):
                f.write(f"{pathway}\n")
        
        print(f"Exported {len(pathways)} pathways to {output_file}")
        

def main():
    parser = argparse.ArgumentParser(
        description='Generate publication-ready metabolic pathway comparison heatmaps with manual Y-axis control',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with SVG output
  python metabolic_heatmap.py genus1.tsv genus2.tsv --output analysis.svg
  
  # Export pathway list for manual grouping
  python metabolic_heatmap.py *.tsv --export-pathways pathways.txt
  
  # With manual pathway groups (JSON format)
  python metabolic_heatmap.py *.tsv --groups groups.json --scheme plasma
  
  # Custom pathway order and styling
  python metabolic_heatmap.py *.tsv --pathway-order M00001 M00002 M00003 --figsize 15 20
        """
    )
    
    parser.add_argument('files', nargs='+', help='TSV files to process (one per genus)')
    parser.add_argument('--names', nargs='+', help='Custom names for genera')
    parser.add_argument('--threshold', type=float, default=75.0, 
                       help='Completeness threshold percentage (default: 75)')
    parser.add_argument('--order', nargs='+', help='Custom order for genera on X-axis')
    parser.add_argument('--pathway-order', nargs='+', help='Custom order for pathways on Y-axis')
    parser.add_argument('--groups', help='Pathway groups file (JSON or text format)', default="c:\\Users\\mlazar\\Desktop\\heatmapper_script\\kegg_module_groups_from_image.json")
    parser.add_argument('--export-pathways', help='Export pathway list to file for manual grouping')
    parser.add_argument('--output', '-o', default='metabolic_heatmap.svg', 
                       help='Output filename (default: metabolic_heatmap.svg)')
    parser.add_argument('--figsize', nargs=2, type=int, default=[12, 16], 
                       help='Figure size in inches (width height)')
    parser.add_argument('--no-grouping', action='store_true', 
                       help='Disable grouping by pathway class')
    parser.add_argument('--no-separators', action='store_true', 
                       help='Hide group separator lines')
    parser.add_argument('--no-labels', action='store_true', 
                       help='Hide group labels')
    parser.add_argument('--dpi', type=int, default=300, 
                       help='Output resolution (default: 300)')
    parser.add_argument('--font-scale', type=float, default=1.0,
                       help='Scale factor for all fonts (default: 1.0)')
    
    args = parser.parse_args()
    
    # Create generator
    generator = MetabolicHeatmapGenerator()
    
    # Load data
    generator.load_files(args.files, args.names)
    
    # Export pathways list if requested
    if args.export_pathways:
        generator.export_pathway_list(args.export_pathways)
        return
    
    # Load pathway groups if provided
    if args.groups:
        generator.load_pathway_groups_from_file(args.groups)
    
    # Set parameters
    generator.set_completeness_threshold(args.threshold)
    if args.order:
        generator.set_genus_order(args.order)
    if args.pathway_order:
        generator.set_pathway_order(args.pathway_order)
    
    # Calculate font sizes based on scale
    font_size = {
        'title': int(16 * args.font_scale),
        'xlabel': int(14 * args.font_scale), 
        'ylabel': int(14 * args.font_scale),
        'xticks': int(10 * args.font_scale),
        'yticks': int(8 * args.font_scale),
        'colorbar': int(10 * args.font_scale),
        'group_labels': int(9 * args.font_scale)
    }
    
      # Create heatmap
    try:
        fig = generator.create_heatmap(
            figsize=tuple(args.figsize),
            group_by_class=not args.no_grouping,
            show_group_separators=not args.no_separators,
            show_group_labels=not args.no_labels,
            save_path=args.output,
            dpi=args.dpi,
            font_size=font_size
        )
        
        print("Heatmap generation completed successfully!")
        plt.show()
        
    except Exception as e:
        print(f"Error creating heatmap: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()