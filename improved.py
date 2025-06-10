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
        self.group_labels = None    # Custom group labels        # Enhanced color schemes
        self.color_schemes = {
            'blue_pink_violet': ['#F8BBD9', '#4FC3F7', '#2196F3', '#673AB7', '#4A148C'],
        }
        
    def _create_colormap(self, scheme_name: str = 'blue_pink_violet') -> LinearSegmentedColormap:
        """Create custom colormap from predefined schemes"""
        if scheme_name not in self.color_schemes:
            scheme_name = 'blue_pink_violet'
        
        colors = self.color_schemes[scheme_name]
        return LinearSegmentedColormap.from_list(f'{scheme_name}_custom', colors, N=256)
        
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
                print(f"✅ Loaded {len(df)} pathways for {genus_name}")
                
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
                    if 'labels' in data:
                        self.group_labels = data['labels']
                    else:
                        self.pathway_groups = data
                        
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
            for group_name, pathway_list in self.pathway_groups.items():
                if pathway_id in pathway_list or pathway_name in pathway_list:
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
        
    def _truncate_name(self, name: str, max_length: int = 30) -> str:
        """Intelligently truncate pathway names for display"""
        if len(name) <= max_length:
            return name
            
        # Try to truncate at meaningful separators
        separators = [' => ', ', ', '; ', ' / ', ' - ']
        for sep in separators:
            if sep in name:
                parts = name.split(sep)
                if len(parts[0]) <= max_length:
                    return parts[0].strip()
        
        # Truncate at word boundary
        words = name.split()
        truncated = ''
        for word in words:
            if len(truncated + ' ' + word) <= max_length - 3:
                truncated += (' ' + word) if truncated else word
            else:
                break
        
        return (truncated + '...') if truncated else name[:max_length-3] + '...'
        
    def prepare_data(self, group_by_class: bool = True) -> Tuple[pd.DataFrame, List[str], Dict[str, List[int]]]:
        """Prepare data matrix for heatmap with manual grouping"""
        if not self.data:
            raise ValueError("No data loaded. Please load TSV files first.")
            
        # Collect all pathways that meet threshold
        pathway_info = {}
        pathway_groups = defaultdict(list)
        
        for genus_name, df in self.data.items():
            for _, row in df.iterrows():
                completeness = float(row.get('completeness', 0))
                pathway_id = row.get('module_accession', '')
                pathway_name = row.get('pathway_name', pathway_id)
                
                if completeness >= self.completeness_threshold:
                    if pathway_id not in pathway_info:
                        pathway_info[pathway_id] = {
                            'name': pathway_name,
                            'short_name': self._truncate_name(pathway_name),
                            'group': self._get_pathway_group(pathway_id, pathway_name),
                            'max_completeness': completeness
                        }
                        pathway_groups[pathway_info[pathway_id]['group']].append(pathway_id)
                    else:
                        # Update max completeness
                        pathway_info[pathway_id]['max_completeness'] = max(
                            pathway_info[pathway_id]['max_completeness'], 
                            completeness
                        )
        
        # Sort pathways within groups by completeness (descending)
        for group in pathway_groups:
            pathway_groups[group].sort(
                key=lambda x: pathway_info[x]['max_completeness'], 
                reverse=True
            )
        
        # Create ordered pathway list
        if self.pathway_order:
            # Use custom order
            pathways_list = []
            used_pathways = set()
            
            for pathway_ref in self.pathway_order:
                # Match by ID or name
                matched_pathway = None
                for pathway_id in pathway_info:
                    if (pathway_id == pathway_ref or 
                        pathway_info[pathway_id]['name'] == pathway_ref or
                        pathway_info[pathway_id]['short_name'] == pathway_ref):
                        matched_pathway = pathway_id
                        break
                
                if matched_pathway and matched_pathway not in used_pathways:
                    pathways_list.append(matched_pathway)
                    used_pathways.add(matched_pathway)
            
            # Add remaining pathways
            remaining = set(pathway_info.keys()) - used_pathways
            if group_by_class and self.pathway_groups:
                # Add by group order
                group_order = list(self.pathway_groups.keys()) if self.pathway_groups else sorted(pathway_groups.keys())
                for group in group_order:
                    for pathway_id in pathway_groups[group]:
                        if pathway_id in remaining:
                            pathways_list.append(pathway_id)
                            remaining.remove(pathway_id)
            
            # Add any remaining pathways
            pathways_list.extend(sorted(remaining, key=lambda x: pathway_info[x]['short_name']))
            
        else:
            # Default ordering by groups
            pathways_list = []
            if group_by_class:
                group_order = list(self.pathway_groups.keys()) if self.pathway_groups else sorted(pathway_groups.keys())
                for group in group_order:
                    pathways_list.extend(pathway_groups[group])
            else:
                pathways_list = sorted(pathway_info.keys(), 
                    key=lambda x: pathway_info[x]['short_name'])
        
        # Set genus order
        genera_list = self.genus_order or sorted(self.data.keys())
        
        # Create data matrix
        matrix_data = []
        pathway_labels = []
        
        for pathway_id in pathways_list:
            row_data = []
            for genus in genera_list:
                if genus in self.data:
                    genus_df = self.data[genus]
                    pathway_row = genus_df[genus_df['module_accession'] == pathway_id]
                    if not pathway_row.empty:
                        completeness = float(pathway_row.iloc[0]['completeness'])
                        row_data.append(completeness if completeness >= self.completeness_threshold else np.nan)
                    else:
                        row_data.append(np.nan)
                else:
                    row_data.append(np.nan)
            
            matrix_data.append(row_data)
            pathway_labels.append(pathway_info[pathway_id]['short_name'])
        
        # Create DataFrame
        heatmap_df = pd.DataFrame(matrix_data, index=pathway_labels, columns=genera_list)
        
        # Create group positions for separators
        group_positions = {}
        if group_by_class:
            current_group = None
            for i, pathway_id in enumerate(pathways_list):
                pathway_group = pathway_info[pathway_id]['group']
                if pathway_group != current_group:
                    if current_group is not None:
                        if current_group not in group_positions:
                            group_positions[current_group] = []
                        group_positions[pathway_group] = [i]
                    else:
                        group_positions[pathway_group] = [i]
                    current_group = pathway_group
        
        return heatmap_df, pathway_labels, group_positions
        
    def create_heatmap(self, 
                      figsize: Tuple[int, int] = (12, 16),
                      color_scheme: str = 'blue_pink_violet',
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
        
        # Set up the plot with extra space for group labels
        fig_width, fig_height = figsize
        if show_group_labels and group_by_class:
            fig_width += 2  # Extra space for group labels
            
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor='white')
        
        # Create heatmap
        cmap = self._create_colormap(color_scheme)
        
        # Create mask for missing values
        mask = heatmap_df.isna()
          # Plot heatmap
        im = ax.imshow(heatmap_df.values, 
                      cmap=cmap, 
                      aspect='auto',
                      vmin=self.completeness_threshold, 
                      vmax=100,
                      interpolation='nearest')
        
        # Handle missing values with light gray
        if mask.any().any():
            ax.imshow(np.where(mask.values, 1, np.nan), 
                     cmap=ListedColormap(['#f5f5f5']), 
                     aspect='auto',
                     vmin=0, vmax=1,
                     interpolation='nearest')
        
        # Add visible white grid lines between all cells
        # Vertical lines between columns
        for i in range(len(heatmap_df.columns) - 1):
            ax.axvline(x=i + 0.5, color='white', linewidth=1.5, alpha=1.0)
        
        # Horizontal lines between rows
        for i in range(len(heatmap_df.index) - 1):
            ax.axhline(y=i + 0.5, color='white', linewidth=1.5, alpha=1.0)
        
        # Set ticks and labels
        ax.set_xticks(range(len(heatmap_df.columns)))
        ax.set_xticklabels(heatmap_df.columns, rotation=45, ha='right', fontsize=font_size['xticks'])
        
        ax.set_yticks(range(len(heatmap_df.index)))
        ax.set_yticklabels(heatmap_df.index, fontsize=font_size['yticks'])
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
        
        title = f'Metabolic Pathway Completeness Comparison\n(Threshold: ≥{self.completeness_threshold}%)'
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
            
            print(f"💾 Heatmap saved to {save_path} (optimized for Inkscape)")
        
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
        
        print(f"📝 Exported {len(pathways)} pathways to {output_file}")
        
    def print_summary(self):
        """Print comprehensive summary statistics"""
        if not self.data:
            print("❌ No data loaded")
            return
            
        print(f"\n{'='*60}")
        print(" METABOLIC PATHWAY ANALYSIS SUMMARY")
        print(f"{'='*60}")
        print(f" Number of genera: {len(self.data)}")
        print(f" Completeness threshold: {self.completeness_threshold}%")
        print(f" Available color schemes: {list(self.color_schemes.keys())}")
        
        if self.pathway_groups:
            print(f"📋 Manual pathway groups: {len(self.pathway_groups)}")
            for group, pathways in self.pathway_groups.items():
                display_label = self.group_labels.get(group, group) if self.group_labels else group
                print(f"  • {display_label}: {len(pathways)} pathways")
        
        print("\nPer-genus statistics:")
        all_pathways = set()
        
        for genus_name, df in self.data.items():
            total_pathways = len(df)
            above_threshold = sum(1 for _, row in df.iterrows() 
                                if float(row.get('completeness', 0)) >= self.completeness_threshold)
            avg_completeness = df['completeness'].mean()
            
            print(f"  • {genus_name}: {total_pathways} total, {above_threshold} above threshold "
                  f"(avg: {avg_completeness:.1f}%)")
            
            all_pathways.update(df['module_accession'].tolist())
            
        print(f"\n Total unique pathways: {len(all_pathways)}")
        
        if self.genus_order:
            print(f" Genus order: {' → '.join(self.genus_order)}")
        
        if self.pathway_order:
            print(f" Custom pathway order: {len(self.pathway_order)} pathways specified")

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
    parser.add_argument('--groups', help='Pathway groups file (JSON or text format)')
    parser.add_argument('--export-pathways', help='Export pathway list to file for manual grouping')
    parser.add_argument('--output', '-o', default='metabolic_heatmap.svg', 
                       help='Output filename (default: metabolic_heatmap.svg)')
    parser.add_argument('--figsize', nargs=2, type=int, default=[12, 16], 
                       help='Figure size in inches (width height)')
    parser.add_argument('--scheme', default='blue_pink_violet', 
                       choices=['blue_pink_violet', 'viridis', 'plasma', 'cool_warm', 'nature', 'magma'],
                       help='Color scheme (default: blue_pink_violet)')
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
    
    # Print summary
    generator.print_summary()
    
    # Create heatmap
    try:
        fig = generator.create_heatmap(
            figsize=tuple(args.figsize),
            color_scheme=args.scheme,
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