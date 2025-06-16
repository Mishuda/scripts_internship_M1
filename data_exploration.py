import pandas as pd
import glob
import os
import traceback

import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Read all TSV files into a list of DataFrames
tsv_files = glob.glob("C:/Users/mlazar/Desktop/heatmapper_script/*.tsv")
print(f"Found {len(tsv_files)} TSV files: {[os.path.basename(f) for f in tsv_files]}\n")

dfs = []

for f in tsv_files:
    df = pd.read_csv(f, sep="\t")
    
    # Extract filename without extension for sample name
    sample_name = os.path.splitext(os.path.basename(f))[0]
    df["sample"] = sample_name
    
    # Truncate pathway_name to text before first comma (if it exists)
    if 'pathway_name' in df.columns:
        df['pathway_name'] = df['pathway_name'].apply(lambda x: str(x).split(',')[0].strip())
    
    dfs.append(df)
    print(f"Successfully loaded: {f}")

# Combine all dataframes
combined_df = pd.concat(dfs, ignore_index=True)

# Create a unique identifier from module accession and pathway name
combined_df['module_info'] = combined_df['module'] + ': ' + combined_df['pathway_name']

# Create the pivot table for the heatmap
# Rows: module_info, Columns: sample names, Values: completeness
pivot_df = combined_df.pivot_table(
    index='module_info',
    columns='sample',
    values='completeness'
)

# Sort rows by mean completeness values (optional)
pivot_df = pivot_df.reindex(pivot_df.mean(axis=1).sort_values(ascending=False).index)

# Create the heatmap
plt.figure(figsize=(12, len(pivot_df) * 0.3))  # Adjust height based on number of modules

# Create heatmap with color scale from 0 to 100%
heatmap = sns.heatmap(
    pivot_df, 
    annot=True,            # Show values in cells
    fmt=".1f",             # Format numbers to 1 decimal place
    cmap="YlGnBu",         # Blue-green-yellow colormap
    linewidths=0.5,        # Add grid lines
    vmin=0,                # Minimum value (0%)
    vmax=100,              # Maximum value (100%)
    cbar_kws={'label': 'Completeness (%)'}  # Add label to color bar
)

plt.title("Module Completeness Across Samples")
plt.ylabel("Modules")
plt.xlabel("Samples")

# Rotate x-axis labels if needed
plt.xticks(rotation=45, ha='right')

# Adjust layout to prevent clipping of labels
plt.tight_layout()

# Save the figure
plt.savefig("module_completeness_heatmap.png", dpi=300, bbox_inches='tight')
plt.show()

# Save the processed data
pivot_df.to_csv("heatmap_data.csv")
print("\nHeatmap created and saved as 'module_completeness_heatmap.png'")
print("Heatmap data saved as 'heatmap_data.csv'")