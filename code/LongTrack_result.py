import sys
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Load the data from the file
results_path = sys.argv[1]
file_path = results_path+'/results_readdistribution_actualreads_confidencescores'
data = pd.read_csv(file_path, sep='\t', index_col=0)

cmap = sns.color_palette(["#DBDBDB", "#A9D18E"])  # grey for 0, light green for 1

# Create the heatmap
plt.figure(figsize=(6, 3))
ax = sns.heatmap(data, cmap=cmap, cbar=False, linewidths=0.5, linecolor='white')

legend_patches = [mpatches.Patch(color="#DBDBDB", label="absence"), mpatches.Patch(color="#A9D18E", label="presence")]
plt.legend(handles=legend_patches, title='Strain', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

ax.set_xlabel('')
ax.set_ylabel('Strain')

ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()

# Save the heatmap to a file
plt.savefig(results_path+'/Strain_tracking_results.png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

