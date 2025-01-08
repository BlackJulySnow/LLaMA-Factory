import json
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import numpy as np

# Read data and count fragments
fragment_counts = []
file_path = "custom/result/ChemDual_pretrained_ChemDual_fragment_test.jsonl"

with open(file_path, "r") as f:
    for line in f:
        data = json.loads(line.strip())
        if data and 'predict' in data:
            count = len(data['predict'].split('.'))
            if count > 15:  # Cap at 15
                count = 15
            fragment_counts.append(count)

# Calculate basic statistics
avg_fragments = np.mean(fragment_counts)
median_fragments = np.median(fragment_counts)
min_fragments = np.min(fragment_counts)
max_fragments = np.max(fragment_counts)

# Create frequency statistics
counter = Counter(fragment_counts)
frequencies = dict(sorted(counter.items()))

# Create figure
plt.figure(figsize=(12, 8))

# Plot histogram and density curve
sns.histplot(data=fragment_counts, bins=15, color='skyblue', stat='density', label='Frequency Distribution')
sns.kdeplot(data=fragment_counts, color='red', linewidth=2, label='Density Curve')

# Add title and labels
# plt.title('Distribution of Molecular Fragments')
plt.xlabel('Number of Fragments', fontsize=24)
plt.ylabel('Density', fontsize=24)
plt.legend(fontsize=14)

# Set x-axis range from 0 to 15
plt.xlim(0, 15)
# Set x-axis ticks to integers
plt.xticks(range(0, 16), fontsize=20)
plt.yticks(fontsize=20)

# Add statistical information
stats_text = f'Statistics:\n' \
             f'Mean: {avg_fragments:.2f}\n' \
             f'Median: {median_fragments:.2f}\n' \
             f'Min: {min_fragments}\n' \
             f'Max: {max_fragments}'
plt.text(0.95, 0.95, stats_text,
         transform=plt.gca().transAxes,
         verticalalignment='top',
         horizontalalignment='right',
         fontsize=24,
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Adjust layout
plt.tight_layout()

# Save the figure
plt.savefig('distribution_fragments.pdf', format='pdf', bbox_inches='tight', dpi=300)

# Show plot
plt.show()