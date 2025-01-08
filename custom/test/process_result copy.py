import json
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

# Set font sizes for readability
plt.rcParams.update({'font.size': 14})  # General font size
sns.set_context("talk", font_scale=1.2)  # Seaborn context for larger font

# Read data and calculate molecular weights
molecular_weights = []
file_path = "custom/result/ChemDual_pretrained_ChemDual_recombination_test.jsonl"

with open(file_path, "r") as f:
    for line in f:
        data = json.loads(line.strip())
        if data and 'predict' in data:
            try:
                mol = Chem.MolFromSmiles(data['predict'])
                if mol:  # Check if molecule is valid
                    mw = Descriptors.ExactMolWt(mol)
                    if mw > 1000:
                        continue  # Ignore molecules with molecular weight > 1000
                    molecular_weights.append(mw)
            except:
                continue
# Calculate basic statistics for molecular weights
avg_weight = np.mean(molecular_weights)
median_weight = np.median(molecular_weights)
min_weight = np.min(molecular_weights)
max_weight = np.max(molecular_weights)
std_weight = np.std(molecular_weights)

# Create figure for molecular weights
plt.figure(figsize=(12, 8))

# Plot histogram and density curve for molecular weights
sns.histplot(data=molecular_weights, bins=30, color='skyblue', stat='density', label='Frequency Distribution')
sns.kdeplot(data=molecular_weights, color='red', linewidth=2, label='Density Curve')

# Add title and labels with consistent font sizes
plt.xlabel('Molecular Weight (g/mol)', fontsize=24)
plt.ylabel('Density', fontsize=24)
plt.legend(fontsize=14)  # Consistent legend font size

# Set x-axis range and tick font size
plt.xlim(0, max_weight + 50)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Add statistical information with consistent font size
stats_text = f'Statistics:\n' \
             f'Mean: {avg_weight:.2f}\n' \
             f'Median: {median_weight:.2f}\n' \
             f'Std: {std_weight:.2f}\n' \
             f'Min: {min_weight:.2f}\n' \
             f'Max: {max_weight:.2f}'
plt.text(0.95, 0.95, stats_text,
         transform=plt.gca().transAxes,
         verticalalignment='top',
         horizontalalignment='right',
         fontsize=24,  # Font size for text box
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Adjust layout and show plot for molecular weights
plt.tight_layout()
plt.show()

