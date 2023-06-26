#############Figure 2##############
############# iModulon Gene Weight ################
import matplotlib.pyplot as plt
from pymodulon.plotting import plot_regulon_histogram
fig, ax = plt.subplots(figsize=(3, 3))
# Call the plot_regulon_histogram function
plot_regulon_histogram(ica_data, 'RP-Operon-2', ax=ax)
# Set axes labels and title
ax.set_xlabel('Gene Expression Level')
ax.set_ylabel('Frequency')
ax.set_title('')
# Save the figure to the specified directory
plt.savefig("RP-Operon-2_Regulon.png", dpi=300, bbox_inches='tight')

#Let's make plots for the Activity
import pandas as pd
import matplotlib.pyplot as plt
color_scale = {
    'Brick 1': '#FFC300',   # Yellow
    'Brick 2': '#FF5733',   # Orange
    'Brick 3': '#C70039',   # Red
    'Brick 4': '#900C3F',   # Maroon
    'Brick 5': '#581845',   # Dark Purple
    'Brick 6': '#FF4500',   # Orange-Red
    'Brick 7': '#FF8C00',   # Dark Orange
    'Brick 8': '#FFA500',   # Orange
    'Brick 9': '#FFD700',   # Gold
    'Brick 10': '#FFFF00',  # Yellow
    'Brick 11': '#ADFF2F',  # Green-Yellow
    'Brick 12': '#32CD32',  # Lime Green
    'Brick 13': '#00FF00',  # Green
    'Brick 14': '#00FF7F',  # Spring Green
    'Brick 15': '#00FFFF',  # Cyan
    'Brick 16': '#00BFFF',  # Deep Sky Blue
    'Brick 17': '#0000FF',  # Blue
    'Brick 18': '#8A2BE2',  # Blue Violet
    'Brick 19': '#9932CC',  # Dark Orchid
    'Brick 20': '#BA55D3',  # Medium Orchid
    'Brick 21': '#800080',   # Purple
    'Brick 22': '#FF00FF',   # Magenta
    'Brick 23': '#FF1493',   # Deep Pink
    'Brick 24': '#FF69B4',   # Hot Pink
    'Brick 25': '#FFC0CB',   # Pink
}

# Read the CSV file
data = pd.read_csv('Figure2_Input.csv', delimiter=',')
# Specify the iModulon you want to plot
regulon = 'T-box(Met)'
# Filter the data for the selected iModulon
selected_data = data[data['iModulon'] == regulon]
# Transpose the data for plotting
selected_data = selected_data.set_index('iModulon').T
# Create the bar plot
fig, ax = plt.subplots(figsize=(10, 6))
selected_data.plot(kind='bar', ax=ax, legend=False, width=0.6)
# Add transparent color blocks
#Control
ax.axvspan(0, 0.5, facecolor='gray', alpha=0.2)  # Example block from x-axis position 2.5 to 4.5
#Salts
ax.axvspan(0.5, 3.5, facecolor=color_scale['Brick 1'], alpha=0.2)
#Carbohydrates
ax.axvspan(3.5, 6.5, facecolor=color_scale['Brick 2'], alpha=0.2)
#pH
ax.axvspan(6.5, 8.5, facecolor=color_scale['Brick 4'], alpha=0.2)
#Amino Acid Supplement
ax.axvspan(8.5, 13.5, facecolor=color_scale['Brick 6'], alpha=0.2)
#Amino Acid Removal
ax.axvspan(13.5, 18.5, facecolor=color_scale['Brick 8'], alpha=0.2)
#Vitamin Supplementation
ax.axvspan(18.5, 23.5, facecolor=color_scale['Brick 11'], alpha=0.2)
#Vitamin Removal
ax.axvspan(23.5, 24.5, facecolor=color_scale['Brick 12'], alpha=0.2)
#Bile Salt
ax.axvspan(24.5, 27.5, facecolor=color_scale['Brick 14'], alpha=0.2)
#Antibiotics
ax.axvspan(27.5, 30.5, facecolor=color_scale['Brick 15'], alpha=0.2)
#Gut Microbial-derived Metabolites
ax.axvspan(30.5, 35.5, facecolor=color_scale['Brick 16'], alpha=0.2)
#Nucleotide Supplementation
ax.axvspan(35.5, 36.5, facecolor=color_scale['Brick 17'], alpha=0.2)
#Nucleotide Removal
ax.axvspan(36.5, 38.5, facecolor=color_scale['Brick 19'], alpha=0.2)
#Human Food Medium
ax.axvspan(38.5, 40.5, facecolor=color_scale['Brick 22'], alpha=0.2)
#Media & Temp
ax.axvspan(40.5, 44.5, facecolor=color_scale['Brick 23'], alpha=0.2)
#Cocultures
ax.axvspan(44.5, 48.5, facecolor=color_scale['Brick 24'], alpha=0.2)
#Fe Supplementation/Removal
ax.axvspan(48.5, 52.5, facecolor=color_scale['Brick 25'], alpha=0.2)
# Customize the plot
plt.xlabel('')
plt.ylabel('')
plt.title(''.format(regulon))
plt.xticks(rotation=45)
plt.tight_layout()

# Save the plot to a file
plt.savefig("RP-Operon-2_Regulon.png", dpi=300, bbox_inches='tight')

# Display the plot
plt.show()
