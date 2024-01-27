import sys
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import seaborn as sns
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.colors as mcolors
from PIL import Image
import os
protein = "Calmodulin"
dnaseq = "ATGGCGGATCAACTCACGGAGGAACAAATCGCTGAATTTAAAGAAGCTTTTTCCTTATTTGACAAAGACGGTGATGGCACCATTACGACAAAAGAACTTGGTACCGTTATGCGTAGTTTGGGGCAAAACCCAACGGAGGCGGAATTACAGGACATGATTAACGAAGTGGATGCAGATGGTAATGGTACAATCGATTTTCCAGAGTTTCTCACGATGATGGCGCGTAAGATGAAAGATACGGATAGCGAGGAGGAGATCCGCGAAGCGTTCCGGGTATTTGATAAAGATGGGAACGGCTACATTTCTGCTGCTGAACTGCGCCACGTTATGACGAATTTGGGCGAAAAACTTACCGACGAGGAAGTGGACGAAATGATTCGGGAAGCGGATATTGACGGCGATGGGCAAGTGAACTACGAGGAGTTTGTACAGATGATGACTGCTAAA".upper()




    


top_lum = 500#max_value = max(cluster_colors.values())


#########SpA/TF3##########
#file_path = f"C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/2023_07_26_spa_tf3/{protein}_alignment"

#########P85a###########
#file_path = f"C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/2023_11_09_p85a_copy"
#dnaseq = "TCTGCAGAAGGTTATCAGTATCGCGCACTGTATGATTACAAAAAAGAACGCGAAGAAGATATTGATCTGCATCTGGGTGATATTCTGACCGTTAATAAAGGTAGCCTGGTTGCACTGGGTTTTTCTGATGGTCAGGAAGCACGTCCGGAAGAAATTGGTTGGCTGAATGGCTATAATGAAACCACCGGTGAACGTGGTGATTTTCCGGGTACGTATGTGGAATATATTGGTCGCAAAAAAATCAGCCCTCCGACCCCGAAACCTCGTCCTCCGCGTCCGCTGCCGGTTGCTCCGGGTAGCAGCAAAACCGAAGCAGATGTTGAACAGCAGGCACTGACCCTGCCGGATCTGGCAGAACAGTTTGCACCTCCGGATATTGCACCTCCGCTGCTGATTAAACTGGTTGAAGCCATTGAAAAAAAAGGTCTGGAATGCAGCACCCTGTATCGTACCCAGAGCAGCAGCAATCTGGCCGAACTGCGTCAGCTGCTGGATTGTGATACCCCGAGCGTTGATCTGGAAATGATTGATGTTCATGTTCTGGCCGATGCCTTTAAACGTTATCTGCTGGATCTGCCGAATCCGGTTATTCCGGCAGCAGTTTATAGCGAAATGATTAGCCTGGCACCGGAAGTTCAGAGCAGCGAAGAATATATTCAGCTGCTGAAAAAACTGATTCGTAGCCCGAGCATTCCGCATCAGTATTGGCTGACCCTGCAGTATCTGCTGAAACATTTCTTTAAACTGAGCCAGACCAGCAGCAAAAATCTGCTGAATGCACGTGTTCTGAGCGAAATTTTTAGCCCGATGCTGTTTCGTTTTAGCGCAGCAAGCAGCGATAACACCGAAAATCTGATTAAAGTGATTGAAATTCTGATTAGCACCGAATGGAATGAACGTCAGCCTGCACCGGCACTGCCTCCGAAACCTCCGAAACCGACCACCGTTGCAAATAATGGCATGAATAACAATATGAGCCTGCAGGATGCAGAATGGTATTGGGGTGATATTAGCCGTGAAGAAGTGAACGAAAAACTGCGTGATACCGCAGATGGCACCTTTCTGGTTCGTGATGCAAGCACCAAAATGCACGGTGATTATACCCTGACCCTGCGTAAAGGTGGTAACAATAAACTGATTAAAATCTTTCATCGTGATGGCAAATATGGTTTTAGCGATCCGCTGACCTTTAGCAGCGTTGTGGAACTGATTAATCATTATCGCAATGAAAGCCTGGCACAGTATAATCCGAAACTGGATGTGAAACTGCTGTATCCGGTTAGCAAATATCAGCAGGATCAGGTGGTGAAAGAAGATAATATTGAAGCCGTGGGCAAAAAACTGCATGAATATAATACCCAGTTTCAGGAAAAAAGCCGCGAATACGATCGCCTGTATGAAGAATATACCCGTACCAGCCAAGAGATTCAGATGAAACGTACCGCAATTGAAGCCTTTAATGAAACCATCAAAATCTTTGAAGAACAGTGCCAGACCCAGGAACGTTATAGCAAAGAATATATTGAAAAATTCAAACGCGAAGGCAATGAAAAAGAAATTCAGCGCATTATGCACAATTATGATAAACTGAAAAGCCGCATTAGCGAAATTATTGATAGCCGTCGTCGTCTGGAAGAAGATCTGAAAAAACAGGCAGCCGAATATCGCGAAATTGATAAACGCATGAATAGCATTAAACCGGATCTGATTCAGCTGCGTAAAACCCGTGATCAGTATCTGATGTGGCTGACCCAGAAAGGTGTTCGTCAGAAAAAACTGAATGAATGGCTGGGCAATGAAAATACCGAAGATCAGTATAGCCTGGTGGAAGATGATGAAGATCTGCCGCATCATGATGAAAAAACCTGGAATGTTGGTAGCAGCAATCGTAATAAAGCCGAAAATCTGCTGCGTGGTAAACGTGATGGTACATTTCTGGTGCGTGAAAGCAGCAAACAGGGTTGTTATGCCTGTAGCGTTGTTGTTGATGGCGAAGTTAAACATTGCGTGATTAACAAAACCGCAACAGGTTATGGTTTTGCCGAACCGTATAATCTGTATAGCAGCCTGAAAGAACTGGTTCTGCATTATCAGCATACCAGCCTGGTTCAGCATAATGATAGCCTGAATGTTACCCTGGCATATCCGGTTTATGCACAGCAGCGTCGT".upper()
  
########CALM/Mem##########
file_path = f"C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/2023_08_16_Calm_mem_2/{protein}"

#Calculate the standard curve
def nMol_calc(x):
    #SpA_TF3
    #result = 1.1498*x - 761.17
    #calmmem 
    result =1.4792*x - 2254.1
    #p85a
    #result = 0.8477*x - 678.47
    #p85a_1000RPM
    #result = 1.8577*x
    return result






df = pd.read_csv(f'{file_path}/NGS_filtered.csv')
df["plate"] = df['platePos'].str.extract(r'(\d+),').astype(int)
#df = df[df["plate"]!=6]##########################JUST FOR P85a!!!##################
#Remove low-expression constructs
#df = df[df["Full_Cell"]>10000]

df=df[df["preinduction_auc"]>0.4]


#df = df[df[]]

df["Cytosol"]= nMol_calc(df["Cytosol"])/150*250
df["Full_Cell"]= nMol_calc(df["Full_Cell"])
#df = df[df["end_position"]-df["start_position"] < 600]
#df = df[df["end_position"]-df["start_position"] >75]

#df = df[df["Full_Cell"] > 1000]


#############Find AA seq################
inseq = []
for i in df.index:
    inseq.append(dnaseq[int(df["start_position"][i]):int(df["end_position"][i])])

df["seqBP"] = inseq


#######aa-translator
def dna_to_protein(dna_sequence):
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    
    protein_sequence = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if codon in codon_table:
            amino_acid = codon_table[codon]
            protein_sequence += amino_acid
        else:
            print("Invalid codon: " + codon)
            return
    return protein_sequence

aaseq = dna_to_protein(dnaseq)


aainseq = []
for i in df.index:
    #print(int(Input["start_position"][i]),int(Input["end_position"][i]))
    aainseq.append(aaseq[round(int(df["start_position"][i])/3):round(int(df["end_position"][i])/3)])


df["seqAA"] = aainseq

####################################

##########Calculate protein weights based on Promegas molecular weights for amino acids########
def calculate_protein_mass(sequence):
    aa_weights = {
        'A': 89.09,  'R': 174.20, 'N': 132.12, 'D': 133.10,
        'C': 121.16, 'E': 147.13, 'Q': 146.15, 'G': 75.07,
        'H': 155.16, 'I': 131.17, 'L': 131.17, 'K': 146.19,
        'M': 149.21, 'F': 165.19, 'P': 115.13, 'S': 105.09,
        'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
    }
    total_weight = sum(aa_weights.get(aa, 0) for aa in sequence)
    total_weight -= (len(sequence) - 1) * 18.015
    return total_weight / 1000  # convert to kDa



df['protein_kda'] = df['seqAA'].apply(calculate_protein_mass)
###########


#######################################
#Calculate the protein in g/L
df["mg_l_cyt"] = df["Cytosol"]*10**(-9)*df["protein_kda"]*(10**6)
df["mg_l_fc"] = df["Full_Cell"]*10**(-9)*df["protein_kda"]*(10**6)
#######################################



# Remove data where Cytosol/Full_Cell > 1.1
#df = df[df['Cytosol'] / df['Full_Cell'] <= 2]
#Remove data below 0 (after calculations)

df = df[df['Cytosol']>0] 
df = df[df['Full_Cell']>0]

n_clusters = 40

kmeans = KMeans(n_clusters=n_clusters)
df['cluster'] = kmeans.fit_predict(df[['start_position', 'end_position']])

# Calculate the y-axis values
#df['y_values'] = df['end_position'] - df['start_position']
df['sizes'] = df['end_position'] - df['start_position']
#numbers = [1*n_dist for n_dist in range(1, df.shape[0]+1)]
#df = df.sort_values(by="sizes")

df['y_values'] = df['sizes']



# Calculate the average ratio of df['Cytosol']/df['Full_Cell'] for each cluster
cluster_colors = {}
cluster_sizes = {}
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]
    avg_ratio = np.mean(cluster_df['Cytosol'] / cluster_df['Full_Cell'])
    #print(avg_ratio)
    cluster_colors[cluster] = avg_ratio
    cluster_sizes[cluster] = len(cluster_df)



# Use a color palette from seaborn
color_palette = sns.color_palette("coolwarm", n_clusters)
colors = [color_palette[int(np.interp(cluster_colors[i], (0, 1), (0, n_clusters-1)))] for i in range(n_clusters)]

# Create a line plot
fig, ax = plt.subplots()

cluster_data = []

# Calculate data for each cluster
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]

    avg_start = np.mean(cluster_df['start_position'])
    avg_end = np.mean(cluster_df['end_position'])
    avg_y = np.mean((cluster_df['end_position'] - cluster_df['start_position']) / 2 + cluster_df['start_position'])

    cluster_data.append({
        'avg_start': avg_start,
        'avg_end': avg_end,
        'avg_y': avg_y,
        'color': colors[cluster]
    })

# Sort the cluster data based on avg_y
cluster_data_sorted = sorted(cluster_data, key=lambda x: x['avg_y'])
cluster_data = []

# Calculate data for each cluster
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]

    avg_start = np.mean(cluster_df['start_position'])
    avg_end = np.mean(cluster_df['end_position'])
    avg_y = np.mean((cluster_df['end_position'] - cluster_df['start_position']) / 2 + cluster_df['start_position'])

    cluster_data.append({
        'avg_start': avg_start,
        'avg_end': avg_end,
        'avg_y': avg_y,
        'color': colors[cluster]
    })

# Sort the cluster data based on avg_y
cluster_data_sorted = sorted(cluster_data, key=lambda x: x['avg_y'])


# Now plot with equal spacing
spacing = 1  # Define the desired spacing between lines
for index, data in enumerate(cluster_data_sorted):
    ax.plot([data['avg_start'], data['avg_end']], [index * spacing, index * spacing], color=data['color'], marker='o', linestyle='-')
    



# Create a colorbar
sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(0, 1))
cbar = plt.colorbar(sm)
cbar.set_label('Solubility')

plt.plot([1*3,(20*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([20*3,(45*3)],[-2,-2],c="green", linewidth=5)
plt.plot([(45*3),62*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(62*3),93*3],[-2,-2],c="yellow", linewidth=5)
plt.plot([(93*3),111*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(111*3),130*3],[-2,-2],c="green", linewidth=5)
plt.plot([(130*3),149*3],[-2,-2],c="purple", linewidth=5)
    
# Label axes
plt.xlabel('BP positions')
#plt.ylabel('Length')
plt.title(f"{protein} Solubility")
# Display the plot
plt.savefig(f"{file_path}/2new_solubility", dpi=600)
plt.show()

#######################################
# Calculate the Full cell luminescence
cluster_colors = {}
cluster_sizes = {}
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]
    FC = np.mean(cluster_df['mg_l_fc'])
    cluster_colors[cluster] = FC
    cluster_sizes[cluster] = len(cluster_df)

# Use a color palette from seaborn
color_palette = sns.color_palette("coolwarm", n_clusters)
colors = [color_palette[int(np.interp(cluster_colors[i], (0, top_lum), (0, n_clusters-1)))] for i in range(n_clusters)]

# Create a line plot
fig, ax = plt.subplots()

cluster_data = []

# Calculate data for each cluster
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]

    avg_start = np.mean(cluster_df['start_position'])
    avg_end = np.mean(cluster_df['end_position'])
    avg_y = np.mean((cluster_df['end_position'] - cluster_df['start_position']) / 2 + cluster_df['start_position'])

    cluster_data.append({
        'avg_start': avg_start,
        'avg_end': avg_end,
        'avg_y': avg_y,
        'color': colors[cluster]
    })

# Sort the cluster data based on avg_y
cluster_data_sorted = sorted(cluster_data, key=lambda x: x['avg_y'])
cluster_data = []

# Calculate data for each cluster
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]

    avg_start = np.mean(cluster_df['start_position'])
    avg_end = np.mean(cluster_df['end_position'])
    avg_y = np.mean((cluster_df['end_position'] - cluster_df['start_position']) / 2 + cluster_df['start_position'])

    cluster_data.append({
        'avg_start': avg_start,
        'avg_end': avg_end,
        'avg_y': avg_y,
        'color': colors[cluster]
    })

# Sort the cluster data based on avg_y
cluster_data_sorted = sorted(cluster_data, key=lambda x: x['avg_y'])

# Now plot with equal spacing
spacing = 1  # Define the desired spacing between lines
for index, data in enumerate(cluster_data_sorted):
    ax.plot([data['avg_start'], data['avg_end']], [index * spacing, index * spacing], color=data['color'], marker='o', linestyle='-')

# Create a colorbar
sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(0,top_lum))
cbar = plt.colorbar(sm)
cbar.set_label('mg/L')

plt.plot([1*3,(20*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([20*3,(45*3)],[-2,-2],c="green", linewidth=5)
plt.plot([(45*3),62*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(62*3),93*3],[-2,-2],c="yellow", linewidth=5)
plt.plot([(93*3),111*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(111*3),130*3],[-2,-2],c="green", linewidth=5)
plt.plot([(130*3),149*3],[-2,-2],c="purple", linewidth=5)
    
# Label axes
plt.xlabel('BP positions')
#plt.ylabel('Length')
plt.title(f"{protein} Full cell (mg/L)")
plt.savefig(f"{file_path}/2new_full_cell", dpi=600)
# Display the plot
plt.show()



#################################################
# Calculate the Cytosolic
cluster_colors = {}
cluster_sizes = {}
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]
    cytosolic = np.mean(cluster_df['mg_l_cyt']) #Cytosol in mg/L
    cluster_colors[cluster] = cytosolic
    cluster_sizes[cluster] = len(cluster_df)

# Use a color palette from seaborn
color_palette = sns.color_palette("coolwarm", n_clusters)
colors = [color_palette[int(np.interp(cluster_colors[i], (0, top_lum), (0, n_clusters-1)))] for i in range(n_clusters)]

# Create a line plot
fig, ax = plt.subplots()

cluster_data = []

# Calculate data for each cluster
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]

    avg_start = np.mean(cluster_df['start_position'])
    avg_end = np.mean(cluster_df['end_position'])
    avg_y = np.mean((cluster_df['end_position'] - cluster_df['start_position']) / 2 + cluster_df['start_position'])

    cluster_data.append({
        'avg_start': avg_start,
        'avg_end': avg_end,
        'avg_y': avg_y,
        'color': colors[cluster]
    })

# Sort the cluster data based on avg_y
cluster_data_sorted = sorted(cluster_data, key=lambda x: x['avg_y'])
cluster_data = []

# Calculate data for each cluster
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]

    avg_start = np.mean(cluster_df['start_position'])
    avg_end = np.mean(cluster_df['end_position'])
    avg_y = np.mean((cluster_df['end_position'] - cluster_df['start_position']) / 2 + cluster_df['start_position'])

    cluster_data.append({
        'avg_start': avg_start,
        'avg_end': avg_end,
        'avg_y': avg_y,
        'color': colors[cluster]
    })

# Sort the cluster data based on avg_y
cluster_data_sorted = sorted(cluster_data, key=lambda x: x['avg_y'])

# Now plot with equal spacing
spacing = 1  # Define the desired spacing between lines
for index, data in enumerate(cluster_data_sorted):
    ax.plot([data['avg_start'], data['avg_end']], [index * spacing, index * spacing], color=data['color'], marker='o', linestyle='-')

# Create a colorbar
sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(0,top_lum))
cbar = plt.colorbar(sm)
cbar.set_label('Cytosolic (mg/L)')

plt.plot([1*3,(20*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([20*3,(45*3)],[-2,-2],c="green", linewidth=5)
plt.plot([(45*3),62*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(62*3),93*3],[-2,-2],c="yellow", linewidth=5)
plt.plot([(93*3),111*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(111*3),130*3],[-2,-2],c="green", linewidth=5)
plt.plot([(130*3),149*3],[-2,-2],c="purple", linewidth=5)

# Label axes
plt.xlabel('BP positions')
#plt.ylabel('Length')
plt.title(f"{protein} Cytosolic (mg/L)")
# Display the plot
plt.savefig(f"{file_path}/2new_cytosolic", dpi=600)
plt.show()

#######################################
# Calculate the postinduction_growth
cluster_colors = {}
cluster_sizes = {}
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]
    FC = np.mean(cluster_df['postinduction_auc'])
    cluster_colors[cluster] = FC
    cluster_sizes[cluster] = len(cluster_df)
top_lum = max_value = max(cluster_colors.values())
# Use a color palette from seaborn
colors = [color_palette[int(np.interp(cluster_colors[i], (0, top_lum), (0, n_clusters-1)))] for i in range(n_clusters)]

# Create a line plot
fig, ax = plt.subplots()

cluster_data = []

# Calculate data for each cluster
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]

    avg_start = np.mean(cluster_df['start_position'])
    avg_end = np.mean(cluster_df['end_position'])
    avg_y = np.mean((cluster_df['end_position'] - cluster_df['start_position']) / 2 + cluster_df['start_position'])

    cluster_data.append({
        'avg_start': avg_start,
        'avg_end': avg_end,
        'avg_y': avg_y,
        'color': colors[cluster]
    })

# Sort the cluster data based on avg_y
cluster_data_sorted = sorted(cluster_data, key=lambda x: x['avg_y'])
cluster_data = []

# Calculate data for each cluster
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]

    avg_start = np.mean(cluster_df['start_position'])
    avg_end = np.mean(cluster_df['end_position'])
    avg_y = np.mean((cluster_df['end_position'] - cluster_df['start_position']) / 2 + cluster_df['start_position'])

    cluster_data.append({
        'avg_start': avg_start,
        'avg_end': avg_end,
        'avg_y': avg_y,
        'color': colors[cluster]
    })

# Sort the cluster data based on avg_y
cluster_data_sorted = sorted(cluster_data, key=lambda x: x['avg_y'])

# Now plot with equal spacing
spacing = 1  # Define the desired spacing between lines
for index, data in enumerate(cluster_data_sorted):
    ax.plot([data['avg_start'], data['avg_end']], [index * spacing, index * spacing], color=data['color'], marker='o', linestyle='-')

# Create a colorbar
sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(0,top_lum))
cbar = plt.colorbar(sm)
cbar.set_label('Postinduction OD600 (AUC)')

plt.plot([1*3,(20*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([20*3,(45*3)],[-2,-2],c="green", linewidth=5)
plt.plot([(45*3),62*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(62*3),93*3],[-2,-2],c="yellow", linewidth=5)
plt.plot([(93*3),111*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(111*3),130*3],[-2,-2],c="green", linewidth=5)
plt.plot([(130*3),149*3],[-2,-2],c="purple", linewidth=5)
    

# Label axes
plt.xlabel('BP positions')
#plt.ylabel('Length')
plt.title(f"{protein} OD600 AUC postiduction")
plt.savefig(f"{file_path}/2new_AUC", dpi=600)
# Display the plot
plt.show()

####################
# Create a scatter plot of cluster size vs average Cytosol/Full_Cell ratio
fig, ax = plt.subplots()
ax.scatter(cluster_sizes.values(), cluster_colors.values(), alpha=0.5)
ax.set_xlabel('Cluster Size')
ax.set_ylabel('Average Ratio of Cytosol/Full_Cell')
plt.title(protein)
plt.show()


dfsave = f"{file_path}/final_overview.xlsx"
df.to_excel(dfsave)



#############Combining figures#############
image1_path = os.path.join(file_path, '2new_cytosolic.png')
image2_path = os.path.join(file_path, '2new_full_cell.png')
image3_path = os.path.join(file_path, '2new_AUC.png')

# Load the PNG images
image1 = Image.open(image1_path)
image2 = Image.open(image2_path)
image3 = Image.open(image3_path)

# Check the sizes of the images to ensure compatibility
width1, height1 = image1.size
width2, height2 = image2.size
width3, height3 = image3.size

# Ensure all images have the same height
if height1 != height2 or height1 != height3:
    raise ValueError("Image heights do not match")

# Set the truncation widths for each image (individual for both sides)
truncate_width1_left = 250
truncate_width1_right = 900

truncate_width2_left = 450
truncate_width2_right = 250

truncate_width3_left = 250
truncate_width3_right = 0

# Calculate the combined width and height
combined_width = (width1 - truncate_width1_left - truncate_width1_right) + \
                 (width2 - truncate_width2_left - truncate_width2_right) + \
                 (width3 - truncate_width3_left - truncate_width3_right)
combined_height = height1

# Create a new image with the combined width and same height
combined_image = Image.new("RGBA", (combined_width, combined_height))

# Paste the images side by side onto the new image with truncation
combined_image.paste(image1.crop((truncate_width1_left, 0, width1 - truncate_width1_right, height1)), (truncate_width1_left, 0))
combined_image.paste(image2.crop((truncate_width2_left, 0, width2 - truncate_width2_right, height2)),
                     (width1 - truncate_width1_right, 0))
combined_image.paste(image3.crop((truncate_width3_left, 0, width3 - truncate_width3_right, height3)),
                     (width1 + width2 - truncate_width1_right - truncate_width2_right-truncate_width2_left-truncate_width1_left+300, 0))

# Save the combined image
combined_image.save(os.path.join(file_path, 'combined_image.png'))


'''
TF3: 
    
protein = "TF3"
dnaseq = "ATGGAAACACCTGCATGGCCGCGCGTCCCACGGCCAGAAACAGCCGTCGCACGCACTCTGTTGTTGGGCTGGGTATTCGCTCAAGTGGCCGGTGCTAGTGGGACTACAAACACAGTTGCAGCATACAACTTGACTTGGAAAAGCACAAATTTTAAGACTATTTTGGAATGGGAACCGAAGCCGGTGAATCAGGTGTACACCGTCCAGATTTCAACAAAGAGCGGTGATTGGAAGAGCAAGTGCTTCTATACAACCGACACGGAGTGCGATTTAACGGACGAGATTGTCAAGGACGTGAAGCAGACTTACTTGGCCCGCGTGTTTAGCTACCCTGCCGGTAATGTCGAATCAACGGGCAGTGCAGGCGAACCTCTTTACGAAAACTCTCCTGAGTTCACTCCATATCTTGAGACAAACCTTGGGCAACCTACAATTCAGTCATTTGAACAAGTCGGTACTAAGGTGAATGTCACAGTTGAGGATGAACGGACCCTTGTGCGTCGCAATAACACTTTTCTCTCGCTCCGTGATGTTTTCGGTAAGGACCTCATCTATACGCTGTACTATTGGAAATCATCGTCATCTGGTAAAAAGACAGCTAAGACAAATACGAATGAATTCTTAATTGATGTGGATAAAGGGGAGAACTATTGCTTCAGTGTGCAAGCAGTGATCCCTTCACGGACTGTGAATCGCAAATCAACAGACAGTCCAGTCGAGTGTATGGGTCAAGAGAAGGGGGAATTCCGGGAAATTTTCTATATCATCGGCGCGGTTGTTTTTGTGGTCATTATTCTTGTCATTATCCTGGCCATCAGTCTTCACAAATGCCGGAAAGCTGGCGTTGGGCAAAGCTGGAAAGAAAACTCGCCATTGAACGTTAGC".upper()
    
    
    
plt.plot([0,(42*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([(42*3),138*3],[-2,-2],c="green", linewidth=5)
plt.plot([(138*3),(241*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([(241*3),(888)],[-2,-2],c="green", linewidth=5)
'''

'''
SpA: 
    
protein = "SpA"
dnaseq = "ATGAAGAAAAAGAATATTTACTCTATCCGTAAACTTGGCGTGGGGATTGCAAGCGTAACACTCGGTACCCTTCTCATCTCAGGCGGTGTAACACCAGCCGCAAACGCGGCTCAGCATGACGAAGCGCAACAAAATGCTTTCTACCAGGTGCTTAACATGCCAAACCTCAACGCCGACCAACGTAACGGCTTCATCCAATCTTTAAAAGACGATCCGTCCCAGTCCGCAAACGTCCTCGGGGAAGCTCAAAAACTGAATGATTCACAGGCCCCAAAAGCCGACGCGCAACAAAACAAGTTCAATAAAGACCAACAATCTGCTTTTTACGAAATCCTTAACATGCCTAACCTGAACGAGGAGCAGCGGAATGGGTTTATTCAATCGCTGAAAGACGACCCTAGCCAATCAACCAATGTATTAGGTGAAGCGAAAAAATTAAACGAATCACAAGCACCGAAAGCCGATAATAACTTCAATAAGGAGCAACAGAACGCCTTTTATGAGATTCTCAATATGCCTAATCTGAATGAGGAACAACGGAATGGTTTTATTCAGTCGTTAAAGGATGATCCGTCACAATCCGCCAACCTGCTGGCAGAAGCTAAGAAACTCAATGAGAGTCAAGCGCCTAAAGCCGACAATAAGTTTAACAAAGAACAGCAAAATGCTTTTTATGAGATCCTGCATTTGCCAAATCTGAACGAAGAGCAACGCAATGGGTTCATCCAGAGTCTGAAGGACGACCCTTCGCAGAGCGCAAATTTGTTAGCTGAAGCCAAAAAGTTGAATGACGCGCAAGCCCCGAAGGCTGACAACAAGTTCAATAAGGAGCAACAAAACGCTTTTTATGAAATTCTCCACCTCCCAAACCTTACTGAGGAGCAGCGTAACGGGTTTATTCAGTCACTGAAGGACGATCCGAGCGTATCTAAGGAAATTCTTGCCGAGGCGAAAAAATTAAATGACGCACAAGCGCCGAAAGAAGAAGACAATAATAAACCTGGCAAAGAAGATAACAATAAGCCTGGGAAAGAAGAT".upper()

    
plt.plot([0,(37*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([(37*3),92*3],[-2,-2],c="green", linewidth=5)
plt.plot([(92*3),(153*3)],[-2,-2],c="yellow", linewidth=5)
plt.plot([(153*3),211*3],[-2,-2],c="green", linewidth=5)
plt.plot([(211*3),(270*3)],[-2,-2],c="yellow", linewidth=5)
plt.plot([(270*3),327*3],[-2,-2],c="green", linewidth=5)
plt.plot([(327*3),1041],[-2,-2],c="purple", linewidth=5)
'''

'''
Calm: 

    
protein = "Calmodulin"
dnaseq = "ATGGCGGATCAACTCACGGAGGAACAAATCGCTGAATTTAAAGAAGCTTTTTCCTTATTTGACAAAGACGGTGATGGCACCATTACGACAAAAGAACTTGGTACCGTTATGCGTAGTTTGGGGCAAAACCCAACGGAGGCGGAATTACAGGACATGATTAACGAAGTGGATGCAGATGGTAATGGTACAATCGATTTTCCAGAGTTTCTCACGATGATGGCGCGTAAGATGAAAGATACGGATAGCGAGGAGGAGATCCGCGAAGCGTTCCGGGTATTTGATAAAGATGGGAACGGCTACATTTCTGCTGCTGAACTGCGCCACGTTATGACGAATTTGGGCGAAAAACTTACCGACGAGGAAGTGGACGAAATGATTCGGGAAGCGGATATTGACGGCGATGGGCAAGTGAACTACGAGGAGTTTGTACAGATGATGACTGCTAAA".upper()

plt.plot([1*3,(20*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([20*3,(45*3)],[-2,-2],c="green", linewidth=5)
plt.plot([(45*3),62*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(62*3),93*3],[-2,-2],c="yellow", linewidth=5)
plt.plot([(93*3),111*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(111*3),130*3],[-2,-2],c="green", linewidth=5)
plt.plot([(130*3),149*3],[-2,-2],c="purple", linewidth=5)
'''


'''
Mem: 
    
protein = "Membrane"
dnaseq = "atggcggatagcaacggcaccattaccgtggaagaactgaaaaaactgctggaacagtggaacctggtgattggctttctgtttctgacctggatttgcctgctgcagtttgcgtatgcgaaccgcaaccgctttctgtatattattaaactgatttttctgtggctgctgtggccggtgaccctggcgtgctttgtgctggcggcggtgtatcgcattaactggattaccggcggcattgcgattgcgatggcgtgcctggtgggcctgatgtggctgagctattttattgcgagctttcgcctgtttgcgcgcacccgcagcatgtggagctttaacccggaaaccaacattctgctgaacgtgccgctgcatggcaccattctgacccgcccgctgctggaaagcgaactggtgattggcgcggtgattctgcgcggccatctgcgcattgcgggccatcatctgggccgctgcgatattaaagatctgccgaaagaaattaccgtggcgaccagccgcaccctgagctattataaactgggcgcgagccagcgcgtggcgggcgatagcggctttgcggcgtatagccgctatcgcattggcaactataaactgaacaccgatcatagcagcagcagcgataacattgcgctgctggtgcag".upper()
    
    
plt.plot([0,(116*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([(116*3),204*3],[-2,-2],c="green", linewidth=5)
'''

'''
#p85a
#dnaseq = "TCTGCAGAAGGTTATCAGTATCGCGCACTGTATGATTACAAAAAAGAACGCGAAGAAGATATTGATCTGCATCTGGGTGATATTCTGACCGTTAATAAAGGTAGCCTGGTTGCACTGGGTTTTTCTGATGGTCAGGAAGCACGTCCGGAAGAAATTGGTTGGCTGAATGGCTATAATGAAACCACCGGTGAACGTGGTGATTTTCCGGGTACGTATGTGGAATATATTGGTCGCAAAAAAATCAGCCCTCCGACCCCGAAACCTCGTCCTCCGCGTCCGCTGCCGGTTGCTCCGGGTAGCAGCAAAACCGAAGCAGATGTTGAACAGCAGGCACTGACCCTGCCGGATCTGGCAGAACAGTTTGCACCTCCGGATATTGCACCTCCGCTGCTGATTAAACTGGTTGAAGCCATTGAAAAAAAAGGTCTGGAATGCAGCACCCTGTATCGTACCCAGAGCAGCAGCAATCTGGCCGAACTGCGTCAGCTGCTGGATTGTGATACCCCGAGCGTTGATCTGGAAATGATTGATGTTCATGTTCTGGCCGATGCCTTTAAACGTTATCTGCTGGATCTGCCGAATCCGGTTATTCCGGCAGCAGTTTATAGCGAAATGATTAGCCTGGCACCGGAAGTTCAGAGCAGCGAAGAATATATTCAGCTGCTGAAAAAACTGATTCGTAGCCCGAGCATTCCGCATCAGTATTGGCTGACCCTGCAGTATCTGCTGAAACATTTCTTTAAACTGAGCCAGACCAGCAGCAAAAATCTGCTGAATGCACGTGTTCTGAGCGAAATTTTTAGCCCGATGCTGTTTCGTTTTAGCGCAGCAAGCAGCGATAACACCGAAAATCTGATTAAAGTGATTGAAATTCTGATTAGCACCGAATGGAATGAACGTCAGCCTGCACCGGCACTGCCTCCGAAACCTCCGAAACCGACCACCGTTGCAAATAATGGCATGAATAACAATATGAGCCTGCAGGATGCAGAATGGTATTGGGGTGATATTAGCCGTGAAGAAGTGAACGAAAAACTGCGTGATACCGCAGATGGCACCTTTCTGGTTCGTGATGCAAGCACCAAAATGCACGGTGATTATACCCTGACCCTGCGTAAAGGTGGTAACAATAAACTGATTAAAATCTTTCATCGTGATGGCAAATATGGTTTTAGCGATCCGCTGACCTTTAGCAGCGTTGTGGAACTGATTAATCATTATCGCAATGAAAGCCTGGCACAGTATAATCCGAAACTGGATGTGAAACTGCTGTATCCGGTTAGCAAATATCAGCAGGATCAGGTGGTGAAAGAAGATAATATTGAAGCCGTGGGCAAAAAACTGCATGAATATAATACCCAGTTTCAGGAAAAAAGCCGCGAATACGATCGCCTGTATGAAGAATATACCCGTACCAGCCAAGAGATTCAGATGAAACGTACCGCAATTGAAGCCTTTAATGAAACCATCAAAATCTTTGAAGAACAGTGCCAGACCCAGGAACGTTATAGCAAAGAATATATTGAAAAATTCAAACGCGAAGGCAATGAAAAAGAAATTCAGCGCATTATGCACAATTATGATAAACTGAAAAGCCGCATTAGCGAAATTATTGATAGCCGTCGTCGTCTGGAAGAAGATCTGAAAAAACAGGCAGCCGAATATCGCGAAATTGATAAACGCATGAATAGCATTAAACCGGATCTGATTCAGCTGCGTAAAACCCGTGATCAGTATCTGATGTGGCTGACCCAGAAAGGTGTTCGTCAGAAAAAACTGAATGAATGGCTGGGCAATGAAAATACCGAAGATCAGTATAGCCTGGTGGAAGATGATGAAGATCTGCCGCATCATGATGAAAAAACCTGGAATGTTGGTAGCAGCAATCGTAATAAAGCCGAAAATCTGCTGCGTGGTAAACGTGATGGTACATTTCTGGTGCGTGAAAGCAGCAAACAGGGTTGTTATGCCTGTAGCGTTGTTGTTGATGGCGAAGTTAAACATTGCGTGATTAACAAAACCGCAACAGGTTATGGTTTTGCCGAACCGTATAATCTGTATAGCAGCCTGAAAGAACTGGTTCTGCATTATCAGCATACCAGCCTGGTTCAGCATAATGATAGCCTGAATGTTACCCTGGCATATCCGGTTTATGCACAGCAGCGTCGT".upper()


plt.plot([(0*3),96*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(115*3),301*3],[-2,-2],c="green", linewidth=5)
plt.plot([(326*3),435*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(441*3),588*3],[-2,-2],c="yellow", linewidth=5)
plt.plot([(616*3),724*3],[-2,-2],c="green", linewidth=5)
'''