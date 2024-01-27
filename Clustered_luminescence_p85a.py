import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import seaborn as sns

protein = "P85a"

P85a = "TCTGCAGAAGGTTATCAGTATCGCGCACTGTATGATTACAAAAAAGAACGCGAAGAAGATATTGATCTGCATCTGGGTGATATTCTGACCGTTAATAAAGGTAGCCTGGTTGCACTGGGTTTTTCTGATGGTCAGGAAGCACGTCCGGAAGAAATTGGTTGGCTGAATGGCTATAATGAAACCACCGGTGAACGTGGTGATTTTCCGGGTACGTATGTGGAATATATTGGTCGCAAAAAAATCAGCCCTCCGACCCCGAAACCTCGTCCTCCGCGTCCGCTGCCGGTTGCTCCGGGTAGCAGCAAAACCGAAGCAGATGTTGAACAGCAGGCACTGACCCTGCCGGATCTGGCAGAACAGTTTGCACCTCCGGATATTGCACCTCCGCTGCTGATTAAACTGGTTGAAGCCATTGAAAAAAAAGGTCTGGAATGCAGCACCCTGTATCGTACCCAGAGCAGCAGCAATCTGGCCGAACTGCGTCAGCTGCTGGATTGTGATACCCCGAGCGTTGATCTGGAAATGATTGATGTTCATGTTCTGGCCGATGCCTTTAAACGTTATCTGCTGGATCTGCCGAATCCGGTTATTCCGGCAGCAGTTTATAGCGAAATGATTAGCCTGGCACCGGAAGTTCAGAGCAGCGAAGAATATATTCAGCTGCTGAAAAAACTGATTCGTAGCCCGAGCATTCCGCATCAGTATTGGCTGACCCTGCAGTATCTGCTGAAACATTTCTTTAAACTGAGCCAGACCAGCAGCAAAAATCTGCTGAATGCACGTGTTCTGAGCGAAATTTTTAGCCCGATGCTGTTTCGTTTTAGCGCAGCAAGCAGCGATAACACCGAAAATCTGATTAAAGTGATTGAAATTCTGATTAGCACCGAATGGAATGAACGTCAGCCTGCACCGGCACTGCCTCCGAAACCTCCGAAACCGACCACCGTTGCAAATAATGGCATGAATAACAATATGAGCCTGCAGGATGCAGAATGGTATTGGGGTGATATTAGCCGTGAAGAAGTGAACGAAAAACTGCGTGATACCGCAGATGGCACCTTTCTGGTTCGTGATGCAAGCACCAAAATGCACGGTGATTATACCCTGACCCTGCGTAAAGGTGGTAACAATAAACTGATTAAAATCTTTCATCGTGATGGCAAATATGGTTTTAGCGATCCGCTGACCTTTAGCAGCGTTGTGGAACTGATTAATCATTATCGCAATGAAAGCCTGGCACAGTATAATCCGAAACTGGATGTGAAACTGCTGTATCCGGTTAGCAAATATCAGCAGGATCAGGTGGTGAAAGAAGATAATATTGAAGCCGTGGGCAAAAAACTGCATGAATATAATACCCAGTTTCAGGAAAAAAGCCGCGAATACGATCGCCTGTATGAAGAATATACCCGTACCAGCCAAGAGATTCAGATGAAACGTACCGCAATTGAAGCCTTTAATGAAACCATCAAAATCTTTGAAGAACAGTGCCAGACCCAGGAACGTTATAGCAAAGAATATATTGAAAAATTCAAACGCGAAGGCAATGAAAAAGAAATTCAGCGCATTATGCACAATTATGATAAACTGAAAAGCCGCATTAGCGAAATTATTGATAGCCGTCGTCGTCTGGAAGAAGATCTGAAAAAACAGGCAGCCGAATATCGCGAAATTGATAAACGCATGAATAGCATTAAACCGGATCTGATTCAGCTGCGTAAAACCCGTGATCAGTATCTGATGTGGCTGACCCAGAAAGGTGTTCGTCAGAAAAAACTGAATGAATGGCTGGGCAATGAAAATACCGAAGATCAGTATAGCCTGGTGGAAGATGATGAAGATCTGCCGCATCATGATGAAAAAACCTGGAATGTTGGTAGCAGCAATCGTAATAAAGCCGAAAATCTGCTGCGTGGTAAACGTGATGGTACATTTCTGGTGCGTGAAAGCAGCAAACAGGGTTGTTATGCCTGTAGCGTTGTTGTTGATGGCGAAGTTAAACATTGCGTGATTAACAAAACCGCAACAGGTTATGGTTTTGCCGAACCGTATAATCTGTATAGCAGCCTGAAAGAACTGGTTCTGCATTATCAGCATACCAGCCTGGTTCAGCATAATGATAGCCTGAATGTTACCCTGGCATATCCGGTTTATGCACAGCAGCGTCGT"


#########P85a#######
#file_path = f"C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/GG768_P85a_2/P85a2/20221101_1527_MN25343_ALL475_b5980e97/guppy_output/pass"
file_path = f"C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/2023_11_09_p85a_copy"


df = pd.read_csv(f'{file_path}/NGS_filtered.csv')

df["Cytosol"]= (df["Cytosol"]/150*250)
df=df[df["preinduction_auc"]>0.4]

#df = df[df["Full_Cell"] > 5000]
#df = df[df["end_position"]-df["start_position"] < 600]
#df = df[df["end_position"]-df["start_position"] >75]



#############Find AA seq################
inseq = []
for i in df.index:
    inseq.append(P85a[int(df["start_position"][i]):int(df["end_position"][i])])

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

aaseq = dna_to_protein(P85a)


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
############


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

#######################################
#######################################
# Calculate the average ratio of df['Cytosol']/df['Full_Cell'] for each cluster
cluster_colors = {}
cluster_sizes = {}
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]
    avg_ratio = np.mean(cluster_df['Cytosol'] / cluster_df['Full_Cell'])
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
    ax.plot([data['avg_start'], data['avg_end']], [index * spacing+2, index * spacing+2], color=data['color'], marker='o', linestyle='-')

# Create a colorbar
sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(0, 1))
cbar = plt.colorbar(sm)
cbar.set_label('Solubility')


plt.plot([(0*3),96*3],[0.04,0.04],c="purple", linewidth=5)
plt.plot([(115*3),301*3],[0.04,0.04],c="green", linewidth=5)
plt.plot([(326*3),435*3],[0.04,0.04],c="purple", linewidth=5)
plt.plot([(441*3),588*3],[0.04,0.04],c="yellow", linewidth=5)
plt.plot([(616*3),724*3],[0.04,0.04],c="green", linewidth=5)


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
    FC = np.mean(cluster_df['Full_Cell']*cluster_df["protein_kda"]/10**6)
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
cbar.set_label('Full cell (nM*KdA/10^6)')


plt.plot([(0*3),96*3],[0.04,0.04],c="purple", linewidth=5)
plt.plot([(115*3),301*3],[0.04,0.04],c="green", linewidth=5)
plt.plot([(326*3),435*3],[0.04,0.04],c="purple", linewidth=5)
plt.plot([(441*3),588*3],[0.04,0.04],c="yellow", linewidth=5)
plt.plot([(616*3),724*3],[0.04,0.04],c="green", linewidth=5)

# Label axes
plt.xlabel('BP positions')
#plt.ylabel('Length')
plt.title(f"{protein} Full cell")
plt.savefig(f"{file_path}/2new_full_cell", dpi=600)
# Display the plot
plt.show()



#################################################
# Calculate the Cytosolic
cluster_colors = {}
cluster_sizes = {}
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]
    cytosolic = np.mean(cluster_df['Cytosol']*cluster_df["protein_kda"]/10**6) 
    cluster_colors[cluster] = cytosolic
    cluster_sizes[cluster] = len(cluster_df)

# Use a color palette from seaborn, here the actual colors in the plot is defined
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
sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(0, top_lum))
cbar = plt.colorbar(sm)
cbar.set_label('Cytosolic (nM*KdA/10^6)')

 
plt.plot([(0*3),96*3],[0.04,0.04],c="purple", linewidth=5)
plt.plot([(115*3),301*3],[0.04,0.04],c="green", linewidth=5)
plt.plot([(326*3),435*3],[0.04,0.04],c="purple", linewidth=5)
plt.plot([(441*3),588*3],[0.04,0.04],c="yellow", linewidth=5)
plt.plot([(616*3),724*3],[0.04,0.04],c="green", linewidth=5)

# Label axes
plt.xlabel('BP positions')
#plt.ylabel('Length')
plt.title(f"{protein} Cytosolic")
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


plt.plot([(0*3),96*3],[0.04,0.04],c="purple", linewidth=5)
plt.plot([(115*3),301*3],[0.04,0.04],c="green", linewidth=5)
plt.plot([(326*3),435*3],[0.04,0.04],c="purple", linewidth=5)
plt.plot([(441*3),588*3],[0.04,0.04],c="yellow", linewidth=5)
plt.plot([(616*3),724*3],[0.04,0.04],c="green", linewidth=5)

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



'''
TF3: 
    
plt.plot([0,(42*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([(42*3),138*3],[-2,-2],c="green", linewidth=5)
plt.plot([(138*3),(241*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([(241*3),(888)],[-2,-2],c="green", linewidth=5)
'''

'''
SpA: 
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
plt.plot([1*3,(43*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([44*3,(79*3)],[-2,-2],c="green", linewidth=5)
plt.plot([(81*3),116*3],[-2,-2],c="purple", linewidth=5)
plt.plot([(117*3),149*3],[-2,-2],c="green", linewidth=5)
'''


'''
Mem: 
plt.plot([0,(116*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([(116*3),204*3],[-2,-2],c="green", linewidth=5)
'''

