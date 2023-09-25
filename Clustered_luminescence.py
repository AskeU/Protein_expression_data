import sys
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import seaborn as sns

protein = "Membrane"

#########SpA/TF3##########

#file_path = f"C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/2023_07_26_spa_tf3/{protein}_alignment"

########CALM/Mem##########
file_path = f"C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/2023_08_16_Calm_mem_2/{protein}"

#########P85a#######
#file_path = f"C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/GG768_P85a_2/P85a2/20221101_1527_MN25343_ALL475_b5980e97/guppy_output/Alignment"

#Calculate the standard curve
def nMol_calc(x):
    #SpA_TF3
    result = 0.2203 * x - 145.82
    #calmmem 
    #result = 0.2834*x - 431.82
    return result




df = pd.read_csv(f'{file_path}/NGS_filtered.csv')
df["plate"] = df['platePos'].str.extract(r'(\d+),').astype(int)

#Remove low-expression constructs
df = df[df["Full_Cell"]>10000]


#df = df[df[]]

df["Cytosol"]= nMol_calc(df["Cytosol"])#/130*230
df["Full_Cell"]= nMol_calc(df["Full_Cell"])
df = df[df["end_position"]-df["start_position"] < 600]
df = df[df["end_position"]-df["start_position"] >75]

#df = df[df["Full_Cell"] > 1000]


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


plt.plot([0,(116*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([(116*3),204*3],[-2,-2],c="green", linewidth=5)

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
    FC = np.mean(cluster_df['Full_Cell'])
    cluster_colors[cluster] = FC
    cluster_sizes[cluster] = len(cluster_df)
top_lum = max_value = max(cluster_colors.values())
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
cbar.set_label('Full cell (nM)')


plt.plot([0,(116*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([(116*3),204*3],[-2,-2],c="green", linewidth=5)
# Label axes
plt.xlabel('BP positions')
#plt.ylabel('Length')
plt.title(f"{protein} Full cell (nM)")
plt.savefig(f"{file_path}/2new_full_cell", dpi=600)
# Display the plot
plt.show()



#################################################
# Calculate the Cytosolic
cluster_colors = {}
cluster_sizes = {}
for cluster in range(n_clusters):
    cluster_df = df[df['cluster'] == cluster]
    cytosolic = np.mean(cluster_df['Cytosol']) 
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
sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(0, top_lum))
cbar = plt.colorbar(sm)
cbar.set_label('Cytosolic (nM)')

 
plt.plot([0,(116*3)],[-2,-2],c="purple", linewidth=5)
plt.plot([(116*3),204*3],[-2,-2],c="green", linewidth=5)
# Label axes
plt.xlabel('BP positions')
#plt.ylabel('Length')
plt.title(f"{protein} Cytosolic (nM)")
# Display the plot
plt.savefig(f"{file_path}/2new_cytosolic", dpi=600)
plt.show()
####################
# Create a scatter plot of cluster size vs average Cytosol/Full_Cell ratio
fig, ax = plt.subplots()
ax.scatter(cluster_sizes.values(), cluster_colors.values(), alpha=0.5)
ax.set_xlabel('Cluster Size')
ax.set_ylabel('Average Ratio of Cytosol/Full_Cell')
plt.title(protein)
plt.show()

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

