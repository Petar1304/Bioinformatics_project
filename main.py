import pandas as pd
import os
import io
import hdbscan
import numpy as np
import matplotlib.pyplot as plt

cnv_path = os.path.join(os.getcwd(), 'cnv.txt')
vcf_path = os.path.join(os.getcwd(), 'Strelka.vcf')

# cnv_path = '/sbgenomics/project-files/cnv.txt'
# vcf_path = '/sbgenomics/project-files/Strelka.vcf'
cnv_df = pd.read_csv(cnv_path, delimiter='\t')

filtered_regions = cnv_df[(cnv_df['cn'] == 2) & ((cnv_df['cn1'].isna()) | (cnv_df['cn1'] == 1)) & ((cnv_df['cn2'].isna()) | (cnv_df['cn2'] == 1))]
filtered_regions = filtered_regions[['chromosome', 'start', 'end', 'cn', 'cn1', 'cn2']]

# reading vcf file (found at: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744)
# PyVCF biblioteka ne radi
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

vcf_df = read_vcf(vcf_path)

# variant allele frequency
def vaf(ref, alt):
    return alt / (alt + ref)

# find value(s) of allelic depth for given base
def get_allelic_depth_str(base):
    return values[format.index(base + 'U')]

# convert value(s) from string to int, and sum if necessary
def convert(val):
    if val.find(',') != -1:
        res = 0
        val_list = val.split(',')
        for v in val_list:
            res += int(v)
        return res
    else:
        return int(val)

# looping over whole dataset

vafs = np.array([])
bases = np.array([])

for i, row in vcf_df.iterrows():
    for j, region in filtered_regions.iterrows():

        drop = True
        
        if (row['CHROM'] == region['chromosome'] and region['start'] <= row['POS'] <= region['end']):
            ref = row['REF']
            alt = row['ALT']
            format = row['FORMAT'].split(':')
            values = row['TUMOR'].split(':')
            if len(ref) != 1:
                ref_val = 0
                for b in ref:
                    ref_val_str = get_allelic_depth_str(b)
                    ref_val = convert(ref_val_str)
            else:
                ref_val = convert(get_allelic_depth_str(ref))
            
            if len(alt) != 1:
                alt_val = 0
                for b in alt:
                    alt_val_str = get_allelic_depth_str(b)
                    alt_val += convert(alt_val_str)
            else:
                alt_val = convert(get_allelic_depth_str(alt))

            vaf_val = vaf(ref_val, alt_val)
            bases = np.append(bases, ref_val + alt_val)
            vafs = np.append(vafs, vaf_val)

            if vaf_val <= 0.6:
                remove = False
            
    if drop:
        vcf_df.drop(i, inplace=True)

# Klasterizacija
# train_data = np.vstack([vafs, bases])
train_data = np.hstack([vafs.reshape(-1, 1), bases.reshape(-1, 1)])
clusterer = hdbscan.HDBSCAN(min_cluster_size=5, gen_min_span_tree=True)
cluster_labels = clusterer.fit_predict(train_data)

plt.figure()
plt.grid()
plt.title('Klasterizovani podaci')
plt.xlabel('VAF')
plt.ylabel('Broj baza')
# plt.scatter((train_data[0]).T, (train_data[1]).T, c=cluster_labels, cmap='viridis', s=50, alpha=0.7)
plt.scatter((train_data[:, 0]), (train_data[:, 1]), c=cluster_labels, cmap='viridis', s=50, alpha=0.7)
plt.show()

# Nalazenje centara klastera
cluster_centers = np.array([])
for cluster in np.unique(cluster_labels):
    cluster_members = train_data[cluster_labels == cluster]
    persistence_values = clusterer.probabilities_[clusterer.labels_ == cluster]
    representative_point = cluster_members[np.argmax(persistence_values)]
    cluster_centers = np.append(cluster_centers, representative_point)
cluster_centers = cluster_centers.reshape(-1, 2)


plt.figure()
plt.scatter(train_data[:, 0], train_data[:, 1], c=clusterer.labels_, cmap='viridis', s=50, alpha=0.7, label='Data Points')
plt.scatter(cluster_centers[:, 0], cluster_centers[:, 1], c='red', marker='X', s=100, label='Cluster Centers')
plt.title('HDBSCAN Clustering Result with Cluster Centers')
plt.xlabel('VAF (Variant Allele Frequency)')
plt.ylabel('Bases')
plt.legend()
plt.show()
