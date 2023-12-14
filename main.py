import pandas as pd
import os

cnv_path = os.path.join(os.getcwd(), 'cnv.txt')
vcf_path = os.path.join(os.getcwd(), 'Strelka.vcf')

cnv = pd.read_csv(cnv_path, delimiter='\t')
print(cnv)