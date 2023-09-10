import glob
import pandas as pd
from shutil import copyfile
import gzip
import subprocess

# wgs_reference =
dragen_ref = 
samtools_path = 
path_script = 

df = pd.read_csv(
    '/home/jsims/viperInversions/result_ids.csv',
    sep=",")
print(df)
print('AppSession,BND,INV')
for index, row in df.iterrows():
    app = str(row['result_id'])
    print(
        'XXXXXXXXXXXXXXXXXXXXXXXX/.ResourceById/AppSession/' + app + '/AppResults.*/Files/*.sv.vcf.gz')

    sv_file = glob.glob(
        'XXXXXXXXXXXXXXXXXXX/.ResourceById/AppSession/' + app + '/AppResults.*/Files/*.sv.vcf.gz')
    # convertInversion script was failing with vcf.gz file, therefore vcf.gz files were copied and gunzipped to separate breakends from inversions
    copyfile(sv_file[0],
             '/home/jsims/Inversions/copybin/sv_' + app + '.vcf.gz')

    # gunzip them
    cmd = "gunzip " + '/home/jsims/Inversions/copybin/sv_' + app + '.vcf.gz'
    run_gz = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE).communicate()

    gunzipped_sv_file = '/home/jsims/Inversions/copybin/sv_' + app + '.vcf'
    subprocess.call(
        "python " + path_script + ' ' + samtools_path + ' ' + dragen_ref + ' ' + gunzipped_sv_file + '>' + app + '_after_conversion.vcf',
        shell=True)

    print('script complete')
    # create 2 vcf files ones for the inversions and one for the breakends only for the PASS variants
    subprocess.call(
        "bcftools view -f PASS -i 'SVTYPE=\"BND\"' " + app + "_after_conversion.vcf" + ">" + app + "_after_conversion_bnd.vcf",
        shell=True)
    subprocess.call(
        "bcftools view -f PASS -i 'SVTYPE=\"INV\"' " + app + "_after_conversion.vcf" + ">" + app + "_after_conversion_inv.vcf",
        shell=True)

    # count_bnd = subprocess.Popen('bcftools stats ' + app + '_after_conversion_bnd.vcf |grep "number of records:"',shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    # count_inv = subprocess.Popen('bcftools stats ' + app + '_after_conversion_inv.vcf |grep "number of records:"',shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    # no_bnd = (count_bnd[0].decode().split('\t')[3][0:-1])
    # no_inv = (count_inv[0].decode().split('\t')[3][0:-1])

    # print(app + ',' + no_bnd + ',' + no_inv )
