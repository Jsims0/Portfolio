from Bio.Seq import Seq
from Bio import SeqIO
import regex

REFERENCE = '/XXXXXXXXXXXXXXXXXX/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/' \
            'WholeGenomeFasta/genome.fa'

SITE = 'GTTTAAAACTGTGCG'

record = SeqIO.index(REFERENCE, 'fasta')
chr_20 = record['chr20'].seq
reg = regex.search(str(chr_20), SITE)

#print(reg)
