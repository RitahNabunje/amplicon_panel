import gzip
import sys
import subprocess
from Bio import SeqIO
import time
import csv

fastafile = open("schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa", 'r')

chrms =['SM_V7_1']
for chrm in chrms:
	for record in SeqIO.parse(fastafile, "fasta"):
		if record.id == chrm:
			chrmSeq = record.seq
		
	with open(chrm+".tsv") as ampfile: #open corresponding file with amp regions
		reader = csv.DictReader(ampfile, delimiter='\t')
		for i, region in enumerate(reader): #each region is a line in the file
			regionID = chrm+'_r'+str(i)
			start = int(region['primer1_start'])
			end = int(region['primer2_end'])
			regionSeq= chrmSeq[start:end]
			target = str(int(region['target_snp'])-start+1)+','+region['gap']
			print(regionID, target)
