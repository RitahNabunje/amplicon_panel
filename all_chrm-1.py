import vcf
import pandas as pd
import numpy as np
import time

start_time = time.time()

#parameters
#I think its worth testing with a big number here to check the logic!
a_size = 300  # length of target amplicon in bases
p_dist = 100  #length in bases to fit a primer (i.e 80-100 bases before target amplicon)
#vcf_file = "100k_SM_chrm1.vcf"
#vcf_file = "1k_SM_chrm1.vcf"
vcf_file = "50k_SM_chrm1.vcf"
#vcf_file = "5_sample.vcf"

vcf_reader = vcf.Reader(open(vcf_file, 'r'))
snp_list=[]
for record in vcf_reader:
        record = (record.CHROM, record.POS) #extract snp chrom and Position info from vcf
        snp_list.append(record) #save snp info to list
extract_time = time.time()
print("extracting snp records : " + str(round(extract_time - start_time, 3)) + " s")

snp_dict = dict() #create empty snp dictionary

#convert snp_list a dict and save it to snp_dict
for chrm, pos in snp_list:
        snp_dict.setdefault(chrm, []).append(pos)
#print (snp_dict)
#print("snp dict : done")

#represented chrms
chrms=[]
for chrm in snp_dict.keys():
	chrms.append(chrm)
#print(type(chrms), chrms)

#number of snps represented in the vcf record
#for chrm,pos in snp_dict.items():
 #       snps = len(list(filter(None, pos)))

#make df for each chrm/contig represented in the vcf
dfs=[]
for chrm, pos in snp_dict.items():
	df = pd.DataFrame.from_dict((list(filter(None, pos))))
	df.columns=['snp_pos']
	df.loc[-1]=[0] #make 0 the first snp position such that the first snp is captured
	df.index = df.index+1
	df=df.sort_index().reset_index(drop=True)
	df['diff'] = df['snp_pos'].diff()
	dfs.append(df)
#at some point you will also need to know the chromosome lengths to
#make sure you check the final window.

#chrms_length= read file with contig lengths
diff_time = time.time()
print(" making diff df: " + str(round(diff_time - extract_time, 3)) + " s")

#extracting snps with p_dist between them for a primer
df1s=[]
for df in dfs:
	snp1 =[]
	snp2=[]
	l_index = 0
	h_index = 1
	last_index = len(df)-1
	pos1 = df.loc[l_index,'snp_pos']
	pos2 = df.loc[h_index,'snp_pos']
	diff = pos2 - pos1
	while h_index < last_index:
        	if diff >=p_dist:
                	snp1.append(pos1)
                	snp2.append(pos2)
        	l_index = l_index + 1
        	h_index = h_index + 1
        	pos1 = df.loc[l_index,'snp_pos']
        	pos2 = df.loc[h_index,'snp_pos']
        	diff = pos2 - pos1
	df1 = pd.DataFrame({'snp1':snp1, 'snp2':snp2})
	df1s.append(df1)
#print(df1s[1])
spacing_time = time.time()
print("extracting snps with enough primer space : "+ str(round(spacing_time - diff_time, 3)) + " s")

#my idea here is that rather than count SNP between each adjacent plausible 'primer sites'
#you might as well find plausible amplicons and then count SNPs in them
# so this code replaces the stuff commented out below
#loop through contigs
amplicons=[]
for i, df1 in enumerate(df1s):
#loop through primer sites
#find primer pairs that look plausibly close
	#start and end of SNP-free regions that could become primer sites
	primer1_starts = []
	primer1_ends = []
	target_snp =[]
	primer2_starts = []
	primer2_ends = []
	#size of gap between end of primer 1 region and start of primer 2 region
	gap = []
	#number of SNPs in this gap
	num_snps = []
	#index_list = list(df1.index) replace in the next loop!
	#loop through all first primer positions
	for primer1_index in list(df1.index):
		#end of possible placements for primer 1
		start_primer1 = (df1.loc[primer1_index, 'snp1'] +1)
		end_primer1 = (df1.loc[primer1_index,'snp2'] -1)
		primer2_index = primer1_index+1
		#loop through primer 2 positions
		#
		while primer2_index < list(df1.index)[-1]:
			start_primer2 = (df1.loc[primer2_index,'snp1'] +1)
			end_primer2 = (df1.loc[primer2_index, 'snp2']-1)
			if (start_primer2 - end_primer1 > a_size - (2*p_dist) ):
			#if the start of primer 2 region is too far, skip to the end of this loop
				primer2_index = list(df1.index)[-1]
			else:
			#otherwise, store this as  a possible amplicon
				primer1_starts.append(start_primer1)
				primer1_ends.append(end_primer1)
				target_snp.append(df1.loc[primer1_index,'snp2'])
				primer2_starts.append(start_primer2)
				primer2_ends.append(end_primer2)
				gap.append(start_primer2 - (end_primer1+1) )
				n=[n for n in dfs[i].loc[:,'snp_pos'] if df1.loc[primer1_index,'snp2'] <= n <=  df1.loc[primer2_index,'snp1']]
				num_snps.append(len(n))
			primer2_index=primer2_index+1
	amplicon_contig = pd.DataFrame({'primer1_start':primer1_starts,'primer1_end':primer1_ends,'target_snp':target_snp,'primer2_start':primer2_starts,'primer2_end':primer2_ends,'gap':gap,'num_snps':num_snps},columns =['primer1_start', 'primer1_end','target_snp', 'primer2_start','primer2_end','gap','num_snps'])
	amplicons.append(amplicon_contig)
#print(amplicons[1])
amplicon_time= time.time()
print("Amplicon saving : " + str(round(amplicon_time - spacing_time, 3)) + " s")

print("All finished! Now saving amplicons to file")
#save each chrm's amplicon regions to a file
for chrm, amplicon in zip(chrms, amplicons):
	amplicon.to_csv(f'{chrm}.tsv', index=False, header=True, sep='\t')


#note at this stage I just know where I *could* put a primer..
#you will want to extract smaller regions for primer3: but I think its good to have the maximum
#in case you need to relax your rules a bit for primer3

#now identify regions to pass to primer 3:
#I think this will be a_size - gap at the end of primer 1 region and start of primer 2 region

#NB : I think it should be possible to add information to the existing 'amplicons' data table from this point on, rather than
#creating lots of new data tables.

print("Runtime: " + str(round(time.time()- start_time, 3)) + " s")
print("That'\s All Folks!")