import vcf
import pandas as pd
#import numpy as np
import time

start_time = time.time()

#parameters
#I think its worth testing with a big number here to check the logic!
a_size = 300 # length of target amplicon in bases
p_dist = 100  #length in bases to fit a primer (i.e 80-100 bases before target amplicon)
#vcf_file = "100k_SM_chrm1.vcf"
vcf_file = "1k_SM_chrm1.vcf"
#vcf_file = "50k_SM_chrm1.vcf"


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

#represented chrms in the vcf
chrms=[]
for chrm in snp_dict.keys():
	chrms.append(chrm)

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
print(dfs[0])
diff_time = time.time()
print(" Calculating diff btn snps: " + str(round(diff_time - extract_time, 3)) + " s")

#at some point you will also need to know the chromosome lengths to make sure you check the final window.
#chrms_length= read file with contig lengths

#extracting snps with p_dist between them for a primer
df1s=[]
for df in dfs:
    snp1 = []
    snp2 = []
    snp_count = []
    l_index = 0
    h_index = 1
    last_diffindex = -1
    last_index = len(df)-1
    pos1 = df.loc[l_index,'snp_pos']
    pos2 = df.loc[h_index,'snp_pos']
    diff = pos2 - pos1
    while h_index < last_index:
        if diff >=p_dist:
            if last_diffindex == -1:
                snp_count_left= -999
            snp_count_left = (l_index - last_diffindex+1)
            snp_count.append(snp_count_left)
            last_diffindex = h_index
            snp1.append(pos1)
            snp2.append(pos2)
        l_index = l_index + 1
        h_index = h_index + 1
        pos1 = df.loc[l_index,'snp_pos']
        pos2 = df.loc[h_index,'snp_pos']
        diff = pos2 - pos1
    df1 = pd.DataFrame({'snp1':snp1, 'snp2':snp2, 'snp_count':snp_count})
    df1s.append(df1)
print(df1s[0])
spacing_time = time.time()
print("spacing snps for a primer: "+ str(round(spacing_time - diff_time, 3)) + " s")

#my idea here is that rather than count SNP between each adjacent plausible 'primer sites'
#you might as well find plausible amplicons and then count SNPs in them
# so this code replaces the stuff commented out below
#loop through contigs

ampliconDfs = []
#df1s is a list of dataframes each for a different chromosome
for i, df1 in enumerate(df1s):
    # loop through primer sites
    # find primer pairs that look plausibly close
    # start and end of SNP-free regions that could become primer sites
    index_list = list(df1.index)
    ampliconsDf = pd.DataFrame(columns=['primer1_start', 'primer1_end', 'target_snp', 'primer2_start', 'primer2_end', 'gap', 'num_snps'])
    # loop through all first primer positions
    for primer1_index in index_list:
        # end of possible placements for primer 1
        start_primer1 = (df1.loc[primer1_index, 'snp1'] + 1)
        end_primer1 = (df1.loc[primer1_index, 'snp2'] - 1)
        primer2_index = primer1_index + 1
        # loop through primer 2 positions
        while primer2_index < index_list[-1]:

            start_primer2 = (df1.loc[primer2_index, 'snp1'] + 1)
            end_primer2 = (df1.loc[primer2_index, 'snp2'] - 1)
            if (start_primer2 - end_primer1 > a_size - (2 * p_dist)):
                # if the start of primer 2 region is too far, skip to the end of this loop
                primer2_index = index_list[-1]
            else:
                # otherwise, store this as  a possible amplicon
                target = df1.loc[primer1_index, 'snp2']
                gap = start_primer2 - (end_primer1 + 1)
                #n = [n for n in dfs[i].loc[:, 'snp_pos'] if
                    # df1.loc[primer1_index, 'snp2'] <= n <= df1.loc[primer2_index, 'snp1']]
                #num_snps = len(n)
                amplicon_info = {'primer1_start':start_primer1, 'primer1_end':end_primer1, 'target_snp':target, 'primer2_start':start_primer2, 'primer2_end':end_primer2, 'gap':gap, 'num_snps':num_snps}
                ampliconsDf = ampliconsDf.append(amplicon_info, ignore_index=True)
            primer2_index = primer2_index + 1
        ampliconDfs.append(ampliconsDf)
#print(ampliconDfs[0])
amplicon_time= time.time()
print("Amplicon saving : " + str(round(amplicon_time - spacing_time, 3)) + " s")

#save each chrm's amplicon regions to a file
for chrm, amplicon in zip(chrms, ampliconDfs):
	amplicon.to_csv(f'{chrm}_edited_amplicons.tsv', index=False, header=True, sep='\t', chunksize= 1000000)
save_file_time = time.time()
print("Saving to file: "+str(round(save_file_time - amplicon_time, 3)) + " s")

print("Runtime: " + str(round(time.time()- start_time, 3)) + " s")
print("That'\s All Folks!")