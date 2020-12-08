import vcf
import pandas as pd
import numpy as np
import time

start_time = time.time()

#parameters
#I think its worth testing with a big number here to check the logic!
a_size = 300 # length of target amplicon in bases
p_dist = 100  #length in bases to fit a primer (i.e 80-100 bases before target amplicon)
#vcf_file = "100k_SM_chrm1.vcf"
#vcf_file = "1k_SM_chrm1.vcf"
#vcf_file = "50k_SM_chrm1.vcf"
vcf_file = "/Users/rnabunje/Projects/SchistoAmpliconPanel/data/5_sample.vcf"
#vcf_file = "/Users/rnabunje/Projects/SchistoAmpliconPanel/data/SM_chrm1.vcf"
#vcf_file = "CRELLEN_FULLFILTER.198.vcf"

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

spaced_dfs = []
for chrm, pos in snp_dict.items():
    chrm_df = pd.DataFrame.from_dict((list(filter(None, pos))))
    chrm_df.columns=['snp_pos']
    chrm_df.loc[-1]=[0] #make 0 the first snp position such that the first snp is captured
    chrm_df.index = chrm_df.index+1
    chrm_df=chrm_df.sort_index().reset_index(drop=True)
    l_index = 0
    h_index = 1
    snp1 = []
    snp2 = []
    snp_count = [-999]
    last_index = len(chrm_df)-1
    last_diffindex = -1
    s1 = chrm_df.loc[l_index,'snp_pos']
    s2 = chrm_df.loc[h_index,'snp_pos']
    diff = s2 - s1
    while h_index < last_index:
        if diff >= p_dist:
            if last_diffindex == -1:
                snp_count_left = -999
            else:
                snp_count_left = (l_index - last_diffindex + 1)
                snp_count.append(snp_count_left)
            last_diffindex = h_index
            snp1.append(s1)
            snp2.append(s2)
        l_index = l_index + 1
        h_index = h_index + 1
        s1 = chrm_df.loc[l_index, 'snp_pos']
        s2 = chrm_df.loc[h_index, 'snp_pos']
        diff = s2 - s1
    spaced_df = pd.DataFrame({'snp1':snp1, 'snp2':snp2, 'snp_count':snp_count})
    spaced_dfs.append(spaced_df)
#print(spaced_dfs[1])
spacing_time = time.time()
print ("spacing snps for a primer: "+ str(round(spacing_time - extract_time, 3)) + " s")

ampliconDfs = []
for i, spaced_df in enumerate(spaced_dfs):
    index_list = list(spaced_df.index)
    ampliconDf = pd.DataFrame(columns=['primer1_start', 'primer1_end', 'target_snp', 'primer2_start', 'primer2_end', 'gap', 'num_snps'])
        # loop through all first primer positions
    for primer1_index in index_list:
            # end of possible placements for primer 1
        start_primer1 = (spaced_df.loc[primer1_index, 'snp1'] + 1)
        end_primer1 = (spaced_df.loc[primer1_index, 'snp2'] - 1)
        primer2_index = primer1_index + 1
        A_snps = 0
        # loop through primer 2 positions
        while primer2_index < index_list[-1]:
            A_snps = A_snps + spaced_df.loc[primer2_index, 'snp_count']
            start_primer2 = (spaced_df.loc[primer2_index, 'snp1'] + 1)
            end_primer2 = (spaced_df.loc[primer2_index, 'snp2'] - 1)
            if (start_primer2 - end_primer1 > a_size - (2 * p_dist)):
            # if the start of primer 2 region is too far, skip to the end of this loop
                primer2_index = index_list[-1]
            else:
            # otherwise, store this as  a possible amplicon
                target = spaced_df.loc[primer1_index, 'snp2']
                gap = start_primer2 - (end_primer1 + 1)
                amplicon_info = {'primer1_start':start_primer1, 'primer1_end':end_primer1, 'target_snp':target, 'primer2_start':start_primer2, 'primer2_end':end_primer2, 'gap':gap, 'num_snps':A_snps}
                ampliconDf = ampliconDf.append(amplicon_info, ignore_index=True)
            primer2_index = primer2_index + 1
    ampliconDfs.append(ampliconDf)
#print(ampliconDfs[1])
amplicon_time = time.time()
print("Amplicon saving : " + str(round(amplicon_time - spacing_time, 3)) + " s")

#save each chrm's amplicon regions to a file
for chrm, ampliconDf in zip(chrms, ampliconDfs):
	ampliconDf.to_csv(f'{chrm}_7_12.tsv', index=False, header=True, sep='\t', chunksize= 1000000)
saving_time = time.time()
print("Saving the amplicon info to a file: "+str(round(saving_time - amplicon_time, 3))+"s")

print("Runtime: " + str(round(time.time()- start_time, 3)) + " s")
print("That'\s All Folks!")