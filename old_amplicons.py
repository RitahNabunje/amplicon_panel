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
	index_list = list(df1.index)
	
	#loop through all first primer positions
	for primer1_index in index_list:
		#end of possible placements for primer 1
		start_primer1 = (df1.loc[primer1_index, 'snp1'] +1)
		end_primer1 = (df1.loc[primer1_index,'snp2'] -1)
		primer2_index = primer1_index+1
		#loop through primer 2 positions
		#
		while primer2_index < index_list[-1]:
			start_primer2 = (df1.loc[primer2_index,'snp1'] +1)
			end_primer2 = (df1.loc[primer2_index, 'snp2']-1)
			if (start_primer2 - end_primer1 > a_size - (2*p_dist) ):
			#if the start of primer 2 region is too far, skip to the end of this loop
				primer2_index = index_list[-1]
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