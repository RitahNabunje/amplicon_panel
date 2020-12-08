import subprocess
from Bio import SeqIO
import time
import csv
import pandas as pd

start_time = time.time()

a_size = 300
primerOpt = 26
primerMin = 24
primerMax = 32
primerNS = 1
product_size = "%i-%i" % ((a_size - 50), a_size)
min_Tm = 58
max_Tm = 62
maxDiff_Tm = 5

chrms = ['SM_V7_1']
# chrms = ['SM_V7_1', 'SM_V7_ZW', 'SM_V7_3', 'SM_V7_2', 'SM_V7_4']
fastafile = "/Users/rnabunje/Projects/SchistoAmpliconPanel/data/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa"

#blast parameters
evalue = str(30000)
fmt = str(6) # no header
wordsize = str(7)
output_file = "blastout_x812"
DB = "genomeDB"

# format DB for the BLAST runs
make_DB = ['makeblastdb',
           '-in', fastafile,
           '-dbtype', 'nucl',
           '-out', DB,
           '-parse_seqids']
subprocess.call(make_DB)

# make primers for the amplicons, run the primers against the fastafile, pick only the specific primers
for chrm in chrms:
	# parse fastafile for each chrom's seq
	for record in SeqIO.parse(fastafile, "fasta"):
		if record.id == chrm:
			chrmSeq = record.seq
			print(record.id)
			# open corresponding file with amplicon
			#with open(chrm+"_new.tsv") as ampfile:
			with open("bits_chrm1.tsv") as ampfile:
				reader = csv.DictReader(ampfile, delimiter='\t')
				# each Amplicon is a line in the file
				for i, amplicon in enumerate(reader):
					ampliconID = chrm+'_A'+str(i)
					start = int(amplicon['primer1_start'])-1
					end = int(amplicon['primer2_end'])
					ampliconSeq = chrmSeq[start:end]
					target = str(int(amplicon['target_snp'])-start+1)+','+amplicon['gap']
					#target is relative to the ampliconSeq, includes the first snp and gap which includes the remaining gaps

					#define contents of boulderfile as a dict
					boulder = {
						"SEQUENCE_ID": ampliconID,
						"SEQUENCE_TEMPLATE": ampliconSeq,
						"SEQUENCE_TARGET": target,
						"PRIMER_TASK": "generic",
						"PRIMER_OPT_SIZE": primerOpt,
						"PRIMER_MIN_SIZE": primerMin,
						"PRIMER_MAX_SIZE": primerMax,
						"PRIMER_MAX_NS_ACCEPTED": primerNS,
						"PRIMER_PRODUCT_SIZE_RANGE": product_size,
						"PRIMER_MIN_TM": min_Tm,
						"PRIMER_MAX_TM": max_Tm,
						"PRIMER_PAIR_MAX_DIFF_TM": maxDiff_Tm,
						"P3_FILE_FLAG": 1,
						"SEQUENCE_INTERNAL_EXCLUDED_REGION": target,
						"PRIMER_EXPLAIN_FLAG": 1,
					}
					# each amplicon has a separate boulder file
					boulderfile = ampliconID+".boulderio"
					# write the boulder contents to the amplicon's file
					with open(boulderfile, "w") as boulderf:
						for param in boulder.keys():
							# add "=" btn parameter and value, then "newline" after each parameter
							boulderf.write(str(param+"="+(str(boulder[param])+"\n")))
						boulderf.write("="+"\n")
					#Run primer3_core via commandline (outside script)
					primer3out = subprocess.check_output(["primer3_core", boulderfile])
					primeroutdict = {}
					primerinfolist = primer3out.decode().strip().split("\n")
					for line in primerinfolist:
						# if line.strip() == "=": next
						entry = line.strip().split("=")
						primeroutdict[entry[0]] = entry[1]
					if primeroutdict.get('PRIMER_PAIR_NUM_RETURNED=') != 0:
						file = open(ampliconID+"_primers.fa", "w")
						
						file.close()
                    #run blast with query as the amplicon primers
                    querySeq = ampliconID+"_primers.fa"
					query_line = ['blastn',
                                  '-query', querySeq,
                                  '-word_size', wordsize,
                                  '-evalue', evalue,
                                  '-max_hsps', '1',
                                  '-out', blast_output,
                                  '-outfmt', fmt,
                                  '-db', DB]
					subprocess.call(query_line)
                    hits_df = pd.read_csv(blast_output, delimiter='\t')
                    hits_df.columns = ["query", "subject", "identity", "alignment_length", "mismatches",
                                           "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit_score"]

                    print(hits_df.head())
                    identical_df = hits_df[hits_df['identity'] == 100.0]

					print(ampliconID, ": done")
					subprocess.call(["rm", ampliconID+".for"])
					subprocess.call(["rm", ampliconID+".rev"])
					subprocess.call(["rm", ampliconID+".boulderio"])
			ampfile.close()
print("Runtime: " + str(round(time.time() - start_time, 3)) + " s")
