import subprocess
from Bio import SeqIO
import time
import csv
import pandas as pd

start_time = time.time()

#constrait for length of the search space 
#track the read length, target position, primer lengths

a_size = 300
primerOpt = 26
primerMin = 24
primerMax = 32
primerNS = 1
product_size = "%i-%i" % ((a_size-50), a_size)
min_Tm = 58
max_Tm = 62
maxDiff_Tm = 5

chrms = ['SM_V7_1']
#chrms = ['SM_V7_1', 'SM_V7_ZW', 'SM_V7_3', 'SM_V7_2', 'SM_V7_4']

fastafile = "/Users/rnabunje/Projects/SchistoAmpliconPanel/data/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa"
DB = "genomeDB"

left = []
right = []
ID = []
left_primer = []
right_primer = []
left_ID = []
right_ID = []

#format DB for the BLAST runs

make_DB = ['makeblastdb',
           '-in', fastafile,
           '-dbtype', 'nucl',
           '-out', DB,
           '-parse_seqids']
#subprocess.call(make_DB)

for chrm in chrms:
	# parse fastafile for each chrom's seq
	for record in SeqIO.parse(fastafile, "fasta"):
		if record.id == chrm:
			chrmSeq = record.seq
			print(record.id)
			# open corresponding file with amp regions
			with open(chrm+"_7_12.tsv") as ampfile:
			#with open("bits_chrm1.tsv") as ampfile:
				reader = csv.DictReader(ampfile, delimiter='\t')
				# each Amplicon region is a line in the file
				for i, region in enumerate(reader):
					regionID = chrm+'_A'+str(i)
					start = int(region['primer1_start'])-1
					end = int(region['primer2_end'])
					regionSeq = chrmSeq[start:end]
					target = str(int(region['target_snp'])-start+1)+','+region['gap']
					#target is relative to the regionSeq, includes the first snp and gap which includes the remaining gaps

					#define contents of boulderfile as a dict
					boulder = {
						"SEQUENCE_ID": regionID,
						"SEQUENCE_TEMPLATE": regionSeq,
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
					# each region has a separate boulder file
					boulderfile = regionID+".boulderio"
					# write the boulder contents to the region's file
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
					ID.append(primeroutdict.get('SEQUENCE_ID'))
					left.append(primeroutdict.get('PRIMER_LEFT_0_SEQUENCE'))
					right.append(primeroutdict.get('PRIMER_RIGHT_0_SEQUENCE'))
					left_ID.append(primeroutdict.get('SEQUENCE_ID')+"_PRIMER_LEFT_0")
					right_ID.append(primeroutdict.get('SEQUENCE_ID')+"_PRIMER_RIGHT_0")
					left_primer.append(primeroutdict.get('PRIMER_LEFT_0_SEQUENCE'))
					right_primer.append(primeroutdict.get('PRIMER_RIGHT_0_SEQUENCE'))
					print(regionID, ": done")
					subprocess.call(["rm", regionID+".for"])
					subprocess.call(["rm", regionID+".rev"])
					subprocess.call(["rm", regionID+".boulderio"])
					ampfile.close()
primersdf = pd.DataFrame({'ID': ID, 'Left': left, 'Right': right}, columns=["ID", "Left", "Right"])
print(primersdf)

#write forward primers to file
#lfile = open("for_primers.fa", "w")
#for i in range(len(left_ID)):
		#	if str(left_primer[i]) != 'None':
		#lfile.write(">" + str(left_ID[i]) + "\n" + str(left_primer[i]) + "\n")
#file.close()

#write reverse primers to file
#rfile = open("rev_primers.fa", "w")
#for i in range(len(right_ID)):#
#	if str(right_primer[i]) != 'None':
#		rfile.write(">" + str(right_ID[i]) + "\n" + str(right_primer[i]) + "\n")
#rfile.close()


print("Runtime: " + str(round(time.time() - start_time, 3)) + " s")
