import sys
import subprocess
from Bio import SeqIO
import pandas as pd
import time

start_time = time.time()

fastafile = open("SM_chrm1.fa", 'r')
seqid = 'target1'

min_Tm = 58
max_Tm = 62
maxDiff_Tm = 5

for record in SeqIO.parse(fastafile, 'fasta'):
	sequence = str(record.seq)
seq1 = sequence[(1-1):(15400)]

#constrait for length of the search space 
#track the read length, target posn, primer lengths


def make_boulderio(id, seq):
	boulder= {
	"SEQUENCE_ID":id,
	"SEQUENCE_TEMPLATE":seq,
	"SEQUENCE_TARGET":"2389,1",
	"PRIMER_TASK":"generic",
	"PRIMER_OPT_SIZE":19,
	"PRIMER_MIN_SIZE":17,
	"PRIMER_MAX_SIZE":21,
	"PRIMER_MAX_NS_ACCEPTED":1,
	"PRIMER_MIN_TM": min_Tm,
	"PRIMER_MAX_TM": max_Tm,
	"PRIMER_PAIR_MAX_DIFF_TM": maxDiff_Tm,
	"PRIMER_PRODUCT_SIZE_RANGE":"250-300",
	"P3_FILE_FLAG":1,
	"SEQUENCE_INTERNAL_EXCLUDED_REGION":"2389,1",
	"PRIMER_EXPLAIN_FLAG":1,
	}
	boulderfile = id+".boulderio"
	with open(boulderfile, "w") as boulderf:
		for param in boulder.keys():
			boulderf.write(str(param+"="+(str(boulder[param])+"\n")))
		boulderf.write("="+"\n")
	return str(boulderfile)


def run_primer3(boulderfile):
	primer3out = subprocess.check_output(["primer3_core", boulderfile])
	return primer3out
	#print(primer3out)

def make_primerout_dict(primer3out):
	ID = []
	Left_primer = []
	Right_primer = []
	pair_penalty = []
	primeroutdict = {}
	primerinfolist = primer3out.decode().strip().split("\n")
	#print(type(primerinfolist))
	for line in primerinfolist:
		#if line.strip() == "=": next
		entry = line.strip().split("=")
		primeroutdict[entry[0]] = entry[1]
	ID.append(primeroutdict.get('SEQUENCE_ID'))
	Left_primer.append(primeroutdict.get('PRIMER_LEFT_0_SEQUENCE'))
	Right_primer.append(primeroutdict.get('PRIMER_RIGHT_0_SEQUENCE'))
	pair_penalty.append(primeroutdict.get('PRIMER_PAIR_0_PENALTY'))
	primersdf = pd.DataFrame({'ID':ID, 'Left_primer':Left_primer, 'Right_primer':Right_primer, 'pair_penalty':pair_penalty}, columns=["ID", "Left_primer", "Right_primer", "pair_penalty"])
	#write df to a file
	print(primersdf)
	return primeroutdict

#def check_primers()

def main():
	#outfile = create_output_file()
	boulderfile = make_boulderio(seqid, seq1)
	primer3out = run_primer3(boulderfile)
	primerdictout = make_primerout_dict(primer3out)
	#candidate_primers = check_primers(primerdictout)
	#cross = dimer_check(primerdict)

fastafile.close()

if __name__ == '__main__':
    try:
        main()

    except KeyboardInterrupt:
        print('Interrupted')

print("Runtime: " + str(round(time.time()- start_time, 3)) + " s")