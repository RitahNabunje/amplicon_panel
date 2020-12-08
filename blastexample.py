#from Bio.Blast.Applications import NcbiblastpCommandline
import subprocess
import time
#import csv
import pandas as pd
start_time = time.time()

#parameters
queryseq = "/Users/rnabunje/Projects/SchistoAmpliconPanel/data/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa"
#DBfile = "/Users/rnabunje/Projects/SchistoAmpliconPanel/data/5_sample_primers.fa"
DBfile = "/Users/rnabunje/Projects/SchistoAmpliconPanel/primers.fa"
#evalue = str(30000)
#evalue = str(1e-30)
#output_fmt = str(7) #with explanations, no. of hits, fields
output_fmt = str(6) # no header
wordsize = str(7)
# 'word_size', '7',
#'-evalue', evalue,
max_target_seqs = str(50000)
penalty = str(-1)
reward = str(1)
#max_hits = '1'
output_file = "blast_out_3_s12"

#Formatting  Database
make_DB = ['makeblastdb',
           '-in', DBfile,
           '-dbtype', 'nucl',
           '-out', 'primersDB',
           '-parse_seqids']
subprocess.call(make_DB)

#'-task', 'blastn-short'
#Querying database (removed the out file '-out','blast_out')
query_line = ['blastn',
              '-query', queryseq,
              '-word_size', wordsize,
              '-max_hsps', '1',
              '-max_target_seqs', max_target_seqs,
              '-out', output_file,
              '-outfmt', output_fmt,
              '-db', 'primersDB']
subprocess.call(query_line)
#blastout = subprocess.check_output(query_line)
#blastoutinfo = blastout.decode().strip().split("\n")

#read output file
hits_df = pd.read_csv(output_file, delimiter='\t')
hits_df.columns = ["query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",
                   "q.start", "q.end", "s.start", "s.end", "evalue", "bit_score"]
print(hits_df.head())
print(hits_df.iloc[:, 1].value_counts(ascending=True))

primers = hits_df.subject.unique()

for primer in primers:
    primerdf = hits_df[hits_df['subject'] == primer]
    id_df = primerdf[primerdf['identity'] == 100.0]
    identicaldf = id_df[id_df['evalue'] < math.exp(-10)]
    print(primer, len(id_df), len(identicaldf), len(primerdf))

#those with hits whose evalue is less than e-10
sig_df = hits_df[hits_df['evalue'] < math.exp(-10)]
print(sig_df['subject'].value_counts(ascending=True))

#hits_df['identity'].value_counts(normalize=True)

print("Runtime: " + str(round(time.time() - start_time, 3)) + " s")
print("That's All Folks!")
