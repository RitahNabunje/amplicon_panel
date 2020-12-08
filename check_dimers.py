import subprocess
import pandas as pd
import primer3


regionID="trial_check"
regionSeq="AGTATTTTAAGAGTATTTGCTTTCTCCTGGATATTCCAGAATACTACTGAACATAGTATTGATTCCTGATATTTCTGAAACTCAGTGAAAGGAGATGCGCTGGAAGTACCAGTGATCTTTGTCGGGAATCCATATGTCGTTTGTTCAGATGTGTGCATTACATCTTACTTATCTAAGCAAGGATGTGTATTATAAATTATAAAAAATGGGACATCTGTGTCAATCAGAGAAGATCAACGACTATGCGAAATAAATTTAGGGGTCTTATCACAGATGAAGCATGTAAAAAGTGACCACATAAATAACCCAAGAATAAGGAGTGCTGTAAATTGTTACCAAGAAATAAAACCTAAGACGATGTGTATAACTAAATAAGATAGATGGTTTGTTAGAACACCGCATCGATCGTTCTACAAAACGAATCTTGAAAAATACTTCAAACTACTATGTAGTGAATAGTAGTGTGGTGTTTTTGTACAGAAAAAGCGCAAGAAGAGGTTTACAAAACAAGAGGCAATGAATTTGGGCTTGTTATTCAAATCAAGGAAATGGGCTAATGATAATCTCATGGAGAGTGACGACAGTCGAGCTTTTCAATGTCTCCATCATGATTTACGGTCAATACGTCATTAACAGTGACGAAGAAAGTTTATTGTTTCTATGAGTTAACAAATGGTAATGAATATCAATATCTGATGTTGAGAAAATAAAATAGAGAGGTCAATTTTGTTTTGATTCACTCAATCAACTAACTGAGATCTATTGTTTGAGATATAGGCAACTCATTACCACACCTCTTTGTGAATATGATGTCTAAATTAGCTAAAATAGTTAGTAGACGTTTTGGCTTATGTTTATCGAAACTTGTTGGTAAAATATTGATGAATATTTAGTGAAATTAGAAGTGGACCGAGTTCTTCCTTATG"
target="502,45"
product_size="250-300"
primerOpt = 19
primerMin = 17
primerMax = 21
primerNS = 1
min_Tm = 58
max_Tm = 62
maxDiff_Tm = 5
ID = []
Left_primer = []
Right_primer = []

boulder= {
    "SEQUENCE_ID":regionID,
    "SEQUENCE_TEMPLATE":regionSeq,
    "SEQUENCE_TARGET":target,
    "PRIMER_TASK":"generic",
    "PRIMER_OPT_SIZE":primerOpt,
    "PRIMER_MIN_SIZE":primerMin,
    "PRIMER_MAX_SIZE":primerMax,
    #"PRIMER_MIN_TM":min_Tm,
    #"PRIMER_MAX_TM":max_Tm,
    #"PRIMER_PAIR_MAX_DIFF_TM":maxDiff_Tm,
    "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT":1
    "PRIMER_MAX_NS_ACCEPTED":primerNS,
    "PRIMER_PRODUCT_SIZE_RANGE":product_size,
    "P3_FILE_FLAG":1,
    "SEQUENCE_INTERNAL_EXCLUDED_REGION":target,
    "PRIMER_EXPLAIN_FLAG":1,
    }
boulderfile = regionID+".boulderio"
# write the boulder contents to the region's file
with open(boulderfile, "w") as boulderf:
    for param in boulder.keys():
    # add "=" btn parameter and value, then "newline" after each parameter
        boulderf.write(str(param+"="+(str(boulder[param])+"\n")))
    boulderf.write("="+"\n")
    #Run primer3_core via commandline (outside script)
primer3out = subprocess.check_output(["primer3_core", boulderfile])
#primer3out = subprocess.check_output(["primer3_core", "-output=trial_check.out", boulderfile])

primer3outdict = {}
primerinfolist = primer3out.decode().strip().split("\n")
for line in primerinfolist:
    entry = line.strip().split("=")
    primer3outdict[entry[0]] = entry[1]
left_primer = primer3outdict.get('PRIMER_LEFT_0_SEQUENCE')
right_primer = primer3outdict.get('PRIMER_RIGHT_0_SEQUENCE')

boulder2 = {
    "SEQUENCE_ID":regionID,
    "SEQUENCE_PRIMER":left_primer,
    "SEQUENCE_PRIMER_REVCOMP":right_primer,
    "PRIMER_PICK_RIGHT_PRIMER":0,
    "PRIMER_PICK_LEFT_PRIMER":0,
    "PRIMER_PICK_ANYWAY":1,
    }
boulderfile2 = regionID+"2.boulderio"
# write the boulder contents to the region's file
with open(boulderfile2, "w") as boulderf2:
    for param in boulder2.keys():
        boulderf2.write(str(param+"="+(str(boulder2[param])+"\n")))
    boulderf2.write("="+"\n")
    #Run primer3_core via commandline (outside script)
check_primersout = subprocess.check_output(["primer3_core", boulderfile2])

print(check_primersout)
check_primersoutdict = {}
check_primersinfo = check_primersout.decode().strip().split("\n")
for line in check_primersinfo:
    check_entry = line.strip().split("=")
    check_primersoutdict[check_entry[0]] = check_entry[1]
    print(line)

#ID.append(primeroutdict.get('SEQUENCE_ID'))
#Left_primer.append(primeroutdict.get('PRIMER_LEFT_0_SEQUENCE'))
#Right_primer.append(primeroutdict.get('PRIMER_RIGHT_0_SEQUENCE'))

#primersdf = pd.DataFrame({'ID':ID, 'Left_primer':Left_primer, 'Right_primer':Right_primer}, columns=["ID", "Left_primer", "Right_primer"])
#print(primersdf)
