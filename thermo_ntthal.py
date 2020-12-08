import csv

#check melting tm and deltaG values ("dG" and "t" in ntthal output)
max_Tm = 62 # max primer Tm
Tm_thresh = max_Tm -10 #Primer3 suggests secondary structures melting temperatures should
                        # be at least 10°C below the primer melting temperature.
dG_thresh = -10.000 #in kcal/mol more negative values suggest stronger binding

#simulate: all forward - all forward, r-r, f-r

# first try with primers saved in a csv
# the full_run.py provides primersdf, plan to work with a df, not a csv
with open(primersout.txt) as primerfile:
    primers = csv.DictReader(primerfile, delimiter='\t')
    for i, primer in enumerate(primers):
        #forward with all forward
        primer1 = primer["Left_primer"]
        primer2 = primer["Right_Primer"]

        crossD_check = subprocess.check_output(["ntthal -s1", primer1, "-s2", primer2])
        if cross_check(dG) <= dG_thresh & cross_check(t) <= Tm_thresh:

        # forward with all reverse
        primer1 = primer["Right_primer"]
        primer2 =
        cross_check = subprocess.check_output(["ntthal -s1", primer1, "-s2", primer2])
        if cross_check(dG) <= dG_thresh & cross_check(t) <= Tm_thresh:

        #Reverse with all reverse
        primer1 =
        primer2 =

