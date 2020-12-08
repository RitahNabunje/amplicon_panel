import primer3
import time

start_time = time.time()

regionID = "SM_V7_1_r0"
left = "CCCCCAAACACAACATCGTC"
right = "ACATGACGTTTTGGGACAGT"

cross = primer3.calcHeterodimer(left, right)
print(cross.dg)

self1 = primer3.calcHomodimer(left)
print(self1.dg)

self2= primer3.calcHomodimer(right)
print(self2.dg)

print("Runtime: " + str(round(time.time()- start_time, 3)) + " s")