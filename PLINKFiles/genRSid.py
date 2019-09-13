import random
import string

#generate random rsIDs for all snps
l=8818790 #number of snps in .bim 
for i in range (l):
	value=''.join(random.choices(string.digits, k=14)) #mke random combinations of 14 digits
	print("rs"+ "" + value)
