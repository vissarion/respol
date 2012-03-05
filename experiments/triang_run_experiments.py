import os, sys


for d in range(7,10):
	for n in range(100,1000,50):
		cmd = "../triang " + str(d) + " " + str(n)   
		#print cmd
		os.system(cmd)
		#print
