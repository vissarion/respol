import os, sys

if len(sys.argv) == 2:
	path = sys.argv[1]
	#print "running experiments for " + path + "..."	
else:
	print "no arguments given"
	exit(1)

cmd = "ls -rS " + path + "  > files"   
os.system(cmd)
fres = open("files", 'r')
input_files=[]
for line in fres:
  input_files.append(line)

for filename in input_files:
	#print str(filename)
	#for i in range(10):  
	cmd = "../res_enum_d < " + path + "/" + str(filename)
	os.system(cmd)
	#cmd = "polymake test_bb.polymake"
	#os.system(cmd)
	cmd = "polymake test_cdd.polymake"
	os.system(cmd)
	cmd = "polymake test_lrs.polymake"  
	os.system(cmd)
	#os.system('echo')
	#print cmd
	
	#print
	#print
