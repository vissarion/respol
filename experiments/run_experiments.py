import os, sys

if len(sys.argv) == 2:
	path = sys.argv[1]
	#print "running experiments for " + path + "..."	
else:
	print "no arguments given"
	exit(1)

cmd = "ls -Sr " + path + "  > files"   
os.system(cmd)
fres = open("files", 'r')
input_files=[]
for line in fres:
  input_files.append(line)

for filename in input_files:
	cmd = "../res_enum_d < " + path + "/" + str(filename)   
	#print cmd
	print str(filename),
	os.system(cmd)
	print
	#print
