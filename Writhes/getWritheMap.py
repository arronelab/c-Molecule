import numpy as np
import sys
import os
import shutil

#read in the initial file
curveFile=sys.argv[1]
#open that file into numpy array
coords = np.loadtxt(curveFile)
#make a list of subsets of this curve
size=int(np.size(coords)/3)
print(size)
# we need at least 5 points for a reasonable writhe make
workingDirectory=os.getcwd()
os.mkdir(workingDirectory+"/tmp")
f3 = open(sys.argv[2]+"DIMap.dat", 'w')
for i in range(5,size):
    for j in range(size-i):
        subsets=coords[j:i+j]
        # write it to file
        subfile=("tmp/subcurve{:d}.dat".format(i-4))
        np.savetxt(subfile,subsets[:])
        fullloc3 =workingDirectory+"/DIwr "+subfile
        f3.write(str(j)+" "+str(i+j)+" "+os.popen(fullloc3).read()+"\n")
        os.remove(subfile)
shutil.rmtree(workingDirectory+"/tmp")
f3.close()
