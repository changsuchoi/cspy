import os
import glob


ch="Calib-*.fits"
inlist=glob.glob(ch)

newch="sa"
complist = glob.glob(newch+ch)


newfiles=[]
for i in inlist :
	if newch+i not in complist : 
		print (i)
		newfiles.append(i)
