import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_Mouse.97.py > RunMouseCDHIT.97.log 2>&1 & 

d = datetime.date.today()
percentID = 0.97
finalOutputPath = "MouseCDHit/CDHITcompare.Mouse." + str(d) 
finalSifFilePath = "MouseCDHit/CDHITcompare.Mouse.cytoscape." + str(percentID * 100) + "percent." + str(d)
fastaPath = "/home/alexs/FastaFiles/uniProtIDsPerPaper.2012-06-12.tsv"
pmidListFile = "MouseCDHit/MousePMIDs.txt"
sifFile = "MouseCDHit/all.Mouse.sif"


print "running cd-hit"

sp_tools.run_cd_hit(pmidListFile, fastaPath, finalOutputPath, finalSifFilePath, percentID, sifFile)
