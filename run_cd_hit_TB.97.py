import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_TB.97.py > RunTBCDHIT.97.log 2>&1 & 

d = datetime.date.today()
percentID = 0.97
finalOutputPath = "TBCDHit/CDHITcompare.TB." + str(d) 
finalSifFilePath = "TBCDHit/CDHITcompare.TB.cytoscape." + str(percentID * 100) + "percent." + str(d)
fastaPath = "/home/alexs/FastaFiles/uniProtIDsPerPaper.2012-06-12.tsv"
pmidListFile = "TBCDHit/TBPMIDs.txt"
sifFile = "TBCDHit/all.TB.sif"


print "running cd-hit"

sp_tools.run_cd_hit(pmidListFile, fastaPath, finalOutputPath, finalSifFilePath, percentID, sifFile)
