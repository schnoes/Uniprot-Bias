import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_SacC.97.py > RunSacCCDHIT.97.log 2>&1 & 

d = datetime.date.today()
percentID = 0.97
finalOutputPath = "SacCCDHit/CDHITcompare.SacC." + str(d) 
finalSifFilePath = "SacCCDHit/CDHITcompare.SacC.cytoscape." + str(percentID * 100) + "percent." + str(d)
fastaPath = "/home/alexs/FastaFiles/uniProtIDsPerPaper.2012-06-12.tsv"
pmidListFile = "SacCCDHit/SacCPMIDs.txt"
sifFile = "SacCCDHit/all.SacC.sif"


print "running cd-hit"

sp_tools.run_cd_hit(pmidListFile, fastaPath, finalOutputPath, finalSifFilePath, percentID, sifFile)
