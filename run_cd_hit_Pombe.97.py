import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_Pombe.97.py > RunPombeCDHIT.97.log 2>&1 & 

d = datetime.date.today()
percentID = 0.97
finalOutputPath = "PombeCDHit/CDHITcompare.Pombe." + str(d) 
finalSifFilePath = "PombeCDHit/CDHITcompare.Pombe.cytoscape." + str(percentID * 100) + "percent." + str(d)
fastaPath = "/home/alexs/FastaFiles/uniProtIDsPerPaper.2012-06-12.tsv"
pmidListFile = "PombeCDHit/PombePMIDs.txt"
sifFile = "PombeCDHit/all.Pombe.sif"


print "running cd-hit"

sp_tools.run_cd_hit(pmidListFile, fastaPath, finalOutputPath, finalSifFilePath, percentID, sifFile)
