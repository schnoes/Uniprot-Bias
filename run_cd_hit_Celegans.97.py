import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_Celegans.97.py > RunCelegansCDHIT.97.log 2>&1 & 

d = datetime.date.today()
percentID = 0.97
finalOutputPath = "CelegansCDHit/CDHITcompare.Celegans." + str(d) 
finalSifFilePath = "CelegansCDHit/CDHITcompare.Celegans.cytoscape." + str(percentID * 100) + "percent." + str(d)
fastaPath = "/home/alexs/FastaFiles/uniProtIDsPerPaper.2012-06-12.tsv"
pmidListFile = "CelegansCDHit/CelegansPMIDs.txt"
sifFile = "CelegansCDHit/all.Celegans.sif"


print "running cd-hit"

sp_tools.run_cd_hit(pmidListFile, fastaPath, finalOutputPath, finalSifFilePath, percentID, sifFile)
