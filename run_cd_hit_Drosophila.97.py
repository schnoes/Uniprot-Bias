import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_Drosophila.97.py > RunDrosophilaCDHIT.97.log 2>&1 & 

d = datetime.date.today()
percentID = 0.97
finalOutputPath = "DrosophilaCDHit/CDHITcompare.Drosophila." + str(d) 
finalSifFilePath = "DrosophilaCDHit/CDHITcompare.Drosophila.cytoscape." + str(percentID * 100) + "percent." + str(d)
fastaPath = "/home/alexs/FastaFiles/uniProtIDsPerPaper.2012-06-12.tsv"
pmidListFile = "DrosophilaCDHit/DrosophilaPMIDs.txt"
sifFile = "DrosophilaCDHit/all.Drosophila.sif"


print "running cd-hit"

sp_tools.run_cd_hit(pmidListFile, fastaPath, finalOutputPath, finalSifFilePath, percentID, sifFile)
