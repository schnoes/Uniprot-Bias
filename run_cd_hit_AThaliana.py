import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_AThaliana.97.py > RunAThalianaCDHIT.97.log 2>&1 & 

d = datetime.date.today()
percentID = 0.97
finalOutputPath = "AThalianaCDHit/CDHITcompare.AThaliana." + str(d) 
finalSifFilePath = "AThalianaCDHit/CDHITcompare.AThaliana.cytoscape." + str(percentID * 100) + "percent." + str(d)
fastaPath = "/home/alexs/FastaFiles/uniProtIDsPerPaper.2012-06-12.tsv"
pmidListFile = "AThalianaCDHit/AThalianaPMIDs.txt"
sifFile = "HumanCDHit/allHumanPaper.sif"


print "running cd-hit"

sp_tools.run_cd_hit(pmidListFile, fastaPath, finalOutputPath, finalSifFilePath, percentID, sifFile)
