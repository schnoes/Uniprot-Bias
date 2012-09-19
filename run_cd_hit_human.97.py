import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_human.97.py > RunCDHIT.97.log 2>&1 

d = datetime.date.today()
percentID = 0.97
finalOutputPath = "HumanCDHit/CDHITcompare.Human." +  str(d) 
finalSifFilePath = "HumanCDHit/CDHITcompare.cytoscape.Human." + str(percentID * 100) + "percent." + str(d)
fastaPath = "/home/alexs/FastaFiles/uniProtIDsPerPaper.2012-06-12.tsv"
pmidListFile = "HumanCDHit/humanPMIDs.txt"
sifFile = "HumanCDHit/allHumanPaper.sif"

print "running cd-hit"

sp_tools.run_cd_hit(pmidListFile, fastaPath, finalOutputPath, finalSifFilePath, percentID, sifFile)
