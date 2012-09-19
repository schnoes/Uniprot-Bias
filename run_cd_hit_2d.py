import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python ./run_cd_hit_2d.py > RunCDHIT2D.log 2>&1 

d = datetime.date.today()
finalOutputPath = "CDHIT2DCompare." + str(d) 
fastaPath = "/home/alexs/FastaFiles/uniProtIDsPerPaper.2012-06-07.tsv"
pmidListFile = "HumanCDHit/humanPMIDs.txt"

print "running cd-hit-2d"

sp_tools.run_cd_hit_2d(pmidListFile, fastaPath, finalOutputPath)
