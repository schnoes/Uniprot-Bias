import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_TB-all.py > RunAllTBCDHIT.log 2>&1 

d = datetime.date.today()
finalOutputPath = "TBCDHit/AllTogetherCDHit/CDHITcompare.TB-all." + str(d) + ".txt"
finalSifFilePath = "TBCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.TB-all." + str(d)
fastaPath = "/home/alexs/FastaFiles/TB-all.fasta"
attrOutFilePath = "TBCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.TB-all." + str(d) + ".attr"
percentID = 0.97

print "running cd-hit"
sp_tools.run_all_cd_hit(fastaPath, attrOutFilePath, finalOutputPath, finalSifFilePath, percentID)
