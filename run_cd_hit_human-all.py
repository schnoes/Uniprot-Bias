import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_human-all.py > RunAllCDHIT.log 2>&1 

d = datetime.date.today()
finalOutputPath = "HumanCDHit/AllTogetherCDHit/CDHITcompare.human-all." + str(d) + ".txt"
finalSifFilePath = "HumanCDHit/AllTogetherCDHit/CDHITcompare.cytoscape." + str(d)
fastaPath = "/home/alexs/FastaFiles/human-all.fasta"
attrOutFilePath = "HumanCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.all-human." + str(d) + ".attr"
percentID = 0.97

print "running cd-hit"
sp_tools.run_all_cd_hit(fastaPath, attrOutFilePath, finalOutputPath, finalSifFilePath, percentID)
