import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_Pombe-all.py > RunAllPombeCDHIT.log 2>&1 

d = datetime.date.today()
finalOutputPath = "PombeCDHit/AllTogetherCDHit/CDHITcompare.Pombe-all." + str(d) + ".txt"
finalSifFilePath = "PombeCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.Pombe-all." + str(d)
fastaPath = "/home/alexs/FastaFiles/Pombe-all.fasta"
attrOutFilePath = "PombeCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.Pombe-all." + str(d) + ".attr"
percentID = 0.97

print "running cd-hit"
sp_tools.run_all_cd_hit(fastaPath, attrOutFilePath, finalOutputPath, finalSifFilePath, percentID)
