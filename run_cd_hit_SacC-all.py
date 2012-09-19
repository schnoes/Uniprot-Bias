import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_SacC-all.py > RunAllSacCCDHIT.log 2>&1 

d = datetime.date.today()
finalOutputPath = "SacCCDHit/AllTogetherCDHit/CDHITcompare.SacC-all." + str(d) + ".txt"
finalSifFilePath = "SacCCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.SacC-all." + str(d)
fastaPath = "/home/alexs/FastaFiles/SacC-all.fasta"
attrOutFilePath = "SacCCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.SacC-all." + str(d) + ".attr"
percentID = 0.97

print "running cd-hit"
sp_tools.run_all_cd_hit(fastaPath, attrOutFilePath, finalOutputPath, finalSifFilePath, percentID)
