import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_AThaliana-all.py > RunAllAThalianaCDHIT.log 2>&1 

d = datetime.date.today()
finalOutputPath = "AThalianaCDHit/AllTogetherCDHit/CDHITcompare.AThaliana-all." + str(d) + ".txt"
finalSifFilePath = "AThalianaCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.AThaliana-all." + str(d)
fastaPath = "/home/alexs/FastaFiles/AThaliana-all.fasta"
attrOutFilePath = "AThalianaCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.AThaliana-all." + str(d) + ".attr"
percentID = 0.97

print "running cd-hit"
sp_tools.run_all_cd_hit(fastaPath, attrOutFilePath, finalOutputPath, finalSifFilePath, percentID)
