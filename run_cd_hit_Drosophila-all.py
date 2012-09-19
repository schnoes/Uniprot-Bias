import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_Drosophila-all.py > RunAllDrosophilaCDHIT.log 2>&1 

d = datetime.date.today()
finalOutputPath = "DrosophilaCDHit/AllTogetherCDHit/CDHITcompare.Drosophila-all." + str(d) + ".txt"
finalSifFilePath = "DrosophilaCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.Drosophila-all." + str(d)
fastaPath = "/home/alexs/FastaFiles/Drosophila-all.fasta"
attrOutFilePath = "DrosophilaCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.Drosophila-all." + str(d) + ".attr"
percentID = 0.97

print "running cd-hit"
sp_tools.run_all_cd_hit(fastaPath, attrOutFilePath, finalOutputPath, finalSifFilePath, percentID)
