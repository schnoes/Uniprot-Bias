import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_Celegans-all.py > RunAllCelegansCDHIT.log 2>&1 

d = datetime.date.today()
finalOutputPath = "CelegansCDHit/AllTogetherCDHit/CDHITcompare.Celegans-all." + str(d) + ".txt"
finalSifFilePath = "CelegansCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.Celegans-all." + str(d)
fastaPath = "/home/alexs/FastaFiles/Celegans-all.fasta"
attrOutFilePath = "CelegansCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.Celegans-all." + str(d) + ".attr"
percentID = 0.97

print "running cd-hit"
sp_tools.run_all_cd_hit(fastaPath, attrOutFilePath, finalOutputPath, finalSifFilePath, percentID)
