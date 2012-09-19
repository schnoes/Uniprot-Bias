import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python run_cd_hit_Mouse-all.py > RunAllMouseCDHIT.log 2>&1 

d = datetime.date.today()
finalOutputPath = "MouseCDHit/AllTogetherCDHit/CDHITcompare.Mouse-all." + str(d) + ".txt"
finalSifFilePath = "MouseCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.Mouse-all." + str(d)
fastaPath = "/home/alexs/FastaFiles/Mouse-all.fasta"
attrOutFilePath = "MouseCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.Mouse-all." + str(d) + ".attr"
percentID = 0.97

print "running cd-hit"
sp_tools.run_all_cd_hit(fastaPath, attrOutFilePath, finalOutputPath, finalSifFilePath, percentID)
