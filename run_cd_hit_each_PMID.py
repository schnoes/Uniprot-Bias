import sp_tools
from time import clock
import datetime
import sys
import subprocess

# Usage:
# python run_cd_hit_each_PMID.py > RunEACHPMIDCDHit.log 2>&1 


d = datetime.date.today()
fastaList = "/home/alexs/FastaFiles/eachPMIDFasta.list"
finalOutputPath = "TBCDHit/AllTogetherCDHit/CDHITcompare.TB-all." + str(d) + ".txt"
finalSifFilePath = "TBCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.TB-all." + str(d)
fastaPath = "/home/alexs/FastaFiles/TB-all.fasta"
attrOutFilePath = "TBCDHit/AllTogetherCDHit/CDHITcompare.cytoscape.TB-all." + str(d) + ".attr"
percentID = 0.97

print "running cd-hit"
sp_tools.run_all_cd_hit(fastaPath, attrOutFilePath, finalOutputPath, finalSifFilePath, percentID)


fastaListFile = open(fastaList, "r")

#print outFile
cmd = "cd-hit -i %s -o %s -c %f -d 0 -g 1"% (combinedFasta, outFile, percentID)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = process.communicate()
#outFile.write(out)
logFile = outpath + "." + id + "vs." + comp +"." + str(percentID) + ".log"
logFileHandle = open(logFile, 'w')
logFileHandle.write(out)
logFileHandle.close()
