import cPickle
import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python ./print_CountExpEvCodesLO-GOA.py > countExpEvCodesLO-GOA.log 2>&1 

d = datetime.date.today()
outpath = "allECCodeCount-LO." + str(d) + ".tsv"


#load up the pre-pickled data
print "load up the pre-pickled data"
papersExpLO_handle = open('Uniprot-Bias/goa_exp_papers_lo.pik', 'rb') #leaves only
papersExpLO_dict = cPickle.load(papersExpLO_handle)


#count all the annotations & the annotations per Ev Code
print "Count up all the experimental evidence codes for all of GOA. Output is allECCodeCount-LO.<date>.tsv"
print clock()
sys.stdout.flush()
allECCode_dict = sp_tools.count_all_annotations_per_ec(papersExpLO_dict, outpath)

print "all done"
print clock()
