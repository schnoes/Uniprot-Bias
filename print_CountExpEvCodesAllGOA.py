import cPickle
import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python ./print_CountExpEvCodesAllGOA.py > countExpEvCodesAllGOA.log 2>&1 

d = datetime.date.today()
outpath = "allECCodeCount." + str(d) + ".tsv"


#load up the pre-pickled data
print "load up the pre-pickled data"
papersExp_handle = open('Uniprot-Bias/goa_exp_papers.pik', 'rb')
papersExp_dict = cPickle.load(papersExp_handle)


#count all the annotations & the annotations per Ev Code
print "Count up all the experimental evidence codes for all of GOA. Output is allECCodeCount.<date>.tsv"
print clock()
sys.stdout.flush()
allECCode_dict = sp_tools.count_all_annotations_per_ec(papersExp_dict, outpath)

print "all done"
print clock()
