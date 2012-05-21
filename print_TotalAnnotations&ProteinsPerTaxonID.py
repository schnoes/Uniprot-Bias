import cPickle
import sp_tools
from time import clock
import sys
import datetime

#USAGE:
# python ./print_TotalAnnotations\&ProteinsPerTaxonID.py > TotalAnnotationsProteinsPerTaxonID.log 2>&1 

d = datetime.date.today()
finalOutputFile = "countTaxonIDAnnots&Prots." + str(d) + ".tsv"


#load up the pre-pickled data
print "load up the pre-pickled data"
papersTaxExp_handle = open('Uniprot-Bias/goa_taxid_exp_papers.pik', 'rb')
papersTaxExp_dict = cPickle.load(papersTaxExp_handle)

print "count all the taxonID annotations & proteins. Output is: countTaxonIDAnnots&Prots.<date>.tsv "
print clock()
sys.stdout.flush()
sp_tools.count_all_annotations_and_proteins_taxonIDs(papersTaxExp_dict, finalOutputFile)
print clock()
