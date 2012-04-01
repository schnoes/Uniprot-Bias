import cPickle
import sp_tools
from time import clock
import sys

#load up the pre-pickled data
print "load up the pre-pickled data"
papersTaxExp_handle = open('Uniprot-Bias/goa_taxid_exp_papers.pik', 'rb')
papersTaxExp_dict = cPickle.load(papersTaxExp_handle)

print "count all the taxonID annotations & proteins"
print clock()
sys.stdout.flush()
sp_tools.count_all_annotations_and_proteins_taxonIDs(papersTaxExp_dict, "countTaxonIDAnnots&Prots.tsv")
print clock()
