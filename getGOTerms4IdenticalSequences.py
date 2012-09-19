import cPickle
import sp_tools
from time import clock
import sys
import datetime

d = datetime.date.today()
inputFile = "/Users/schnoes/Alex/AFP_CAFA_2011/References Work/AThaliana/CDHITcompare.AThaliana.cytoscape.100.0percent.2012-07-04.all.sif"
outFile = "/Users/schnoes/Alex/AFP_CAFA_2011/References Work/AThaliana/AThaliana.100percentID.GOtermsAndText.txt"

#load up the pre-pickled data
print "load up the pre-pickled data"
papersTaxExp_handle = open('Uniprot-Bias/goa_taxid_exp_papers.pik', 'rb')
papersTaxExp_dict = cPickle.load(papersTaxExp_handle)

#load up the file that lists the identical sequences. Will be a sif file. Identity is identified by 'ss'
print "parse out the identitical relationships"
relationship_dict = sp_tools.parse_sif_file_for_ident_seqs(inputFile)


# print out the identity relationships with GO terms and term text
print "print it out"
sp_tools.print_go_terms_and_text_for_sim_seqs(relationship_dict, papersTaxExp_dict, outFile)

