import cPickle
import sp_tools
from time import clock
import sys
import datetime

#USAGE:
# python ./print_Top50ProtsPMIDperTaxonID.py > Top50ProtsPMIDperTaxonID.log 2>&1 

d = datetime.date.today()
finalOutputFile = "Top50ProtSpeciesWithTop50AnnotPapers." + str(d) + ".tsv"
taxonIDFile = 'Uniprot-Bias/handTax50list.txt'


#load up the pre-pickled data
print "load up the pre-pickled data"
#papers_protsExp_handle = open('Uniprot-Bias/goa_exp_papers_prots.pik', 'rb')
#papers_protsExp_dict = cPickle.load(papers_protsExp_handle)
papersTaxExp_handle = open('Uniprot-Bias/goa_taxid_exp_papers.pik', 'rb')
papersTaxExp_dict = cPickle.load(papersTaxExp_handle)
#papersExp_handle = open('Uniprot-Bias/goa_exp_papers.pik', 'rb')
#papersExp_dict = cPickle.load(papersExp_handle)
ncbi_paper_handle = open('Uniprot-Bias/ncbi_paper_info.pik')
ncbi_paper_dict = cPickle.load(ncbi_paper_handle)


#Print out PMID info per TaxonID
top = 50
pmidTaxonIDs_dict = sp_tools.count_PMID_in_taxonID(papersTaxExp_dict)
#papers_annots2_dict = sp_tools.top_papers_dict(papersExp_dict, papers_protsExp_dict, top=500)
sp_tools.print_papers_from_TaxonID_list(taxonIDFile,  pmidTaxonIDs_dict, ncbi_paper_dict, finalOutputFile, papersTaxExp_dict, True,  top)

