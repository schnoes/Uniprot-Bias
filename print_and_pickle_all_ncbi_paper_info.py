import cPickle
import sp_tools
from time import clock
import sys
import datetime


#usage: python ./print_and_pickle_all_ncbi_paper_info.py > AllNCBIPaperInfo.log 2>&1 

#load up the pre-pickled data
papersExp_handle = open('Uniprot-Bias/goa_exp_papers.pik', 'rb')
papersExp_dict = cPickle.load(papersExp_handle)
papers_protsExp_handle = open('Uniprot-Bias/goa_exp_papers_prots.pik', 'rb')
papers_protsExp_dict = cPickle.load(papers_protsExp_handle)


print "all_NCBI_paper_info_dict: get all the pubmed info we want"
print clock()
sys.stdout.flush()
ncbi_paper_dict = sp_tools.all_NCBI_paper_info_dict(papersExp_dict, papers_protsExp_dict, outpath=None,delim="\t", top=None)

print "cPikleDump: save that ncbi_paper_dict for later"
print clock()
sys.stdout.flush()
sp_tools.cPickleDump(ncbi_paper_dict, "Uniprot-Bias/ncbi_paper_info.pik")

print "Done!"
print clock()
sys.stdout.flush()
