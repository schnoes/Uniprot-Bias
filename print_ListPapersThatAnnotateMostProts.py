import cPickle
import sp_tools
from time import clock
import datetime
import sys

# Usage:
# python ./print_ListPapersThatAnnotateMostProts.py > ListPapersThatAnnotateMostProts.log 2>&1 

d = datetime.date.today()
finalOutputFile = "allExpPaperInfoTop50." + str(d) + ".tsv"


#load up the pre-pickled data
print "load up the pre-pickled data"
papersExp_handle = open('Uniprot-Bias/goa_exp_papers.pik', 'rb')
papersExp_dict = cPickle.load(papersExp_handle)
papers_protsExp_handle = open('Uniprot-Bias/goa_exp_papers_prots.pik', 'rb')
papers_protsExp_dict = cPickle.load(papers_protsExp_handle)
all_tt_count_handle = open("Uniprot-Bias/all_tt_count.pik", 'rb')
all_tt_count = cPickle.load(all_tt_count_handle)

# going for the list of papers That annotate most proteins. Top designates how far down the list we go
top = 50
print "top_papers_dict: Get all the PMID info for the top papers"
print clock()
sys.stdout.flush()
papers_annots2_dict = sp_tools.top_papers_dict(papersExp_dict, papers_protsExp_dict, top=top)

print "term_types_all_papers: Count up all the terms types for each paper"
print clock()
sys.stdout.flush()
all_tt_count = sp_tools.term_types_all_papers(papersExp_dict) #takes a really long time

print "cPikleDump: save that all_tt_count for later"
print clock()
sys.stdout.flush()
sp_tools.cPickleDump(all_tt_count, "Uniprot-Bias/all_tt_count.pik")

print "go_terms_with_ec_per_paper: Create a dict that counts up how many times a specific (GO ID, GO Term Text, EvCode) tuple occurs for each paper"
print clock()
sys.stdout.flush()
go_ec_count = sp_tools.go_terms_with_ec_per_paper(papersExp_dict, top=top) # this takes a bit of time too

print "ev_codes_all_papers: Calculate the number of times a paper gives a certain experimental evidence code."
print clock()
sys.stdout.flush()
allEvCodes_dict = sp_tools.ev_codes_all_papers(papersExp_dict)

print "sort_papers_prots: Sort the dictionary papers_prots according to the number of proteins annotated by a particular paper (PMID)."
print clock()
sys.stdout.flush() 
sortedProtsPerPaper_tuple = sp_tools.sort_papers_prots(papers_protsExp_dict)

print "print_paper_per_prots: print out the results of the the top papers per proteins. Final Ouptfile: allExpPaperInfoTop50.<date>.tsv"
print clock()
sys.stdout.flush()
sp_tools.print_paper_per_prots_go(papers_annots2_dict, all_tt_count, go_ec_count, allEvCodes_dict, 
                         sortedProtsPerPaper_tuple, finalOutputFile, top=top)
print "all done"
print clock()

