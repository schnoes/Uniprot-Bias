import cPickle
import sp_tools

#load up the pre-pickled data
papersExp_handle = open('Uniprot-Bias/goa_exp_papers.pik', 'rb')
papersExp_dict = cPickle.load(papersExp_handle)
papers_protsExp_handle = open('Uniprot-Bias/goa_exp_papers_prots.pik', 'rb')
papers_protsExp_dict = cPickle.load(papers_protsExp_handle)
papersTaxExp_handle = open('Uniprot-Bias/goa_taxid_exp_papers.pik', 'rb')
papersTaxExp_dict = cPickle.load(papersTaxExp_handle)

# going for the top fifty papers
top = 50
papers_annots2_dict = sp_tools.top_papers_dict(papersExp_dict, papers_protsExp_dict, top=top)
all_tt_count = sp_tools.term_types_all_papers(papersExp_dict) #takes a really long time
go_ec_count = sp_tools.go_terms_with_ec_per_paper(papersExp_dict, top=top) # this takes a bit of time too
GO_EC_Count_collect_tuple = sp_tools.sort_go_ec_count(go_ec_count)
allEvCodes_dict = sp_tools.ev_codes_all_papers(papersExp_dict)
sortedProtsPerPaper_tuple = sp_tools.sort_papers_prots(papers_protsExp_dict)
sp_tools.print_paper_per_prots_sorted_go(papers_annots2_dict, all_tt_count, go_ec_count, allEvCodes_dict, 
                         sortedProtsPerPaper_tuple, "allExpPaperInfoTop50_sortedGO.txt", top=top)

