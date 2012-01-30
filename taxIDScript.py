import cPickle
import sp_tools

log_handle = open('logfile', 'w')
papersExp_handle = open('goa-pickles/goa_exp_papers.pik', 'rb')
papersExp_dict = cPickle.load(papersExp_handle)
papers_protsExp_handle = open('goa-pickles/goa_exp_papers_prots.pik', 'rb')
papers_protsExp_dict = cPickle.load(papers_protsExp_handle)

log_handle.write("Make the TaxID files")
#Make the taxon ID pickle info
goa_path = 'gene_association.goa_uniprot'
papersTaxID_dict, papers_prots_dict = sp_tools.go_tax_in_papers_goa(goa_path)
papersTaxIDExp_dict, papers_protsExp_dict = sp_tools.exp_in_papers(papersTaxID_dict, papers_prots_dict)
papersTaxIDExp_handle = open('/home/alexs/goa_taxid_exp_papers.pik', 'wb')
cPickle.dump(papersTaxIDExp_dict, papersTaxIDExp_handle)

log_handle.write("print out the TaxID files")
#Print the taxon ID info
top = 50
papers_annots2_dict = sp_tools.top_papers_dict(papersExp_dict, papers_protsExp_dict, top=top)
sortedProtsPerPaper_tuple = sp_tools.sort_papers_prots(papers_protsExp_dict)
papersTaxonIDs_dict = sp_tools.count_taxonIDs(papersTaxIDExp_dict)
sp_tools.print_papers_taxonIDs(papersTaxonIDs_dict, papers_annots2_dict, sortedProtsPerPaper_tuple, "allExpTaxonIDsTop50.txt", top=top)
log_handle.write("all done")


log_handle.close()
papersExp_handle.close()
papers_protsExp_handle.close()
papersTaxIDExp_handle.close()