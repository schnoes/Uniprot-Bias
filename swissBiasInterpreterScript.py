import cPickle
import sp_tools

#load up the pre-pickled data
papersExp_handle = open('Uniprot-Bias/goa_exp_papers.pik', 'rb')
papersExp_dict = cPickle.load(papersExp_handle)
papers_protsExp_handle = open('Uniprot-Bias/goa_exp_papers_prots.pik', 'rb')
papers_protsExp_dict = cPickle.load(papers_protsExp_handle)
papersTaxExp_handle = open('Uniprot-Bias/goa_taxid_exp_papers.pik', 'rb')
papersTaxExp_dict = cPickle.load(papersTaxExp_handle)



#Print out PMID info per TaxonID
top = 50
pmidTaxonIDs_dict = sp_tools.count_PMID_in_taxonID(papersTaxExp_dict)
papers_annots2_dict = sp_tools.top_papers_dict(papersExp_dict, papers_protsExp_dict, top=500)
taxonIDFile = 'Uniprot-Bias/handTax50list.txt'
sp_tools.print_papers_from_TaxonID_list(taxonIDFile,  pmidTaxonIDs_dict, papers_annots2_dict, "AllPMIDvTaxonIDsTop50.txt", papersExp_dict, True,  top)





#count all the annotations & the annotations per Ev Code
allECCode_dict = sp_tools.count_all_annotations_per_ec(papersExp_dict)

# count all the term type annotations
sum_tt_count = sp_tools.count_all_term_types(papersExp_dict)

# count all the species annotations
all_taxonID_dict = sp_tools.count_all_annotations_taxonIDs(papersTaxExp_dict)

#top go codes
top = 50
ec_go_code_count = sp_tools.count_top_go_terms_per_ecode_all_entries(papersExp_dict, "sortedECGO.txt", top)

# going for the top fifty papers
top = 50
papers_annots2_dict = sp_tools.top_papers_dict(papersExp_dict, papers_protsExp_dict, top=top)
all_tt_count = sp_tools.term_types_all_papers(papersExp_dict) #takes a really long time
go_ec_count = sp_tools.go_terms_with_ec_per_paper(papersExp_dict, top=top) # this takes a bit of time too
allEvCodes_dict = sp_tools.ev_codes_all_papers(papersExp_dict)
sortedProtsPerPaper_tuple = sp_tools.sort_papers_prots(papers_protsExp_dict)
sp_tools.print_paper_per_prots_go(papers_annots2_dict, all_tt_count, go_ec_count, allEvCodes_dict, 
                         sortedProtsPerPaper_tuple, "allExpPaperInfoTop50.txt", top=top)


# going for the top fifty papers
top = 50
papers_annots2_dict = sp_tools.top_papers_dict(papersExp_dict, papers_protsExp_dict, top=top)
all_tt_count = sp_tools.term_types_all_papers(papersExp_dict) #takes a really long time
go_ec_count = sp_tools.go_terms_with_ec_per_paper(papersExp_dict, top=top) # this takes a bit of time too
GO_EC_Count_collect_tuple = sp_tools.sort_go_ec_count(go_ec_count)
allEvCodes_dict = sp_tools.ev_codes_all_papers(papersExp_dict)
sortedProtsPerPaper_tuple = sp_tools.sort_papers_prots(papers_protsExp_dict)
sp_tools.print_paper_per_prots_go(papers_annots2_dict, all_tt_count, go_ec_count, allEvCodes_dict, 
                         sortedProtsPerPaper_tuple, "allExpPaperInfoTop50.txt", top=top)



#Make the taxon ID pickle info
goa_path = 'gene_association.goa_uniprot'
papersTaxID_dict, papers_prots_dict = sp_tools.go_tax_in_papers_goa(goa_path)
papersTaxIDExp_dict, papers_protsExp_dict = sp_tools.exp_in_papers(papersTaxID_dict, papers_prots_dict)
papersTaxIDExp_handle = open('goa_taxid_exp_papers.pik', 'wb')
cPickle.dump(papersTaxIDExp_dict, papersTaxIDExp_handle)


#Print the taxon ID info
top = 50
papers_annots2_dict = sp_tools.top_papers_dict(papersExp_dict, papers_protsExp_dict, top=top)
sortedProtsPerPaper_tuple = sp_tools.sort_papers_prots(papers_protsExp_dict)
papersTaxonIDs_dict = sp_tools.count_taxonIDs(papersTaxIDExp_dict)
sp_tools.print_papers_taxonIDs(papersTaxonIDs_dict, papers_annots2_dict, sortedProtsPerPaper_tuple, "allExpTaxonIDsTop50.txt", top=top)




"""
# Print out the dicts. Can only do this on the command line, not in the interpreter. It gets cranky otherwise.
sp_tools.printDict(papersExp_dict, "papersExp_dict.txt")
sp_tools.printDict(papers_protsExp_dict, "papers_protsExp_dict.txt")

# Want to know how many proteins are annotated by each paper.
# Only print out this number if it is above 100.
count = 0
for key, values in papers_protsExp_dict.iteritems():
    count = int(count) + 1
    if len(values) > 100:
        print len(values)

print "Papers Count: " + str(count)

# Just trying to get a handle on what this dict is made of (its a dict of dicts)
for key, values in papers_protsExp_dict.iteritems():
    print values


# This onlys serves to confirm that both Swiss Prot and Trembl sequences are annotated with GO
# Experimental evidence codes. The numbers calculated below are useless though.
spCount = 0
trCount = 0
for pmid, uniprotID_dicts in papers_protsExp_dict.iteritems():
    for uniprotID, annotationCounts in uniprotID_dicts.iteritems():
        if ':P' in uniprotID:
            spCount = int(spCount) + 1
        else:
            trCount = int(trCount) + 1

print "SwissProt Dup Count: " + str(spCount)
print "Trembl Dup Count: " + str(trCount)
"""

sortedLen = sorted(papers_protsExp_dict.iteritems(), key=lambda x: len(x[1]))
sortedLen2 = sorted([(key, len(value)) for (key, len(value)) in papers_protsExp_dict.items()])
sortedLen2 = sorted([(key, len(value)) for (key, value) in papers_protsExp_dict.items()])
sortedLen2 = sorted([(len(value), key) for (key, value) in papers_protsExp_dict.items()])
import collections
ProtsPerPaper_collect = collections.namedtuple('ProtsPerPaper_collect', 'numProts PMID')
sortedProtsPerPaper_tuple = sorted([ProtsPerPaper_collect(len(value), key) for (key, value) in papers_protsExp_dict.items()])

# evCode_dict = sp_tools.ec_stats(papersExp_dict)
# sp_tools.printDict(evCode_dict, "evCode_dict.txt")

top50Papers_list = sp_tools.top_papers(papersExp_dict, "topPapers50.txt", "\t", 50)
top50Papers_dict = sp_tools.top_papers_dict(papersExp_dict, "topPapers50.txt", "\t", 50)

top50Ontology_list = sp_tools.top_ontology(papersExp_dict, "topOntology50.txt", 50)
top50GOTerms_list = sp_tools.top_go_terms(papersExp_dict, "topGoTerms50.txt", 50)
top50GoTermsPerPaper_dict = sp_tools.go_terms_per_paper(papersExp_dict, "top50GoTermsPerPaper.txt", 50)
top50GoTermsWithEvCodePerPaper_dict = sp_tools.go_terms_per_paper(papersExp_dict, "top50GoTermsWithEvCodePerPaper.txt", 50)
sp_tools.printDict(top50GoTermsWithEvCodePerPaper_dict, "top50GoTermsWithEvCodePerPaper.txt")


protStats_dict = sp_tools.prot_stats(papersExp_dict)
sp_tools.printDict(protStats_dict, "protStats_dict.txt")

all_tt_count = sp_tools.term_types_all_papers(papersExp_dict) #takes a really long time

#go_ids, exp_go_ids = sp_tools.go_freq_in_papers(data1)



