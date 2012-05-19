
from Bio import SwissProt as SP
from Bio import Entrez, Medline
#import matplotlib.pyplot as pyplot
from GO import go_utils as gu
import MySQLdb
import collections
import getpass
import os.path
import os
import cPickle
import datetime

# edits: AMS Jan 4, 2012

#edit IF 2/17/2012
# Set of GO experimental evidence codes
EEC = set(['EXP','IDA','IPI','IMP','IGI','IEP'])

#Email address used for NCBI efetch tools
MY_EMAIL = 'schnoes@gmail.com'

#cPickle
#When using the pickeled data...
#    papers corresponds to any pickle that has 'papers' in the name
#    papers_prots corresponds to any pickle that has 'papers_prots' in the name
#    must load pickles first to use this code. e.g.
#        test1 = open('goa-pickles/goa_exp_papers.pik', 'rb')
#        data1 = cPickle.load(test1)

###########################################################
# mysqlConnect
###########################################################
def mysqlConnect():
    """This is where you set all the mysql login, db info"""
    configFilePath = os.path.join(os.path.expanduser("~" + getpass.getuser()), '.my.cnf')
    return MySQLdb.connect(db="GeneOntology", read_default_file=configFilePath, host="mysql-dev.cgl.ucsf.edu", 
                           port=13308)
     


###########################################################
# cPickleDump
###########################################################
def cPickleDump(toDump_dict, pickleFileName):
    """Dump the input dict to the input filename"""
    dump_handle = open(pickleFileName, 'wb')
    cPickle.dump(toDump_dict, dump_handle)
    
    
    
def exp_in_papers(papers,papers_prots):
    # Pulls out Annotations with experimental evidence codes only
    # exp_papers_prots: key: pubmed_id; value: protein count dictionary
    # A protein count dictionary: key: UniProtID; value: count of that protein.
    # The protein count is the number of different experimentally-evidenced GO codes this
    # protein receives for a given paper "p".
    #
    # exp_papers: key: pubmed_id; value: a list of "go_rec" records.
    # go_rec record is a dictionary. Keys (values): 'sp_id' (swissprot id); 
    # 'go_id': (GO ID); 'go_ec': (GO Evidence Code).
    # 
    # Can be used with SP & GOA data
    # Note: This creates a non-redundant set of proteins. It does not dubble count proteins!
    
    exp_papers = {}
    exp_papers_prots = {}
    for p in papers:
        for go_rec in papers[p]:
            if go_rec['go_ec'] in EEC:
                exp_papers.setdefault(p,[]).append(go_rec)
                if p not in exp_papers_prots:
                    exp_papers_prots[p] = {go_rec['sp_id']: 1}
                else:
                    exp_papers_prots[p][go_rec['sp_id']] = \
                    exp_papers_prots[p].get(go_rec['sp_id'],0)+1
                
    return exp_papers, exp_papers_prots

def not_exp_in_papers(papers,papers_prots):
    # Annotations with no experimental evidence codes
    # nexp_papers_prots: key: pubmed_id; value: protein count dictionary
    # A protein count dictionary: key: SwissProtID; value: count of that protein.
    # The protein count is the number of different experimentally-evidenced GO codes this
    # protein receives for a given paper "p".
    #
    # nexp_papers: key: pubmed_id; value: a list of "go_rec" records.
    # go_rec record is a dictionary. Keys (values): 'sp_id' (swissprot id); 
    # 'go_id': (GO ID); 'go_ec': (GO Evidence Code).
    #
    # Can be used with SP & GOA data
    
    nexp_papers = {}
    nexp_papers_prots = {}
    for p in papers:
        for go_rec in papers[p]:
            if go_rec['go_ec'] not in EEC:
                nexp_papers.setdefault(p,[]).append(go_rec)
                if p not in nexp_papers_prots:
                    nexp_papers_prots[p] = {go_rec['sp_id']: 1}
                else:
                    nexp_papers_prots[p][go_rec['sp_id']] = \
                    nexp_papers_prots[p].get(go_rec['sp_id'],0)+1
                
    return nexp_papers, nexp_papers_prots

def plot_papers_prots(papers_prots, exp_papers_prots, do_plot=True):
    # hist: key: number of proteins; value: number of papers describing that number of proteins
    #
    # Can be used with SP & GOA data
    
    hist  = {}
    exp_hist  = {}
    for p in papers_prots:
        n_prots = len(papers_prots[p].keys())
        hist[n_prots] = hist.get(n_prots,0) + 1
    for p in exp_papers_prots:
        n_prots = len(exp_papers_prots[p].keys())
        exp_hist[n_prots] = hist.get(n_prots,0) + 1
    # Make a list in  the format: [(n_papers, n_proteins), ...]
    hist_table = [(hist[i],i) for i in hist]
    hist_table.sort()
    hist_table.reverse()
    exp_hist_table = [(hist[i],i) for i in hist]
    exp_hist_table.sort()
    exp_hist_table.reverse()
    if do_plot:
        pyplot.figure()
        pyplot.loglog()
        pyplot.xlabel('papers')
        pyplot.ylabel('proteins')
        pyplot.plot([x[0] for x in hist_table],[y[1] for y in hist_table],'ob')
        pyplot.plot([x[0] for x in exp_hist_table],[y[1] for y in exp_hist_table],'xr')
        pyplot.show()
        
    return None
 
def ec_stats(papers):
    """Creat a dict that holds the types and number of uses of each evidence code associate 
    with that particular paper."""
    # FN Can be used with SP & GOA data
    
    ec_count = {}
    #ec_count[PMID] = { {Ev Code : # of times it was used} } 
    for p in papers:
        #Papers[PMID] = [a list of dicts: each dict containing 3 entries: swissProt ID entry (key='sp_id'), 
        # go_id entry (key = 'go_id'), GO evidence code (key = 'go_ec')]    
        ec_count[p] = {}
        for annot in papers[p]:
            ec = annot['go_ec']
            ec_count[p][ec] = ec_count[p].get(ec,0) + 1
    return ec_count

def prot_stats(papers):
    #
    # Can be used with SP & GOA data
    
    prot_count = {}
    for p in papers:
        prot_count[p] = {}
        for annot in papers[p]:
            prot = annot['sp_id']
            prot_count[p][prot] = prot_count[p].get(prot,0) + 1
    return prot_count

def top_go_terms(papers,outpath=None,top=20):
    """Determines the top GO terms annotated in the analysis set and 1) puts it in 
    the output list top_go and 2) writes it out to a tab delim file 'outpath'
    
    Note: this function is currently identical to top_ontology()"""    #
    # Can be used with SP & GOA data
    
    go_count = {}
    for p in papers:
        for rec in papers[p]:
            go_id = rec['go_id']
            go_count[go_id] = go_count.get(go_id,0) + 1
    # top_go = [(# of uses, GO ID)]
    top_go = [(i[1],i[0]) for i in go_count.items()]
    top_go.sort()
    if outpath:
        go_con = mysqlConnect()
        go_cur = go_con.cursor()
        f = open(outpath,"w")
        for i in top_go[-top:]:
            name = gu.go_acc_to_name(i[1],go_cur)
            f.write("%d\t%s\t%s\n" % (i[0], i[1], name))
        go_hist = {}
        for i in top_go:
            go_hist[i[0]] = go_hist.get(i[0],0) + 1
        go_hist_list = [(h[1],h[0]) for h in go_hist.items()] 
        go_hist_list.sort()
        fhist = open("hist_%s" % outpath, "w")
        for h in go_hist_list:
#            print h
            fhist.write("%d\t%d\n" % h)
        
        f.close()
        fhist.close()
        go_con.close()
    return top_go
        
def go_terms_per_paper(papers,outpath=None,top=20):
    #
    # Can be used with SP & GOA data
    
    go_count = {}
    
    go_con = mysqlConnect()
    go_cur = go_con.cursor()
    
    for p in papers:
        for rec in papers[p]:
            go_id = rec['go_id']
            try:
                name = gu.go_acc_to_name(go_id,go_cur)
            except IndexError: #sometimes the GO ID given is actually a synonym
                try:
                    name = gu.go_acc_to_synonym_name(go_id, go_cur)
                except IndexError: #sometimes it just doesn't work
                    print "problem with GO ID", go_id
                    name = ''
            gokey = (go_id, name)
            # go_count[PMID] = {{(GO ID, GO Term Text) : # times paper gives this annotaion}}
            if p in go_count:
                go_count[p][gokey] = go_count[p].get(gokey,0) + 1
            else:
                go_count[p] = {gokey: 1}
    go_con.close()
    return go_count

###########################################################
# go_terms_with_ec_per_paper
###########################################################
def go_terms_with_ec_per_paper(papers,outpath=None,top=20):
    """Create a dict that counts up how many times a specific (GO ID, GO Term Text, EvCode) 
    tuple occurs for each paper"""
    #
    # Can be used with SP & GOA data
    
    go_ec_count = {}
    
    go_con = mysqlConnect()
    go_cur = go_con.cursor()
    for p in papers:
        for rec in papers[p]:
            go_id = rec['go_id']
            go_ec = rec['go_ec']
            try:
                name = gu.go_acc_to_name(go_id,go_cur)
            except IndexError: #sometimes the GO ID given is actually a synonym
                try:
                    name = gu.go_acc_to_synonym_name(go_id, go_cur)
                except IndexError: #sometimes it just doesn't work
                    print "problem with GO ID", go_id
                    name = ''
            gokey = (go_id, name, go_ec)
            # go_ec_count[PMID] = {{(GO ID, GO Term Text, Ev Code) : # times paper gives this annotaion}}
            if p in go_ec_count:
                go_ec_count[p][gokey] = go_ec_count[p].get(gokey,0) + 1
            else:
                go_ec_count[p] = {gokey: 1}
    go_con.close()
    return go_ec_count


###########################################################
# print_paper_per_prots_go
###########################################################
def print_paper_per_prots_go(papers_annots2_dict, all_tt_count, go_ec_count,
                             allEvCodes_dict, sortedProtsPerPaper_tuple, outpath, top=20):
    """Prints out all information that we have for each paper(PMID).
    papers_annots2_dict: dict of the top X papers, with title, year and journal name
    all_tt_count: this is a dict that gives how many term types each paper annotates
    go_ec_count: this is a dict that gives how many times a paper gives a specific (go ID, go Name, Ev code) annotation
    allEvCodes_dict: this is a dict that gives how many times a paper supports a given experimental ev Code (EEC global).
    sortedProtsPerPaper_tuple: this is a sorted named tuple (largest first) of all the papers, sorted by the number of proteins the paper annotates
    outpath: where the data will be printed out.
    top: the number of papers we want to print out. Sorted by the number of proteins the paper annotates.
    """
    outFile = open(outpath, 'w')
    #print out header line
    outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                  ('Num Prots', 'Num Annots', 'PMID', 'Title', 'Year', 'Journal', 'MFO Annot', 'BPO Annot', 'CCO Annot', 
                   'Num EXP', 'Num IDA', 'Num IEP', 'Num IGI', 'Num IMP', 'Num IPI',   'GO ID', 'GO Name', 'Ev Code', 'Num Used'))
    out_list = []
    type_list = ['molecular_function', 'biological_process', 'cellular_component']
    for sortedPaper in sortedProtsPerPaper_tuple[0:top]:
        # add the number of proteins & the PMID
        out_list.append(str(sortedPaper.numProts)) #[0]
        out_list.append(str(sortedPaper.PMID)) # will be [2]
        # get paper info
        out_list.extend(papers_annots2_dict[sortedPaper.PMID][1:]) #add title, year & journal [3, 4, 5]
        out_list.insert(1, str(papers_annots2_dict[sortedPaper.PMID][0])) # add number of annotations, [1]
        # get GO type info, [6, 7, 8]
        for t in type_list:
            try:
                out_list.append(str(all_tt_count[sortedPaper.PMID][t]))
            except:
                out_list.append('0')
        #get ev code data [9, 10, 11, 12, 13, 14]
        for ec in sorted(EEC):
            try:
                out_list.append(str(allEvCodes_dict[(sortedPaper.PMID, ec)]))
            except:
                out_list.append('0')
        # get GO ID info
        goInfo_dict = go_ec_count[sortedPaper.PMID]
        for key, value in goInfo_dict.iteritems(): # [15, 16, 17, 18....]
            out_list.extend(list(key))
            out_list.append(str(value))
        outFile.write('\t'.join((out_list)))
        outFile.write('\n')
        out_list = []
        
    outFile.close()  
    return


###########################################################
# print_paper_per_prots_sorted_go
###########################################################
def print_paper_per_prots_sorted_go(papers_annots2_dict, all_tt_count, go_ec_count,
                             allEvCodes_dict, sortedProtsPerPaper_tuple, outpath, top=20):
    """Prints out all information that we have for each paper(PMID).
    papers_annots2_dict: dict of the top X papers, with title, year and journal name
    all_tt_count: this is a dict that gives how many term types each paper annotates
    go_ec_count: this is a dict that gives how many times a paper gives a specific (go ID, go Name, Ev code) annotation
    allEvCodes_dict: this is a dict that gives how many times a paper supports a given experimental ev Code (EEC global).
    sortedProtsPerPaper_tuple: this is a sorted named tuple (largest first) of all the papers, sorted by the number of proteins the paper annotates
    outpath: where the data will be printed out.
    top: the number of papers we want to print out. Sorted by the number of proteins the paper annotates.
    """
    outFile = open(outpath, 'w')
    #print out header line
    outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                  ('Num Prots', 'Num Annots', 'PMID', 'Title', 'Year', 'Journal', 'MFO Annot', 'BPO Annot', 'CCO Annot', 
                   'Num EXP', 'Num IDA', 'Num IEP', 'Num IGI', 'Num IMP', 'Num IPI',   'GO ID', 'GO Name', 'Ev Code', 'Num Used'))
    out_list = []
    type_list = ['molecular_function', 'biological_process', 'cellular_component']
    for sortedPaper in sortedProtsPerPaper_tuple[0:top]:
        # add the number of proteins & the PMID
        out_list.append(str(sortedPaper.numProts)) #[0]
        out_list.append(str(sortedPaper.PMID)) # will be [2]
        # get paper info
        out_list.extend(papers_annots2_dict[sortedPaper.PMID][1:]) #add title, year & journal [3, 4, 5]
        out_list.insert(1, str(papers_annots2_dict[sortedPaper.PMID][0])) # add number of annotations, [1]
        # get GO type info, [6, 7, 8]
        for t in type_list:
            try:
                out_list.append(str(all_tt_count[sortedPaper.PMID][t]))
            except:
                out_list.append('0')
        #get ev code data [9, 10, 11, 12, 13, 14]
        for ec in sorted(EEC):
            try:
                out_list.append(str(allEvCodes_dict[(sortedPaper.PMID, ec)]))
            except:
                out_list.append('0')
        # get GO ID info
        
        go_ec_dict = go_ec_count[sortedPaper.PMID]
        GO_EC_Count_collect_tuple = sort_go_ec_dict(go_ec_dict)
        
        for sortedGO in GO_EC_Count_collect_tuple: # [15, 16, 17, 18....]
            out_list.extend(list(sortedGO.go_tuple))
            out_list.append(str(sortedGO.numAnnots))
        outFile.write('\t'.join((out_list)))
        outFile.write('\n')
        out_list = []
        
    outFile.close()  
    return



###########################################################
# sort_papers_prots
###########################################################
def sort_papers_prots(papers_prots):
    """`Return the sorted tuple (sorted highest to lowest) as a 
    named tuple."""
    ProtsPerPaper_collect = collections.namedtuple('ProtsPerPaper_collect', 'numProts PMID')
    #creates a tuple of (numProts, PMID) sorted on numProts highest to lowest
    sortedProtsPerPaper_tuple = sorted([ProtsPerPaper_collect(len(value), key) for (key, value) in papers_prots.items()],
                                       reverse=True)
    return sortedProtsPerPaper_tuple


###########################################################
# sort_pmidTaxonIDs
###########################################################
def sort_pmidTaxonIDs(pmidTaxonIDs_dict):
    """Sort the dictionary pmidTaxonIDs according to the number of proteins annotated in 
    a particular species (taxonID). Return the sorted tuple (sorted highest to lowest) as a 
    named tuple."""
    #this has to have it keyed beforehand for the right taxonID... oddly enough, PMID and numprots again?
    PMIDPerTaxon_collect = collections.namedtuple('PMIDPerTaxon_collect', 'numProts taxonID_tuple')
    #creates a tuple of (numProts, taxonID_tuple) sorted on numProts highest to lowest
    sortedPMIDPerTaxon_tuple = sorted([PMIDPerTaxon_collect(value, key) for (key, value) in pmidTaxonIDs_dict.items()],
                                       reverse=True)
    return PMIDPerTaxon_collect


###########################################################
# sort_one_TaxonIDvPMIDs
###########################################################
def sort_one_TaxonIDvPMIDs(taxonIDvPMIDs_dict):
    """Sort the dictionary taxonIDvPMIDs_dict according to the number of proteins annotated in 
    a particular species (taxonID). Return the sorted tuple (sorted highest to lowest) as a 
    named tuple."""
    #this has to have it keyed beforehand for the right taxonID... oddly enough, PMID and numprots again?
    PMIDPerTaxon_collect = collections.namedtuple('PMIDPerTaxon_collect', 'numProts PMID')
    #creates a tuple of (numProts, PMID) sorted on numProts highest to lowest
    sortedPMIDPerTaxon_tuple = sorted([PMIDPerTaxon_collect(value, key) for (key, value) in taxonIDvPMIDs_dict.items()],
                                       reverse=True)
    return sortedPMIDPerTaxon_tuple

###########################################################
# sort_go_code
###########################################################
def sort_go_code(go_code_count):
    """Sort the dictionary go_code_count according to the number of times a certain go code is used  
    Return the sorted tuple (sorted highest to lowest) as a 
    named tuple."""
    GOCode_collect = collections.namedtuple('GOCode_collect', 'numAnnots go_code')
    #creates a tuple of (numAnnots, go_code) sorted on numAnnots highest to lowest
    GOCode_collect_tuple = sorted([GOCode_collect(value, key) for (key, value) in go_code_count.items()],
                                       reverse=True)
    return GOCode_collect_tuple


###########################################################
# sort_go_ec_dict
###########################################################
def sort_go_ec_dict(go_ec_dict):
    """Sort the dictionary go_ec_dict (a section of go_ec_count) according to the number of times a certain 
    go code is used. Return the sorted tuple (sorted highest to lowest) as a named tuple.
    """
    GO_EC_Count_collect = collections.namedtuple('GO_EC_Count_collect', 'numAnnots go_tuple')
    #creates a tuple of (numAnnots, go_code) sorted on numAnnots highest to lowest
    GO_EC_Count_collect_tuple = sorted([GO_EC_Count_collect(value, key) for (key, value) in go_ec_count.items()],
                                       reverse=True)
    print GO_EC_Count_collect_tuple
    return GO_EC_Count_collect_tuple



def top_ontology(papers,outpath=None,top=20):
    """Determines the top GO terms annotated in the analysis set and 1) puts it in 
    the output dict top_go and 2) writes it out to a tab delim file 'outpath'
    
    Note: this function is currently identical to top_go_terms()"""
    #
    # Can be used with SP & GOA data
    
    go_count = {}
    for p in papers:
        for rec in papers[p]:
            go_id = rec['go_id']
            go_count[go_id] = go_count.get(go_id,0) + 1
    top_go = [(i[1],i[0]) for i in go_count.items()]
    top_go.sort()
    if outpath:
        go_con = mysqlConnect()
        go_cur = go_con.cursor()
        f = open(outpath,"w")
        for i in top_go[-top:]:
            name = gu.go_acc_to_name(i[1],go_cur)
            f.write("%d\t%s\t%s\n" % (i[0], i[1], name))
        go_hist = {}
        for i in top_go:
            go_hist[i[0]] = go_hist.get(i[0],0) + 1
        go_hist_list = [(h[1],h[0]) for h in go_hist.items()] 
        go_hist_list.sort()
        fhist = open("hist_%s" % outpath, "w")
        for h in go_hist_list:
#            print h
            fhist.write("%d\t%d\n" % h)
        
        f.close()
        fhist.close()
        go_con.close()
    return top_go
 
def top_papers(papers,outpath=None,delim="\t", top=20):
    """This function fetches all the relevent PubMed info for each PMID in 'papers' and 
    1) puts it into a list and 2) outputs it to a file named in outpath."""
    #
    # Can be used with SP & GOA data
    
    papers_annots = [(len(papers[p]), p) for p in papers]
    papers_annots2 = []
        
    papers_annots.sort()
    idlist = [p[1] for p in papers_annots[-top:]]
    Entrez.email = MY_EMAIL
    h = Entrez.efetch(db="pubmed", id=idlist, 
                          rettype="medline", retmode="text")
    medrecs = list(Medline.parse(h))
    titles = [medrec.get("TI","?") for medrec in medrecs]
    years = [medrec.get("DP","?") for medrec in medrecs]
    journals = [medrec.get("JT", "?") for medrec in medrecs]
    for p, title, year, journal in zip(papers_annots[-top:], titles,years, journals):
        papers_annots2.append((p[0],p[1], title, year.split()[0].strip(), journal))
    if outpath:
        fout = open(outpath,"w")
        print >> fout, "num proteins\tpubmed ID\tTitle\tYear\tJournal"
        for p in papers_annots2:
            print >> fout, "%d\t%s\t%s\t%s\t%s" % p
        fout.close()
    #papers_annots2 = [(# all annotations, PMID, Title, Year, Journal)] 
    return papers_annots2

###########################################################
# pickle_journal_info
###########################################################
def pickle_journal_info(papers_prots, outpath=None,delim="\t", top=None):
    """This function fetches all the relevent PubMed info for each PMID in 'papers' 
    (at the limit supplied in 'top') and 1) puts it into a dict."""
    #
    # Can be used with SP & GOA data
    
    journalInfo_dict = {}
    pmid_list = papers_prots.keys()    
    

    idlist = [p[1] for p in papers_annots[negTop:]]
    print idlist
    Entrez.email = MY_EMAIL
    
    h = Entrez.efetch(db="pubmed", id=idlist, 
                          rettype="medline", retmode="text")
    medrecs = list(Medline.parse(h))
    titles = [medrec.get("TI","?") for medrec in medrecs]
    years = [medrec.get("DP","?") for medrec in medrecs]
    journals = [medrec.get("JT", "?") for medrec in medrecs]
    for p, title, year, journal in zip(papers_annots[negTop:], titles,years, journals):
        #papers_annots2_dict[PMID] = [# of total annotations, Title, Year, Journal] 
        papers_annots2_dict[p[1]] = [len(papers[p[1]]), title, year.split()[0].strip(), journal]
    """if outpath:
        fout = open(outpath,"w")
        print >> fout, "num proteins\tpubmed ID\tTitle\tYear\tJournal"
        for p in papers_annots2:
            print >> fout, "%d\t%s\t%s\t%s\t%s" % p
        fout.close()
    """
    return papers_annots2_dict



###########################################################
# top_papers_dict
###########################################################
def top_papers_dict(papers, papers_prots, outpath=None,delim="\t", top=None):
    """This function fetches all the relevent PubMed info for each PMID in 'papers' 
    (at the limit supplied in 'top') and 1) puts it into a dict."""
    #
    # Can be used with SP & GOA data
    
    papers_annots = [(len(papers_prots[p]), p) for p in papers_prots]
    papers_annots2_dict = {}
        
    papers_annots.sort()
    if top is None:
        negTop = 0
    else:
        negTop = -top
    idlist = [p[1] for p in papers_annots[negTop:]]
    print idlist
    Entrez.email = MY_EMAIL
    
    h = Entrez.efetch(db="pubmed", id=idlist, 
                          rettype="medline", retmode="text")
    medrecs = list(Medline.parse(h))
    titles = [medrec.get("TI","?") for medrec in medrecs]
    years = [medrec.get("DP","?") for medrec in medrecs]
    journals = [medrec.get("JT", "?") for medrec in medrecs]
    for p, title, year, journal in zip(papers_annots[negTop:], titles,years, journals):
        #papers_annots2_dict[PMID] = [# of total annotations, Title, Year, Journal] 
        papers_annots2_dict[p[1]] = [len(papers[p[1]]), title, year.split()[0].strip(), journal]
    """if outpath:
        fout = open(outpath,"w")
        print >> fout, "num proteins\tpubmed ID\tTitle\tYear\tJournal"
        for p in papers_annots2:
            print >> fout, "%d\t%s\t%s\t%s\t%s" % p
        fout.close()
    """
    return papers_annots2_dict


###########################################################
# term_types_all_papers
###########################################################
def ev_codes_all_papers(papers):
    """Calculate the number of times a paper gives a certain experimental evidence code.
    Possible evidence codes: EEC global defined above.
    """
    allEvCodes_dict = {}
    for pmid, go_annot_list in papers.iteritems():
        for go_annot in go_annot_list:
            go_ec = go_annot['go_ec']
            if go_ec in EEC:
                #allEvCodes_dict[(PMID, evidence code)] = # times evidence code supported by that paper
                allEvCodes_dict[(pmid, go_ec)] = allEvCodes_dict.get((pmid, go_ec), 0) + 1
    return allEvCodes_dict


###########################################################
# term_types_all_papers
###########################################################
def term_types_all_papers(papers):
    """Count up how many times each paper annotates to a certain GO group:
    Possible term types:
    biological_process
    molecular_function
    cellular_component
    """
    all_tt_count = {}
    for pmid, annot_list in papers.iteritems():
        all_tt_count[pmid] = term_types_in_paper(annot_list)
    # all_tt_count[pmid] = { {term type : term type counts} }
    return all_tt_count



def term_types_in_paper(paper):
    """For each paper, count how often different term types appear.
    Possible term types:
    biological_process
    molecular_function
    cellular_component
    """
    # term types (ontologies) in a given paper
    tt_count = {}
    go_con = mysqlConnect()
    go_cursor = go_con.cursor()
    for prec in paper:
        go_id = prec['go_id']
        term_type = gu.go_acc_to_term_type(go_id, go_cursor)
        # tt_count[term_type] = count of that term
        tt_count[term_type] = tt_count.get(term_type,0) + 1
    go_con.close()
    return tt_count # count of term types
        

###########################################################
# count_all_term_types
###########################################################
def count_all_term_types(papers):
    """For all entries, count how often different term types appear.
    Possible term types:
    biological_process
    molecular_function
    cellular_component
    """
    # term types (ontologies) in a given paper
    sum_tt_count = {}
    go_con = mysqlConnect()
    go_cursor = go_con.cursor()
    for pmid, dict_list in papers.iteritems():
        for go_dict in dict_list:
            go_id = go_dict['go_id']
            term_type = gu.go_acc_to_term_type(go_id, go_cursor)
            # sum_tt_count_dict[term_type] = count of that term over all annotations
            sum_tt_count_dict[term_type] = sum_tt_count_dict.get(term_type,0) + 1
            #print 1
            #print term_type
    go_con.close()
    printDict_one_value(sum_tt_count_dict, "TermTypeAllCount.txt")
    return sum_tt_count_dict # count of term types
        
    
def go_freq_in_papers(papers):
    # Frequency of GO terms in specific papers
    # Key: GO ID; Value: count of that GO ID.
    #
    # Can be used with SP & GOA data  
   
    go_ids = {} # all GO terms
    exp_go_ids = {} # GO terms with experimental evidence codes only
    for p in papers:
        for go_rec in papers[p]:
            # All evidence codes
            go_id = go_rec['go_id']
            go_ids[go_id] = go_ids.get(go_id,0) + 1
            # Experimental evidence codes only
            if go_rec['go_ec'] in EEC:
                exp_go_ids[go_id] = exp_go_ids.get(go_id,0) + 1
         
    return go_ids, exp_go_ids
            
    
    
def go_in_papers(sp_path):
    # Returns: papers: key: pubmed_id; value: list of go_rec records
    # go_rec record is a dictionary. Keys (values): 'sp_id' (swissprot id); 
    # 'go_id': (GO ID); 'go_ec': (GO Evidence Code).
    
    # To be used with SP data, not GOA
    
    papers = {}
    go_ids = {}
    sp_recs = {}
    papers_prots = {}
    sph = open(sp_path)
    for sp_rec in SP.parse(sph):
        cur_go_recs = get_go_evidence_codes(sp_rec)
#        print cur_go_recs
        if not cur_go_recs: 
            continue
        cur_papers = get_papers(sp_rec)
        for paper in cur_papers:
            if paper not in papers_prots:
                papers_prots[paper] = {sp_rec.entry_name: 1}
            else:
                papers_prots[paper][sp_rec.entry_name] = \
                    papers_prots[paper].get(sp_rec.entry_name,0)+1
            for cur_go_rec in cur_go_recs:
                d1 = dict(sp_id=sp_rec.entry_name,
                          go_id=cur_go_rec[0],
                          go_ec=cur_go_rec[1])
                papers.setdefault(paper,[]).append(d1)
    return papers, papers_prots        
    
###########################################################
# printDict
###########################################################
def printDict(generic_dict, fileName):
    """Just does a quick print of generic_dict to file fileName)"""
    outFile = open(fileName, "w")

    for key, value in generic_dict.iteritems():
        outFile.write('Key: ' + str(key) + '\n')
        outFile.write("Value: " + str(value) + '\n')
    outFile.write('\n')
    outFile.close()

###########################################################
# printDict
###########################################################
def printDict_one_value(generic_dict, fileName):
    """Just does a quick print of generic_dict to file fileName)"""
    outFile = open(fileName, "w")

    for key, value in generic_dict.iteritems():
        outFile.write(str(key) + '\t')
        outFile.write(str(value) + '\n')
    outFile.write('\n')
    outFile.close()




def go_in_papers_goa(goa_path):
    """Extract the GO data from the Uniprot GOA download"""
    papers = {}
    papers_prots = {}
    for inline in file(goa_path):
        if inline[0] == '!': continue
        db, db_object_id, db_object_symbol, qualifier, go_id, \
        db_reference, evidence, withit, aspect, \
        db_object_name, synonym, db_object_type, \
        taxon_id, date, assigned_by, \
        annotation_extension, gene_product_form_id = inline.rstrip('\n').split('\t')
        key_id = "%s:%s" % (db, db_object_id)

        if db_reference[:4] != "PMID": # only take the PMIDs, don't care about anything else
            continue
        paper = db_reference.split(':')[1] # Paper = the PMID
        if paper not in papers_prots:
            #papers_prots holds how many proteins each paper annotates
            #papers_prots[PMID] = {{key=Uniprot ID(key_id), value=# of times this paper annotates this protein}}
            #it is possible for one paper to produce more than one GO annotation for the same protien
            #Example: Unpript ID = Q9H1C4 and PMID = 19006693 
            papers_prots[paper] = {key_id: 1}
        else:
            papers_prots[paper][key_id] = \
                 papers_prots[paper].get(key_id,0)+1
    
        #Papers[PMID] = [a list of dicts: each dict containing 3 entries: swissProt ID entry (key='sp_id'), 
        # go_id entry (key = 'go_id'), GO evidence code (key = 'go_ec')]
        d1 = dict(sp_id=key_id,
                  go_id=go_id,
                  go_ec=evidence)
        papers.setdefault(paper,[]).append(d1)
        
        printDict(papers, 'papers.txt')
        printDict(papers_prots, 'papers_prots.txt')
    return papers, papers_prots
        
###########################################################
# go_tax_in_papers_goa
###########################################################
def go_tax_in_papers_goa(goa_path):
    """Extract the GO & taxon ID data from the Uniprot GOA download"""
    papers = {}
    papers_prots = {}
    for inline in file(goa_path):
        if inline[0] == '!': continue
        db, db_object_id, db_object_symbol, qualifier, go_id, \
        db_reference, evidence, withit, aspect, \
        db_object_name, synonym, db_object_type, \
        taxon_id, date, assigned_by, \
        annotation_extension, gene_product_form_id = inline.rstrip('\n').split('\t')
        key_id = "%s:%s" % (db, db_object_id)

        if db_reference[:4] != "PMID": # only take the PMIDs, don't care about anything else
            continue
        paper = db_reference.split(':')[1] # Paper = the PMID
        if paper not in papers_prots:
            #papers_prots holds how many proteins each paper annotates
            #papers_prots[PMID] = {{key=Uniprot ID(key_id), value=# of times this paper annotates this protein}}
            #it is possible for one paper to produce more than one GO annotation for the same protien
            #Example: Unpript ID = Q9H1C4 and PMID = 19006693 
            papers_prots[paper] = {key_id: 1}
        else:
            papers_prots[paper][key_id] = \
                 papers_prots[paper].get(key_id,0)+1
    
        #Papers[PMID] = [a list of dicts: each dict containing 4 entries: swissProt ID entry (key='sp_id'), 
        # go_id entry (key = 'go_id'), GO evidence code (key = 'go_ec'), NCBI taxon ID (key='taxon_id']
        #print taxon_id
        #taxon_id =  taxon_id[6:] # get rid of leading 'taxon:'
        #Note, in rare cases, taxon_id can contain more than one taxon ID. Currenly not represented in the 
        #top 50 papers, but it should be noted that the code currently doesn't have any way to deal with that.
        #It will just print it out.
        d1 = dict(sp_id=key_id, 
                  go_id=go_id, 
                  go_ec=evidence, 
                  taxon_id=taxon_id)
        papers.setdefault(paper,[]).append(d1)
        
    printDict(papers, 'papers.txt')
    printDict(papers_prots, 'papers_prots.txt')
    return papers, papers_prots
 
        
###########################################################
# count_top_go_terms_per_ecode_all_entries
###########################################################
def count_top_go_terms_per_ecode_all_entries(papers, outpath=None, top=20):
    """Count up how many time a GO term is used for each experimental evidence code and print out the top 'top'
    GO codes per evidence codes.
    """
    ec_go_code_count = {}
    for p in papers:
        for rec in papers[p]:
            go_id = rec['go_id']
            go_ec = rec['go_ec']
            if go_ec not in ec_go_code_count:
                ec_go_code_count[go_ec] = {go_id: 1}
            else:
                ec_go_code_count[go_ec][go_id] = \
                ec_go_code_count[go_ec].get(go_id,0)+1
    
    sorted_collection_dict = {}
    for ec in EEC:
        sorted_collection_dict[ec] = sort_go_code(ec_go_code_count[ec])
    #top_go = [(i[1],i[0]) for i in go_count.items()]
    #top_go.sort()
    if outpath:
        go_con = mysqlConnect()
        go_cur = go_con.cursor()
        f = open(outpath,"w")
        for ec in sorted(EEC):
            f.write('Evidence Code: ' + ec +'\n')
            for sortedGO in sorted_collection_dict[ec][0:top]: 
                name = gu.go_acc_to_name(sortedGO.go_code,go_cur)
                f.write("%d\t%s\t%s\n" % (sortedGO.numAnnots, sortedGO.go_code, name))
        
        
        f.close()
        go_con.close()
    return ec_go_code_count
    
###########################################################
# count_all_annotations_per_ec
###########################################################
def count_all_annotations_per_ec(papers, outpath):
    """Make and print out a dict that counts how many times a particular evidence code is used in annotation.
    Also includes the total count of how many annotations there are total for all evidence codes.
    *** Must use the 'papers' dict created from go_tax_in_papers_goa
    """
    allECCode_dict = {}
    allAnnotCount = 0
    for pmid, data_list in papers.iteritems():
        allAnnotCount = allAnnotCount + len(data_list)
        for go_dict in data_list:
            go_ec = go_dict['go_ec']
            # allECCode_dict [go ec code] = count 
            allECCode_dict[go_ec] = allECCode_dict.get(go_ec, 0) + 1 #how many times is this ev code used?
    # add 'all' count number
    allECCode_dict['all'] = allAnnotCount
    d = datetime.date.today()
    printDict_one_value(allECCode_dict, outpath)
    return allECCode_dict
        



###########################################################
# count_taxonIDs
###########################################################
def count_taxonIDs(papers):
    """Print out the taxon id info plus how many times that paper annotates that species.
    *** Must use the 'papers' dict created from go_tax_in_papers_goa
    """
    taxonID_dict = {}
    papersTaxonIDs_dict = {}
    for pmid, data_list in papers.iteritems():
        for go_dict in data_list:
            taxonID = go_dict['taxon_id']
            taxonID_dict[taxonID] = taxonID_dict.get(taxonID, 0) + 1 #how many times is this specie (taxonID) annotated
        # papersTaxonIDs_dict [PMID] = { { taxonID : count } }
        papersTaxonIDs_dict[pmid] = taxonID_dict  #associate that taxonID dict with the appropriate paper (pmid)
        taxonID_dict = {}
    printDict(papersTaxonIDs_dict, "papersTaxonIDs_dict.txt")
    return papersTaxonIDs_dict
###########################################################
# count_PMID_in_taxonID
###########################################################
def count_PMID_in_taxonID(papers):
    """Print out the taxon id info plus how many times that paper annotates that species.
    *** Must use the 'papers' dict created from go_tax_in_papers_goa
    """
    pmidTaxonIDs_dict = {}
    PMID_dict = {}
    for pmid, data_list in papers.iteritems():
        for go_dict in data_list:
            taxonID = go_dict['taxon_id']
            PMID_dict = pmidTaxonIDs_dict.setdefault(taxonID, {})            
            PMID_dict[pmid] = PMID_dict.get(pmid, 0) + 1 #increase the pmid dict by 1 
            pmidTaxonIDs_dict[taxonID] = PMID_dict  #associate that taxonID dict with the appropriate paper dict
        # pmidTaxonIDs_dict [taxonID] = { { PMID : count } }       
        PMID_dict = {}
    printDict(pmidTaxonIDs_dict, "pmidTaxonIDs_dict.txt")
    return pmidTaxonIDs_dict

###########################################################
# count_all_annotations_taxonIDs
###########################################################
def count_all_annotations_taxonIDs(papers):
    """Count up how many times a particular species is annotated in goa. Print that out to a file.
    *** Must use the 'papers' dict created from go_tax_in_papers_goa
    """
    all_taxonID_dict = {}
    for pmid, dict_list in papers.iteritems():
        for go_dict in dict_list:
            taxonID = go_dict['taxon_id']
            #all_taxonID_dict [ taxon ID ] = count
            all_taxonID_dict[taxonID] = all_taxonID_dict.get(taxonID, 0) + 1 #how many times is this specie (taxonID) annotated
    printDict_one_value(all_taxonID_dict, "AllTaxonIDsCount_dict.txt")
    return all_taxonID_dict


###########################################################
# count_all_annotations_and_proteins_taxonIDs
###########################################################
def count_all_annotations_and_proteins_taxonIDs(papers, outpath):
    """Count up how many times a particular species is annotated in goa and see how many proteins it annotates. 
    Print that out to a file.
    *** Must use the 'papers' dict created from go_tax_in_papers_goa
    """
    all_taxonID_annot_dict = {}
    all_taxonID_prot_dict = {}
    for pmid, dict_list in papers.iteritems():
        for go_dict in dict_list:
            taxonID = go_dict['taxon_id']
            uniprotID = go_dict['sp_id']
            #all_taxonID_annot_dict [ taxon ID ] = count
            #how many times is this species (taxonID) annotated?
            all_taxonID_annot_dict[taxonID] = all_taxonID_annot_dict.get(taxonID, 0) + 1 
            if taxonID in all_taxonID_prot_dict:
                if uniprotID in all_taxonID_prot_dict[taxonID]:
                    pass
                else:
                    all_taxonID_prot_dict[taxonID].append(uniprotID)
            else:
                all_taxonID_prot_dict[taxonID] = [uniprotID]
    
    out_handle = open(outpath, "w")
    
    out_handle.write("Num Prots\tNum Annots\tTaxonID\tSpecies\n")

    go_con = mysqlConnect()
    go_cur = go_con.cursor()
    for taxonID, count in all_taxonID_annot_dict.iteritems():
        speciesName = gu.go_species_name_from_taxonID(taxonID, go_cur)
        out_handle.write(str(len(all_taxonID_prot_dict[taxonID])))
        out_handle.write("\t")
        out_handle.write(str(count))
        out_handle.write("\t")
        out_handle.write(taxonID)
        out_handle.write("\t")
        out_handle.write(speciesName)
        out_handle.write("\n")
    #printDict_one_value(all_taxonID_annot_dict, "AllTaxonIDsCount_dict.txt")
    #return all_taxonID_dict
    out_handle.close()
    go_con.close()


###########################################################
# print_papers_taxonIDs
###########################################################
def print_papers_taxonIDs(papersTaxonIDs_dict, papers_annots2_dict, sortedProtsPerPaper_tuple, outpath, top=20):            
    """Print out what Taxon IDs are annotated by a paper and how many times it annotates this species. Print out only 
    the top 'top'. Also print out the related information about the paper.
    """
    outFile = open(outpath, 'w')
    #print out header line
    outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" % 
                  ('Num Prots', 'Num Annotations', 'PMID', 'Title', 'Year', 'Journal', 'Taxon ID', 'Num Used'))
    out_list = []
    for sortedPaper in sortedProtsPerPaper_tuple[0:top]:
        # add the number of proteins & the PMID
        out_list.append(str(sortedPaper.numProts)) #[0]
        out_list.append(str(sortedPaper.PMID)) # will be [2]
        # get paper info
        out_list.extend(papers_annots2_dict[sortedPaper.PMID][1:]) #add title, year & journal [3, 4, 5]
        out_list.insert(1, str(papers_annots2_dict[sortedPaper.PMID][0])) # add number of annotations, [1]
        taxonID_dict = papersTaxonIDs_dict[sortedPaper.PMID]
        print taxonID_dict
        for taxonID, num in taxonID_dict.iteritems():
            out_list.extend(([str(taxonID), str(num)])) # add taxon id info [6, 7....]
        outFile.write('\t'.join((out_list)))
        outFile.write('\n')
        out_list = []
        
    outFile.close()  
    return
    
###########################################################
# readInTaxon
###########################################################
def readInTaxon(taxonIdFile):
    """Read in a list of taxons from a file and make a list out of them"""
    inFile = open(taxonIdFile, 'r') 
    taxonID_list = []
    discardLine = inFile.readline() #File should have a title line that needs to be discarded
    for line in inFile.readlines():
        line = line.strip()
        taxonID_list.append(line)
    inFile.close()
    return taxonID_list


###########################################################
# print_papers_from_TaxonID_list
###########################################################
def print_papers_from_TaxonID_list(taxonIDFile,  pmidTaxonIDs_dict, papers_annots2_dict,  outpath, papers, oneFile=True, top=20):
    """Take in a list of taxon IDs and then print out the top 'top' of the papers that annotate that species"""
    taxonID_list = readInTaxon(taxonIDFile)
    if os.path.exists(outpath):
        os.remove(outpath)
    for taxon in taxonID_list:
        if oneFile != True: # do not print all the output to one file. 
            outpathNew = outpath + "_" + taxon + ".txt" #print it id'd by the taxonID
            print_papers_for_one_taxonID(taxon,  pmidTaxonIDs_dict, papers_annots2_dict,  outpathNew, papers, top)
        else:
            print_papers_for_one_taxonID(taxon,  pmidTaxonIDs_dict, papers_annots2_dict,  outpath, papers, top)
    


###########################################################
# print_papers_for_one_taxonID
###########################################################
def print_papers_for_one_taxonID(taxonID,  pmidTaxonIDs_dict, papers_annots2_dict, outpath, papers, top=20):            
    """Print out what Taxon IDs are annotated by a paper and how many times it annotates this species. Print out only 
    the top 'top'. Also print out the related information about the paper.
    
    NOTE BUG HERE! NOT COUNTING PROTEINS CORRECTLY!!!!!!!!!!!!!!
    
    """
    outFile = open(outpath, 'a')
    #print out header line
    #things to add... num annotations, name of species
    outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                  ('Num Prots', 'NumAnnots', 'PMID', 'Title', 'Year', 'Journal', 'Taxon ID'))
    
    out_list = []
    sortedPMIDPerTaxon_tuple =  sort_one_TaxonIDvPMIDs(pmidTaxonIDs_dict[taxonID])
    
    for sortedEntry in sortedPMIDPerTaxon_tuple[0:top]:
        out_list.append(str(sortedEntry.numProts))  #[0]
        out_list.append(str(sortedEntry.PMID)) # will be [2]
        if sortedEntry.PMID in papers_annots2_dict:
            out_list.extend((papers_annots2_dict[sortedEntry.PMID][1:])) #add title, year & journal [3, 4, 5]
        else:
            """Entrez.email = MY_EMAIL
            h = Entrez.efetch(db="pubmed", id=sortedEntry.PMID, 
                          rettype="medline", retmode="text")
            medrecs = list(Medline.parse(h))
            title = [medrec.get("TI","?") for medrec in medrecs]
            year = [medrec.get("DP","?") for medrec in medrecs]
            journal = [medrec.get("JT", "?") for medrec in medrecs]
            out_list.extend(str([title, year , journal]))
            """
            out_list.extend([""])
        
        tmp_list =  papers[sortedEntry.PMID]
        annots = len(tmp_list)   
        out_list.insert(1, str(annots)) # add number of annotations, [1]
        out_list.append(str(taxonID))  #[6]
        outFile.write('\t'.join((out_list)))
        outFile.write('\n')
        out_list = []
        
    outFile.close()  
    return
###########################################################
# sum_num_prots_annot_over_taxID
###########################################################
#def sum_num_prots_annot_over_taxID(sortedProtsPerPaper_tuple, papers):
    """Give continual sum over all the number of additional proteins annotated by each additional paper in each species"""
"""    taxID_round_dict = {}
    for tup in sortedProtsPerPaper_tuple:
        
    
    
    round = 1
    for tup in sortedProtsPerPaper_tuple:
        prot_list = papers[tup.PMID]
        if round == 1:
            
        for dict_in_list in prot_list:
            
 """           


def get_go_evidence_codes(sp_rec):
    # isolate go evidence codes from sp_rec
    go_ids = []
    for xref in sp_rec.cross_references:
#        print xref
        if xref[0] == 'GO' and len(xref) == 4:
            go_ids.append((xref[1], xref[3][:3]))
    return go_ids

def get_papers(sp_rec):
    # get all pubs association with this sp_rec
    refs = {}
    for ref in sp_rec.references:
        for refid in ref.references:
            if refid[0] == 'PubMed':
                refs[refid[1]] = ref.title
    return refs
    
from sets import Set
def redundant_annotations(go_papers_dict):
    go_con, go_cur = gu.open_go()
    ancestors_found = {}
    to_remove = {}
    gpd_leaves_only = {}
    for pmid in go_papers_dict.keys(): #[:1000]:
        if len(go_papers_dict[pmid]) < 2:
            continue
        for i in range(len(go_papers_dict[pmid])):
            for j in range(len(go_papers_dict[pmid])):
                if j >= i: continue
                go_id_1 = go_papers_dict[pmid][i]['go_id']
                go_id_2 = go_papers_dict[pmid][j]['go_id']
                if go_id_1 == go_id_2: continue
                sp_id_1 = go_papers_dict[pmid][i]['sp_id']
                sp_id_2 = go_papers_dict[pmid][j]['sp_id']
                if sp_id_1 != sp_id_2: continue
                if gu.is_ancestor(go_id_1, go_id_2, go_cur):
                    ancestors_found[(pmid,go_id_1,go_id_2)] = \
                        ancestors_found.get((pmid,go_id_1,go_id_2),0) + 1
                    to_remove.setdefault(pmid,Set([])).add(j)
                    #ancestors_found.setdefault(pmid,[]).append(go_id_1,go_id_2)
                elif gu.is_ancestor(go_id_2, go_id_1, go_cur):
                    ancestors_found[(pmid,go_id_2,go_id_1)] = \
                        ancestors_found.get((pmid,go_id_2,go_id_1),0) + 1
                    to_remove.setdefault(pmid,Set([])).add(i)
                    #ancestors_found.setdefault(pmid,[]).append(go_id_2,go_id_1)
    go_con.close()
    for pmid in go_papers_dict:
        if pmid not in to_remove:
            gpd_leaves_only[pmid] = go_papers_dict[pmid]
            continue
        else:
            gpd_leaves_only[pmid] = []
            for i in range(len(go_papers_dict[pmid])):
                if i not in to_remove[pmid]:
                    gpd_leaves_only[pmid].append(go_papers_dict[pmid][i])
            
            

    return ancestors_found, to_remove,gpd_leaves_only
