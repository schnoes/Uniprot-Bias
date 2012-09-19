import cPickle
import sp_tools
from time import clock
import sys
import datetime

d = datetime.date.today()
goa_path = "/home/alexs/gene_association.goa_uniprot-Sept162012"
#pfamList_path = 

print clock()
sys.stdout.flush()

#uniprotGOA_dict = sp_tools.go_in_papers_goa_CAFA(goa_path)
#print clock()
#sys.stdout.flush()
#sp_tools.cPickleDump(uniprotGOA_dict, "/home/alexs/uniprotGOAEEC-CAFA.pik")
#print clock()
#sys.stdout.flush()
#pfamGOA_dict = sp_tools.parse_pfam_domain_ec_CAFA(unprotGOA_dict, pfamList_path)
#sp_tools.print_pfamGOA_CAFA(pfamGOA_dict)


#load up the pre-pickled data
print "load up the pre-pickled data"
uniprotGOA_dict_handle = open('/home/alexs/uniprotGOAEEC-CAFA.pik', 'rb')
uniprotGOA_dict = cPickle.load(uniprotGOA_dict_handle)
print clock()
sys.stdout.flush()

print "RNase_PH_PF01138_UniprotIDs"
pfam_path = "/home/alexs/RNase_PH_PF01138_UniprotIDs.txt"
outpath = "/home/alexs/RNase_PH_PF01138_GOAannotations.txt"
sp_tools.print_PFAM_GOAs_CAFA(uniprotGOA_dict, pfam_path, outpath)
print clock()
sys.stdout.flush()

print "RNase_PH_C_PF03725_full_length_UniprotIDs"
pfam_path = "/home/alexs/RNase_PH_C_PF03725_full_length_UniprotIDs.txt"
outpath = "/home/alexs/RNase_PH_C_PF03725_full_length_GOAannotations.txt"
sp_tools.print_PFAM_GOAs_CAFA(uniprotGOA_dict, pfam_path, outpath)
print clock()
sys.stdout.flush()

print "PNPase_PF03726_full_length_UniprotIDs"
pfam_path = "/home/alexs/PNPase_PF03726_full_length_UniprotIDs.txt"
outpath = "/home/alexs/PNPase_PF03726_full_length_GOAannotations.txt"
sp_tools.print_PFAM_GOAs_CAFA(uniprotGOA_dict, pfam_path, outpath)
print clock()
sys.stdout.flush()

print "KH_1_PF00013_full_length_UniprotIDs"
pfam_path = "/home/alexs/KH_1_PF00013_full_length_UniprotIDs.txt"
outpath = "/home/alexs/KH_1_PF00013_full_length_GOAannotations.txt"
sp_tools.print_PFAM_GOAs_CAFA(uniprotGOA_dict, pfam_path, outpath)
print clock()
sys.stdout.flush()

print "S1_PF00575_full_length_UniprotIDs"
pfam_path = "/home/alexs/S1_PF00575_full_length_UniprotIDs.txt"
outpath = "/home/alexs/S1_PF00575_full_length_GOAannotations.txt"
sp_tools.print_PFAM_GOAs_CAFA(uniprotGOA_dict, pfam_path, outpath)
print clock()
sys.stdout.flush()
