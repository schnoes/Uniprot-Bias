#!/usr/bin/env python
# ======================================================================
# WSDbfetch REST using urllib2
#
# See:
# http://www.ebi.ac.uk/Tools/webservices/services/dbfetch_rest
# http://www.ebi.ac.uk/Tools/webservices/tutorials/python
# ======================================================================
# Load libraries
import sys
import urllib2

usage = """
  %prog [[[dbName] entryId] format]

dbName    database name ['uniprot']
entryId   entry identifier(s) ['wap_rat']
format    data format ['default']
"""

# Defaults
dbName = 'uniprot'
#entryId = 'wap_rat'
#entryId = 'Q8WWK9+OR+Q8WWL7+OR+Q8WWQ8+OR+Q8WWS9+OR+Q8WWX4+OR+Q8WWX9+OR+Q8WX92+OR+Q8WXA3+OR+Q8WXA9+OR+Q8WXB4+OR+Q8WXE9+OR+Q8WXF1+OR+Q8WXF8+OR+Q8WXG6+OR+Q8WXH0+OR+Q8WXI9+OR+Q8WXS3+OR+Q8WXT5+OR+Q8WY54+OR+Q8WYA6+OR+Q8WYB'
entryId = 'Q8WWK9'
format = 'fasta'

# Process command-line
if len(sys.argv) > 1 and sys.argv[1] == '-h':
    print usage
    exit(0)
if len(sys.argv) > 1:
    dbName = sys.argv[1]
if len(sys.argv) > 2:
    entryId = sys.argv[2]
if len(sys.argv) > 3:
    format = sys.argv[3]

# Construct URL
baseUrl = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch'
url = baseUrl + '/' + dbName + '/' + entryId
if format != None:
    url += '/' + format

# Get the entry
fh = urllib2.urlopen(url)
result = fh.read()
fh.close()

# Print the entry
print result,
