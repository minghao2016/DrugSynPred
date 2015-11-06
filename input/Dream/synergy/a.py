import urllib2
import xml.etree.ElementTree as ET

with open("chembl_id_file.txt", "r") as CHEMBL:
    i=0
    for chembl_id in CHEMBL:
        try:
            url="http://cts.fiehnlab.ucdavis.edu/service/convert/ChEMBL/PubChem%20CID/"+chembl_id
            #url="https://www.ebi.ac.uk/chembl/api/data/substructure/"+smi
            response_xml = urllib2.urlopen(url).read()
            print response_xml
            
        except:
            print "aa"