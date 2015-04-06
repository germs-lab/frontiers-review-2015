import sys 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

handle = open(sys.argv[1])
iterator = SeqIO.parse(handle, "Genbank")
records = list(SeqIO.parse(sys.argv[1], 'genbank'))
for x in records:
    taxonomy = x.annotations['taxonomy']
    for feat in x.features:
        if feat.type == "CDS":
            #gene_id = feat.qualifiers['gene']
            if feat.qualifiers.has_key("locus_tag"):
                id = feat.qualifiers['locus_tag']
            else:
                id = "NCBI unnamed protein"
            if feat.qualifiers.has_key('translation'):
                sequence = feat.qualifiers['translation']
                print id[0] + '\t' + '\t'.join(taxonomy)
            else:
                continue



