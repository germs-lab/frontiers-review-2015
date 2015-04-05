import sys 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord




def translate_record(record):
    """Returns a new SeqRecord with translated sequence."""
    return SeqRecord(seq = record.seq.translate(), id =  record.id + "_translated", description = record.description)

genome=SeqIO.read(sys.argv[1], 'genbank')

proteins = []
for record in list(SeqIO.parse(sys.argv[1], 'genbank')):
    for feat in record.features:
        if feat.type == "CDS":
            #gene_id = feat.qualifiers['gene']
            if feat.qualifiers.has_key("protein_id"):
                id = feat.qualifiers['protein_id']
            else:
                id = "NCBI unnamed protein"
            if feat.qualifiers.has_key('translation'):
                sequence = feat.qualifiers['translation']
                print '>' + id[0]
                print sequence[0]
            else:
                continue
