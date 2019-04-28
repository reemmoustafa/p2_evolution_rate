from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
def make_prot_rec(nuc_rec):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(seq = nuc_rec.seq.translate(cds=True),
                     id = "trans_" + nuc_rec.id,
                     description = "translation of CDS, using default table")
