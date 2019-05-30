#https://www.biostars.org/p/327003/"""
#from Bio import AlignIO
#from Bio import SeqIO
#alignments = AlignIO.parse("x.fa", "fasta")
#with open("example.faa", "w") as handle:
 #   count = SeqIO.write(alignments, handle, "fasta")
from Bio import SeqIO
from Bio import AlignIO
#with open("example.faa", "w") as handle:
 #   for record in AlignIO.parse("xA.fa", "fasta"):
  #      count = SeqIO.write(record, handle, "phylip")
count = AlignIO.convert("x.fa", "fasta", "example.phy", "phylip-relaxed")
print("Converted %i alignments" % count)
