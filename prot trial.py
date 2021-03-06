
from Bio import SeqIO


nuc_inFile = open(r'C:\Users\Mina Saleh\PycharmProjects\p2_evolution_rate\Human_ADAMTS6_orthologues_2.fa','r')
prot_inFile = open(r'C:\Users\Mina Saleh\PycharmProjects\p2_evolution_rate\Alg_protSeq_Human_ADAMTS6_orthologues_2.fa','r')

for n_record in SeqIO.parse(nuc_inFile, "fasta"):
    unalg_nuc_seq = n_record.seq

for p_record in SeqIO.parse(prot_inFile, "fasta"):
    alg_protein = p_record.seq

alg_nuc_seq: str = ''
print(alg_nuc_seq)
i = 0

for aa in alg_protein:
    if aa == '-':
        alg_nuc_seq+= '---'
    else:
        codon = unalg_nuc_seq[i:i + 3]
        i =i+3
        alg_nuc_seq+= codon

print(alg_nuc_seq)
print(len(alg_nuc_seq))

fname = 1
with open(fpath) as file:
    fname_alg_nuc = "alg_nuc_seq_" + fname  # fname_prot: variable of the alg_protein file name (translation step output)
    # translate neucleotide seqrecords to alg_protein seqrecords and store it in an output file

    SeqIO.write(str(alg_nuc_seq), fname_alg_nuc, "fasta")