
""""
*
* The development repository for this tool is here:
      https://github.com/reemmoustafa/p2_evolution_rate.git
* copy rights for: Reem Mostafa, Nehal Ghonim, Mina S. A. Saleh, Manar M. Rashad, Abeer Shalaby
* Under supervision of: Dr. Rosina
* This algorithm is part of CIT-656:Programming for Bioinformatics Course (WINTER 2019)
* Bioinformatics Diploma, Nile University, Cairo, Egypt
* To run this script , please type python evolution_rate.py in terminal
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO


def make_prot_rec(nuc_rec):
    """function that returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(seq=nuc_rec.seq.translate(cds=True),
                     id=nuc_rec.id, description=""
                     )

def algprot_to_algdna(mus_alg_p, unalg_nuc_origf):
    """function module that takes aligned proteins and convert them back to aligned DNA sequences with fasta file generation"""
    str_to_write =""
    with open(mus_alg_p) as algd_p_file:  # opening the aligned protein file
        with open(unalg_nuc_origf) as file:  # opening the unaligned neucelotide sequence file
            unalg_nuc_seqs = list(SeqIO.parse(file, "fasta"))  # unalg_nuc_seqs vriable
            # list of seqrecords of unaligned neuclotide sequences
            p_alignment = AlignIO.read(algd_p_file, "fasta")  # variable for the aligned
            #str_to_write += str((len(p_alignment)))
            # protein seuqnces
            # print(len(p_alignment))
            for unalg_nuc in range(len(unalg_nuc_seqs)):
                # for loop on the unalg_nuc_seqs, this for loop will be work in range
                # of the length of the aunlaigned nuc sequences , it will extract its id
                # and sequence in two variables unalg_nuc_id & unalg_nuc_seq. then it
                # will contain nested for loop for the rest of program logic
                i = 0  # counter = 0 , reset will always happen when loop starts
                alg_nuc_seq = ""  # alg_nuc_seq variable that will carry the aligned , set as empty
                # string when the for loop starts neucleotide sequence
                unalg_nuc_id = unalg_nuc_seqs[unalg_nuc].id
                unalg_nuc_seq = unalg_nuc_seqs[unalg_nuc].seq
                for prot_index in range(len(p_alignment)):  # A nested for loop that
                    # will loop on the proteins contained in p_alignment variable.
                    if p_alignment[prot_index].id == unalg_nuc_id:
                        # if condition that will compare the protein id to the unaligned
                        # nuceleotide sequence id , if they are equal then the protein id
                        # and its sequence will be stored in two variables alg_protein_id
                        # and alg_protein_seq. Also a for loop will be executed for every
                        # amino acid the protein sequence and convert it the corresponding
                        # triple codon
                        alg_protein_id = p_alignment[prot_index].id
                        alg_protein_seq = p_alignment[prot_index].seq
                        for aa in alg_protein_seq:
                            # for loop on amino acids present in a protein sequence for
                            # prtoein back-translation into dna, where a gap - will translated
                            # into --- and an aminoacid will be translated into
                            # triple neucelotides inside the alg_nuc_seq that will carry the converted
                            # aligned neucloetide sequence.
                            if aa == '-':
                                alg_nuc_seq += '---'
                            else:
                                codon = unalg_nuc_seq[i:i + 3]  # the codon will be 3 positions the the
                                # unalg_nuc_seq in order to be 3 neucleotide bases
                                i = i + 3  # counter i will be incremented by 3 bases
                                alg_nuc_seq += codon
                        # print(len(unalg_nuc_seq))#print(len(alg_protein_seq)) #print(len(alg_nuc_seq))
                        # print(unalg_nuc_seq) #print(alg_protein_seq) #print(alg_nuc_seq)

                        str_to_write += '>' + unalg_nuc_id + "\n" + str(alg_nuc_seq) + "\n"
            return (str_to_write)


# A welcome message that briefly explains the aim of the script
print(
    'Welcome to the "Determining the rate of evolution of alg_protein-coding sequences" script.'
    + '\nThe purpose of the program is to determine the rate of evolution of alg_protein coding sequences'
      '\nby calculating the dn/ds ratio of a certain gene among different species presented in a FASTA file.')
print('*' * 100)

# subtask 0: ask user for input fasta file
# file validation for the entered path and file extension
# this subtask wasn't in the project plan but it will make the program more interactive
#############################################################################

fpath_flag = False  # fpath_flag variable used for file path validation
while fpath_flag is False:
    # ask for file path as an input parameter from user
    fpath = input("Please enter the full path of your FASTA file that contains CDS sequences:  ")
    fname = os.path.basename(fpath)  # parameter that extracts file name from file path
    # validate the existence of the file path with success message or failure message
    try:
        assert os.path.exists(fpath)
        print('File ' + str(fpath) + ' was found. ')
        if fpath.endswith('.fasta') or fpath.endswith('.FASTA')or fpath.endswith('.fa')or fpath.endswith('.FA'):
            print('File extension is correct')
            fpath_flag = True  # flag is true if user entered correct file path
        else:
            print('File extension is incorrect, Please try again with a correct file extension')
            fpath_flag = False
    except AssertionError as e:
        print('You have entered the following file path: ' + str(
            fpath) + '.' + ' This file does not exist.' + '\nPlease try again with correct file path')
        fpath_flag = False  # fpath_flag is false if file path is incorrect

# Sub-task 1: parse input file

with open(fpath) as file:
    # with as method for proper handling of large files regardless of the Operating system styling file:
    # a parameter that will be used as a handle
    # subtask 2 + 3: A- validation of coding sequences only.
    #                B- Convert the coding sequences to alg_protein sequences by translation
    try:
        unalg_nuc_seqs = SeqIO.parse(file, "fasta")
        fname_prot = "protSeq_" + fname  # fname_prot: variable of the alg_protein file name (translation step output)
        # translate neucleotide seqrecords to alg_protein seqrecords and store it in an output file
        proteins = (make_prot_rec(nuc_rec) for nuc_rec in SeqIO.parse(file, "fasta"))
        SeqIO.write(proteins, fname_prot, "fasta")
    except Exception as cds_error:
        print('Error in Sequence id ' + str(record.id) + ':\n\t' + str(cds_error) + '\n\t'
              + str(gene_seq))  # print of sequence id + exact error + coding sequence
        print('\tUnfortunately, The Program will terminate now.')  # exit message for the user
        exit()  # program will crash and exit as required
    print("All the sequences successfully passed filters for ORF integrity.File contains complete CDS only."
          "\nFile is accepted and translated to proteins")

# subtask 4: Align the alg_protein sequences using MUSCLE program

muscle_exe = input("Please enter the full path of muscle program:  ")
fname_prot_musc_in = fname_prot  # variable: input file for muscle
fname_prot_musc_out = "Alg_" + fname_prot_musc_in  # variable: output file from muscle
muscle_cline = MuscleCommandline(muscle_exe, input=fname_prot_musc_in, out=fname_prot_musc_out)
# variable for muscle commandline
print(muscle_cline)  # print statement of the commandline variable
stdout, stderr = muscle_cline()  # stdout, stderr runs muscle command variable

#subtask 5: converting protein back to dna (Protein back translation to DNA)

prot_covert_dna = algprot_to_algdna(fname_prot_musc_out,fpath) #coverting protein
#to dna through calling fuction prot_covert_dna and storing its output in
# prot_covert_dna variable
no_of_species = prot_covert_dna.count(">")#get number of species in file
for line in prot_covert_dna: #for loop to calculate sequence length
    if not line.startswith('>'):
        seq_len = len(line)
        break

f_phylip = "Alg_NucSeq_" +fname[:-3]+'.phy' #phylip file creation
f = open(f_phylip, 'w+')
f.write(str(no_of_species)+"      "+str(seq_len)+'\n') #add phylip file header
print(prot_covert_dna.count('\n'))
for species.split('\n') in prot_covert_dna:
    #if species.startswith('>'):
        #species_name = species
    print (species)

f.close()

