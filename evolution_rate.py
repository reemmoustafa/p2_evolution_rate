
""""
*
* The development repository for this tool is here:
      https://github.com/reemmoustafa/p2_evolution_rate.git
* copy rights for: Reem Mostafa, Nehal Ghonim, Mina S. A. Saleh, Manar M. Rashad, Abeer Shalaby
* Under supervision of: Dr. Rosina
* This algorithm is part of CIT-656:Programming for Bioinformatics Course (WINTER 2019)
* Bioinformatics Diploma, Nile University, Cairo, Egypt
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
# from io import StringIO
# from Bio import AlignIO


def make_prot_rec(nuc_rec):
    """function that returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(seq=nuc_rec.seq.translate(cds=True),
                     id=nuc_rec.id, description=""
                     )


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
# codon_table = CodonTable.ambiguous_dna_by_id[1]
with open(fpath) as file:  # with as method for proper handling of large files
    # regardless of the Operating system styling
    # file: a parameter that will be used as a handle
    # subtask 1 + 2: continue with file parsing + validation of coding sequences only
    print("Validating the input file for the presence of complete coding sequences only")
    for record in SeqIO.parse(file, "fasta"):  # for loop to parse the input file
        # record: variable of type SeqRecord object
        gene_seq = record.seq  # gene_seq: is a variable of type seq object
        try:
            gene_seq.translate(cds=True)
        # parameter cds - Boolean, indicates this is a complete CDS.
        # If True, this checks the sequence starts with a valid alternative start codon
        # (which will be translated as methionine, M),
        # that the sequence length is a multiple of three, and that there is a single in frame stop codon at the end.
        # If these tests fail, an exception is raised.
        except Exception as cds_error:
            print('Error in Sequence id '+str(record.id)+':\n\t'+str(cds_error)+'\n\t'
                  + str(gene_seq))  # print of sequence id + exact error + coding sequence
            print('\tUnfortunately, The Program will terminate now.')  # exit message for the user
            exit()  # program will crash and exit as required
    print("File contains complete CDS only. File is accepted")


# subtask 3: Convert the coding sequences to alg_protein sequences by translation
with open(fpath) as file:
    fname_prot = "protSeq_" + fname  # fname_prot: variable of the alg_protein file name (translation step output)
    # translate neucleotide seqrecords to alg_protein seqrecords and store it in an output file
    proteins = (make_prot_rec(nuc_rec) for nuc_rec in SeqIO.parse(file, "fasta"))
    SeqIO.write(proteins, fname_prot, "fasta")


# subtask 4: Align the alg_protein sequences using MUSCLE program
# fn_p_muscle = "Alg_"+fname_prot
muscle_exe = "muscle3.8.31_i86win32.exe"
# muscle_exe: variable containing the path of muscle program
fname_prot_musc_in = fname_prot  # variable: input file for muscle
fname_prot_musc_out = "Alg_" + fname_prot_musc_in  # variable: output file from muscle
muscle_cline = MuscleCommandline(muscle_exe, input=fname_prot_musc_in, out=fname_prot_musc_out)
# variable for muscle commandline
print(muscle_cline)  # print statement of the commandline variable
stdout, stderr = muscle_cline()  # stdout, stderr runs muscle command variable
# align = AlignIO.read(StringIO(stdout), "fasta")
# print(align)
