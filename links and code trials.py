#http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:SeqIO-translate
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
def make_prot_rec(nuc_rec):
    """function that returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(seq = nuc_rec.seq.translate(cds=True),
                     id = "trans_" + nuc_rec.id,
                     description = "translation of CDS, using default table")
# A welcome message that briefly explains the aim of the script
print(
    'Welcome to the "Determining the rate of evolution of protein-coding sequences" script.' + '\nThe purpose of the program is to determine the rate of evolution of protein coding sequences\nby calculating the dn/ds ratio of a certain gene among different species presented in a FASTA file.')
print('*' * 100)
# subtask 0: ask user for input fasta file + file validation for the entered
# path and file extension(this subtask
# wasn't in the project plan but it will make the program more interactive)
#############################################################################
fp_flag = 'false'  # fpath_flag variable used for file path validation
while fp_flag is 'false':
    # ask for file path as an input parameter from user
    fp = input("Please enter the full path of your FASTA file that contains CDS sequences:  ")
    fn = os.path.basename(fp)  # parameter that extracts file name from file path
    # validate the existence of the file path with success message or failure message
    try:
        assert os.path.exists(fp)
        print('File ' + str(fp) + ' was found. ')
        if fp.endswith('.fasta') or fp.endswith('.FASTA')or fp.endswith('.fa')or fp.endswith('.FA'):
            print('File extension is correct')
            fp_flag = 'true'  # flag is true if user entered correct file path
        else:
            print('File extension is incorrect, Please try again with a correct file extension')
            fp_flag = 'false'
    except AssertionError as e:
        print('You have entered the following file path: ' + str(
            fp) + '.' + ' This file does not exist.' + '\nPlease try again with correct file path')
        fp_flag = 'false'  # fpath_flag is false if file path is incorrect

# Sub-task 1: parse input file
#codon_table = CodonTable.ambiguous_dna_by_id[1]
with open(fp) as file:  #with as method for proper handling of large files
    # regardless of the Operating system styling
    # file: a parameter that will be used as a handle
#subtask 1 + 2: continue with file parsing + validation of coding sequences only
    print("Validating the input file for the presence of complete coding sequences only")
    for record in SeqIO.parse(file, "fasta"): #for loop to parse the input file
        # record: variable of type SeqRecord object
        gene_seq = record.seq #gene_seq: is a variable of type seq object
        try:
            gene_seq.translate(cds=True)
            print(record.id , 'is ok!')
        #parameter cds - Boolean, indicates this is a complete CDS.
        # If True, this checks the sequence starts with a valid alternative start
        # codon (which will be translated as methionine, M), that the sequence
        # length is a multiple of three, and that there is a single in frame stop
        # codon at the end. If these tests fail, an exception is raised.
        except Exception as cds_error:
           print(cds_error)