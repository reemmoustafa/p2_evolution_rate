import os
from Bio import SeqIO
from Bio.Data import CodonTable
# A welcome message that briefly explains the aim of the script
print(
    'Welcome to the "Determining the rate of evolution of protein-coding sequences" script.' + '\nThe purpose of the program is to determine the rate of evolution of protein coding sequences\nby calculating the dn/ds ratio of a certain gene among different species presented in a FASTA file.')
print('*' * 100)
# subtask 0: ask user for input fasta file + file validation for the entered
# path and file extension(this subtask
# wasn't in the project plan but it will make the program more interactive)
#############################################################################
fp_flag = 'false'  # fp_flag variable used for file path validation
while fp_flag is 'false':
    # ask for file path as an input parameter from user
    fp = input("Please enter the full path of your FASTA file that contains CDS sequences:  ")
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
        fp_flag = 'false'  # fp_flag is false if file path is incorrect

# Sub-task 1: parse input file
#codon_table = CodonTable.ambiguous_dna_by_id[1]
with open(fp) as file:  #with as method for proper handling of large files
    # 'rU' The argument 'rU' means open for reading using universal readline mode
    # regardless of the Operating system styling
    # file is a parameter that will be used as a handle
    for record in SeqIO.parse(file, "fasta"): #for loop to parse the input file
        # , record is a seqrecord object
        gene_seq = record.seq #gene_seq is a seq object
        print("Validating the input file for the presence of complete coding sequences only")
        try:
            gene_seq.translate(cds=True)
            print("File contains complete CDS only. File is accepted")
        except Exception as cds_error:
            print('Error in Sequence id '+str(record.id)+':\n\t'+str(cds_error)+'\n\t'
                  +str(gene_seq))
            print('\tUnfortunately, The Program will terminate')
            break

