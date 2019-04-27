import os

#A welcome message that briefly explains the aim of the script
print('Welcome to the "Determining the rate of evolution of protein-coding sequences" script.' + '\n The purpose of the program is to determine the rate of evolution of protein coding sequences by calculating the dn/ds ratio of a certain gene among different species presented in a FASTA file.\n')


fp_flag = 'false' #fp_flag variable used for file path validation
while fp_flag is 'false':
    # ask for file path as an input parameter from user
    fp = input("Please enter the full path of your FASTA file that contains CDS sequences:  ")
    # validate the existence of the file path with success message or failure message
    try:
        assert os.path.exists(fp)
        #fp_flag = true
        print('File ' + str(fp) + ' was found. ')
        fp_flag = 'true' #flag is true if user entered correct file path
    except AssertionError as e:
        print('You have entered the following file path: ' + str(
        fp) + '.'+ ' This file either does not exist.' + '\n Please try again with correct file path')
        fp_flag = 'false' #fp_flag is false if file path is incorrect

     #Subtask 1: parse input file
