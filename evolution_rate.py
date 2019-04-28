import os

#A welcome message that briefly explains the aim of the script
print('Welcome to the "Determining the rate of evolution of protein-coding sequences" script.' + '\nThe purpose of the program is to determine the rate of evolution of protein coding sequences\nby calculating the dn/ds ratio of a certain gene among different species presented in a FASTA file.')
print('*'*100)
fp_flag = 'false' #fp_flag variable used for file path validation
while fp_flag is 'false':
    # ask for file path as an input parameter from user
    fp = input("Please enter the full path of your FASTA file that contains CDS sequences:  ")
    # validate the existence of the file path with success message or failure message
    try:
        assert os.path.exists(fp)
        print('File ' + str(fp) + ' was found. ')
        if fp.endswith('.fasta')or fp.endswith('.FASTA'):
            print('File extension is correct')
            fp_flag = 'true' #flag is true if user entered correct file path
        else:
            print('File extension is incorrect, Please try again with a correct file extension')
            fp_flag = 'false'
    except AssertionError as e:
        print('You have entered the following file path: ' + str(
        fp) + '.'+ ' This file does not exist.' + '\nPlease try again with correct file path')
        fp_flag = 'false' #fp_flag is false if file path is incorrect

     #Subtask 1: parse input file
