# ------------------------------------------------
#                 STARTER MODULE
# ------------------------------------------------
# Basic functions to check files and prepare the program's run.

# MODULE IMPORTATION
import os
import re
from subprocess import Popen, PIPE
from subprocess import call
from time import sleep

from Bio import Seq
from Bio import SeqIO

# FUNCTION DEFINITION
def organizer():
    """This function creates the project folders (if they don't already exist),
    and returns a numeric id for the program run.
    """
    
    max_files = []
    for dir in ("data","logs","results","./results/blastp","./results/muscle","./results/prosite"):
        # Create the directory if it doesn't exist.
        if not os.path.exists(dir):
            os.mkdir(dir)
        # List all the elements on the diferent folders to get the highest number id of the files on them.
        # Except results, because we are interested on its subdirectories.
        if dir != "results":
            ls = Popen(['ls',dir], stdout=PIPE, stderr=PIPE)
            ls_error,ls_results = ls.stderr.read(),ls.stdout.read()
            ls_error,ls_results = ls_error.decode("utf-8"),ls_results.decode("utf-8")
            ls.stderr.close()
            ls.stdout.close()
            if not ls_error:
                if not ls_results:
                    max_files.append(0)
                else:
                    if dir == "./results/muscle":
                        numbers = list(map(int,re.findall(r"(\d{1,})\.[ta]",ls_results)))
                        max_files.append(max(numbers))
                    else:
                        numbers = list(map(int,re.findall(r"\d{1,}",ls_results)))
                        max_files.append(max(numbers))
            else:
                print("The following error has ocurred: %s" % ls_error)
                exit()
    # Number ID for the process: next integrer to the highest id stored.
    run_id = max(max_files) + 1
    
    return run_id


def querys(fasta,output):
    """fasta: name of the file containing the query sequence(s) (str).
    output: name of the output file (str). 

    This function checks the following:
    - If the given file exists
    - If the given file is on FASTA format
    - If the given file contains protein sequences

    If the file meets the requirements, the user can choose the query sequences 
    to use (if there's more than one). At the end, a file containing the query 
    sequences to use is created. Also, the total number of querys is returned.
    """

    for x in range(0,4):
        print(">> Checking query file"+"."*x, end="\r")
        sleep(0.2)
    print()
    
    # Check if the file exists and give the oportunity of trying again before exiting.
    if not os.path.exists(fasta):
        print("\nERROR. The file %s doesn't exist or isn't found on this directory.\n" % fasta)
        selection = input("Try again: T | Exit: any other key\n")
        if selection.upper() == "T":
            new = input("> Name of FASTA file containing the query(s): ")
            querys(new,output)
        else:
            exit()
    else:
        # Check FASTA format (search for lines that start by >)
        with open(fasta,"r") as f:
            n = 0
            seqs = []
            print("\nFASTA format will be checked. ")
            print("Do you wish to check if there are strange characters on the sequences to prevent further errors also?")
            print("Yes: Y | No: any other key")
            aacheck = input("> ")
            for line in f.readlines():
                if line[0] == ">":
                    n += 1
                    seqs.append(line[1:-1])
                else:
                    if n == 0:
                        print("\nERROR. The query file is not on FASTA format.\n")
                        exit()
                    else:
                        # Check that there are no strange characters in the sequence
                        if aacheck.upper() == "Y":
                            for letter in line:
                                if letter.upper() not in "\nACDEFGHIKLMNPQRSTVWYUO":
                                    print("ERROR. The query file contains characters that don't correspond with protein sequences")
                                    exit()
                    
        print("\nQuery file is correct.")
        print("You have provided a FASTA file with %i sequences.\n" % n)

        # If more than one query on the file, give the option to select which of them to use.
        # Get a list of the numbers of the querys to use.
        if n > 1:
            print("Do you wish to use all of these querys?")
            print("Yes: Y | No: N | Exit: any other key")
            selection = input("> ")
            if selection.upper() == "N":
                print("\nSequences contained on the file:")
                for i in range(0,n):
                    print("%i -> %s" % (i, seqs[i]))

                print("\nPlease write the numbers of the querys you want to use (0 - %i), separed by commas. To exit press Q." % n) 
                use = input("> ")
                use = use.split(",")
                
                # Check if the given numbers are correct. If not, give the option to try again.
                correct = False
                while correct == False:
                    error = 0
                    if use[0].upper() == "Q":
                        exit()
                    elif use[0] == "":
                        print("\nYou haven't selected any querys. Please try again (0 - %i) or press Q to exit." % n)
                        use = input("> ")
                        use = use.split(",")
                        correct = False
                    else:
                        for x in use:
                            try:
                                int(x)
                                if int(x) >= n or int(x) < 0:
                                    error +=1
                            except:
                                error += 1
                            
                        if error > 0:
                            print("\nYou have introduced characters that are not numbers or avaiable options. Please try again (0 - %i) or press Q to exit." % n)
                            use = input("> ")
                            use = use.split(",")
                            correct = False
                        else:
                            correct = True
                
                print("\nYou have selected the querys %s." % ', '.join(use))
                use = list(map(int, use))
                nquerys = len(use)
    
            elif selection.upper() == "Y":
                use = list(range(0,n))
                nquerys = len(use)
            else:
                exit()
        else:
            use = [0]
            nquerys = len(use)

        # Parse the original file, and copy the desired sequences (those on the list "use") in new fasta file.
        s,query = 0,False
        with open(fasta,"r") as f, open(output,"w") as qtemp:
            for line in f.readlines():
                if line[0] == ">":
                    if s in use:
                        qtemp.write("%s" % line)
                        query = True
                    else:
                        query = False    
                    s += 1
                elif query == True:
                    qtemp.write("%s" % line)
        return nquerys

                
                  
def genbanks(files):
    """files: list of GenBank files (list).

    This function checks if all the files provided exist.
    """
    
    for x in range(0,4):
        print(">> Checking GenBank files"+"."*x, end="\r")
        sleep(0.2)
    print()

    # Check if all the given files exist.
    errors = []
    for file in files:
        if not os.path.exists(file):
            errors.append(file)
    
    # Delete repeated names of files.
    i,j = 0,0
    while i < len(files):
        while j < len(files):
            if (files[i] == files[j]) and (i != j):
                repeat = files.pop(j)
                print("\n(File %s was repeated, so one of its instances was ommited.)" % repeat)
            else:
                j+=1
        i+=1

    # If any of them is not found, give the chance to try again.
    if errors != []:
        print("\nERROR. The following GenBank files don't exist or are not found on this directory: %s.\n" % ', '.join(errors))
        selection = input("Try again: T | Exit: any other key\n")
        if selection.upper() == "T":
            new = input("\n> Name(s) of GenBank file(s). If more than one, please use commas to separe them: ")
            new = new.split(",")
            new_genbanks = []
            for item in new:
                if item.strip(" ") != "":
                    new_genbanks.append(item.strip(" "))
            if new_genbanks == []:
                print("\nERROR. You haven't provided any GenBank files.\n")
                exit()
            else:
                genbanks(new_genbanks)
        else:
            exit()
    else:
        n = len(files)
        print("\nSUCCESS. %i GenBank files have been provided.\n" % n)
    

def multifaster(genbanks,output):
    """genbanks: list of GenBank files (list).
    output: name of output file (str).

    This function creates a multifasta file with all the protein sequences on the provided genbanks.
    """

    # With Bio.SeqIO, parse GenBank file and store coding sequences on multifasta file.
    for file in genbanks:
        with open(file, "r") as input_handle, open(output,"w") as gtemp:
            for record in SeqIO.parse(input_handle, "genbank"):
                for feature in record.features:
                    if feature.type == 'CDS':
                        gtemp.write(">"+feature.qualifiers['locus_tag'][0]+"\n"+feature.qualifiers['translation'][0]+"\n")


