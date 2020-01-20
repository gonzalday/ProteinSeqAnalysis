# ------------------------------------------------
#                 MESSAGES MODULE
# ------------------------------------------------
# Functions to print program's help and other messages for the user.    
# MODULE IMPORTATION
from time import sleep
    
# FUNCTION DEFINITION    
def arguments():   
    print("""
O---o   This program recieves the following arguments:
 O-o     
  O     [-h | --help]  Shows the program's manual
 o-O    [-q file]      Name of the query file
o---O   [-g file]      Name(s) of the subject GenBank file(s) 
O---o   
 O-o    Please remember: the query file has to be on FASTA
  O     format. If you want to use more than one query,
 o-O    please submit all of them in a multi-FASTA file. 
o---O   
O---o   All the files to use have to be on the same
 O-o    directory as this program. Otherwise, you'll have
  O     to write the absolute route to them.
 o-O    
o---O   EXAMPLE: python main.py -q query.fasta -g genbank.gbff.txt  

Thanks for using ProteinSeqAnalysis!
    """)

def help():
    print("""
                                ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                ~ ProteinSeqAnalysis HELP ~
                                ~~~~~~~~~~~~~~~~~~~~~~~~~~~

SYNOPSIS:
~~~~~~~~~

This program makes recieves one query protein sequence (or several),
and one or more GenBank files. Extracts all the protein sequences from
the GenBank file(s) and uses them as a subject on a Blastp, with the
given query sequence(s).

With the sequences of the hits obtained on the Blastp, makes an alignment
using Muscle (or more than one, if there are several querys), and its 
corresponding phylogenetic tree.

After this, the program performs a search on query and hit sequences using
the Prosite database, looking for protein domains of interest.


USAGE:
~~~~~~

O---o   This program recieves the following arguments:
 O-o     
  O     [-h | --help]  Shows the program's manual
 o-O    [-q file]      Name of the query file
o---O   [-g file]      Name(s) of the subject GenBank file(s) 
O---o   
 O-o    Please remember: the query file has to be on FASTA
  O     format. If you want to use more than one query,
 o-O    please submit all of them in a multi-FASTA file. 
o---O   
O---o   All the files to use have to be on the same
 O-o    directory as this program. Otherwise, you'll have
  O     to write the absolute route to them.
 o-O    
o---O   EXAMPLE: python main.py -q query.fasta -g genbank.gbff.txt 


REQUIREMENTS:
~~~~~~~~~~~~~

This program requires to have Muscle and Blast installed. Also, the 
Prosite database files have to be on the program folder. Please do
NOT remove them. Another requirements are to have Biopython
and Matplitlib python modules installed.


RESULT, OUTPUTS AND SOME INDICATIONS FUNCTIONING:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The program will show you on screen all the important results, or
ask you if you want to see them. Anyway, all of them will be stored on
the program folder, in the directories created with this purpose.

Any time you run the program, you'll be given a numeric ID. All of the
files created on this run will be identified by this number, so it's easier
for you to find them. If you can't remember the numeric ID after running
the program, you can look at the LOG files. On them, the exact date and
time of start of every run is stored, so you can relate it to its ID.

The folders with the files resulting of the run of this program are:

- "logs" containing all the log files.

- "data" containing the fasta files with the query sequences and the
   protein sequences retieved from the genbanks.

- "results" containing all the results, in different subdirectories:
    - "Blastp": results on a tsv and Blastp query and hit sequences on fasta format.
    - "Muscle": alignment results on fasta format and tree results on newick format.
    (if more than one query, there will be an alignment and tree file for
    each query, identified by a number that refers to their order on the
    query fasta file, stored on the "data" directory)
    - "Prosite": prosite domain search on txt files.



This is a project for the course on Programming for Bioinformatics,
Biotechnology degree, UPM.

~~~~~~~~ Program created by Raquel González Alday, 2020 ~~~~~~~~~~~

    """)


def header():
    print("\nWELCOME TO ProteinSeqAnalysis! Thanks for using this program.\n")

def byebye():
    print("""

>> Your protein analysis has been succesfully completed.
   You can check the data, results and log files on the
   corresponding directories.

Thank you for using this program. See you!
    """)

    for x in [" \\ "," / "," - "]*5:
            print("\t"+x+"Exiting"+x, end="\r")
            sleep(0.1)
