# ------------------------------------------------
#                 MUSCLE MODULE
# ------------------------------------------------
# Functions to check muscle installation and to make alignments and trees.

# MODULE IMPORTATION
import re
from subprocess import Popen, PIPE
import sys
from time import sleep

# FUNCTION DEFINITION
def check():
    """Function to check it muscle is installed"""

    # Try to call muscle to check if it is installed
    try:
        test = Popen(['muscle'],stdout=PIPE,stderr=PIPE)
        test.stderr.close()
        test.stdout.close()
        return True
    except:
        return False

def align(subject,alignment):
    """subject: name of the FASTA file containing the sequences to align (str).
    alignment: name of the output file to save the obtained alignment (str).

    This function uses muscle to align the given sequences, and stores the result.
    """

    for x in range(0,4):
        print(">> Starting muscle alignment"+"."*x, end="\r")
        sleep(0.2)
    print()

    # Make alignment
    alignment_process = Popen(['muscle','-in',subject,'-out',alignment], stderr=PIPE)
    alignment_error = alignment_process.stderr.read()
    alignment_error = alignment_error.decode("utf-8")
    alignment_process.stderr.close()

    # Check if there has been any errors. 
    if re.search(r"\*\*\* ERROR \*\*\*",alignment_error) != None:
        return alignment_error
    else:
        #If not, print standard error output, because it has interesting info.
        print(alignment_error)
        return None
    

def tree(alignment,tree):
    """alignment: name of the file containing an alignemnt of sequences (str).
    tree: name of the outpt file to save the obtained tree (str).

    This function uses muscle to create a neighbor-joining tree using the given alignment
    and stores the result.
    """

    for x in range(0,4):
        print(">> Starting muscle tree build"+"."*x, end="\r")
        sleep(0.2)
    print()

    # Make phylogenetic tree with muscle.
    tree_process = Popen(['muscle','-maketree','-in',alignment,'-out',tree,'-cluster','neighborjoining'], stderr=PIPE)
    tree_error = tree_process.stderr.read()
    tree_error = tree_error.decode("utf-8")
    tree_process.stderr.close()

    # Check if there has been any errors.
    if re.search(r"\*\*\* ERROR \*\*\*",tree_error) != None:
        return tree_error
    else:
        #If not, print standard error output, because it has interesting info.
        print(tree_error)
        # Show phylogenetic tree.
        print("Tree obtained, on Newick format (one-line):\n")
        with open(tree,"r") as tree:
            for line in tree:
                print(line.replace("\n",""),end="")
        print()
        return None
