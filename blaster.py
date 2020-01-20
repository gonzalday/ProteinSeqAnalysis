# ------------------------------------------------
#                 BLASTER MODULE
# ------------------------------------------------
# Functions to check blastp installation and to make a blastp.

# MODULE IMPORTATION
from subprocess import call
from subprocess import Popen, PIPE
from time import sleep

from Bio import SeqIO

# FUNCTION DEFINITION
def check():
    """Function to check if blastp is installed."""

    # Try to open blastp help to check if it's installed.
    try:
        test = Popen(['blastp','-h'],stdout=PIPE,stderr=PIPE)
        test.stderr.close()
        test.stdout.close()
        return True
    except:
        return False

def blastp(query,subject,output):
    """query: name of FASTA file containing the query sequeces (str).
    subject: name of FASTA file containing the subject sequeces (str).
    output: name of tsv file to store the results of the blast (str).

    This function makes a blastp. The output is stored in a file, and if any error
    during the blast, it is returned.
    """

    for x in range(0,4):
        print(">> Starting blast"+"."*x, end="\r")
        sleep(0.2)
    print()

    # Call Blastp.
    blast_process = Popen(['blastp','-query',query,'-subject',subject,'-outfmt','6 qseqid sseqid qcovs pident evalue'], stdout=PIPE, stderr=PIPE)
    blast_error,blast_results = blast_process.stderr.read(),blast_process.stdout.read()
    blast_error,blast_results = blast_error.decode("utf-8"),blast_results.decode("utf-8")

    blast_process.stderr.close()
    blast_process.stdout.close()

    # Make output file with table header.
    header = "#query_id\tsubject_id\tquery_cov\tp_ident\tevalue\n"
    with open(output,"w") as file:
        file.write(header)
        file.write(blast_results)

    if not blast_error:
        # Give to the user the choice to print the results on screen or not.
        print("\n>> BLASTP succesfully performed. If you want to see a graphical representation of your Blastp results, go to \"results/blastp\" folder and look at the png image(s).")
        print("\nDo you want to see the results table? It's also stored on \"results/blastp\" folder.")
        print("Yes: Y | No: any other key")
        selection = input("> ")
        if selection.upper() == "Y":
            print("\n> BLASTP results summary. For more info, check \"results\" folder\n") 
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(" query id ~~ subject id ~~ query cov ~~ % ident ~~ evalue ")
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
            with open(output,"r") as results:
                for line in results:
                    data = line.split("\t")
                    if line[0]!="#":
                        print(' ~~ '.join(data[0:5]))
            input("\n\nPRESS ENTER TO CONTINUE ")
        return None 
    else: 
        print("\nThe following error ocurred while performing the blast: %s\n" % blast_error)
        return blast_error


def hitfile(blastp_result,query_fasta,subject_fasta,filename):
    """blastp_result: name of tsv file with the blastp results (str).
    query_fasta: name of the FASTA file containing the query sequences (str).
    subject_fasta: name of the FASTA file containing all the subject sequences used on blastp (str).
    filename: name of output file (str).

    This function creates a FASTA file with the sequences of the querys (with sufix _QUERYSEQ
    after their name to recognise them later) and their blastp hits, keeping the order of the
    blastp output. After each query sequence, its hits.
    """

    with open(blastp_result,"r") as hits, open(filename,"w") as output:
        qnames = []
        # Parse the tsv. If not the header, look into the line.
        for line in hits:
            data = line.split("\t")
            if line[0]!="#":
                if data[0] not in qnames:
                    # If the query seq has still not been stored, store it.
                    qnames.append(data[0])
                    with open(query_fasta, "r") as querys:
                        for record in SeqIO.parse(querys, "fasta"):
                            if record.id == data[0]:
                                output.write(">%s_QUERYSEQ\n" % record.id)
                                output.write("%s\n" % str(record.seq))
                # Store sequences of all subject sequences that are a hit.
                with open(subject_fasta, "r") as subjects:
                    for record in SeqIO.parse(subjects, "fasta"):
                        if record.id == data[1]:
                            output.write(">%s\n" % record.id)
                            output.write("%s\n" %str(record.seq))
                    