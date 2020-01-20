# Program created by Raquel González Alday - 2020

# ------------------------------------------------
#                   MODULE IMPORT
# ------------------------------------------------
import os
import sys
from datetime import datetime
from subprocess import call

import blaster
import graph
import messages
import muscle
import prosite
import starter


# ------------------------------------------------
#          PROGRAM INIZIALITATION
# ------------------------------------------------
call('clear')

# Check the given arguments.
if len(sys.argv) == 1:
    print("\n>> You have provided no arguments.\n")
    messages.arguments()
    exit()
elif len(sys.argv) == 2 and sys.argv[1] in ["-h","--help"]:
    messages.help()
    exit()
elif len(sys.argv) >= 5:
    query,genbank = False,False
    query_file,genbank_files = [],[]
    for arg in sys.argv[1:]:
        if arg == "-q":
            query,genbank = True,False
        elif arg == "-g":
            query,genbank = False,True
        else:
            if query == True:
                query_file.append(arg)
            elif genbank == True:
                genbank_files.append(arg)
    if query_file == []:
        print("\n>> You haven't provided a query file. Remember:\n")
        messages.arguments()
        exit()
    elif genbank_files == []:
        print("\n>> You haven't provided any GenBank files. Remember:\n")
        messages.arguments()
        exit()
    elif len(query_file) > 1:
        print("\n>> You have provided more than one query file.")
        print("Please submit only one query file. Remember:\n")
        messages.arguments()
        exit()
else:
    print("\n>> You have provided not enough arguments. You have to submit")
    print("at least 4 arguments (-q query_file -g genbank_file). Remember:\n")
    messages.arguments()
    exit()

# Check if the needed programs and databases are avaiable.
blastp_inst = blaster.check()
if blastp_inst == False:
    print("\nERROR. Blastp is not installed. Pease install it before using this program.\n")
    exit()

muscle_inst = muscle.check()
if muscle_inst == False:
    print("\nERROR. Muscle is not installed. Pease install it before using this program.\n")
    exit()

prosite_inst = prosite.check()
if prosite_inst != []:
    print("\nERROR. Prosite file(s) missing: %s. Please place them on the program folder again.\n" % ', '.join(prosite_inst))
    exit()

messages.header()

# Prepare directories, log file and data and result filenames.
run_id = starter.organizer()

log_filename = "./logs/log" + str(run_id) + ".txt"
q_filename = "./data/" + str(run_id) + "_query.fasta"
gb_filename = "./data/" + str(run_id) + "_gbsequences.fasta"

now = datetime.now()
now = now.strftime("%m/%d %H:%M:%S")

# Create log file.
log = open(log_filename,"w")
log.write("PROCESS ID: %d. Time of start (day/month h:m:s): %s.\n" % (run_id,now))

# Check wether given GenBank files exist.
starter.genbanks(genbank_files)
# Check query file and save number of querys to use.
nquerys = starter.querys(query_file[0],q_filename)
log.write("> Query sequences succesfully read and stored on directory data.\n")
# Create multifasta with all protein sequences from GenBank files.
starter.multifaster(genbank_files,gb_filename)
log.write("> GenBank files succesfully read and stored on directory data.\n")


# ------------------------------------------------
#           FUNCTIONING PROGRAM
# ------------------------------------------------
call('clear')
print("\nPROCESS NUMBER ID: %d\n" % run_id)
print(">> Succesfully created and stored data files on \"data\" directory.\n")

# Perform Blastp with given query and proteins from GenBank files.
blast_result_filename = "./results/blastp/blastp" + str(run_id) + ".result.tsv"
blast_error = blaster.blastp(q_filename,gb_filename,blast_result_filename)
# Check if any errors happened to write them on the log file.
if blast_error != None:
    log.write("> The following error ocurred while performing the blastp: %s\n" % blast_error)
    exit()
else:
    log.write("> Blastp succesfully performed.\n")

# Create a file with query sequences and their blast hits.
seqs_filename = "./results/blastp/blastp" + str(run_id) + ".seqs.fasta"
blaster.hitfile(blast_result_filename,q_filename,gb_filename,seqs_filename)

# Create graphs with Blastp results.
graph.blastp(blast_result_filename,nquerys,run_id)

call('clear')
print("\nPROCESS NUMBER ID: %d\n" % run_id)
print(">> Blastp performed and stored results and sequences on \"results/blastp\" directory.\n")

# Make alignment using muscle, with query + blast hits. One alignment for each query.
for i in range(0,nquerys):
    # If more than one query, make different temporal files with the needed sequences.
    alignment_filename = "./results/muscle/muscle" + str(run_id) + ".alignment" + str(i) + ".fasta"
    temp_seqs = "temp_file_seqs" + str(run_id) + "." + str(i)
    j,query = 0,False
    with open(seqs_filename,"r") as seqs, open(temp_seqs,"w") as temp:
        for line in seqs:
            if (line[0] == ">") and ("QUERYSEQ" in line):
                if  j == i:
                    query = True
                    j += 1
                    temp.write(line)
                else:
                    query = False
                    j += 1
            if query == True:
                temp.write(line)     
            
    alignment_error = muscle.align(temp_seqs,alignment_filename)
    # Check if any errors happened to write them on the log file. Remove temporal files.
    if alignment_error != None:
        print("\n> Alignment %d of %d failed. An error ocurred:" % (i+1,nquerys))
        print(alignment_error)
        log.write("> Alignment %d of %d failed. An error ocurred:\n" % (i+1,nquerys))
        log.write(alignment_error)
        os.remove(temp_seqs)
        exit()
    else:
        os.remove(temp_seqs)
        print("\n> Alignment %d of %d completed.\n" % (i+1,nquerys))
        log.write("> Alignment %d of %d completed.\n" % (i+1,nquerys))


input("\n\nPRESS ENTER TO CONTINUE ")
call('clear')
print("\nPROCESS NUMBER ID: %d\n" % run_id)
print(">> Alignment(s) performed and stored on \"results/muscle\" directory.\n")

# Create neighbour-joining tree(s) using muscle, given the previous alignmet(s).
for i in range(0,nquerys):
    alignment_filename = "./results/muscle/muscle" + str(run_id) + ".alignment" + str(i) + ".fasta"
    tree_filename = "./results/muscle/muscle" + str(run_id) + ".tree" + str(i) + ".phy"            
    
    tree_error = muscle.tree(alignment_filename,tree_filename)
    # Check if any errors happened to write them on the log file.
    if tree_error != None:
        print("\n> Tree %d of %d failed. An error ocurred:" % (i+1,nquerys))
        print(tree_error)
        log.write("> Tree %d of %d failed. An error ocurred:\n" % (i+1,nquerys))
        log.write(tree_error)
        exit()
    else:
        print("\n> Tree %d of %d completed.\n" % (i+1,nquerys))
        log.write("> Tree %d of %d completed.\n" % (i+1,nquerys))


input("\n\nPRESS ENTER TO CONTINUE ")
call('clear')
print("\nPROCESS NUMBER ID: %d\n" % run_id)
print(">> Tree(s) built and stored on \"results/muscle\" directory.\n")

# Prosite pattern search on querys and blast hits.
prohits_filename = "./results/prosite/prosite" + str(run_id) + ".hits.txt"
prosite.patfinder(seqs_filename,prohits_filename)
log.write("> Prosite domain search finished.\n")
log.write(">> Finished analysis.")

call('clear')
print("\nPROCESS NUMBER ID: %d\n" % run_id)
print(">> Finished prosite search, results stored on \"results/prosite\" directory.\n")

messages.byebye()

# Close log file before exiting.
log.close()
exit()
