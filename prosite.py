# ------------------------------------------------
#                 PROSITE MODULE
# ------------------------------------------------
# Functions to check Prosite database and search domain patterns.

# MODULE IMPORTATION
import os
import re
from subprocess import call
from time import sleep

from Bio import SeqIO
from Bio.ExPASy import Prosite,Prodoc

# FUNCTION DEFINITION
def check():
	"""Function to check if prosite. dat and prosite.doc are on the program folder."""

	# Search for the files on the program folder.
	errors = []
	for file in ["prosite.dat","prosite.doc"]:
		if not os.path.exists(file):
			errors.append(file)
	return errors


def patfinder(fasta,output):
	"""fasta: name of the FASTA file with the sequences to search domains on (str).
	output: name of the file where the results will be stored (str).

	This function searchs for domains of the prosite database on the given sequences,
	and shows some information about them. 
	"""

	for x in range(0,4):
		print(">> Starting Prosite pattern search "+"."*x, end="\r")
		sleep(0.2)
	print()

	# Count number of proteins to analyze.
	nseqs = 0
	with open(fasta,"r") as file:
		for line in file:
			if line[0] == ">":
				nseqs += 1
	
	print("\n> %d proteins to analyze." % nseqs)
	print("\nThe results will be shown sequence by sequence. Each query sequence (name ended by _QUERYSEQ) will be followed by its blastp hits.")
	print("On the results folder, you'll find also a txt file with all the results, including also the exact sequence of the domains on the proteins and their position.")
	print("You'll be given also the option to see further information of the found domains of your choice.")
	input("\nPRESS ENTER TO CONTINUE\n")

	initial = [".","-","<",">","x","X","{","}","(",")"]
	final = ["","","^","$",".",".","[^","]","{","}"]

	out_file = open(output,"w")

	j = 0
	with open(fasta,"r") as seqs_handle:
		# Parse prosite.dat
		for seq_record in SeqIO.parse(seqs_handle, "fasta"):
			call('clear')
			j += 1
			seq,seq_id = seq_record.seq,seq_record.id
			print("\n---------------------------------------")
			print(">> (%d of %d) Prosite domains on sequence %s" % (j,nseqs,seq_id))
			print("-----------------------------------------")
			out_file.write("\n\n>> Prosite domains on sequence %s" % seq_id)
			with open("prosite.dat","r") as handle:
				pat_records = Prosite.parse(handle)
				total = 0
				results = []
				for record in pat_records:
					
					# Some patterns are empty. If not, convert them to regular expresions.
					if record.pattern != "":
						pattern = record.pattern
						for i in range(0,len(initial)):
							pattern = pattern.replace(initial[i],final[i])

						# Search domains.
						matches = re.finditer(pattern, str(seq))
						hit = False 
						domains,pos = [],[]
						for m in matches:
							domains.append(m.group())
							pos.append(m.start())
							hit = True

						# Show found domains.
						if hit == True:
							total += 1
							print("\n> Found %d hits for domain %s." % (len(domains),record.name))
							out_file.write("\n> Found %d hits for domain %s:\n" % (len(domains),record.name))
							out_file.write("Pos\tHit sequence\n")
							out_file.write("---\t------------\n")
							for i in range(0,len(domains)):
								out_file.write("%s\t%s\n" % (pos[i],domains[i]))

							results.append(record.accession)
							print("Domain accesion id: %s" % record.accession)
							print("Description: %s" % record.description)
							print("Pattern: %s" % record.pattern)
							out_file.write("Domain accesion id: %s\n" % record.accession)
							out_file.write("Description: %s\n" % record.description)
							out_file.write("Pattern: %s\n" % record.pattern)
							
				if total == 0:
					print("No domains found for this protein.")
					out_file.write("No domains found for this protein.")
					print("\n----------------------------")
					input("PRESS ENTER TO CONTINUE\n")
				else:
					print("\nTotal: %d different domains found.\n" % total)
					out_file.write("\n\nTotal: %d different domains found.\n\n" % total)
					print("\n---------------------------------------------------------")
					print("If you want further information of these domains, press Y.")
					print("Press ENTER or any other key to continue.")
					selection = input("> ")
					if selection.upper() == "Y":
						# Parse prosite.doc for further information.
						with open("prosite.doc","r") as doc_handle:
							doc_records = Prodoc.parse(doc_handle)
							for doc_record in doc_records:
								for x in results:
									if x in str(doc_record.prosite_refs):
										print("> %s domain." % x)
										print(doc_record.text)
							print("\n----------------------------")
							input("PRESS ENTER TO CONTINUE\n")

	out_file.close()	
