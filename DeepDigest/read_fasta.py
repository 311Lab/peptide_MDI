import os
import sys
import re


# read fasta file
def read_fasta(file, regular):
    # check the file path
    if not os.path.exists(file):
        print("Error: The input file does not exist. Please check again.")
        sys.exit(1)
    
    # record protein names and sequences
    with open(file) as fasta:
        records = fasta.read()
    
    # FASTA file must start with character '>'
    if re.search('>', records) == None:
        print("Error: The input file seems not in FASTA format!")
        sys.exit(1)
    # check the regular expression
    elif re.search(regular, records) == None:
        print("Error: Cannot parse the fasta file by the regular expression.")
        sys.exit(1)
    
    # fasta list
    records = re.split('(>.*?)\\n', records)[1:]
    length = len(records)
    fasta_list = []
    for ind in range(0, length, 2):
        name = re.split(regular, records[ind])[1]
        sequence = records[ind + 1].replace('\n', '')
        fasta_list += [[name, sequence]]
    
    return fasta_list