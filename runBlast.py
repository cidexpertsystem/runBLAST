# runBlast.py
# Cindy Reichel
# 2018-06-15
#
# usage: runBlastV1.py [-h] -i in-file -o out-file
#
# arguments:
#  -h, --help   show this help message and exit
#  -i in-file   input file in fasta format
#  -o out-file  output file for blast results

import argparse
import sqlite3
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# parse command line arguments using argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', metavar='in-file', required=True, type=argparse.FileType('rt'), help="input file in fasta format")
parser.add_argument('-o', metavar='out-file', required=True, type=argparse.FileType('wt'), help="output file for blast results")
results = parser.parse_args()
print("Opening input file '" + results.i.name + "'")

# Create SQLite database
print("Opening SQLite database '" + results.o.name + "'")

# Connecting to the database file
conn = sqlite3.connect(results.o.name)
c = conn.cursor()

# Creating a new SQLite tables
c.execute('CREATE TABLE IF NOT EXISTS query_sequences (query TEXT, query_letters TEXT)')
c.execute('CREATE TABLE IF NOT EXISTS blast_hits (query TEXT, title TEXT, length INT)')
c.execute('CREATE TABLE IF NOT EXISTS hsps (query TEXT, title TEXT, score FLOAT, bits FLOAT, expect FLOAT, coverage FLOAT, identity FLOAT)')

# run BLAST on fasta sequences from input file
out_handle = open("my_blast.xml", "w")
for seq_record in SeqIO.parse(results.i, "fasta"):
    # run blast on this sequence
    print("Running blast on " + seq_record.id)
    result_handle = NCBIWWW.qblast("blastn", "nt", seq_record.format("fasta"))

    # save results to XML file
    out_handle.write(result_handle.read())
    result_handle.close()
out_handle.close()

# Parse xml file into blast_record, save to SQLite database
blast_record = NCBIXML.parse(open('my_blast.xml', 'r'))
for record in blast_record:
        hit_count = 0
        c.execute('''INSERT INTO query_sequences(query, query_letters) \
        VALUES(?,?)''', (record.query, record.query_letters))
        for alignment in record.alignments:
            hit_count+=1
            c.execute('''INSERT INTO blast_hits(query, title, length) \
            VALUES(?,?,?)''', (record.query, alignment.title, alignment.length))
            for hsp in alignment.hsps:
                # calculate percent identity
                coverage = hsp.align_length / record.query_length
                identity = hsp.identities/ hsp.align_length
                c.execute('''INSERT INTO hsps(query, title, score, bits, expect, coverage, identity) \
                VALUES(?,?,?,?,?,?,?)''', (record.query, alignment.title, hsp.score, hsp.bits, hsp.expect, coverage, identity))
        print("Found " + str(hit_count) + " hits for " + record.query)

# Committing changes and closing the connection to the database file
conn.commit()

# Close SQLite connection
conn.close()
