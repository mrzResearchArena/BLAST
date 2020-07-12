File = '/home/mrz/MyDrive/Education/Bioinformatics/sequences.txt'

from Bio import SeqIO

C= 1
for record in SeqIO.parse(File, 'fasta'):
    openFile = open(str(C) + '.fasta', 'w')
    SeqIO.write(record, openFile, 'fasta')
    C += 1
#end-for

