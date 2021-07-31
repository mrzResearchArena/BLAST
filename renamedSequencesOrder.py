File = '/home/rafsanjani/Downloads/TS88.fa'

from Bio import SeqIO

C= 1

print('Original-Sequence-Name, Renamed-Sequence, Corresponding-PSSM')
for record in SeqIO.parse(File, 'fasta'):
    print('{}, {}.fasta, {}.fasta.pssm'.format(record.id, C, C))
    C += 1
#end-for
