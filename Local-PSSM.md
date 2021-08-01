
# Local PSSM:


### Step 1: Create Local Database:


```console
user@machine:~$ /home/user/ncbi-blast-2.12.0+/bin/makeblastdb -in protein.fa -dbtype prot -out Pluto -parse_seqids
```

&nbsp;

After that, we will get below files.
```
Pluto.pdb
Pluto.phr
Pluto.pin
Pluto.pog
Pluto.pos
Pluto.pot
Pluto.psq
Pluto.ptf
Pluto.pto
```

&nbsp;
&nbsp;

### Step 2: Split Multile FASTA File into Single FASTA Files:

```python
File = '/home/user/Bioinformatics/multiSequences.fa'

from Bio import SeqIO # Install (If you don't have it.): pip install biopython

C= 1
for record in SeqIO.parse(File, 'fasta'):
    openFile = open(str(C) + '.fasta', 'w')
    SeqIO.write(record, openFile, 'fasta')
    C += 1
#end-for
```

&nbsp;
