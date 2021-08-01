
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

###### Notes:
1. I renamed the origial name of FASTA sequence as it is helpful for tracking the implementation.
2. I used sequential numerical order rather than the original sequence name.
3. Renaming the sequence is optional.
4. We can also use Colab for the splitting Multiple FASTA sequence into single sequences [[Update Implementation](https://github.com/mrzResearchArena/BLAST/blob/master/Split-FASTA-using-BioPython-Colab.ipynb)].


