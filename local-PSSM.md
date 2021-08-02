
# Local PSSM:

### Step 0: Download `UniRef` Protein Database:

We can download `UniRef` (e.g., UniRef50, UniRef90, UniRef100) protein database from the given website (https://www.uniprot.org/downloads). After successfully downloading, we have to uncompress the file.

&nbsp;

```console
user@machine:~$ gunzip unirefX.fasta.gz (e.g., X = {50, 90, 100})
```


&nbsp;
&nbsp;

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
4. We will find the updated FASTA splitting procedure from given [URL](https://github.com/mrzResearchArena/BLAST/blob/master/splitFASTA.py).
5. We can also use Colab for the splitting Multiple FASTA sequence into single sequences [[Update Implementation](https://github.com/mrzResearchArena/BLAST/blob/master/Split-FASTA-using-BioPython-Colab.ipynb)].

&nbsp;
&nbsp;


### Step 3: Implementation/Generate Local PSSMs:

```python
import glob
import os

### Parameters:
iteration = 3       # If we increase the number iteration, then we will get the good quality of PSSM.
evalue    = 0.001   # E-value
###

C = 0
for file in sorted(glob.glob('*.fasta')):
    os.system('/home/user/ncbi-blast-2.12.0+/bin/psiblast -query {} -db Pluto -num_iterations={} -evalue={} -out psiblastout.txt -out_ascii_pssm {}.pssm'.format(file, iteration, evalue, file))
#end-for
```

&nbsp;

###### Notes:
1. We will find the updated FASTA splitting procedure from given [URL](https://github.com/mrzResearchArena/BLAST/blob/master/local-PSSM.py).
