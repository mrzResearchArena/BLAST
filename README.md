# BLAST: Basic Local Alignment Search Tool

I describe the procedure for the PSSM generation from the FASTA sequences. It is asynchronous parallel processing that can process up to n-sequence at a time. I spend a considerable amount of time on the PSSM generation purpose. It is definitely a hard and tedious procedure, but I make it easy so that other researchers can use it efficiently. People can use it for PSSM generation; unfortunately, I did not check the benchmark yet.

&nbsp;
&nbsp;


### Step 0: Graphical Representation: How can we generate PSSM using asynchronous parallel processing?

<img src="https://github.com/mrzResearchArena/BLAST/blob/master/asyn-PSSM.jpeg" class="center" title="asyn-PSSM" width="850" height="450" />

&nbsp;
&nbsp;

### Step 1: Download the BLAST Tool:

Please find the latest version of BLAST tool from the given website (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). Then download one of the files as per your operating system (OS) requirement. As I am a Linux OS user, that is why I downloaded "ncbi-blast-...-linux.tar.gz". Please don't worry about the version; it usually changes over time.

&nbsp;

```console
user@machine:~$ wget 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.1+-x64-linux.tar.gz' ### Fetch from website
user@machine:~$ tar -xvzf ncbi-blast-2.10.1+-x64-linux.tar.gz                                                           ### Extract the tool after the download
```

&nbsp;
&nbsp;


### Step 2: Download the Non-redundant (NR) Proteins Database:

We can download the `nr` database from official website (https://ftp.ncbi.nlm.nih.gov/blast/db/), and the downloading processes are given below.

&nbsp;

#### Option-1: Downloading Process:

```console
user@machine:~$ wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz'
user@machine:~$ wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz.md5'
```

#### Option-2: Downloading Process:
```console
user@machine:~$ /home/user/ncbi-blast-2.10.1+/bin/update_blastdb.pl --decompress nr [*]
```


###### Notes:
1. The database had 39 segments, and the initial file size was approximately 100GB; we would get around 450GB after extraction (Until 2020).
2. The database is now 54 segments (Last Update August 1, 2021).
3. The segments change frequently.
4. We can get the `update_blastdb.pl` file from BLAST tool.

&nbsp;
&nbsp;


### Step 3: Update the Non-redundant (NR) Proteins Database (Optional):
When the `nr` database will be old, no need to download (or upgrade) rather than update the previous one.

&nbsp;

#### Update Process:
```console
user@machine:~$ /home/user/ncbi-blast-2.10.1+/bin/update_blastdb.pl --decompress nr [*]
```

###### Notes:
1. We can get the `update_blastdb.pl` file from BLAST tool.
2. Plese run the script from the `nr` directory (or folder), otherwise it won't work.


&nbsp;
&nbsp;

### Step 4: Extract the Non-redundant (NR) Proteins Database:

#### Extract  `*.tar.gz` Files:
```bash
n=38   ### If the number of the segment is n, then we will use n-1.

for i in $(seq 0 1 $n); do
    if [ $i -lt 10 ]; then
        tar -xvzf nr.0$i.tar.gz
    else
        tar -xvzf nr.$i.tar.gz
    fi
done
```

###### Notes:
1. A question can arise why I used 38 in the loop? The answer is, I got the 39 segments in the `nr` directory (or folder).
2. Please make sure that how many segments you have, then update the value of `n`. We can find it from `nr.pal` in `nr` directory (or folder).
3. Plese run the script from the `nr` directory (or folder), otherwise it won't work.
4. We will find the update decompress procedure from given [URL](https://github.com/mrzResearchArena/BLAST/blob/master/decompress-NR.sh).

&nbsp;
&nbsp;

### Step 5: Split Multile FASTA File into Single FASTA Files:

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

### Step 6: Tracking the Original Sequence Name (Optional):

```python
File = '/home/rafsanjani/Downloads/TS88.fa'

from Bio import SeqIO

C= 1

print('Original-Sequence-Name, Renamed-Sequence, Corresponding-PSSM')
for record in SeqIO.parse(File, 'fasta'):
    print('{}, {}.fasta, {}.fasta.pssm'.format(record.id, C, C))
    C += 1
#end-for
```

&nbsp;

###### Notes:
1. As I rename the sequences, you can track the original sequences.
2. You will find the update procedure from given [URL](https://github.com/mrzResearchArena/BLAST/blob/master/renamedSequencesOrder.py).

&nbsp;
&nbsp;

### Step 7: Implementation/Generate PSSMs:
```python
###
database = '/home/learning/mrzResearchArena/NR/nr'   # Please, set path where "nr" database directory is located.
PSSM = '/home/learning/mrzResearchArena/PSSM'        # Please, set path where PSSM directory is located.
core = 8                                             # multiprocessing.cpu_count()
###

###
import multiprocessing
import time
import glob
import os
os.chdir(PSSM)
###

###
def runPSIBLAST(file):
    try:
        os.system('/home/learning/ncbi-blast-2.10.1+/bin/psiblast -query {} -db {} -out {}.out -num_iterations 3 -out_ascii_pssm {}.pssm -inclusion_ethresh 0.001 -comp_based_stats 0 -num_threads 1'.format(file, database, file, file))
    except:
        print('PSI-BLAST is error for the sequence {}!'.format(file))
        return '{}, is error.'.format(file)

    return '{}, is done.'.format(file)
#end-def
###

###
begin   = time.time()
pool    = multiprocessing.Pool(processes=core)
results = [ pool.apply_async(runPSIBLAST, args=(file,)) for file in glob.glob('*.fasta') ] # for x in range(1, 10)
###

###
outputs = [result.get() for result in results]
end = time.time()
###

###
print(sorted(outputs))
print()
print('Time elapsed: {} seconds.'.format(end - begin))
###
```

&nbsp;

###### Notes:
1. You will find the update procedure from given [URL](https://github.com/mrzResearchArena/BLAST/blob/master/asynParallel.py).

&nbsp;
&nbsp;


**Acknowledgement:** I would like to thank you to Professor [Iman Dehzangi](https://scholar.google.com/citations?user=RkamSRYAAAAJ&hl=en), who helped me initially for the PSSM generation.
