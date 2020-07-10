## BLAST: Basic Local Alignment Search Tool

I describe the procedure for the PSSM generation from the FASTA sequences.


&nbsp;
&nbsp;

#### Step 1: Download the Non-redundant (NR) Proteins Database:
```console
user@machine:~$ wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz'
user@machine:~$ wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz.md5'

Note: The database has 39 segments, and the initial file size is ~100GB but ~450GB after the extraction. The segments and the file size changes frequently.
```

&nbsp;
&nbsp;

#### Step 2: Download the BLAST Tool [[Website](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)]:

```console
user@machine:~$ wget 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.1+-x64-linux.tar.gz'

COUTION: The download link may not remain same every time, it will change after the version upgrade. Please make sure your OS and BLAST version from the given website.
```

&nbsp;
&nbsp;

#### Step 3: Extract the Non-redundant (NR) Proteins Database:
```console
user@machine:~$ gunzip nr.*.tar.gz
```

```bash
for i in $(seq 0 1 38); do   ### If segment is n then, the loop goes upto n-1.
    if [ $i -lt 10 ]; then
        tar -xvf nr.0$i.tar
    else
        tar -xvf nr.$i.tar
    fi
done
```

&nbsp;
&nbsp;

#### Step 4: Generate PSSM
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
        os.system('psiblast -query {} -db {} -out {}.out -num_iterations 3 -out_ascii_pssm {}.pssm -inclusion_ethresh 0.001 -comp_based_stats 0 -num_threads 1'.format(file, database, file, file))
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
