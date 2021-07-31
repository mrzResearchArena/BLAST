### Load Essential Files:
database = '/home/learning/mrzResearchArena/NR/nr'   # Please, set path where "nr" database directory is located.
PSSM = '/home/learning/mrzResearchArena/PSSM'        # Please, set path where PSSM directory is located.
###

### Parameters:
core = 8            # multiprocessing.cpu_count(). Please don't use the maximum core.
iteration = 3       # If we increase the number iteration, then we will get the good quality of PSSM.
evalue = 0.001      # E-value
###

### Load Essential Modules:
import multiprocessing
import time
import glob
import os

os.chdir(PSSM) # PSSM file will generate in PSSM directory.
###


### Generate PSSM:
def runPSIBLAST(file):
    try:
        os.system('/home/learning/ncbi-blast-2.10.1+/bin/psiblast -query {} -db {} -out {}.out -num_iterations {} -out_ascii_pssm {}.pssm -inclusion_ethresh {} -comp_based_stats 0 -num_threads 1'.format(file, database, file, iteration, file, evalue))
    except:
        print('PSI-BLAST is error for the sequence {}!'.format(file))
        return '{}, is error.'.format(file)

    return '{}, is done.'.format(file)
#end-def


### Run Procedure:
begin   = time.time()
pool    = multiprocessing.Pool(processes=core)
results = [ pool.apply_async(runPSIBLAST, args=(file,)) for file in glob.glob('*.fasta') ] # for x in range(1, 10)

outputs = [result.get() for result in results]
end = time.time()
###


### Verdict:
print(sorted(outputs))
print()
print('Time elapsed: {} seconds.'.format(end - begin))
###
