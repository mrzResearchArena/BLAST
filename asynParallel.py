database = '/home/learning/mrzResearchArena/NR/nr'   # Please, set path where "nr" database directory is located.
PSSM = '/home/learning/mrzResearchArena/PSSM'        # Please, set path where PSSM directory is located.
core = 8       										 # numberThread=3000 --> Berkeley Lab # multiprocessing.cpu_count()

###


import multiprocessing
import time
import glob
import os

os.chdir(PSSM)

###


def runPSIBLAST(file):
    try:
        os.system('/home/learning/ncbi-blast-2.10.1+/bin/psiblast -query {} -db {} -out {}.out -num_iterations 3 -out_ascii_pssm {}.pssm -inclusion_ethresh 0.001 -comp_based_stats 0 -num_threads 1'.format(file, database, file, file))
    except:
        print('PSI-BLAST is error for the sequence {}!'.format(file))
        return '{}, is error.'.format(file)

    return '{}, is done.'.format(file)

###


begin   = time.time()
pool    = multiprocessing.Pool(processes=core)
results = [ pool.apply_async(runPSIBLAST, args=(file,)) for file in glob.glob('*.fasta') ] # for x in range(1, 10)

###


outputs = [result.get() for result in results]
end = time.time()

###


print(sorted(outputs))
print()
print('Time elapsed: {} seconds.'.format(end - begin))

###

