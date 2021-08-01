import glob
import os

### Parameters:
iteration = 3       # If we increase the number iteration, then we will get the good quality of PSSM.
evalue    = 0.001   # E-value
###

C = 0
for file in sorted(glob.glob('*.fasta')):
    # print(file)
    if C==10: break
    os.system('/home/rafsanjani/Dropbox/CompBio/PSI-BLAST/ncbi-blast-2.12.0+/bin/psiblast -query {} -db Pluto -num_iterations={} -evalue={} -out psiblastout.txt -out_ascii_pssm {}-.PSSM -save_each_pssm'.format(file, iteration, evalue, file))
    C = C+1
#end-for
