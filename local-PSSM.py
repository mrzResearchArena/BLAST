import glob
import os

### Parameters:
iteration = 10       # If we increase the number iteration, then we will get the good quality of PSSM.
evalue    = 0.001    # E-value
###


for file in sorted(glob.glob('*.fasta')):
    os.system('/home/rafsanjani/Dropbox/CompBio/PSI-BLAST/ncbi-blast-2.12.0+/bin/psiblast -query {} -db Pluto -num_iterations={} -evalue={} -out psiblastout.txt -out_ascii_pssm {}.pssm'.format(file, iteration, evalue, file))
#end-for
