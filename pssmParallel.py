
# In[46]:

numberThread=2500
import os
import glob
from multiprocessing import Process
os.chdir('/home/mrz/Desktop/Bk/PSSM/') #Please, set path where PSSM directory is located.

Files = os.listdir()

print(Files)


# In[47]:


def run(file):
    try:
        os.system('psiblast -query {} -db ../../nr/nr -out {}.out -num_iterations 3 -out_ascii_pssm {}.pssm -inclusion_ethresh 0.001 -num_threads {}'.format(file, file, file, numberThread))
    except:
        print('File Error is {}!'.format(file))
    
    
P = []

def beginParallel(): [p.start() for p in P]    
def endParallel():   [p.join() for p in P]    

def parallel():
    for file in glob.glob('*.fasta'): # Run only *.fasta
        P.append(Process(target=run, args=(file,)))
        
    
    beginParallel()
    endParallel()


#%%
import time

begin = time.time()
parallel()
print(time.time()-begin)

