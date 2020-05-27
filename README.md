# BLAST: Basic Local Alignment Search Tool

##### Step 1: Download the Non-redundant (NR) Proteins Database:
```console
user@machine:~$ wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz'
user@machine:~$ wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz.md5'

Note: The database has 39 segments, initially ~100GB, but ~450GB after extraction.
```

##### Step 2: Extract the Non-redundant (NR) Proteins Database:
```bash
for i in $(seq 0 1 38); do   ### If segment is n then, the loop goes upto n-1.
    if [ $i -lt 10 ]; then
        tar -xvf nr.0$i.tar
    else
        tar -xvf nr.$i.tar
    fi
done
```
