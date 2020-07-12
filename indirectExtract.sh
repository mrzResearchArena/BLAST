#gunzip nr.*.tar.gz
#tar -xvf nr.00.tar
#%%bash
#%%bash
for i in $(seq 0 1 38); do
    if [ $i -lt 10 ]; then
        tar -xvf nr.0$i.tar
    else
        tar -xvf nr.$i.tar
    fi
done
