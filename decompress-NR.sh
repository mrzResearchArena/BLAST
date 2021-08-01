n=38   ### If the number of the segment is n, then we will use n-1.

for i in $(seq 0 1 $n); do
    if [ $i -lt 10 ]; then
        tar -xvzf nr.0$i.tar.gz
    else
        tar -xvzf nr.$i.tar.gz
    fi
done
