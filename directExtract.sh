for i in $(seq 1 1 3); do
    if [ $i -lt 10 ]; then
        tar -xvzf nr.0$i.tar.gz
    else
        tar -xvzf nr.$i.tar.gz
    fi
done
