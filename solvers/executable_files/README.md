### Run

$instance: input cnf

$threads: number of threads

$time-limit: maximum running time



treengeling:
```
time timeout $time-limit ./treengeling $instance -t $threads
```

p-mcomsps:
```
time ./painless-mcomsps $instance -c=(($threads-2)) -shr-sleep=750000 -t=$time-limit
```

pakis(kissat):
```
./pakis(kissat) $instance $threads $time-limit
```

pakis(kissat-MAB):
```
./pakis(kissat-MAB) $instance $threads $time-limit
```

lstech_maple:
```
./lstech_maple $instance
```

lstech+LA:
```
./lstech+LA $instance $threads $time-limit
```

lstech+LA+s:
```
time ./lstech+LA+s -c=(($threads-1)) -shr-sleep=500000 -shr-lit=1500 $instance -t=$time-limit 
```

lstech+RS:
```
./lstech+RS $instance $threads $time-limit
```

lstech+RS+s:
```
time ./lstech+RS+s -c=(($threads-1)) -shr-sleep=500000 -shr-lit=1500 $instance -t=$time-limit -initshuffle
```

kissat-MAB:
```
./kissat-MAB $instance
```

kissat-MAB-SAT:
```
./kissat-MAB --sat $instance
```

kissat+LA:
```
./kissat+LA $instance $threads $time-limit
```

kissat+LA+s:
```
time ./kissat+LA+s -c=(($threads-1)) -shr-sleep=500000 -shr-lit=1500 $instance -t=$time-limit 
```

kissat+RS:
```
./kissat+RS $instance $threads $time-limit
```

kissat+RS+s:
```
time ./kissat+RS+s -c=(($threads-1)) -shr-sleep=500000 -shr-lit=1500 $instance -t=$time-limit -initshuffle
```

kissat+RS+d:
```
./kissat+RS+d $instance $threads $time-limit
```

kissat+RS+d+s:
```
time ./kissat+RS+s -c=(($threads-1)) -shr-sleep=500000 -shr-lit=1500 $instance -t=$time-limit -initshuffle
```






