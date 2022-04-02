### Command for competitors

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

lstech_maple:
```
./lstech_maple $instance
```

kissat-MAB:
```
./kissat $instance
```

kissat-MAB-SAT:
```
./kissat --sat $instance
```


## The command lines of how to build and run for other solvers can be found in the subdirectories of their own.



