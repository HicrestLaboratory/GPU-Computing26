# Laboratory 2

Cluster usage, ssh/scp connection, SLURM usage and runtime measurement.

## Exercise 1

log on the Baldo log-in node and run the “hostname” command. Save them on a file and copy the result back to the local machine.

```
kaijiefan@thinkdeeply:~$ ssh kaijie.fan@baldo.disi.unitn.it
kaijie.fan@baldo.disi.unitn.it's password: 
Last login: Tue Dec  2 11:16:15 2025 from 10.196.222.230

#####################################
# Login node for SSH connection of DISI Cluster 
#
# Wiki at https://wiki.cluster.disi.unitn.it
#
# Please note that edu01 node has been designed for educational purposes: 
# hyper threading is enabled and it should not be used for benchmarking
#
#####################################

[kaijie.fan@baldo ~]$ hostname
baldo
[kaijie.fan@baldo ~]$ echo “Hi from $(hostname)” > result.txt
[kaijie.fan@baldo ~]$ exit
logout
Connection to baldo.disi.unitn.it closed.
kaijiefan@thinkdeeply:~$ scp -i kaijie.fan@baldo.disi.unitn.it :~/result.txt ./local_path
```

## Exercise 2

Run "hostname" on the compute node edu01

```
[kaijie.fan@baldo ~]$ srun --nodes=1 --ntasks=1 --cpus-per-task=1 --gres=gpu:0 --partition=short --pty bash 
[kaijie.fan@gpu00 ~]$ hostname
gpu00
[kaijie.fan@gpu00 ~]$ exit
exit
[kaijie.fan@baldo ~]$

```


## Exercise 3

Write and compile a hello world program directly on the Baldo node and run it on the edu01 compute node.

```
kaijiefan@thinkdeeply:~$ ssh kaijie.fan@baldo.disi.unitn.it
[…]
[kaijie.fan@baldo ~]$ vi helloworld.c
[kaijie.fan@baldo ~]$ gcc helloworld.c -o helloworld
[kaijie.fan@baldo ~]$ srun --nodes=1 --ntasks=1 --cpus-per-task=1 --gres=gpu:0 --partition=short --pty bash 
[kaijie.fan@gpu00 ~]$ ./helloworld
Hello world!
[kaijie.fan@gpu00 ~]$ exit
[kaijie.fan@baldo ~]$

```

## Exercise 4

Compile and run the following code on the DISI cluster:

```
#include <stdio.h>
#include <cblas.h>

int main() {
    int m = 2, n = 2, k = 2;
    double A[4] = {1.0, 2.0, 3.0, 4.0};
    double B[4] = {5.0, 6.0, 7.0, 8.0};
    double C[4] = {0.0, 0.0, 0.0, 0.0};

    // Perform C = A * B using BLAS
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0, A, k, B, n, 0.0, C, n);

    // Print result
    printf("Result matrix C:\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", C[i * n + j]);
        }
        printf("\n");
    }
    return 0;
}
```

Since it require the OpenBLAS library, you need to load the module and link the library:

```
[kaijie.fan@baldo ~]$  module load OpenBLAS/
[kaijie.fan@baldo ~]$ gcc -o program program.c -L/opt/shares/openfoam/software/OpenBLAS/0.3.23-GCC-12.3.0/lib -I/opt/shares/openfoam/software/OpenBLAS/0.3.23-GCC-12.3.0/include/ -lopenblas
[kaijie.fan@baldo ~]$
```

Once you done, write a makefile to compile the source automatically.

Do you run on the log-in node or on the compute node? Write a sbatch script to run on the compute node ;)

## Exercise 5

Benchmark the gemm runtime
```
#include <stdio.h>
#include <cblas.h>

int main(void) {
    int m = 2, n = 2, k = 2;
    double A[4] = {1.0, 2.0, 3.0, 4.0};
    double B[4] = {5.0, 6.0, 7.0, 8.0};
    double C[4] = {0.0, 0.0, 0.0, 0.0};

     [ … timer start … ]
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0, A, k, B, n, 0.0, C, n);
     [ … timer stop … ]

     [ … compute mean and print out the result … ]
    return 0;
}

```
Discuss whether the runtime is reliable. 


