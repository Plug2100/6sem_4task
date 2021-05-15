#include <mpi.h>
#include <iostream>
#include <complex>
#include <assert.h>
#include "omp.h"
#include <fstream>
#include <cmath>
#include "time.h"
#include <cmath>
#include "sys/time.h"
#include <stdio.h>
#include <stdlib.h>
using namespace std;
typedef complex<double> complexd;




complexd* generate_condition(unsigned long long seg_size, int rank, int size){
    double module = 0;
    unsigned int seed = time(NULL) + rank;
    complexd *VV;
    VV = (complexd*) malloc(sizeof(complexd) * seg_size);
    for (long long unsigned  i = 0; i < seg_size; i++){
        VV[i] = complexd(rand_r(&seed)%100 + 1.0, rand_r(&seed)%100 + 1.0);
        module += abs(VV[i] * VV[i]);
    }
    int rc;
    double new_m;
    MPI_Status stat;
    if(rank != 0){
        module += 1;
        rc = MPI_Send(&module, 1, MPI_DOUBLE, 0, 999, MPI_COMM_WORLD);
        MPI_Recv(&module, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &stat);
    }
    else{
        for(int i = 1; i < size; i++){
            MPI_Recv(&new_m, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 999, MPI_COMM_WORLD, &stat);
            module += new_m;
        }
        module = sqrt(module);
        for(int i = 1; i < size; i++){
            rc = MPI_Send(&module, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
        }
    }
    for (long long unsigned j = 0; j < seg_size; j++) {
        VV[j] /= module;
    }
    return VV;
}




int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	complexd *V;
	
	complexd *y;
    int rank;
    int size;
    unsigned n;
    n = atoi(argv[1]);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    unsigned long long index = 1LLU << n;
    unsigned long long seg_size = index / size;
	V = generate_condition(seg_size, rank, size);
	unsigned long long N = 1;
	for(int i = 0; i < n; i++){
		N *= 2;
	}
	complexd w;
	w = 2.718281828459;
	float hlp = 2*3.141592653589793/N;

	w = pow(w, complexd(0, hlp));

	y = (complexd*) malloc(sizeof(complexd) * N);
	for(int i = 0; i < N; i++){
		y[i] = 0;
	}
	double begin = MPI_Wtime();

	for(int i = 0; i < N; i++){
		for(int j = rank * seg_size; j < (rank+1)*seg_size; j++){
			y[i] += V[j - rank * seg_size] * pow(w, j*i);
		}
	}
	int rc;
	complexd buf;
	MPI_Status stat;

	for(int i = 0; i < size; i++){
		if(rank == i){
			for (int j = 0; j < seg_size; ++j){
				V[j] = y[rank * seg_size + j];
			}
			for(int j = 0; j < size - 1; j++){
				for(int f = 0; f < seg_size; f++){
					MPI_Recv(&buf, 2, MPI_COMPLEX, MPI_ANY_SOURCE, i * seg_size + f, MPI_COMM_WORLD, &stat);
					V[f] += buf;
				}
			}
		}
		else{
			for (int j = 0; j < seg_size; ++j){
				buf = complexd(y[i * seg_size + j]);
				rc = MPI_Send(&buf, 2, MPI_COMPLEX, i, i * seg_size + j, MPI_COMM_WORLD);
			}
		}
	}
	double stop = MPI_Wtime();
	if(rank == 0)
	{
	//	cout << stop - begin << endl;
	}	
	float ans = 0;
	float need;
	for(int i = 0; i < seg_size; i++){
		ans +=  abs((V[i] / sqrt(N)) * (V[i] / sqrt(N)));
	}
	if(rank == 0){
		for(int i = 0; i < size - 1; i++){
			MPI_Recv(&need, 1, MPI_FLOAT, MPI_ANY_SOURCE, 999999, MPI_COMM_WORLD, &stat);
			ans += need;
		}
		cout << ans << endl;
	}
	else{
		rc = MPI_Send(&ans, 1, MPI_FLOAT, 0, 999999, MPI_COMM_WORLD);
	}
	MPI_Finalize();
    delete[] V;
    delete[] y;
}