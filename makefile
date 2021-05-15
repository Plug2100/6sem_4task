test: main.cpp
	mpic++ -g -o main main.cpp -lm
	mpirun -np 1 main 2
	mpirun -np 2 main 2
	mpirun -np 4 main 2
	mpirun -np 1 main 4
	mpirun -np 2 main 4
	mpirun -np 4 main 4
	mpirun -np 1 main 8
	mpirun -np 2 main 8
	mpirun -np 4 main 8

