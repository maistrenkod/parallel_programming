#include <iostream>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv){
    int lengh = 0, rank = 0;
    double res_sum = 0, part_sum = 0;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &lengh);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status Stat = {};

    int N = atoi(argv[1]);
    int h = N / lengh;
    int l_st = 1 + rank * h;
    int r_st = h + l_st;

    for (int i = l_st; i < r_st; ++i){
        res_sum += 1.0 / i;
    }

    if (rank == 0){
        for (int i = 1; i < lengh; ++i){
            MPI_Recv(&part_sum, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &Stat);
            res_sum += part_sum;
        }
        
        std::cout << res_sum;
    }
    else{
        MPI_Send(&res_sum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}