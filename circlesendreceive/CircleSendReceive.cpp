#include <iostream>
#include <mpi.h>


int main(int argc, char* argv[]) {
	int processNumber;
	int processRank;
	MPI_Status stat;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
	MPI_Comm_size(MPI_COMM_WORLD, &processNumber);

	if (processRank == 0) {
		int info = 0;
		MPI_Send(&info, 1, MPI_INT, processRank + 1, 1, MPI_COMM_WORLD);
		std::cout << "Send info " << info << ". Proccess " << processRank << std::endl;
		MPI_Recv(&info, 1, MPI_INT, processNumber - 1, 1, MPI_COMM_WORLD, &stat);
		std::cout << "Receive info " << info << ". Proccess " << processRank << std::endl;
	}
	else {
		int info;
		MPI_Recv(&info, 1, MPI_INT, processRank - 1, 1, MPI_COMM_WORLD, &stat);
		info += 1;
		if (processRank == processNumber - 1)
		{ 
			MPI_Send(&info, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); 
		}
		else 
		{ 
			MPI_Send(&info, 1, MPI_INT, processRank + 1, 1, MPI_COMM_WORLD); 
		}
	
		std::cout << "Receive info " << info - 1 << ", send info " << info << ". Proccess " << processRank << std::endl;

	}

	MPI_Finalize();
	return 0;
}