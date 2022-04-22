#include<stdio.h>
#include<iostream>
#include<mpi.h>
#include<ctime>
#include<vector>
#include<fstream>

using std::chrono::milliseconds;

int N_x = 1000, M_t = 1000;

//линейное
double tau = 0.01;
double h = 0.01;

double func(double t, double x) {
	return 0;
}

double phi(double x) {
	return x;
}

double psi(double t) {
	return t;
}

template <typename T>
T** new_matrix(int height, int width) {
	T** matrix = new T * [height];
	for (int i = 0; i < height; ++i) {
		matrix[i] = new T[width];
	}
	return matrix;
}

template <typename T>
void delete_matrix(T** matrix, int height) {
	for (int i = 0; i < height; ++i) {
		delete[] matrix[i];
	}
	delete[] matrix;
}

//центральна€ трехточечна€ схема на нижней границе
double centre(double left, double right, double f) {
	return 0.5 * (left + right) + tau * (f + (left - right) / (2 * h));
}

//уголок на правой границе
double angle(double left, double central, double f) {
	return central + tau * (f + (left - central) / h);
}

//крест
double cross(double left, double right, double bottom, double f) {
	return bottom + 2 * tau * (f + (left - right) / (2 * h));
}

void send(double** value, MPI_Status& status, int ProcNum, int ProcRank, int N_PerProc, int n) {
	double send_right, send_left;

	if (ProcRank == 0) {
		send_right = value[n][N_PerProc - 1];
		MPI_Send(&send_right, sizeof(send_right), MPI_BYTE, 1, 1, MPI_COMM_WORLD);
	}
	else if (ProcRank == ProcNum - 1) {
		send_left = value[n][0];
		MPI_Send(&send_left, sizeof(send_left), MPI_BYTE, ProcRank - 1, 1, MPI_COMM_WORLD);
	}
	else {
		send_right = value[n][N_PerProc - 1];
		MPI_Send(&send_right, sizeof(send_right), MPI_BYTE, ProcRank + 1, 1, MPI_COMM_WORLD);
		send_left = value[n][0];
		MPI_Send(&send_left, sizeof(send_left), MPI_BYTE, ProcRank - 1, 1, MPI_COMM_WORLD);
	}
}

std::pair<double, double> recive(double** value, MPI_Status& status, int ProcNum, int ProcRank, int N_PerProc, int n) {
	double recv_right, recv_left;

	if (ProcRank == 0) {
		MPI_Recv(&recv_right, sizeof(recv_right), MPI_BYTE, 1, 1, MPI_COMM_WORLD, &status);
		recv_left = 0;
	}
	else if (ProcRank == ProcNum - 1) {
		MPI_Recv(&recv_left, sizeof(recv_left), MPI_BYTE, ProcRank - 1, 1, MPI_COMM_WORLD, &status);
		recv_right = 0;
	}
	else {
		MPI_Recv(&recv_right, sizeof(recv_right), MPI_BYTE, ProcRank + 1, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&recv_left, sizeof(recv_left), MPI_BYTE, ProcRank - 1, 1, MPI_COMM_WORLD, &status);
	}

	return std::pair<double, double>(recv_left, recv_right);
}

const void curTime() {
	auto millisec_since_epoch = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
	//time_t t;
	//time(&t);
	//std::cout << "Time: " << ctime(&t);
	std::cout << "Time: " << millisec_since_epoch << endl;
}

int main(int argc, char* argv[]) {

	int ProcNum, ProcRank;
	MPI_Status stat;

	int N_ForProc, N_Add;
	int X1;

	int t_d = 0, t_e = 0, t_s = 0;

	t_s = std::clock();

	MPI_Init(&argc, &argv);
	auto start_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

	if (ProcRank == 0) {
		curTime();
		std::cout << "From process " << ProcRank << ": start" << std::endl;
	}

	N_ForProc = N_x / ProcNum;
	N_Add = N_x % ProcNum;

	if (ProcRank < N_Add) {
		N_ForProc++;
		X1 = N_ForProc * ProcRank;
	}
	else {
		X1 = N_Add + N_ForProc * ProcRank;
	}

	double** value = new_matrix<double>(M_t, N_ForProc);

	for (int n = 0; n < N_ForProc; ++n) {
		value[0][n] = phi((n + X1) * h);
	}

	send(value, stat, ProcNum, ProcRank, N_ForProc, 0);

	for (int n = 1; n < N_ForProc - 1; ++n) {
		value[1][n] = centre(value[0][n - 1], value[0][n + 1], func(tau, (X1 + n) * h));
	}

	std::pair<double, double> value_recv = recive(value, stat, ProcNum, ProcRank, N_ForProc, 0);

	if (ProcRank == 0) {
		value[1][0] = psi(tau);
	}
	else {
		value[1][0] = centre(value_recv.first, value[0][1], func(tau, X1 * h));
	}

	if (ProcRank == ProcNum - 1) {
		value[1][N_ForProc - 1] = angle(value[0][N_ForProc - 2], value[0][N_ForProc - 1], func(tau, h * (X1 + N_ForProc - 1)));
	}
	else {
		value[1][N_ForProc - 1] = centre(value[0][N_ForProc - 2], value_recv.second, func(tau, h * (X1 + N_ForProc - 1)));
	}

	for (int k = 2; k < M_t; ++k) {

		send(value, stat, ProcNum, ProcRank, N_ForProc, k - 1);

		for (int n = 1; n < N_ForProc - 1; ++n) {
			value[k][n] = cross(value[k - 1][n - 1], value[k - 2][n], value[k - 1][n + 1], func(k * tau, h * (X1 + n)));
		}

		std::pair<double, double> value_recv = recive(value, stat, ProcNum, ProcRank, N_ForProc, k - 1);

		if (ProcRank == 0) {
			value[k][0] = psi(k * tau);
		}
		else {
			value[k][0] = cross(value_recv.first, value[k - 2][0], value[k - 1][1], func(k * tau, X1));
		}
 
		if (ProcRank == ProcNum - 1) {
			value[k][N_ForProc - 1] = angle(value[k - 1][N_ForProc - 2], value[k - 1][N_ForProc - 1], func(k * tau, (X1 + N_ForProc - 1) * h));
		}
		else {
			value[k][N_ForProc - 1] = cross(value[k - 1][N_ForProc - 2], value[k - 2][N_ForProc - 1], value_recv.second, func(k * tau, (X1 + N_ForProc - 1) * h));
		}
	}


	if (ProcRank != 0) {
		double* rev_value = new double[N_ForProc * M_t];

		for (int iy = 0; iy < M_t; ++iy){
			for (int ix = 0; ix < N_ForProc; ++ix) {
				rev_value[iy * N_ForProc + ix] = value[iy][ix];
			}
		}

		MPI_Send(&N_ForProc, sizeof(N_ForProc), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&X1, sizeof(X1), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
		MPI_Send(rev_value, N_ForProc * M_t * sizeof(double), MPI_BYTE, 0, 1, MPI_COMM_WORLD);

		curTime();
		std::cout << "From process " << ProcRank << ": sent to process 0 succesfully " << std::endl;

		delete[] rev_value;
	}
	else {
		double** fin_value = new_matrix<double>(M_t, N_x);
		for (int iy = 0; iy < M_t; ++iy) {
			for (int ix = 0; ix < N_ForProc; ++ix) {
				fin_value[iy][ix] = value[iy][ix];
			}
		}

		for (int i = 1; i < ProcNum; ++i) {
			int rec_N_PerProcess;
			MPI_Recv(&rec_N_PerProcess, sizeof(rec_N_PerProcess), MPI_BYTE, i, 1, MPI_COMM_WORLD, &stat);

			int proc_FirstX;
			MPI_Recv(&proc_FirstX, sizeof(proc_FirstX), MPI_BYTE, i, 1, MPI_COMM_WORLD, &stat);

			double* rec_value = new double[rec_N_PerProcess * M_t];
			MPI_Recv(rec_value, rec_N_PerProcess * M_t * sizeof(double), MPI_BYTE, i, 1, MPI_COMM_WORLD, &stat);

			curTime();
			std::cout << "From process " << ProcRank << ": data recived from process  " << i << " succesfully" << std::endl;

			for (int iy = 0; iy < M_t; ++iy) {
				for (int ix = 0; ix < rec_N_PerProcess; ++ix) {
					fin_value[iy][ix + proc_FirstX] = rec_value[iy * rec_N_PerProcess + ix];
				}
			}

		}

		t_e = std::clock();
		t_d = t_e - t_s;

		curTime();
		std::cout << "From process " << ProcRank << ": completed. Time " << double(t_d) / CLOCKS_PER_SEC << " seconds" << std::endl;

		std::ofstream out;
		char file_name[] = "output.txt";
		out.open(file_name);
		for (int iy = 0; iy < M_t; ++iy) {
			for (int ix = 0; ix < N_x; ++ix) {
				out << fin_value[iy][ix] << " ";
			}
			out << std::endl;
		}
		out.close();
		
		curTime();
		std::cout << "From process " << ProcRank << ": data saved in file: " << file_name << std::endl;

		delete_matrix<double>(fin_value, M_t);
	}

	delete_matrix<double>(value, M_t);

	auto stop_time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
	std::cout << stop_time - start_time << std::endl;

	MPI_Finalize();
	return 0;
}