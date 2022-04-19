#include<stdio.h>
#include<iostream>
#include<fstream>
#include <windows.h>

int main(int argc, char* argv[]) {

	int N_x = 1000, M_t = 1000;
	int NORM = 10;

	std::ifstream fin;

	fin.open("output.txt");
	double data;
	HDC hdc = GetDC(GetConsoleWindow());

	for (int y = 0; y < M_t; ++y) {
		for (int x = 0; x < N_x; ++x) {
			fin >> data;

			SetPixel(hdc, y, x, RGB(255 * (data / NORM), 255 * (data / NORM), 255 * (data / NORM)));
		}
	}

	fin.close();

	while (true);
	return 0;
}