#include "fftw3.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>


// ausführen mit: ./fft_spektrum i.txt if.txt


void readData(const int N, double* inR, double& xmax, const char* fname);
void writeData(const fftw_complex* const f, const int N, const double L, const char* const fname);

using namespace std;


int main(int argc, char** argv) {

	if(argc != 3) {
		cout << "Usage: " << argv[0] << " input_file \t output_file" << endl;
		exit(1);
	}

	char *in_file  = argv[1];
	char *out_file = argv[2];

	const int N = 16384;
	double xmax;


	// Allocate memory
	fftw_complex* f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
	double* inR  = (double*) malloc(sizeof(double)*N);

  	// Create plan
	fftw_plan FW  = fftw_plan_dft_r2c_1d(N, inR, f, FFTW_ESTIMATE);

	// Read input data
	readData(N, inR, xmax,in_file);
	double L = xmax;
	cout << L << endl;

	// Calculate FFT
	fftw_execute(FW);

	// write output file
	writeData(f, N, L, out_file);

	// Clean up
	fftw_destroy_plan(FW);
	fftw_free(f);
	free(inR);

	return 0;
}


void readData(const int N, double* inR, double& xmax, const char* fname) {
  
	// einlesen auf der Konsole:
	/*ifstream in(...);
	 * in >> a; */

	
	ifstream input(fname);

	for(int i = 0; i < N; i++) {
		input >> xmax >> inR[i];
	}

	input.close();
}


void writeData(const fftw_complex* const f, const int N, const double L, const char* const fname) {

	ofstream out(fname);
	const double dk = 2*M_PI/L;
	double pk;

	for(int i = 0; i <= N/2; i++) {
		pk = sqrt(f[i][0]*f[i][0] + f[i][1]*f[i][1])/N;
		out << i*dk << "\t" << pk << "\t" << f[i][0] << "\t" << f[i][1] << endl;
	}

	out.close();
}



