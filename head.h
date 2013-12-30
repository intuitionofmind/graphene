#include <iostream>
#include <fstream>
#include<iomanip>
#include <memory.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <malloc.h>
#include <ctime>
#include <algorithm>
#include <omp.h>
#include <fftw3.h>  //FFTW

using namespace std;
const double Pi = 3.141592653589793238462643383279;

//parameters
const int N_R1 = 7;
const int N_R2 = 7;
const int N_R3 = 2;
const int N_tau = 100;
const int N_traj = 500;
const int N_md = 10;
const int N_sample = 30;
const int N_interval = 40;
const int N_iter = 10000;
const int N_threads = 2;
const int REAL = 0;
const int IMAG = 1;
const double beta = 10.0;
const double t = 1.0;
const double U = 0.1;
const double lam = 0.0;
const double mu = 0.0;
const double dt = 0.002;
const double Accuracy = 1e-14;
const double Omega0 = 2*Pi/double(N_tau);
const double Omega1 = 2*Pi/double(N_R1); //for fourier transform
const double Omega2 = 2*Pi/double(N_R2);
const double Phase = Pi/double(N_tau);

const int TOT = N_R1*N_R2*N_tau*N_R3;
const int SPACE = N_R1*N_R2;
const double a_t = beta/N_tau;
const double t_hat = a_t*t;
const double U_hat = a_t*U;
const double lam_hat = a_t*lam;
const double mu_hat = a_t*mu;
const double U_hat_sq = sqrt(U_hat);
const double Nor = sqrt(double(N_R1*N_R2));
const double NOR = sqrt(double(N_R1*N_R2*N_R3));

//class
class graphene //graphene
{
private:
    double p[TOT];

public:
    void output();
    void sample();
    int input();
    void set_unit();
    void set_one(int j);
    void set_zero();
    void generate(double sigma);
    double & take(int i);
    double & value(int x, int y, int z, int t);

    graphene operator =(graphene X); //operator overloading
    graphene operator +(graphene X);
    graphene operator -(graphene X);
    graphene operator *(double x);
    graphene operator *(graphene X);
    graphene operator /(double x);
};

class complex
{
private:
	double Re;
	double Im;

public:
	void set_zero();
	double & take_Re();
	double & take_Im();
	void generate();

	complex operator =(complex x);
	complex operator +(complex x);
	complex operator -(complex x);
	complex operator *(double x);
	complex operator *(complex x);
	complex operator /(double x);
};

//main.cpp
double rand_gauss();
double dot(graphene X, graphene Y);
int search(int x, int y, int z, int t);
void info(time_t start, time_t end);

//ditribution.cpp
int CG(graphene Phi, graphene K, graphene & X, graphene (*A)(graphene, graphene));
graphene T(graphene X, graphene K);
graphene M(graphene X, graphene K);
graphene MT(graphene X, graphene K);
graphene MMT(graphene X, graphene K);
int force(graphene & X1, graphene & X2, graphene & Phi1, graphene & Phi2, graphene & F, graphene K);
double hamiltonian(graphene & X1, graphene & X2, graphene & Phi1, graphene & Phi2, graphene X, graphene K);
int trajectory(graphene & X1, graphene & X2, graphene & Phi1, graphene & Phi2, graphene & K);

//acceleration
complex g0(int k0);
double g1(int k1, int k2);
double g2(int k1, int k2);
double S_check();
void permutation(complex *A, int N);
int fft(complex *A, int N);
int fft2d(complex A[N_R1][N_R2]);
int fft3d(complex A[N_R1][N_R2][N_tau]);
int fft3dT(complex A[N_R1][N_R2][N_tau]);
int ifft(complex *A, int N);
int ifft2d(complex A[N_R1][N_R2]);
int ifft3d(complex A[N_R1][N_R2][N_tau]);
int ifft3dT(complex A[N_R1][N_R2][N_tau]);
graphene T0(graphene X);
graphene M0I(graphene X);
graphene M0IT(graphene X);
graphene M_tilde(graphene X, graphene K);
graphene MT_tilde(graphene X, graphene K);
graphene MMT_tilde(graphene X, graphene K);
int force_tilde(graphene & F, graphene K);
double hamiltonian_tilde(graphene X, graphene K);
int trajectory_tilde(graphene & K);

//check
double F_function();
int deal_sample();

//energyband
int ffts(graphene & X, graphene & Y);
int iffts(graphene & X, graphene & Y);
int twopoint(graphene K);

//spin_cor
int fft_spin_cor(double* X, double* Y);
int spin_cor(graphene K, double* M_inverse);
