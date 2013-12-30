//main.cpp
int main(int argc, char** argv);

//operation.cpp
int MPI_Assign(double* X, double* Y, int numPorcs, int myID);
int MPI_Add(double* X, double* Y, double* Z, int numPorcs, int myID);
int MPI_Subtract(double* X, double* Y, double* Z, int numPorcs, int myID);
int MPI_Multiply(double* X, double a, double* Y, int numPorcs, int myID);
double MPI_Dot(double* X, double* Y, int numProcs, int myID);
int MPI_Setunit(double* X, int numProcs, int myID);
int MPI_Setzero(double* X, int numProcs, int myID);
int MPI_Setone(double *X, int index, int numProcs, int myID);

//extra.cpp
int MPI_Setone(double *X, int index, int numProcs, int myID);
int MPI_Generate(double* X, double sigma, int numProcs, int myID);
int MPI_Output(double* X, int myID);
int MPI_Input(double* X);
int Search(int x, int y, int z, int t);
int MPI_Search(int x, int y, int z, int t, int myTau);
void Info(double start, double end, int acc, int all);

//distribution.cpp
int CG(double* Phi, double* X, double* K, int numProcs, int myID);
int T(double* X, double* Y, double* K, int numProcs, int myID);
int M(double* X, double* Y, double* K, int numProcs, int myID);
int MT(double* X, double* Y, double* K, int numProcs, int myID);
int MMT(double* X, double* Y, double* K, int numProcs, int myID);
int MMT(double* X, double* Y, double* K, int numProcs, int myID);
double Hamiltonian(double* X1, double* X2, double* Phi1, double* Phi2, double* P, double* K, int numProcs, int myID);
int Trajectory(double* X1, double* X2, double* Phi1, double* Phi2, double* K, int numProcs, int myID);

//twopoint.cpp
int FFT_Space(double* X, double* Y);
int IFFT_Space(double* X, double* Y);
int Twopoint(double* K, int numProcs, int myID);
