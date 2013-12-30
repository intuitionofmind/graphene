/*
 * =====================================================================================
 *
 *       Filename:  distribution.cpp
 *
 *    Description:  project_Graphene
 *
 *        Version:  2.2
 *        Created:  04/09/2013 02:28:24 PM
 *       Revision:  none
 *       Compiler:  gcc, mpicxx
 *
 *         Author:  Wei Zheng (), intuitionofmind@gmail.com
 *   Organization:  School of Physics, Peking University
 *
 * =====================================================================================
 */
#include"head.h"
#include"parameters.h"
#include"declaration.h"

int CG(double* Phi, double* X, double* K, int numProcs, int myID)
{
        double* R = new double[TOT];
        double* P = new double[TOT];
        double* Ap = new double[TOT];
        double* temBuf = new double[TOT];
        int n = 0;
        MPI_Setzero(X, numProcs, myID);
        MPI_Assign(Phi, R, numProcs, myID);
        MPI_Assign(R, P, numProcs, myID);
        double rnrn = MPI_Dot(R, R, numProcs, myID);
        while((rnrn>Accuracy) && (n<N_iter))
        {
            MMT(P, Ap, K, numProcs, myID);
            double a = rnrn/MPI_Dot(P, Ap, numProcs, myID);
            MPI_Multiply(P, a, temBuf, numProcs, myID);
            MPI_Add(X, temBuf, X, numProcs, myID);
            MPI_Multiply(Ap, a, temBuf, numProcs, myID);
            MPI_Subtract(R, temBuf, R, numProcs, myID);
            double b = MPI_Dot(R, R, numProcs, myID)/rnrn;
            MPI_Multiply(P, b, temBuf, numProcs, myID);
            MPI_Add(R, temBuf, P, numProcs, myID);
            rnrn = MPI_Dot(R, R, numProcs, myID);
            n++;
            }
        if(myID == ROOT)
        {
            ofstream file_log("evolution.log", ios_base::app);
            file_log<<n<<endl;
            file_log.close();
            }
        delete [] R;
        delete [] P;
        delete [] Ap;
        delete [] temBuf;
        if(n == N_iter)
        {
            return 0;
            }
        else
        {
            return 1;
            }
        }

int T(double* X, double* Y, double* K, int numProcs, int myID)
{
        int myTau = N_tau/numProcs;
        int length = TOT/numProcs;
        int myLag = myID*length;
        double* myArrayX = new double[length];
        double* myArrayK = new double[length];
        double* myArrayY = new double[length];
        for(int i = 0; i<length; i++)
        {
            myArrayX[i] = X[myLag+i];
            myArrayK[i] = K[myLag+i];
            }
        for(int l = 0; l<myTau; l++)
        {
            for(int i = 0; i<N_R1; i++)
            {
                for(int j = 0; j<N_R2; j++)
                {
                    for(int k = 0; k<N_R3; k++)
                    {
                        int index = MPI_Search(i, j, k, l, myTau);
                        myArrayY[index] = -t_hat*(myArrayX[MPI_Search(i, j, (k+1), l, myTau)]+myArrayX[MPI_Search(i+(2*k-1), j, k+1, l, myTau)]+myArrayX[MPI_Search(i, j+(2*k-1), k+1, l, myTau)])+(U_hat/2.0-mu_hat+U_hat_sq*myArrayK[index])*myArrayX[index];
                        }
                    }
                }
            }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allgather(myArrayY, length, MPI_DOUBLE, Y, length, MPI_DOUBLE, MPI_COMM_WORLD); //Gather and broadcast.
        delete [] myArrayX;
        delete [] myArrayK;
        delete [] myArrayY;
        return 1;
        }

int M(double* X, double* Y, double* K, int numProcs, int myID)
{
        MPI_Status status;
        int myTau = N_tau/numProcs;
        int length = TOT/numProcs;
        int myLag = myID*length;
        double* myArrayX = new double[length];
        double* myArrayK = new double[length];
        double* myArrayT = new double[length];
        double* myArrayY = new double[length];
        double* sendBuf = new double[SPACE];  //These buffers are for the boundary values between different processes.
        double* recvBuf = new double[SPACE]; 
        double* temBuf = new double[TOT];
        T(X, temBuf, K, numProcs, myID);
        for(int i = 0; i<length; i++)
        {
            myArrayX[i] = X[myLag+i];
            myArrayK[i] = K[myLag+i];
            myArrayT[i] = temBuf[myLag+i];
            }
        for(int i = 0; i<N_R1; i++)
        {
            for(int j = 0; j<N_R2; j++)
            {
                for(int k = 0; k<N_R3; k++)
                {
                    int index = MPI_Search(i, j, k, 0, myTau);
                    sendBuf[index] = myArrayX[index];
                    }
                }
            }
        int dest = (numProcs+myID-1)%numProcs;
        int source = (numProcs+myID+1)%numProcs;
        int sendTag = myID; //Messages are noted by the sending process. Deadlock may occur if you do not note them clearly.
        int recvTag = source;
        MPI_Sendrecv(sendBuf, SPACE, MPI_DOUBLE, dest, sendTag, recvBuf, SPACE, MPI_DOUBLE, source, recvTag, MPI_COMM_WORLD, &status);
        MPI_Barrier(MPI_COMM_WORLD);
        for(int l = 0; l<myTau; l++)
        {
            for(int i = 0; i<N_R1; i++)
            {
                for(int j = 0; j<N_R2; j++)
                {
                    for(int k = 0; k<N_R3; k++)
                    {
                        int index = MPI_Search(i, j, k, l, myTau);
                        double d, t1, t2 = 0.0;
                        if(l == (myTau-1))
                        {
                            if(myID == (numProcs-1))  //anti-boundary condition
                            {   
                                d = -recvBuf[MPI_Search(i, j, k, 0, myTau)]-myArrayX[index];
                                }
                            else
                            {  
                                d = recvBuf[MPI_Search(i, j, k, 0, myTau)]-myArrayX[index];
                                }
                            }
                        else
                        {
                            d = myArrayX[MPI_Search(i, j, k, l+1, myTau)]-myArrayX[index];
                            }
                        t1 = -t_hat*(myArrayX[MPI_Search(i, j, (k+1), l, myTau)]+myArrayX[MPI_Search(i+(2*k-1), j, k+1, l, myTau)]+myArrayX[MPI_Search(i, j+(2*k-1), k+1, l, myTau)])+(U_hat/2.0-mu_hat+U_hat_sq*myArrayK[index])*myArrayX[index];
                        t2 = -t_hat*(myArrayT[MPI_Search(i, j, (k+1), l, myTau)]+myArrayT[MPI_Search(i+(2*k-1), j, k+1, l, myTau)]+myArrayT[MPI_Search(i, j+(2*k-1), k+1, l, myTau)])+(U_hat/2.0-mu_hat+U_hat_sq*myArrayK[index])*myArrayT[index];
                        myArrayY[index] = d+t1-0.5*t2;
                        }
                    }
                }
            }
        MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(myArrayY, length, MPI_DOUBLE, Y, length, MPI_DOUBLE, MPI_COMM_WORLD);
  delete [] myArrayX;
  delete [] myArrayK;
  delete [] myArrayT;
  delete [] myArrayY;
  delete [] sendBuf;
  delete [] recvBuf;
  delete [] temBuf;
  return 1;
}

int MT(double* X, double* Y, double* K, int numProcs, int myID)
{
        MPI_Status status;
        int myTau = N_tau/numProcs;
        int length = TOT/numProcs;
        int myLag = myID*length;
        double* myArrayX = new double[length];
        double* myArrayK = new double[length];
        double* myArrayT = new double[length];
        double* myArrayY = new double[length];
        double* sendBuf = new double[SPACE];  //These buffers are for the boundary values between different nodes.
        double* recvBuf = new double[SPACE]; 
        double* temBuf = new double[TOT];
        T(X, temBuf, K, numProcs, myID);
        for(int i = 0; i<length; i++)
        {
            myArrayX[i] = X[myLag+i];
            myArrayK[i] = K[myLag+i];
            myArrayT[i] = temBuf[myLag+i];
            }
        for(int i = 0; i<N_R1; i++)
        {
            for(int j = 0; j<N_R2; j++)
            {
                for(int k = 0; k<N_R3; k++)
                {
                    int index1 = MPI_Search(i, j, k, 0, myTau);
                    int index2 = MPI_Search(i, j, k, myTau-1, myTau);
                    sendBuf[index1] = myArrayX[index2];
                    }
                }
            }
        int dest = (numProcs+myID+1)%numProcs;
        int source = (numProcs+myID-1)%numProcs;
        int sendTag = myID;
        int recvTag = source;
        MPI_Sendrecv(sendBuf, SPACE, MPI_DOUBLE, dest, sendTag, recvBuf, SPACE, MPI_DOUBLE, source, recvTag, MPI_COMM_WORLD, &status);
        MPI_Barrier(MPI_COMM_WORLD);
        for(int l = 0; l<myTau; l++)
        {
            for(int i = 0; i<N_R1; i++)
            {
                for(int j = 0; j<N_R2; j++)
                {
                    for(int k = 0; k<N_R3; k++)
                    {
                        int index = MPI_Search(i, j, k, l, myTau);
                        double d, t1, t2 = 0.0;
                        if(!l)
                        {
                            if(myID == ROOT) //anti-boundary condition
                            {
                                d = -recvBuf[MPI_Search(i, j, k, 0, myTau)]-myArrayX[index];
                                }
                            else
                            {
                                d = recvBuf[MPI_Search(i, j, k, 0, myTau)]-myArrayX[index];
                                }
                            }
                        else
                        {
                            d = myArrayX[MPI_Search(i, j, k, l-1, myTau)]-myArrayX[index];
                            }
                        t1 = -t_hat*(myArrayX[MPI_Search(i, j, (k+1), l, myTau)]+myArrayX[MPI_Search(i+(2*k-1), j, k+1, l, myTau)]+myArrayX[MPI_Search(i, j+(2*k-1), k+1, l, myTau)])+(U_hat/2.0-mu_hat+U_hat_sq*myArrayK[index])*myArrayX[index];
                        t2 = -t_hat*(myArrayT[MPI_Search(i, j, (k+1), l, myTau)]+myArrayT[MPI_Search(i+(2*k-1), j, k+1, l, myTau)]+myArrayT[MPI_Search(i, j+(2*k-1), k+1, l, myTau)])+(U_hat/2.0-mu_hat+U_hat_sq*myArrayK[index])*myArrayT[index];
                        myArrayY[index] = d+t1-0.5*t2;
                        }
                    }
                }
            }
        MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(myArrayY, length, MPI_DOUBLE, Y, length, MPI_DOUBLE, MPI_COMM_WORLD);
  delete [] myArrayX;
  delete [] myArrayK;
  delete [] myArrayT;
  delete [] myArrayY;
  delete [] sendBuf;
  delete [] recvBuf;
  delete [] temBuf;
  return 1;
}

int MMT(double* X, double* Y, double* K, int numProcs, int myID)
{
        double* temBuf = new double[TOT];
        MT(X, temBuf, K, numProcs, myID);
        M(temBuf, Y,  K, numProcs, myID);
        delete [] temBuf;
        return 1;
        }

int FermionForce(double* X1, double* X2, double* Phi1, double* Phi2, double* F, double* K, int numProcs, int myID)
{
        int length = TOT/numProcs;
        int myLag = myID*length;
        double* myArray = new double[length];
        double* Y1 = new double[TOT];
        double* Y2 = new double[TOT];
        double* temBuf = new double[TOT];
        int flag;
        flag = CG(Phi1, X1, K, numProcs, myID);
        if(!flag)
        {
            return 0;
            }
        MT(X1, Y1, K, numProcs, myID);
        T(Y1, temBuf, K, numProcs, myID);
        MPI_Subtract(Y1, temBuf, Y1, numProcs, myID);
        flag = CG(Phi2, X2, K, numProcs, myID);
        if(!flag)
        {
            return 0;
            }
        MT(X2, Y2, K, numProcs, myID);
        T(Y2, temBuf, K, numProcs, myID);
        MPI_Subtract(Y2, temBuf, Y2, numProcs, myID);
        for(int i = 0; i<length; i++)
        {
            int index = myLag+i;
            myArray[i] = -K[index]+2*U_hat_sq*(X1[index]*Y1[index]+X2[index]*Y2[index]);
            }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allgather(myArray, length, MPI_DOUBLE, F, length, MPI_DOUBLE, MPI_COMM_WORLD);
        delete [] myArray;
        delete [] Y1;
        delete [] Y2;
        delete [] temBuf;
        return 1;
        }

double Hamiltonian(double* X1, double* X2, double* Phi1, double* Phi2, double* P, double* K, int numProcs, int myID)
{
        double myHam = 0.0;
        double ham = 0.0;
        int numReduce = 1;
        int length = TOT/numProcs;
        int myLag = myID*length;
        int flag;
        flag = CG(Phi1, X1, K, numProcs, myID);
        if(!flag)
        {
            return 0;
            }
        flag = CG(Phi2, X2, K, numProcs, myID);
        if(!flag)
        {
            return 0;
            }
        for(int i = 0; i<length; i++)
        {
            int index = myLag+i;
            myHam = myHam+0.5*(P[index]*P[index]+K[index]*K[index])+Phi1[index]*X1[index]+Phi2[index]*X2[index];
            }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&myHam, &ham, numReduce, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return ham;
        }

int Trajectory(double* X1, double* X2, double* Phi1, double* Phi2, double* K, int numProcs, int myID)
{
        double* Eta = new double[TOT];
        double* P = new double[TOT];
        double* F = new double[TOT];
        double* temK = new double[TOT];
        double* temBuf = new double[TOT];
        double t = sqrt(0.5);
        int flag;
        MPI_Generate(Eta, t, numProcs, myID);
        M(Eta, Phi1, K, numProcs, myID);
        MPI_Generate(Eta, t, numProcs, myID);
        M(Eta, Phi2, K, numProcs, myID);
        MPI_Setzero(X1, numProcs, myID);
        MPI_Setzero(X2, numProcs, myID);
        MPI_Generate(P, 1.0, numProcs, myID);
        MPI_Assign(K, temK, numProcs, myID);
        double hamStart = Hamiltonian(X1, X2, Phi1, Phi2, P, K, numProcs, myID);
        flag = FermionForce(X1, X2, Phi1, Phi2, F, K, numProcs, myID);
        if(!flag)
        {
            MPI_Assign(temK, K, numProcs, myID);
            return 0;
            }
        MPI_Multiply(F, 0.5*dt, temBuf, numProcs, myID);
        MPI_Add(P, temBuf, P, numProcs, myID);
        for(int i = 0; i<(N_md-1); i++)
        {
            MPI_Multiply(P, dt, temBuf, numProcs, myID);
            MPI_Add(K, temBuf, K, numProcs, myID);
            flag = FermionForce(X1, X2, Phi1, Phi2, F, K, numProcs, myID);
            if(!flag)
            {
                MPI_Assign(temK, K, numProcs, myID);
                return 0;
                }
            MPI_Multiply(F, dt, temBuf, numProcs, myID);
            MPI_Add(P, temBuf, P, numProcs, myID);
            }
            MPI_Multiply(P, dt, temBuf, numProcs, myID);
            MPI_Add(K, temBuf, K, numProcs, myID);
            flag = FermionForce(X1, X2, Phi1, Phi2, F, K, numProcs, myID);
            if(!flag)
            {
                MPI_Assign(temK, K, numProcs, myID);
                return 0;
                }
            MPI_Multiply(F, 0.5*dt, temBuf, numProcs, myID);
            MPI_Add(P, temBuf, P, numProcs, myID);
            double hamEnd = Hamiltonian(X1, X2, Phi1, Phi2, P, K, numProcs, myID);
            double r = 0.0;
            if(myID == ROOT)
            {
                random_device rd;
                mt19937 generator(rd());
                uniform_real_distribution<double> dis(0.0, 1.0);
                r = dis(generator);
                }
            MPI_Bcast(&r, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
            double s = exp(hamStart-hamEnd);
            if(myID == ROOT)
            {
                ofstream file_log("evolution.log", ios_base::app);
                file_log<<"Hamiltoian: "<<hamStart<<" "<<hamEnd<<" "<<r<<" "<<s<<endl;
                file_log.close();
                }
            delete [] Eta;
            delete [] P;
            delete [] F;
            delete [] temK;
            delete [] temBuf;
            if(r<s)
            {
                return 1;
                }
           else
           {
                MPI_Assign(temK, K, numProcs, myID);
                return 0;
                }
           }
