/*
 * =====================================================================================
 *
 *       Filename:  operation.cpp
 *
 *    Description:  
 *
 *        Version:  2.2
 *        Created:  07/05/2013 01:01:53 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Zheng (), intuitionofmind@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include"head.h"
#include"parameters.h"
#include"declaration.h"

int MPI_Assign(double* X, double* Y, int numProcs, int myID)
{
    int length = TOT/numProcs;
    int myLag = myID*length;
    double* myArray = new double[length];
    for(int i = 0; i<length; i++)
    {
        myArray[i] = X[myLag+i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(myArray, length, MPI_DOUBLE, Y, length, MPI_DOUBLE, MPI_COMM_WORLD);
    delete [] myArray;
    return 1;
}

int MPI_Add(double* X, double* Y, double* Z, int numProcs, int myID)
{
    int length = TOT/numProcs;
    int myLag = myID*length;
    double* myArray = new double[length];
    for(int i = 0; i<length; i++)
    {
        myArray[i] = Y[myLag+i]+X[myLag+i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(myArray, length, MPI_DOUBLE, Z, length, MPI_DOUBLE, MPI_COMM_WORLD);
    delete [] myArray;
    return 1;
}

int MPI_Subtract(double* X, double* Y, double* Z, int numProcs, int myID)
{
    int length = TOT/numProcs;
    int myLag = myID*length;
    double* myArray = new double[length];
    for(int i = 0; i<length; i++)
    {
        myArray[i] = X[myLag+i]-Y[myLag+i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(myArray, length, MPI_DOUBLE, Z, length, MPI_DOUBLE, MPI_COMM_WORLD);
    delete [] myArray;
    return 1;
}

int MPI_Multiply(double* X, double a, double* Y, int numProcs, int myID)
{
    int length = TOT/numProcs;
    int myLag = myID*length;
    double* myArray = new double[length];
    for(int i = 0; i<length; i++)
    {
        myArray[i] = a*(X[myLag+i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(myArray, length, MPI_DOUBLE, Y, length, MPI_DOUBLE, MPI_COMM_WORLD);
    delete [] myArray;
    return 1;
}

double MPI_Dot(double* X, double* Y, int numProcs, int myID)
{
    double myRes = 0.0;
    double res = 0.0;
    int numReduce = 1;
    int length = TOT/numProcs;
    int myLag = myID*length;
    for(int i = 0; i<length; i++)
    {
        myRes = myRes+(X[myLag+i])*(Y[myLag+i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&myRes, &res, numReduce, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return res;
}

int MPI_Setunit(double* X, int numProcs, int myID)
{
    int length = TOT/numProcs;
    double* myArray = new double[length];
    for(int i = 0; i<length; i++)
    {
        myArray[i] = 1.0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(myArray, length, MPI_DOUBLE, X, length, MPI_DOUBLE, MPI_COMM_WORLD);
    delete [] myArray;
    return 1;
}

int MPI_Setzero(double* X, int numProcs, int myID)
{
    int length = TOT/numProcs;
    double* myArray = new double[length];
    for(int i = 0; i<length; i++)
    {
        myArray[i] = 0.0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(myArray, length, MPI_DOUBLE, X, length, MPI_DOUBLE, MPI_COMM_WORLD);
    delete [] myArray;
    return 1;
}

int MPI_Setone(double *X, int index, int numProcs, int myID)
{
    MPI_Setzero(X, numProcs, myID);
    X[index] = 1.0;
    return 1;
}
