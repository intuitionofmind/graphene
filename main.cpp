/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  project_Graphene
 *
 *        Version:  2.2
 *        Created:  04/09/2013 02:26:50 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Zheng (), intuitionofmind@gmail.com
 *   Organization:  School of Physics, Peking University
 *
 * =====================================================================================
 */

#include"head.h"
#include"parameters.h"
#include"declaration.h"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int numProcs, myID;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myID);
    int flag;
    double* X1 = new double[TOT];
    double* X2 = new double[TOT];
    double* Phi1 = new double[TOT];
    double* Phi2 = new double[TOT];
    double* K = new double[TOT];
    double* M_Inverse = new double[TOT*TOT];
    int acc = 0;
    int all = 0;
    MPI_Generate(K, 1.0, numProcs, myID);
    //	MPI_Input(K);
    double start = MPI_Wtime();
    for(int i = 0; i<N_traj; i++)
    {
        flag = Trajectory(X1, X2, Phi1, Phi2, K, numProcs, myID);
        if(flag)
        {
            acc++;
            all++;
        }
        else
        {
            all++;
        }
        if(myID == ROOT)
        {
            ofstream file_Log("tracjectroy_num.log", ios_base::app);
            file_Log<<acc<<"  "<<all<<endl;
            file_Log.close();
        }
    }
    for(int i = 0; i<N_sample; i++)
    {
        for(int j = 0; j<N_interval; j++)
        {
            flag = Trajectory(X1, X2, Phi1, Phi2, K, numProcs, myID);
            if(flag)
            {
                acc++;
                all++;
            }
            else
            {
                all++;
            }
            if(myID == ROOT)
            {
                ofstream file_Log("tracjectroy_num.log", ios_base::app);
                file_Log<<acc<<"  "<<all<<endl;
                file_Log.close();
            }
        }
        flag = Twopoint(K, numProcs, myID);
    }
    delete [] X1;
    delete [] X2;
    delete [] Phi1;
    delete [] Phi2;
    delete [] M_Inverse;
    double end = MPI_Wtime();
    if(myID == ROOT)
    {
        Info(start, end, acc, all);
    }
    MPI_Output(K, myID);
    MPI_Finalize();
    return 1;
}
