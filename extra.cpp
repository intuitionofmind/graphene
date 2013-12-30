/*
 * =====================================================================================
 *
 *       Filename:  extra.cpp
 *
 *    Description:  
 *
 *        Version:  2.2
 *        Created:  07/04/2013 11:10:27 PM
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

double RandomGauss()
{
        random_device rd;
        mt19937 generator(rd());
        normal_distribution<double> dis(0.0, 1.0);
        return dis(generator);
        }

int MPI_Generate(double* X, double sigma, int numProcs, int myID)
{
        int length = TOT/numProcs;
        double* myArray = new double[length];
        for(int i = 0; i<length; i++)
        {
            myArray[i] = sigma*RandomGauss();
            }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allgather(myArray, length, MPI_DOUBLE, X, length, MPI_DOUBLE, MPI_COMM_WORLD);
        delete [] myArray;
        return 1;
        }

int MPI_Output(double* X, int myID)
{
        if(ROOT == myID)
        {
            ofstream file_auxi("auxiliary", ios_base::out);
            for(int i = 0; i<TOT; i++)
            {
                file_auxi<<setprecision(15)<<X[i]<<endl;
                }
            file_auxi.close();
            }
        return 1;
        }

int MPI_Input(double *X)
{
        ifstream file_auxi("auxiliary");
        if(!file_auxi.is_open())
        {
            return 0;
            }
        for(int i = 0; i<TOT; i++)
        {
            file_auxi>>X[i];
            }
        file_auxi.close();
        return 1;
        }

int Search(int x, int y, int z, int t)
{
        int index;
        x = (x+N_R1)%N_R1;
        y = (y+N_R2)%N_R2;
        z = (z+N_R3)%N_R3;
        t = (t+N_tau)%N_tau;
        index = z+y*N_R3+x*N_R2*N_R3+t*N_R1*N_R2*N_R3;
        return index;
        }

int MPI_Search(int x, int y, int z, int t, int myTau)
{
        int index;
        x = (x+N_R1)%N_R1;
        y = (y+N_R2)%N_R2;
        z = (z+N_R3)%N_R3;
        t = (t+myTau)%myTau;
        index = z+y*N_R3+x*N_R2*N_R3+t*N_R1*N_R2*N_R3;
        return index;
        }

void Info(double start, double end, int acc, int all)
{
        ofstream file_Info("Information");
        file_Info<<"size: "<<N_R1<<"*"<<N_R2<<endl;
        file_Info<<"beta = "<<beta<<endl;
        file_Info<<"N_tau = "<<N_tau<<endl;
        file_Info<<"U = "<<U<<endl;
        file_Info<<"t = "<<t<<endl;
        file_Info<<"dt = "<<dt<<endl;
        file_Info<<"N_traj = "<<N_traj<<endl;
        file_Info<<"N_interval = "<<N_interval<<endl;
        file_Info<<"N_sample = "<<N_sample<<endl;
        file_Info<<"Acc = "<<acc<<"  "<<"All = "<<all<<" Ratio = "<<(double(acc)/double(all))*100<<"%"<<endl;
        file_Info<<"Time = "<<(end-start)/60.0<<" min"<<endl;
        file_Info.close();
        }
