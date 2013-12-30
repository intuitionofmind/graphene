/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  project_graphene
 *
 *        Version:  1.1
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
#include "globals.h"
ofstream file_log_main("main.log");

int main()
{
    static graphene K;
    static graphene X1, X2;
    static graphene Phi1, Phi2;
    time_t seed, start, end;
		double* M_inverse = new double[TOT*TOT];
    seed = time(NULL);
    srand(seed);
    K.generate(1);
    start = time(NULL);
    for(int i = 0; i<N_traj; i++)
    {
        int flag = trajectory(X1, X2, Phi1, Phi2, K);
        while(!flag)
        {
            flag = trajectory(X1, X2, Phi1, Phi2, K);
        }
        file_log_main<<i<<endl;
    }
    for(int i = 0; i<N_sample; i++)
    {
        for(int j = 0; j<N_interval; j++)
        {
            int flag = trajectory(X1, X2, Phi1, Phi2, K);
            while(!flag)
            {
                flag = trajectory(X1, X2, Phi1, Phi2, K);
            }
        }
//				int flag = spin_cor(K, M_inverse);
        int flag = twopoint(K);
        if(!flag)
        {
            file_log_main<<"Failure."<<endl;
            break;
        }
        file_log_main<<(N_traj+i)<<endl;
    }
		delete [] M_inverse;
    end = time(NULL);
    info(start, end);
    return 0;
}

double rand_gauss()
{
    double u,v;
    u = (double)rand()/RAND_MAX;
    v = (double)rand()/RAND_MAX;
    return sqrt(-2*log(u))*cos(2*Pi*v);
}

double dot(graphene X, graphene Y)
{
    double res = 0;
    for(int i = 0; i<TOT; i++)
    {
        res = res+X.take(i)*Y.take(i);
    }
    return res;
}

int search(int x, int y, int z, int t)
{
    int index;
    x = (x+N_R1)%N_R1;
    y = (y+N_R2)%N_R2;
    z = (z+N_R3)%N_R3;
    t = (t+N_tau)%N_tau;
    index = z+y*N_R3+x*N_R2*N_R3+t*N_R1*N_R2*N_R3;
    return index;
}

void info(time_t start, time_t end)
{
    ofstream file_Info("Info");
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

