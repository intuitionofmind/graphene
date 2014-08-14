/*
 * =====================================================================================
 *
 *       Filename:  function_class.cpp
 *
 *    Description:  project_graphene
 *
 *        Version:  1.1
 *        Created:  04/09/2013 02:29:58 PM
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

void graphene::output()
{
    ofstream file_Kesi("Kesi");
    for(int i=0; i<TOT; i++)
    {
        file_Kesi<<setprecision(15)<<p[i]<<endl;
    }
    file_Kesi.close();
}

void graphene::sample()
{
    ofstream file_sample("sample", ios_base::app);
    for(int i=0; i<TOT; i++)
    {
        file_sample<<setprecision(15)<<p[i]<<endl;
    }
    file_sample.close();
}

int graphene::input()
{
    ifstream fin("kesi");
    if(!fin.is_open())
    {
        cerr<<"Failed to open the file!"<<endl;
        return 0;
    }
    for(int i=0; i<TOT; i++)
    {
        fin>>p[i];
    }
    fin.close();
    return 1;
}

void graphene::set_unit()
{
    for(int i=0; i<TOT; i++)
    {
        p[i]=1;
    }
}

void graphene::set_one(int j)
{
    for(int i=0; i<TOT; i++)
    {
        p[i]=0;
    }
    p[j]=1.0;
}

void graphene::set_zero()
{
    for(int i=0; i<TOT; i++)
    {
        p[i]=0;
    }
}

void graphene::generate(double sigma)
{
    for(int i=0; i<TOT; i++)
    {
        p[i]=sigma*rand_gauss();
    }
}

double & graphene::take(int i)
{
    return p[i];
}

double & graphene::value(int x, int y, int z, int t) //search matrix, return the matrix element value
{
    int index;
    x=(x+N_R1)%N_R1;
    y=(y+N_R2)%N_R2;
    z=(z+N_R3)%N_R3;
    t=(t+N_tau)%N_tau;
    index=z+y*N_R3+x*N_R2*N_R3+t*N_R1*N_R2*N_R3;
    return p[index];
}

graphene graphene::operator =(graphene X)
{
    for(int i=0; i<TOT; i++)
    {
        p[i]=X.p[i];
    }
    return *this;
}

graphene graphene::operator +(graphene X)
{
    graphene Y;
    for(int i=0; i<TOT; i++)
    {
        Y.p[i]=p[i]+X.p[i];
    }
    return Y;
}

graphene graphene::operator -(graphene X)
{
    graphene Y;
    for(int i=0; i<TOT; i++)
    {
        Y.p[i]=p[i]-X.p[i];
    }
    return Y;
}

graphene graphene::operator *(double x)
{
    graphene Y;
    for(int i=0; i<TOT; i++)
    {
        Y.p[i]=x*p[i];
    }
    return Y;
}

graphene graphene::operator *(graphene X)
{
    graphene Y;
    for(int i=0; i<TOT; i++)
    {
        Y.p[i]=(X.p[i])*p[i];
    }
    return Y;
}

graphene graphene::operator /(double x)
{
    graphene Y;
    for(int i=0; i<TOT; i++)
    {
        Y.p[i]=p[i]/x;
    }
    return Y;
}

//complex

void complex::set_zero()
{
    Re=0;
    Im=0;
}

double & complex::take_Re()
{
    return Re;
}

double & complex::take_Im()
{
    return Im;
}

void complex::generate()
{
    Re=(double)rand()/RAND_MAX;
    Im=0;//(double)rand()/RAND_MAX;
}

complex complex::operator =(complex x)
{
    Re=x.Re;
    Im=x.Im;
    return *this;
}

complex complex::operator +(complex x)
{
    complex y;
    y.Re=Re+x.Re;
    y.Im=Im+x.Im;
    return y;
}

complex complex::operator -(complex x)
{
    complex y;
    y.Re=Re-x.Re;
    y.Im=Im-x.Im;
    return y;
}

complex complex::operator *(double x)
{
    complex y;
    y.Re=Re*x;
    y.Im=Im*x;
    return y;
}

complex complex::operator *(complex x)
{
    complex y;
    y.Re=Re*x.take_Re()-Im*x.take_Im();
    y.Im=Re*x.take_Im()+Im*x.take_Re();
    return y;
}

complex complex::operator /(double x)
{
    complex y;
    y.Re=Re/x;
    y.Im=Im/x;
    return y;
}
