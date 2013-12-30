/*
 * =====================================================================================
 *
 *       Filename:  check.cpp
 *
 *    Description:  project_Graphene
 *
 *        Version:  2.2
 *        Created:  04/09/2013 02:34:07 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Zheng (), intuitionofmind@gmail.com
 *   Organization:  School of Physics, Peking University 
 *
 * =====================================================================================
 */

#include"head.h"
#include"globals.h"

double F_function()
{
        double energy_hat=0.0;
        double res=0.0;
        for(int k=0; k<N_R3; k++)
        {
            int flag=1;
            for(int j=0; j<N_R2; j++)
            {
                for(int i=0; i<N_R1; i++)
                {
                    if(!k)
                    {
                        flag=-1;
                        }
                    energy_hat=flag*sqrt(g1(i, j)*g1(i, j)+g2(i, j)*g2(i, j));
                    res+=pow((1-U_hat+energy_hat), (N_tau-1))/(1+pow((1-U_hat+energy_hat), (N_tau)));
                    }
                }
            }
        res=res/(N_R1*N_R2);
        return res;
        }

int deal_sample()
{
        ofstream file_check_first_order("check_first_order");
        ofstream file_log_check("check.log");
        double s[TOT*N_sample];
        double sum=0, ave=0, sta_dev=0;
        ifstream fin("sample");
        if(!fin.is_open())
        {
            file_log_check<<"Failed to open the sample file!"<<endl;
            return 0;
            }
        else
        {
            file_log_check<<"Succeeded in opening the sample file."<<endl;
            }
        for(int i=0; i<TOT*N_sample; i++)
        {
            fin>>s[i];
            sum+=s[i];
            }
        fin.close();
        ave=sum/(TOT*N_sample);
        if(!myID)
        {
            file_check_first_order<<"xi_pertubation= "<<-(U_hat_sq*F_function())<<endl;
            file_check_first_order<<"average= "<<ave<<endl;
            }
        sum=0;
        for(int i=0; i<TOT*N_sample; i++)
        {
            sum+=(s[i]-ave)*(s[i]-ave);
            }
        sum=sqrt(sum);
        sta_dev=sum/(TOT*N_sample-1);
        file_check_first_order<<"standard deviation= "<<sta_dev<<endl;
        file_check_first_order<<endl<<"Inforamtion:"<<endl;
        file_check_first_order<<"size: "<<N_R1<<"*"<<N_R2<<endl;
        file_check_first_order<<"beta= "<<beta<<endl;
        file_check_first_order<<"N_tau= "<<N_tau<<endl;
        file_check_first_order<<"U= "<<U<<endl;
        file_check_first_order<<"t= "<<t<<endl;
        file_check_first_order<<"dt= "<<dt<<endl;
        file_check_first_order<<"N_traj= "<<N_traj<<endl;
        file_check_first_order<<"N_interval= "<<N_interval<<endl;
        file_check_first_order<<"N_sample= "<<N_sample<<endl;
        return 1;
        }
