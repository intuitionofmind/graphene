/*
 * =====================================================================================
 *
 *       Filename:  spin_cor.cpp
 *
 *    Description:  project_graphene
 *
 *        Version:  1.0
 *        Created:  04/09/2013 02:34:48 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Zheng (), intuitionofmind@gmail.com
 *   Organization:  School of Physics, Peking University
 *
 * =====================================================================================
 */

#include"head.h"
int sam = 0;
ofstream file_spin_log("spin_cor.log");

int fft_spin_cor(double* X, double* Y) //fft for SPACE(N_R1*N_R2) lattice of graphene
{
		fftw_complex *in, *out;  
		fftw_plan p;
		in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N_R1*N_R2);
		out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N_R1*N_R2);
		p = fftw_plan_dft_2d(N_R1, N_R2, in, out, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

		if((in == NULL) || (out == NULL))
		{
        fftw_destroy_plan(p);
		    fftw_free(in);
		    fftw_free(out);
				return 0;
		}

    for(int i = 0; i<N_R1; i++)
    {
        for(int j = 0; j<N_R2; j++)
        {
						int index = i*N_R2+j;
            in[index][0] = X[index];
            in[index][1] = 0.0;
				}
		}
		fftw_execute(p);
    for(int i=0; i<N_R1; i++)
    {
        for(int j=0; j<N_R2; j++)
        {
             int index = i*N_R2+j;
             Y[index] = out[index][0]/sqrt(N_R1*N_R2); //pay attention to that the FFTW computes an unnormalized DFT!!!
				}
    }

		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return 1;
}

int spin_cor(graphene K, double* M_inverse)
{
    static graphene X, Y, Z;
		double* S_space_zz = new double[N_R1*N_R2];
		double* S_space_pm = new double[N_R1*N_R2];
		double* S_space_mp = new double[N_R1*N_R2];
		double* S_fourier_zz = new double[N_R1*N_R2];
		double* S_fourier_pm = new double[N_R1*N_R2];
		double* S_fourier_mp = new double[N_R1*N_R2];

    ofstream file_spin_space_zz("spin_space_zz.dat", ios_base::app);
    ofstream file_spin_fourier_zz("spin_fourier_zz.dat", ios_base::app);
    ofstream file_spin_fourier_tot("spin_fourier_tot.dat", ios_base::app);
		if(!(file_spin_space_zz && file_spin_fourier_zz && file_spin_fourier_tot))
		{
				file_spin_log<<"Failed to open the files!"<<endl;
				return 0;
		}
    file_spin_log<<" "<<acc<<" / "<<all<<" | "<<sam<<endl;

/*		for(int k = 0; k<N_R3; k++)
		{
				for(int i = 0; i<N_R1; i++)
				{
						for(int j = 0; j<N_R2; j++)
						{
								int index1 = search(i, j, k, 0)
						}
				}
		}*/

    for(int i = 0; i<TOT; i++)
    {
        X.set_one(i);
        int flag = CG(X, K, Y, MMT);
        if(!flag)
        {
            return 0;
        }
        Z = MT(Y, K);
        for(int j = 0; j<TOT; j++)
        {
						int index = i*TOT+j;
            M_inverse[index] = Z.take(j);
        }
    }

    for(int l = 0; l<N_tau; l++)
		{
        for(int i = 0; i<N_R1; i++)
        {
            for(int j = 0; j<N_R2; j++)
            {
						    int index = i*N_R2+j;
								int index_0tau = search(N_R1/2, N_R2/2, 0, l); //centre point
								int index_0taul = search(N_R1/2, N_R2/2, 0, l+1);
								int index_rtau = search(i, j, 0, l);
								int index_rtaul = search(i, j, 0, l+1);
                double s_zz = 2*(M_inverse[index_0taul*TOT+index_rtau])*(M_inverse[index_rtaul*TOT+index_0tau]);
								double s_pm = (M_inverse[index_rtaul*TOT+index_0tau])*(M_inverse[index_rtaul*TOT+index_0tau]);
								double s_mp = (M_inverse[index_0taul*TOT+index_rtau])*(M_inverse[index_0taul*TOT+index_rtau]);

								S_space_zz[index] = s_zz;
								S_space_pm[index] = s_pm;
								S_space_mp[index] = s_mp;

								file_spin_space_zz<<setprecision(15)<<s_zz<<endl;
            }
				}
				fft_spin_cor(S_space_zz, S_fourier_zz);
				fft_spin_cor(S_space_pm, S_fourier_pm);
				fft_spin_cor(S_space_mp, S_fourier_mp);

				for(int i = 0; i<N_R1; i++)
				{
						for(int j = 0; j<N_R1; j++)
						{
								int index = i*N_R1+j;
								file_spin_fourier_zz<<setprecision(15)<<S_fourier_zz[index]<<endl;
								file_spin_fourier_tot<<setprecision(15)<<(S_fourier_zz[index]+0.5*(S_fourier_pm[index]+S_fourier_mp[index]))<<endl;
						}
				}
		}

		delete [] S_space_zz;
		delete [] S_space_pm;
		delete [] S_space_mp;
		delete [] S_fourier_zz;
		delete [] S_fourier_pm;
		delete [] S_fourier_mp;
		file_spin_space_zz.close();
		file_spin_fourier_zz.close();
		file_spin_fourier_tot.close();

		sam++;
		return 1;
}
