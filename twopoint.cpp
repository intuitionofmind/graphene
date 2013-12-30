/*
 * =====================================================================================
 *
 *       Filename:  Twopoint.cpp
 *
 *    Description:  project_Graphene
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
#include"parameters.h"
#include"declaration.h"

int FFT_Space(double* X, double* Y) //fft for SPACE(N_R1*N_R2) lattice of Graphene
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
    for(int k = 0; k<N_R3; k++)
    {
        for(int l = 0; l<N_tau; l++)
        {
            for(int i = 0; i<N_R1; i++)
            {
                for(int j = 0; j<N_R2; j++)
                {
										int indexFFT = i*N_R2+j;
										int indexData = Search(i, j, k, l);
                    in[indexFFT][0] = X[indexData];
                    in[indexFFT][1] = Y[indexData];
                }
            }
						fftw_execute(p);
            for(int i=0; i<N_R1; i++)
            {
                for(int j=0; j<N_R2; j++)
                {
                    int indexFFT = i*N_R2+j;
										int indexData = Search(i, j, k, l);
                    X[indexData] = out[indexFFT][0]/sqrt(N_R1*N_R2); //pay attention to that the FFTW computes an unnormalized DFT!!!
                    Y[indexData] = out[indexFFT][1]/sqrt(N_R1*N_R2);
                }
            }
        }
    }
		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return 1;
}

int IFFT_Space(double* X, double* Y)
{
		fftw_complex *in, *out;  
		fftw_plan p;
		in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N_R1*N_R2);
		out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N_R1*N_R2);
		p = fftw_plan_dft_2d(N_R1, N_R2, in, out, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
		if((in == NULL) || (out == NULL))
		{
        fftw_destroy_plan(p);
		    fftw_free(in);
		    fftw_free(out);
				return 0;
		}
    for(int k = 0; k<N_R3; k++)
    {
        for(int l = 0; l<N_tau; l++)
        {
            for(int i = 0; i<N_R1; i++)
            {
                for(int j = 0; j<N_R2; j++)
                {
					int indexFFT = i*N_R2+j;
					int indexData = Search(i, j, k, l);
                    in[indexFFT][0] = X[indexData];
                    in[indexFFT][1] = Y[indexData];
                }
            }
						fftw_execute(p);
            for(int i=0; i<N_R1; i++)
            {
                for(int j=0; j<N_R2; j++)
                {
  					int indexFFT = i*N_R2+j;
					int indexData = Search(i, j, k, l);
                    X[indexData] = out[indexFFT][0]/sqrt(N_R1*N_R2); //Pay attention to that the FFTW computes an unnormalized DFT.
                    Y[indexData] = out[indexFFT][1]/sqrt(N_R1*N_R2);
                }
            }
        }
    }
		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return 1;
}

int Twopoint(double* K, int numProcs, int myID)
{
		double* X = new double[TOT];
		double* Y = new double[TOT];
		double* XX = new double[TOT];
		double* YY = new double[TOT];
		for(int k = 0; k<N_R3; k++)
    {
        for(int i = 0; i<N_R1; i++)
        {
            for(int j = 0; j<N_R2; j++)
            {
	              int index = Search(i, j, k, 0);
								MPI_Setone(X, index, numProcs, myID);
								MPI_Setzero(Y, numProcs, myID);
                IFFT_Space(X, Y);
                int flag1 = CG(X, XX, K, numProcs, myID);
                int flag2 = CG(Y, YY, K, numProcs, myID);
                if(!(flag1 && flag2))
                {
                    return 0;
                }
                MT(XX, X, K, numProcs, myID);
                MT(YY, Y, K, numProcs, myID);
                FFT_Space(X, Y);
				if(!myID)
				{
					ofstream file_re("twopoint_re.dat", ios_base::app);
					ofstream file_im("twopoint_im.dat", ios_base::app);
                    for(int l = 0; l<N_tau; l++)
                    {
						int indexPrint = Search(i, j, k, l);
						file_re<<setprecision(15)<<X[indexPrint]<<endl;
						file_im<<setprecision(15)<<Y[indexPrint]<<endl;
                    }
						file_re.close();
						file_im.close();
				}
				else
				{
					for(int l = 0; l<N_tau; l++)
					{
						int indexPrint = Search(i, j, k, l);
						indexPrint++;
					}
				}
            }
        }
    }
		delete [] X;
		delete [] Y;
		delete [] XX;
		delete [] YY;
    return 1;
}
