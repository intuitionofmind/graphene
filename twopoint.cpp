#include"head.h"
#include "globals.h"

int ffts(graphene & X, graphene & Y) //fft for SPACE(N_R1*N_R2) lattice of graphene
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
						int index = 0;
            for(int i = 0; i<N_R1; i++)
            {
                for(int j = 0; j<N_R2; j++)
                {
										index = i*N_R2+j;
                    in[index][0] = X.value(i, j, k, l);
                    in[index][1] = Y.value(i, j, k, l);
                }
            }
						fftw_execute(p);
            for(int i=0; i<N_R1; i++)
            {
                for(int j=0; j<N_R2; j++)
                {
                    index = i*N_R2+j;
                    X.value(i, j, k, l) = out[index][0]/sqrt(N_R1*N_R2); //pay attention to that the FFTW computes an unnormalized DFT!!!
                    Y.value(i, j, k, l) = out[index][1]/sqrt(N_R1*N_R2);
                }
            }
        }
    }

		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return 1;
}

int iffts(graphene & X, graphene & Y)
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
						int index = 0;
            for(int i = 0; i<N_R1; i++)
            {
                for(int j = 0; j<N_R2; j++)
                {
										index = i*N_R2+j;
                    in[index][0] = X.value(i, j, k, l);
                    in[index][1] = Y.value(i, j, k, l);
                }
            }
						fftw_execute(p);
            for(int i=0; i<N_R1; i++)
            {
                for(int j=0; j<N_R2; j++)
                {
                    index = i*N_R2+j;
                    X.value(i, j, k, l) = out[index][0]/sqrt(N_R1*N_R2); //pay attention to that the FFTW computes an unnormalized DFT!!!
                    Y.value(i, j, k, l) = out[index][1]/sqrt(N_R1*N_R2);
                }
            }
        }
    }

		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return 1;
}

int twopoint(graphene K)
{
		ofstream file_re("twopoint_re.dat", ios_base::app);
    ofstream file_im("twopoint_im.dat", ios_base::app);
    
    for(int k = 0; k<N_R3; k++)
    {
        for(int i = 0; i<N_R1; i++)
        {
            for(int j = 0; j<N_R2; j++)
            {
                int flag1, flag2;
                graphene X, Y, XX, YY;
                int index = search(i, j, k, 0);
                X.set_one(index);
                Y.set_zero();
                iffts(X, Y);
                flag1 = CG(X, K, XX, MMT);
                flag2 = CG(Y, K, YY, MMT);
                if(!(flag1 && flag2))
                {
                    return 0;
                }
                X = MT(XX, K);
                Y = MT(YY, K);
                ffts(X, Y);
                for(int l = 0; l<N_tau; l++)
                {
                    file_re<<setprecision(15)<<X.value(i, j, k, l)<<endl;
                    file_im<<setprecision(15)<<Y.value(i, j, k, l)<<endl;
                }
            }
        }
    }

		file_re.close();
		file_im.close();
    return 1;
}
