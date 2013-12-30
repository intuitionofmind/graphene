/*
 * =====================================================================================
 *
 *       Filename:  acceleration.cpp
 *
 *    Description:  project_graphene
 *
 *        Version:  1.0
 *        Created:  04/09/2013 02:30:51 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Zheng (), intuitionofmind@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */

#include"head.h"
graphene X1_tilde, X2_tilde;
graphene Phi1_tilde, Phi2_tilde;
ofstream file_log_tilde("tilde.log");

complex g0(int k0)
{
    complex a;
    a.take_Re()=cos(Omega0*k0)-1.0+U_hat;
    a.take_Im()=-sin(Omega0*k0);
    return a;
}

double g1(int k1, int k2)
{
    return t_hat*(1+cos(k1*Omega1)+cos(k2*Omega2));
}

double g2(int k1, int k2)
{
    return t_hat*(sin(k1*Omega1)+sin(k2*Omega2));
}

double S_check()
{
    double s1=0, s2=0.1, e_hat;
    int n=0;
    while(abs(s1-s2)>Accuracy)
    {
        //    cout<<n<<endl
        s1=s2;
        for(int k1=0; k1<N_R1; k1++)
        {
            for(int k2=0; k2<N_R2; k2++)
            {
                for(int a=0; a<2; a++)
                {
                    e_hat=sqrt(g1(k1, k2)*g1(k1, k2)+g2(k1, k2)*g2(k1, k2));
                    s2+=pow((1-U_hat-s1+e_hat), (N_tau-1))/(1+pow((1-U_hat-s1+e_hat), (N_tau)));
                }
            }
        }
        s2=s2*(-U_hat/(N_R1*N_R2));
        n++;
    }
    return s1;
}

void permutation(complex *A, int N)
{
    int l=N;
    int *rev=new int[N];
    int *loc=new int[N];
    complex *B=new complex [N];
    for(int i=0; i<N; i++)
    {
        rev[i]=0;
        loc[i]=i;

    }
    while(l>1)
    {
        for(int i=0; i<N; i++)
        {
            if(loc[i]%2==1)
            {
                rev[i]=rev[i]+(l>>1);
            }
            loc[i]>>=1;
        }
        l>>=1;
    }
    for(int i=1; i<(N-1); i++)
    {
        B[rev[i]]=A[i];
    }
    for(int i=1; i<(N-1); i++)
    {
        A[i]=B[i];
    }
    delete [] rev;
    delete [] loc;
    delete [] B;
}

int fft(complex *A, int N)
{
    complex wm, w, t, u;
    double nor=sqrt(N);
    int m=1;
    permutation(A, N);
    for(int s=0; s<log(double(N))/log(2.0); s++)
    {
        m<<=1;
        wm.take_Re()=cos(2*Pi/double(m));
        wm.take_Im()=sin(2*Pi/double(m));
        for(int k=0; k<=(N-1); k+=m)
        {
            w.take_Re()=1;
            w.take_Im()=0;
            for(int j=0; j<=((m>>1)-1); j++)
            {
                t=(A[k+j+(m>>1)])*w;
                u=A[k+j];
                A[k+j]=u+t;
                A[k+j+(m>>1)]=u-t;
                w=wm*w;
            }
        }
    }
    for(int i=0; i<N; i++)
    {
        A[i]=A[i]/nor;
    }
    return 1;
}

int fft2d(complex A[N_R1][N_R2])
{
    complex *X=new complex[N_R1];
    complex *Y=new complex[N_R2];
    for(int j=0; j<N_R2; j++)
    {
        for(int i=0; i<N_R1; i++)
        {
            X[i]=A[i][j];
        }
        fft(X, N_R1);
        for(int i=0; i<N_R1; i++)
        {
            A[i][j]=X[i];
        }
    }
    for(int i=0; i<N_R1; i++)
    {
        for(int j=0; j<N_R2; j++)
        {
            Y[j]=A[i][j];
        }
        fft(Y, N_R2);
        for(int j=0; j<N_R2; j++)
        {
            A[i][j]=Y[j];
        }
    }
    delete [] X;
    delete [] Y;
    return 1;
}


int fft3d(complex A[N_R1][N_R2][N_tau])
{
    double t;
    complex m; //for anti-boundary condition
    complex *X=new complex[N_R1];
    complex *Y=new complex[N_R2];
    complex *Z=new complex[N_tau];
    for(int l=0; l<N_tau; l++)
    {
        for(int j=0; j<N_R2; j++)
        {
            for(int i=0; i<N_R1; i++)
            {
                X[i]=A[i][j][l];
            }
            fft(X, N_R1);
            for(int i=0; i<N_R1; i++)
            {
                A[i][j][l]=X[i];
            }
        }
    }
    for(int l=0; l<N_tau; l++)
    {
        for(int i=0; i<N_R1; i++)
        {
            for(int j=0; j<N_R2; j++)
            {
                Y[j]=A[i][j][l];
            }
            fft(Y, N_R2);
            for(int j=0; j<N_R2; j++)
            {
                A[i][j][l]=Y[j];
            }
        }
    }
    for(int j=0; j<N_R2; j++)
    {
        for(int i=0; i<N_R1; i++)
        {
            for(int l=0; l<N_tau; l++)
            {
                Z[l]=A[i][j][l];
            }
            fft(Z, N_tau);
            for(int l=0; l<N_tau; l++)
            {
                t=l*Phase;
                m.take_Re()=cos(t);
                m.take_Im()=sin(t);
                A[i][j][l]=Z[l]*m;
            }
        }
    }
    delete [] X;
    delete [] Y;
    delete [] Z;
    return 1;
}

int fft3dT(complex A[N_R1][N_R2][N_tau])
{
    double t;
    complex m; //for anti-boundary condition
    complex *X=new complex[N_R1];
    complex *Y=new complex[N_R2];
    complex *Z=new complex[N_tau];
    for(int l=0; l<N_tau; l++)
    {
        for(int j=0; j<N_R2; j++)
        {
            for(int i=0; i<N_R1; i++)
            {
                X[i]=A[i][j][l];
            }
            ifft(X, N_R1);
            for(int i=0; i<N_R1; i++)
            {
                A[i][j][l]=X[i];
            }
        }
    }
    for(int l=0; l<N_tau; l++)
    {
        for(int i=0; i<N_R1; i++)
        {
            for(int j=0; j<N_R2; j++)
            {
                Y[j]=A[i][j][l];
            }
            ifft(Y, N_R2);
            for(int j=0; j<N_R2; j++)
            {
                A[i][j][l]=Y[j];
            }
        }
    }
    for(int j=0; j<N_R2; j++)
    {
        for(int i=0; i<N_R1; i++)
        {
            for(int l=0; l<N_tau; l++)
            {
                Z[l]=A[i][j][l];
            }
            ifft(Z, N_tau);
            for(int l=0; l<N_tau; l++)
            {
                t=l*Phase;
                m.take_Re()=cos(t);
                m.take_Im()=-sin(t);
                A[i][j][l]=Z[l]*m;
            }
        }
    }
    delete [] X;
    delete [] Y;
    delete [] Z;
    return 1;
}

int ifft(complex *A, int N)
{
    complex wm, w, t, u;
    double nor=sqrt(N);
    int m=1;
    permutation(A, N);
    for(int s=0; s<log(double(N))/log(2.0); s++)
    {
        m<<=1;
        wm.take_Re()=cos(2*Pi/double(m));
        wm.take_Im()=-sin(2*Pi/double(m));
        for(int k=0; k<=(N-1); k+=m)
        {
            w.take_Re()=1;
            w.take_Im()=0;
            for(int j=0; j<=((m>>1)-1); j++)
            {
                t=(A[k+j+(m>>1)])*w;
                u=A[k+j];
                A[k+j]=u+t;
                A[k+j+(m>>1)]=u-t;
                w=wm*w;
            }
        }
    }
    for(int i=0; i<N; i++)
    {
        A[i]=A[i]/nor;
    }
    return 1;
}

int ifft2d(complex A[N_R1][N_R2])
{
    complex *X=new complex[N_R1];
    complex *Y=new complex[N_R2];
    for(int j=0; j<N_R2; j++)
    {
        for(int i=0; i<N_R1; i++)
        {
            X[i]=A[i][j];
        }
        ifft(X, N_R1);
        for(int i=0; i<N_R1; i++)
        {
            A[i][j]=X[i];
        }
    }
    for(int i=0; i<N_R1; i++)
    {
        for(int j=0; j<N_R2; j++)
        {
            Y[j]=A[i][j];
        }
        ifft(Y, N_R2);
        for(int j=0; j<N_R2; j++)
        {
            A[i][j]=Y[j];
        }
    }
    delete [] X;
    delete [] Y;
    return 1;
}

int ifft3d(complex A[N_R1][N_R2][N_tau])
{
    double t;
    complex m; //for anti-boundary condition
    complex *X=new complex[N_R1];
    complex *Y=new complex[N_R2];
    complex *Z=new complex[N_tau];
    for(int l=0; l<N_tau; l++)
    {
        for(int j=0; j<N_R2; j++)
        {
            for(int i=0; i<N_R1; i++)
            {
                X[i]=A[i][j][l];
            }
            ifft(X, N_R1);
            for(int i=0; i<N_R1; i++)
            {
                A[i][j][l]=X[i];
            }
        }
    }
    for(int l=0; l<N_tau; l++)
    {
        for(int i=0; i<N_R1; i++)
        {
            for(int j=0; j<N_R2; j++)
            {
                Y[j]=A[i][j][l];
            }
            ifft(Y, N_R2);
            for(int j=0; j<N_R2; j++)
            {
                A[i][j][l]=Y[j];
            }
        }
    }
    for(int j=0; j<N_R2; j++)
    {
        for(int i=0; i<N_R1; i++)
        {
            for(int l=0; l<N_tau; l++)
            {
                t=l*Phase;
                m.take_Re()=cos(t);
                m.take_Im()=-sin(t);
                Z[l]=A[i][j][l]*m;
            }
            ifft(Z, N_tau);
            for(int l=0; l<N_tau; l++)
            {
                A[i][j][l]=Z[l];
            }
        }
    }
    delete [] X;
    delete [] Y;
    delete [] Z;
    return 1;
}

int ifft3dT(complex A[N_R1][N_R2][N_tau])
{
    double t;
    complex m; //for anti-boundary condition
    complex *X=new complex[N_R1];
    complex *Y=new complex[N_R2];
    complex *Z=new complex[N_tau];
    for(int l=0; l<N_tau; l++)
    {
        for(int j=0; j<N_R2; j++)
        {
            for(int i=0; i<N_R1; i++)
            {
                X[i]=A[i][j][l];
            }
            fft(X, N_R1);
            for(int i=0; i<N_R1; i++)
            {
                A[i][j][l]=X[i];
            }
        }
    }
    for(int l=0; l<N_tau; l++)
    {
        for(int i=0; i<N_R1; i++)
        {
            for(int j=0; j<N_R2; j++)
            {
                Y[j]=A[i][j][l];
            }
            fft(Y, N_R2);
            for(int j=0; j<N_R2; j++)
            {
                A[i][j][l]=Y[j];
            }
        }
    }
    for(int j=0; j<N_R2; j++)
    {
        for(int i=0; i<N_R1; i++)
        {
            for(int l=0; l<N_tau; l++)
            {
                t=l*Phase;
                m.take_Re()=cos(t);
                m.take_Im()=sin(t);
                Z[l]=A[i][j][l]*m;
            }
            fft(Z, N_tau);
            for(int l=0; l<N_tau; l++)
            {
                A[i][j][l]=Z[l];
            }
        }
    }
    delete [] X;
    delete [] Y;
    delete [] Z;
    return 1;
}

graphene T0(graphene X)
{
    graphene Y;
    #pragma omp parallel for num_threads(N_threads)
    for(int l=0; l<N_tau; l++)
    {
        for(int i=0; i<N_R1; i++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int k=0; k<N_R3; k++)
                {
                    Y.value(i, j, k, l)=-t_hat*(X.value(i, j, k+1, l)+X.value(i+(2*k-1), j, k+1, l)+X.value(i, j+(2*k-1), k+1, l))+(U_hat/2.0-mu_hat)*X.value(i, j, k, l);
                }
            }
        }
    }
    return Y;
}

graphene T0T(graphene X)
{
    graphene Y;
    #pragma omp parallel for num_threads(N_threads)
    for(int l=0; l<N_tau; l++)
    {
        for(int i=0; i<N_R1; i++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int k=0; k<N_R3; k++)
                {
                    Y.value(i, j, k, l)=-t_hat*(X.value(i, j, k+1, l)+X.value(i+(2*k-1), j, k+1, l)+X.value(i, j+(2*k-1), k+1, l))+(U_hat/2.0-mu_hat)*X.value(i, j, k, l);
                }
            }
        }
    }
    return Y;
}

graphene M0I(graphene X)
{
    complex den_conj, g;
    complex Temp1[N_R1][N_R2][N_tau];
    complex Temp2[TOT];
    complex Temp3[TOT];
    graphene R;
    for(int k=0; k<N_R3; k++)
    {
        #pragma omp parallel for num_threads(N_threads)
        for(int l=0; l<N_tau; l++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int i=0; i<N_R1; i++)
                {
                    Temp1[i][j][l].take_Re()=X.value(i, j, k, l);
                    Temp1[i][j][l].take_Im()=0.0;
                    //cout<<Temp1[i][j][l].take_Re()<<endl;
                }
            }
        }
        fft3d(Temp1);
        #pragma omp parallel for num_threads(N_threads)
        for(int l=0; l<N_tau; l++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int i=0; i<N_R1; i++)
                {
                    int index=search(i, j, k, l);
                    Temp2[index]=Temp1[i][j][l];
                    //cout<<Temp2[index].take_Re()<<endl;
                }
            }
        }
    }
    for(int k=0; k<N_R3; k++)
    {
        for(int l=0; l<N_tau; l++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int i=0; i<N_R1; i++)
                {
                    int flag;
                    if(k)
                    {
                        flag=1;
                    }
                    else
                    {
                        flag=-1;
                    }
                    den_conj.take_Re()=g0(l).take_Re()*g0(l).take_Re()-g0(l).take_Im()*g0(l).take_Im()-g1(i, j)*g1(i, j)-g2(i, j)*g2(i, j);
                    den_conj.take_Im()=-2.0*g0(l).take_Re()*g0(l).take_Im();
                    int index1=search(i, j, k, l);
                    int index2=search(i, j, k+1, l);
                    g.take_Re()=g1(i, j);
                    g.take_Im()=flag*g2(i, j);
                    Temp3[index1]=g0(l)*Temp2[index1]+g*Temp2[index2];
                    Temp3[index1]=Temp3[index1]*den_conj/(den_conj.take_Re()*den_conj.take_Re()+den_conj.take_Im()*den_conj.take_Im());
                }
            }
        }
    }
    for(int k=0; k<N_R3; k++)
    {
        #pragma omp parallel for num_threads(N_threads)
        for(int l=0; l<N_tau; l++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int i=0; i<N_R1; i++)
                {
                    int index=search(i, j, k, l);
                    Temp1[i][j][l]=Temp3[index];
                }
            }
        }
        ifft3d(Temp1);
        #pragma omp parallel for num_threads(N_threads)
        for(int l=0; l<N_tau; l++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int i=0; i<N_R1; i++)
                {
                    int index=search(i, j, k, l);;
                    Temp3[index]=Temp1[i][j][l];
                }
            }
        }
    }
    for(int i=0; i<TOT; i++)
    {
        R.take(i)=Temp3[i].take_Re();
    }
    return R;
}

graphene M0IT(graphene X)
{
    complex den_conj, g;
    complex Temp1[N_R1][N_R2][N_tau];
    complex Temp2[TOT];
    complex Temp3[TOT];
    graphene R;

    for(int k=0; k<N_R3; k++)
    {
        #pragma omp parallel for num_threads(N_threads)
        for(int l=0; l<N_tau; l++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int i=0; i<N_R1; i++)
                {
                    Temp1[i][j][l].take_Re()=X.value(i, j, k, l);
                    Temp1[i][j][l].take_Im()=0.0;
                }
            }
        }
        fft3dT(Temp1);
        #pragma omp parallel for num_threads(N_threads)
        for(int l=0; l<N_tau; l++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int i=0; i<N_R1; i++)
                {
                    int index=search(i, j, k, l);
                    Temp2[index]=Temp1[i][j][l];
                }
            }
        }
    }
    for(int k=0; k<N_R3; k++)
    {
        #pragma omp parallel for num_threads(N_threads)
        for(int l=0; l<N_tau; l++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int i=0; i<N_R1; i++)
                {
                    int flag;
                    if(k)
                    {
                        flag=-1;
                    }
                    else
                    {
                        flag=1;
                    }
                    den_conj.take_Re()=g0(l).take_Re()*g0(l).take_Re()-g0(l).take_Im()*g0(l).take_Im()-g1(i, j)*g1(i, j)-g2(i, j)*g2(i, j);
                    den_conj.take_Im()=-2.0*g0(l).take_Re()*g0(l).take_Im();
                    int index1=search(i, j, k, l);
                    int index2=search(i, j, k+1, l);
                    g.take_Re()=g1(i, j);
                    g.take_Im()=flag*g2(i, j);
                    Temp3[index1]=g0(l)*Temp2[index1]+g*Temp2[index2];
                    Temp3[index1]=Temp3[index1]*den_conj/(den_conj.take_Re()*den_conj.take_Re()+den_conj.take_Im()*den_conj.take_Im());
                }
            }
        }
    }
    for(int k=0; k<N_R3; k++)
    {
        #pragma omp parallel for num_threads(N_threads)
        for(int l=0; l<N_tau; l++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int i=0; i<N_R1; i++)
                {
                    int index=search(i, j, k, l);
                    Temp1[i][j][l]=Temp3[index];
                }
            }
        }
        ifft3dT(Temp1);
        #pragma omp parallel for num_threads(N_threads)
        for(int l=0; l<N_tau; l++)
        {
            for(int j=0; j<N_R2; j++)
            {
                for(int i=0; i<N_R1; i++)
                {
                    int index=search(i, j, k, l);
                    Temp3[index]=Temp1[i][j][l];
                }
            }
        }
    }
    for(int i=0; i<TOT; i++)
    {
        R.take(i)=Temp3[i].take_Re();
    }
    return R;
}
graphene M_tilde(graphene X, graphene K)
{
    graphene R;
    graphene L;
    /*   graphene Ph1, Ph2, Ph3, Ph4, Ph5;
       Ph1=T0(T0(X));
       Ph2=K*X;
       Ph3=T0(Ph2);
       Ph4=K*(T0(X));
       Ph5=K*Ph2;
       R=X+M0I(Ph1)*0.5-M0I(Ph2)*U_hat_sq+M0I(Ph3)*U_hat_sq*0.5+M0I(Ph4)*U_hat_sq*0.5-M0I(Ph5)*U_hat*0.5; */
    graphene Ph;
    Ph=(K*U_hat_sq)*X;
    R=M0I(Ph)+X;
    return R;
}

graphene MT_tilde(graphene X, graphene K)
{
    graphene R;
    /*  graphene Ph1, Ph2, Ph3, Ph4, Ph5;
      Ph1=T0(T0(X));
      Ph2=K*X;
      Ph3=T0(Ph2);
      Ph4=K*(T0(X));
      Ph5=K*Ph2;
      R=X+M0IT(Ph1)*0.5-M0IT(Ph2)*U_hat_sq+M0IT(Ph3)*U_hat_sq*0.5+M0IT(Ph4)*U_hat_sq*0.5-M0IT(Ph5)*U_hat*0.5; */
    graphene Ph, L;
    Ph=(K*U_hat_sq)*X;
    R=M0IT(Ph)+X;
    return R;
}

graphene MMT_tilde(graphene X, graphene K)
{
    return M_tilde(MT_tilde(X, K), K);
}

int force_tilde(graphene & F, graphene K)
{
    graphene Y1, Y2;
    int flag;
    flag=CG(Phi1_tilde, K, X1_tilde, MMT_tilde);
    if(!flag)
    {
        return 0;
    }
    Y1=MT_tilde(X1_tilde, K);
   // Y1=Y1-T(Y1, K);
    flag=CG(Phi1_tilde, K, X1_tilde, MMT_tilde);
    if(!flag)
    {
        return 0;
    }
    Y2=MT_tilde(X2_tilde, K);
   // Y2=Y2-T(Y2, K);
    for(int i=0; i<TOT; i++)
    {
        F.take(i)=-K.take(i)+2*U_hat_sq*(X1_tilde.take(i)*Y1.take(i)+X2_tilde.take(i)*Y2.take(i));
    }
    return 1;
}

double hamiltonian_tilde(graphene X, graphene K)
{
    double ham;
    int flag;
    flag=CG(Phi1_tilde, K, X1_tilde, MMT_tilde);
    if(!flag)
    {
        return 0;
    }
    flag=CG(Phi2_tilde, K, X2_tilde, MMT_tilde);
    if(!flag)
    {
        return 0;
    }
    ham=dot(X, X)/2.0+dot(K, K)/2.0+dot(Phi1_tilde, X1)+dot(Phi2_tilde, X2);
    return ham;
}

int trajectory_tilde(graphene & K)
{
    graphene Eta, P, Force;
    graphene K_temp;
    double ham1, ham2;
    double r, s;
    double t=sqrt(0.5);
    int flag;
    Eta.generate(t);
    Phi1_tilde=M_tilde(Eta, K);
    Eta.generate(t);
    Phi2_tilde=M_tilde(Eta, K);
    P.generate(1);
    K_temp=K;
    ham1=hamiltonian_tilde(P, K);
    flag=force_tilde(Force, K);
    if(!flag)
    {
        K=K_temp;
        all++;
        return 0;
    }
    P=P+Force*(dt/2.0);
    for(int i=0; i<(N_md-1); i++)
    {
        K=K+P*dt;
        flag=force_tilde(Force, K);
        if(!flag)
        {
            K=K_temp;
            all++;
            return 0;
        }
        P=P+Force*dt;
    }
    K=K+P*dt;
    flag=force_tilde(Force, K);
    if(!flag)
    {
        K=K_temp;
        all++;
        return 0;
    }
    P=P+Force*(dt/2.0);
    ham2=hamiltonian_tilde(P, K);
    r=(double)rand()/RAND_MAX;
    s=exp(ham1-ham2);
    file_log_tilde<<ham1<<" "<<ham2<<" "<<r<<" "<<s;
    cout<<ham1<<" "<<ham2<<" "<<r<<" "<<s;
    if(r<s)
    {
        all++;
        acc++;
        file_log_tilde<<" | "<<acc<<"  "<<all<<endl;
        cout<<" | "<<acc<<"  "<<all<<endl;
        return 1;
    }
    else
    {
        K=K_temp;
        all++;
        file_log_tilde<<" | "<<acc<<"  "<<all<<endl;
        cout<<" | "<<acc<<"  "<<all<<endl;
        return 0;
    }
}

