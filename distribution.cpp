/*
 * =====================================================================================
 *
 *       Filename:  distribution.cpp
 *
 *    Description:  project_graphene
 *
 *        Version:  1.1
 *        Created:  04/09/2013 02:28:24 PM
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
graphene X1, X2;
graphene Phi1, Phi2;
ofstream file_log_dis("distribution.log");
ofstream file_ham("hamiltonian");

//couter
int all = 0;
int acc = 0;

int CG(graphene Phi, graphene K, graphene & X, graphene (*A)(graphene, graphene))
{
    graphene R, P, Ap;
    double a, b;
    double rnrn;
    int n = 0;
    X.set_zero();
    R = Phi;
    P = R;
    rnrn = dot(R, R);
    while((rnrn>Accuracy) && (n<N_iter))
    {
        Ap = A(P, K);
        a = rnrn/dot(P, Ap);
        X = X+P*a;
        R = R-Ap*a;
        b = dot(R, R)/rnrn;
        P = R+P*b;
        rnrn = dot(R, R);
        n++;
//       cout<<n<<" "<<rnrn<<endl;
    }
   file_log_dis<<n<<endl;
//    cout<<n<<" "<<endl;
    if(n == N_iter)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

graphene T(graphene X, graphene K)
{
    graphene Y;
    #pragma omp parallel for num_threads(N_threads)
    for(int l = 0; l<N_tau; l++)
    {
        for(int i = 0; i<N_R1; i++)
        {
            for(int j = 0; j<N_R2; j++)
            {
                for(int k = 0; k<N_R3; k++)
                {
                    Y.value(i, j, k, l) = -t_hat*(X.value(i, j, k+1, l)+X.value(i+(2*k-1), j, k+1, l)+X.value(i, j+(2*k-1), k+1, l))+(U_hat/2.0-mu_hat+U_hat_sq*K.value(i, j, k, l))*X.value(i, j, k, l);
                }
            }
        }
    }
    return Y;
}

graphene M(graphene X, graphene K)
{
    graphene Temp, Y;
    Temp = T(X, K);
    #pragma omp parallel for num_threads(N_threads)
    for(int l = 0; l<N_tau; l++)
    {
        for(int i = 0; i<N_R1; i++)
        {
            for(int j = 0; j<N_R2; j++)
            {
                for(int k = 0; k<N_R3; k++)
                {
                    int flag = 1;
                    double d, t1, t2 = 0.0;
                    if(l == (N_tau-1))
                    {
                        flag = -1;
                    }
                    d = flag*X.value(i, j, k, l+1)-X.value(i, j, k, l);
                    t1 = -t_hat*(X.value(i, j, k+1, l)+X.value(i+(2*k-1), j, k+1, l)+X.value(i, j+(2*k-1), k+1, l))+(U_hat/2.0-mu_hat+U_hat_sq*K.value(i, j, k, l))*X.value(i, j, k, l);
                    t2 = -t_hat*(Temp.value(i, j, k+1, l)+Temp.value(i+(2*k-1), j, k+1, l)+Temp.value(i, j+(2*k-1), k+1, l))+(U_hat/2.0-mu_hat+U_hat_sq*K.value(i, j, k, l))*Temp.value(i, j, k, l);
                    Y.value(i, j, k, l) = d+t1-0.5*t2;
                }
            }
        }
    }
    return Y;
}

graphene MT(graphene X, graphene K)
{
    graphene Temp, Y;
    Temp = T(X, K);
    #pragma omp parallel for num_threads(N_threads)
    for(int l = 0; l<N_tau; l++)
    {
        for(int i = 0; i<N_R1; i++)
        {
            for(int j = 0; j<N_R2; j++)
            {
                for(int k = 0; k<N_R3; k++)
                {
                    int flag = 1;
                    double d, t1, t2 = 0.0;
                    if(!l)
                    {
                        flag = -1;
                    }
                    d = flag*X.value(i, j, k, l-1)-X.value(i, j, k, l);
                    t1 = -t_hat*(X.value(i, j, k+1, l)+X.value(i+(2*k-1), j, k+1, l)+X.value(i, j+(2*k-1), k+1, l))+(U_hat/2.0-mu_hat+U_hat_sq*K.value(i, j, k, l))*X.value(i, j, k, l);
                    t2 = -t_hat*(Temp.value(i, j, k+1, l)+Temp.value(i+(2*k-1), j, k+1, l)+Temp.value(i, j+(2*k-1), k+1, l))+(U_hat/2.0-mu_hat+U_hat_sq*K.value(i, j, k, l))*Temp.value(i, j, k, l);
                    Y.value(i, j, k, l) = d+t1-0.5*t2;
                }
            }
        }
    }
    return Y;
}

graphene MMT(graphene X, graphene K)
{
    return M(MT(X, K), K);
}

int force(graphene & F, graphene K)
{
    graphene Y1, Y2;
    int flag;
    flag = CG(Phi1, K, X1, MMT);
    if(!flag)
    {
        return 0;
    }
    Y1 = MT(X1, K);
    Y1 = Y1-T(Y1, K);
    flag = CG(Phi1, K, X1, MMT);
    if(!flag)
    {
        return 0;
    }
    Y2 = MT(X2, K);
    Y2 = Y2-T(Y2, K);
    for(int i = 0; i<TOT; i++)
    {
        F.take(i) = -K.take(i)+2*U_hat_sq*(X1.take(i)*Y1.take(i)+X2.take(i)*Y2.take(i));
    }
    return 1;
}

double hamiltonian(graphene X, graphene K)
{
    double k, v, ham;
    int flag;
    flag = CG(Phi1, K, X1, MMT);
    if(!flag)
    {
        return 0;
    }
    flag = CG(Phi2, K, X2, MMT);
    if(!flag)
    {
        return 0;
    }
    k = dot(X, X)/2.0;
    v = dot(K, K)/2.0+dot(Phi1, X1)+dot(Phi2, X2);
    file_log_dis<<k<<"  + "<<v<<endl;
    ham = k+v;
		file_ham<<ham<<endl;
    return ham;
}

int trajectory(graphene & K)
{
    graphene Eta, P, Force;
    graphene K_t;
    double ham1, ham2;
    double r, s;
    double t = sqrt(0.5);
    int flag;
    Eta.generate(t);
    Phi1 = M(Eta, K);
    Eta.generate(t);
    Phi2 = M(Eta, K);
    X1.set_zero();
    X2.set_zero();
    P.generate(1.0);
    K_t = K;
    ham1 = hamiltonian(P, K);
    flag = force(Force, K);
    if(!flag)
    {
        K = K_t;
        all++;
        return 0;
    }
    P = P+Force*(dt/2.0);
    for(int i = 0; i<(N_md-1); i++)
    {
        K = K+P*dt;
        flag = force(Force, K);
        if(!flag)
        {
            K = K_t;
            all++;
            return 0;
        }
        P = P+Force*dt;
    }
    K = K+P*dt;
    flag = force(Force, K);
    if(!flag)
    {
        K = K_t;
        all++;
        return 0;
    }
    P = P+Force*(dt/2.0);
    ham2 = hamiltonian(P, K);
    r = (double)rand()/RAND_MAX;
    s = exp(ham1-ham2);
    file_log_dis<<ham1<<" "<<ham2<<" "<<r<<" "<<s;
//    cout<<ham1<<" "<<ham2<<" "<<r<<" "<<s;
    if(r<s)
    {
        all++;
        acc++;
        file_log_dis<<" | "<<acc<<"  "<<all<<endl;
//        cout<<" | "<<acc<<"  "<<all<<endl;
        return 1;
    }
    else
    {
        K = K_t;
        all++;
        file_log_dis<<" | "<<acc<<"  "<<all<<endl;
//        cout<<" | "<<acc<<"  "<<all<<endl;
        return 0;
    }
}
