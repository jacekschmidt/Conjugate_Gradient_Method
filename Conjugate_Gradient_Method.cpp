#include <iostream>
#include "mvector.h"
#include "mmatrix.h"
#include <cmath>
#include "timing.h"

void PoiDisEqnMatrix(MMatrix &A, MVector &b, int n)
{
    MMatrix C(n,n,0);
    MVector d(n,1.0/((n+1)*(n+1)));
    C(0,0)=2;
    for(int i=1; i<C.Rows();i++)
    {
        C(i,i)=2;
        C(i,i-1)=-1;
        C(i-1,i)=-1;
    }
    A=C; b=d;
}

void PoiDisEqn2dMatrix(MMatrix &A, MVector &b, int n)
{
    MMatrix C(n*n,n*n,0);
    MVector d(n*n,1.0/((n+1)*(n+1)));
    for(int i=0; i<C.Rows();i++)
    {
        C(i,i)=4;
        if(i%n!=0)
        {
            C(i,i-1)=-1;
            C(i-1,i)=-1;
        }
        if(i>=n)
        {
            C(i-n,i)=-1;
            C(i,i-n)=-1;
        }
    }
    A=C; b=d;
}

void PoiDisEqnMatrixalt(MMatrix &A, MVector &b, int n, double alpha)
{
    MMatrix C(n,n,0);
    MVector d(n,2.5);
    C(0,0)=2+alpha;
    for(int i=1; i<C.Rows();i++)
    {
        C(i,i)=2+alpha;
        C(i,i-1)=-1;
        C(i-1,i)=-1;
    }
    A=C; b=d;
}

MVector ConjugateGradient(MMatrix A, MVector b, MVector &x, double tolerance)
{
    if (A.Cols()!=x.size() || A.Rows()!=b.size())
    {
        std::cout<<"incompatible arguments"<<std::endl;
        return {0};
    }
    double alpha, beta;
    MVector r=b-A*x, currentr;
    MVector p=r;
    for (int iter=0; iter<b.size(); iter++)
    {
        // ...calculate new values for x and r here...
        alpha=dot(r,r)/dot(p,A*p);
        x=x+alpha*p;
        currentr=r;
        r=r-alpha*(A*p);
        // check if solution is accurate enough
        if (r.L2Norm() < tolerance) 
        {
            std::cout<<iter<<std::endl;
            return x;
        }
        // ...calculate new conjugate vector p here...
        beta=dot(r,r)/dot(currentr,currentr);
        p=r+beta*p;
    }
    std::cout<<b.size()<<std::endl;
    return x;
}

MVector FastConjugateGradient(MMatrix A, MVector b, MVector &x, double tolerance)
{
    if (A.Cols()!=x.size() || A.Rows()!=b.size())
    {
        std::cout<<"incompatible arguments"<<std::endl;
        return {0};
    }
    double alpha, beta;
    MVector r=b-A*x, currentr;
    MVector p=r;
    for (int iter=0; iter<b.size(); iter++)
    {
        // ...calculate new values for x and r here...
        MVector q=A*p;
        alpha=dot(r,r)/dot(p,q);
        x=x+alpha*p;
        currentr=r;
        r=r-alpha*(q);
        // check if solution is accurate enough
        if (r.L2Norm() < tolerance) 
        {
            std::cout<<iter<<std::endl;
            return x;
        }
        // ...calculate new conjugate vector p here...
        beta=dot(r,r)/dot(currentr,currentr);
        p=r+beta*p;
    }
    std::cout<<b.size()<<std::endl;
    return x;
}

MVector Steepestdescent(MMatrix A, MVector b, MVector &x, double tolerance)
{
    if (A.Cols()!=x.size() || A.Rows()!=b.size())
    {
        std::cout<<"incompatible arguments"<<std::endl;
        return {0};
    }
    double alpha, beta;
    MVector r=b-A*x, currentr;
    MVector p=r;
    for (int iter=0; iter<100; iter++)
    {
        std::cout<<(2*p)/std::sqrt(p[0]*p[0]+p[1]*p[1])<<std::endl;
        // ...calculate new values for x and r here...
        alpha=dot(r,p)/dot(p,A*p);
        x=x+alpha*p;
        currentr=r;
        r=r-alpha*(A*p);
        // check if solution is accurate enough
        if (r.L2Norm() < tolerance) 
        {
            std::cout<<iter<<std::endl; 
            return x;
        }
        // ...calculate new conjugate vector p here...
        beta=dot(r,r)/dot(currentr,currentr);
        p=r;
    }
    return x;
}


int main()
{
    /*
     // ...initialise here...
    double alpha, beta, maxIterations = 25, tolerance = 1e-6;
    MMatrix A;
    MVector b, x(maxIterations,0);
    PoiDisEqnMatrix(A, b, maxIterations);
    MVector r=b-A*x, currentr;
    MVector p=r;

    for (int iter=0; iter<maxIterations; iter++)
    {
        // ...calculate new values for x and r here...
        alpha=dot(r,r)/dot(p,A*p);
        x=x+alpha*p;
        currentr=r;
        r=r-alpha*(A*p);
        // check if solution is accurate enough
        if (r.L2Norm() < tolerance) 
        {
            std::cout<<iter<<std::endl; 
            break;
        }
        // ...calculate new conjugate vector p here...
        beta=dot(r,r)/dot(currentr,currentr);
        p=r+beta*p;
    }
    std::cout<<r.L2Norm()<<std::endl;
    std::cout<<x<<std::endl;
    */

    /*
    int n=25;
    MMatrix A;
    MVector b, x(n,0);
    PoiDisEqnMatrix(A, b, n);
    std::cout<<ConjugateGradient(A,b,x,1e-6)<<std::endl;
    */

    /*
    int n=25;
    MMatrix A(2,2,0);
    A(0,0)=2; A(1,1)=2; A(0,1)=1; A(1,0)=1;
    MVector b={9,6}, x={1,1};
    std::cout<<ConjugateGradient(A,b,x,1e-6)<<std::endl;
    */

    /*
    MMatrix A(2,2,0);
    A(0,0)=2; A(1,1)=2; A(0,1)=1; A(1,0)=1;
    MVector b={9,6}, x={1,1};
    double alpha, beta;
    MVector r=b-A*x, currentr;
    MVector p={0,1};
        std::cout<<A*x<<std::endl;
        std::cout<<(2*p)/std::sqrt(p[0]*p[0]+p[1]*p[1])<<std::endl;
        std::cout<<x<<std::endl;
        // ...calculate new values for x and r here...
        alpha=dot(r,p)/dot(p,A*p);
        x=x+alpha*p;
        currentr=r;
        r=r-alpha*(A*p);
        // ...calculate new conjugate vector p here...
        beta=dot(r,r)/dot(currentr,currentr);
        p={2,-1};
        std::cout<<A*x<<std::endl;
        std::cout<<(2*p)/std::sqrt(p[0]*p[0]+p[1]*p[1])<<std::endl;
        std::cout<<x<<std::endl;
        // ...calculate new values for x and r here...
        alpha=dot(r,p)/dot(p,A*p);
        x=x+alpha*p;
        currentr=r;
        r=r-alpha*(A*p);
        // ...calculate new conjugate vector p here...
        beta=dot(r,r)/dot(currentr,currentr);
        */

        /*gradient descent
        {
        MMatrix A(2,2,0);
        A(0,0)=2; A(1,1)=2; A(0,1)=1; A(1,0)=1;
        MVector b={9,6}, x={1,1};        
        Steepestdescent(A, b, x, 1e-2);
        }
        */

        /*investigating x_i
        {
            int n=25;
            MMatrix A;
            MVector b, x(n,0);
            PoiDisEqnMatrix(A, b, n);
            std::cout<<ConjugateGradient(A,b,x,1e-6)<<std::endl;
        }
        */

        /*investigating iteration requirements
        {
            for(int n=1; n<101; n++)
            {
            MMatrix A;
            MVector b, x(n,0);
            PoiDisEqnMatrix(A, b, n);
            double start=Timer();
            ConjugateGradient(A,b,x,1e-6);
            double end=Timer();
            std::cout<<end-start<<std::endl;
            }
        }
        */

        /*effect of alpha
        {
            for(double alpha=0; alpha<2; alpha+=0.5)
            {
                for(int n=1; n<101; n++)
                {
                    MMatrix A;
                    MVector b, x(n,0);
                    PoiDisEqnMatrixalt(A, b, n,alpha);
                    ConjugateGradient(A,b,x,1e-6);
                }
            }
        }
        */

        /*POisson 2d fast
        {
            int n=50;
            for(int n=1;n<51;n++)
            {
                MMatrix A;
                MVector b,x(n*n,0);
                PoiDisEqn2dMatrix(A,b,n);
                double start=Timer();
                FastConjugateGradient(A,b,x,1e-6);
                double end=Timer();
                std::cout<<end-start<<std::endl;
            }
        }
        */
       
        {
            int n=60;
            MMatrix A;
            MVector b, x(n*n,0);
            PoiDisEqn2dMatrix(A, b, n);
            std::cout<<FastConjugateGradient(A,b,x,1e-6)<<std::endl;
        }
    return 0;
}