#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>


//Auxiliar function
double gammln(const double xx);
//////////////////////////////////////////////////////
template <class T>
double qgaus(T &func, const double a, const double b){
    
    static const double x[]={0.1488743389816312,0.4333953941292472,
                           0.6794095682990244,0.8650633666889845,0.9739065285171717};
    static const double w[]={0.2955242247147529,0.2692667193099963,
                           0.2190863625159821,0.1494513491505806,0.0666713443086881};
    double xm=0.5*(b+a);
    double xr=0.5*(b-a);
    double s=0;
 
    for (int j=0;j<5;j++) {
        double dx=xr*x[j];
        s += w[j]*(func(xm+dx)+func(xm-dx));
    }
    return s *= xr;
}

//*************************************************
//*                With Pointers!!
//*************************************************

//* Gauss-Legendre (gauleg) int_x1^x2
//* Numerical Recipies page 183
void gauleg(const double x1, const double x2, double *x, double *w,const int n);

//* Gauss-Laguerre
//* improper integral_0^inf  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=x^alpha exp(-x)
//* Numerical Recipies page 184
void gaulag(double *x, double *w, const double alf, const int n);

//* Gauss-Hermite
//* improper integral_-inf^inf  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=exp(-x^2)
//* Numerical Recipies page 185

void gauher(double *x, double *w, const int n);



//* Gauss-Jacobi
//* improper integral_-1^1  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=(1-x)^alpha (1+x)^beta
//* Numerical Recipies page 186
void gaujac(double *x, double *w, const double alf, const double bet, const int n);

//*************************************************
//*            With vector<double>
//*************************************************

//* Gauss-Legendre (gauleg) int_x1^x2
//* Numerical Recipies page 183
void gauleg(const double x1, const double x2, std::vector<double> &x, std::vector<double> &w);


//* Gauss-Laguerre
//* improper integral_0^inf  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=x^alpha exp(-x)
//* Numerical Recipies page 184
void gaulag(std::vector<double> &x, std::vector<double> &w, const double alf);

//* Gauss-Hermite
//* improper integral_-inf^inf  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=exp(-x^2)
//* Numerical Recipies page 185

void gauher(std::vector<double> &x, std::vector<double> &w);



//* Gauss-Jacobi
//* improper integral_-1^1  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=(1-x)^alpha (1+x)^beta
//* Numerical Recipies page 186
void gaujac(std::vector<double> &x,std::vector<double> &w, const double alf, const double bet);



//! ////////////////////////////////////////////////////////
//! ////////////////////////////////////////////////////////
//!                numerical integration                  //    
//! ////////////////////////////////////////////////////////
//! S////////////////////////////////////////////////////////

//Here the weight function is chosen as 1
template<class T>
inline double integral1D_gauleg(T &f,const double x1,const double x2, const int n){
    
    double* x;
    double* w;
    x=new double [n];
    w=new double [n];
    /////////////////
    try{
        gauleg(x1,x2,x,w,n);
    }catch(const char *errtext){
        std::cerr<<errtext<<std::endl;
    }
    
    double s=0;
    #pragma omp parallel for reduction (+:s)
    for (int j=0;j<n;j++) s+=w[j]*f(x[j]);
    
    ///////////
    delete[] x;
    delete[] w;
    return s;
}

//In the function "f" is not included the weight function.
template<class T>
inline double integral1D_gaulag(T &f,const double alpha, const int n){
    
    double* x;
    double* w;
    x=new double [n];
    w=new double [n];
    /////////////////
    try{
        gaulag(x,w,alpha,n);
    }catch(const char *errtext){
        std::cerr<<errtext<<std::endl;
    }
    
    double s=0;
    #pragma omp parallel for reduction (+:s)
    for (int j=0;j<n;j++) s+=w[j]*f(x[j]);
    
    ///////////
    delete[] x;
    delete[] w;
    return s;
}

//In the function "f" is not included the weight function.
template<class T>
inline double integral1D_gauher(T &f, const int n){
    
    double* x;
    double* w;
    x=new double [n];
    w=new double [n];
    /////////////////
    try{
        gauher(x,w,n);
    }catch(const char *errtext){
        std::cerr<<errtext<<std::endl;
    }
    
    double s=0;
    #pragma omp parallel for reduction (+:s)
    for (int j=0;j<n;j++) s+=w[j]*f(x[j]);
    
    ///////////
    delete[] x;
    delete[] w;
    return s;
}

//In the function "f" is not included the weight function.
template<class T>
inline double integral1D_gaujac(T &f, const double alpha, const double beta, const int n){
    
    double* x;
    double* w;
    x=new double [n];
    w=new double [n];
    /////////////////
    try{
        gaujac(x,w,alpha,beta,n);
    }catch(const char *errtext){
        std::cerr<<errtext<<std::endl;
    }
    
    double s=0;
    #pragma omp parallel for reduction (+:s)
    for (int j=0;j<n;j++) s+=w[j]*f(x[j]);
    
    ///////////
    delete[] x;
    delete[] w;
    return s;
}

//In these methods "n" must be an even number
template<class T>
inline double integralSimpson1D(T &f, double extrInf,double extrSup, int n=1e3){

    double dx=(extrSup-extrInf)/(double)n;
    
    double S_i=0;
    double S_p=0;
    
    #pragma omp parallel for reduction (+:S_p) reduction (+:S_i)
    for(int i=1;i<n;i+=2){

        S_i+=f(extrInf+i*dx);
        if(i!=n-1) S_p+=f(extrInf+(i+1)*dx);
    }

    return dx*(f(extrInf)+f(extrSup)+4.0*S_i+2.0*S_p)/3.0;
}

template<class T>
inline double integralSimpson2D(T &f,double extrInfX,double extrInfY,
                         double extrSupX, double extrSupY, int n_x=1e3, int n_y=1e3){

    double dx=(extrSupX-extrInfX)/(double)n_x;
    double dy=(extrSupY-extrInfY)/(double)n_y;

    double S_i=0;
    double S_p=0;
    
    #pragma omp parallel for reduction (+:S_p) reduction (+:S_i)
    for(int i=1;i<n_y;i+=2){
    
        S_i+=Sx_yj(f,dx,extrInfY+i*dy,extrInfX,extrSupX,n_x);
        if(i!=n_y-1) S_p+=Sx_yj(f,dx,extrInfY+(i+1)*dy,extrInfX,extrSupX,n_x);
    }
    
    return dx*dy*(Sx_yj(f,dx,extrInfY,extrInfX,extrSupX,n_x)+
           Sx_yj(f,dx,extrSupY,extrInfX,extrSupX,n_x)+4.0*S_i+2.0*S_p)/9.0;
}

template<class T>
inline double integralSimpson3D(T &f,double extrInfX, double extrInfY, double extrInfZ,
                         double extrSupX, double extrSupY, double extrSupZ, 
                         int n_x=1e2, int n_y=1e2, int n_z=1e2){

    double dx=(extrSupX-extrInfX)/(double)n_x;
    double dy=(extrSupY-extrInfY)/(double)n_y;
    double dz=(extrSupZ-extrInfZ)/(double)n_z;
     
    double S_i=0;
    double S_p=0;
    #pragma omp parallel for reduction (+:S_p) reduction (+:S_i)
    for(int k=1;k<n_y;k+=2){

        S_i+=Sy_zk(f,dx,dy,extrInfY+k*dz,extrInfX,extrInfY,extrSupX,extrSupY,n_x,n_y);
        if(k!=n_z-1) S_p+=Sy_zk(f,dx,dy,extrInfZ+(k+1)*dz,extrInfX,extrInfY,extrSupX,extrSupY,n_x,n_y);
    }
    
    
    return dx*dy*dz*(Sy_zk(f,dx,dy,extrInfZ,extrInfX,extrInfY,extrSupX,extrSupY,n_x,n_y)+
                     Sy_zk(f,dx,dy,extrSupZ,extrInfX,extrInfY,extrSupX,extrSupY,n_x,n_y)+4.0*S_i+2.0*S_p)/27.0;
    
}




//**************************************************
//*         3D integration Gauss-Legendre     
//**************************************************
//* 3D integration!
//* (x,y,z)

class Integral3DGaussLegendre{
    
    private:
        unsigned N_x,N_y,N_z;
        //* Integration limits
        double x_0,x_f;
        double y_0,y_f;
        double z_0,z_f;

        //* Gauss-Legendre weights and zeros for x, y and z
        double *w_x,*w_y,*w_z;
        double *zeros_x,*zeros_y,*zeros_z;


        double get_zero_x(const unsigned) const;
        double get_zero_y(const unsigned) const;
        double get_zero_z(const unsigned) const;

        double integral_factor() const;
    
    public:
        
        Integral3DGaussLegendre(const double x_0,const double x_f,const double y_0,const double y_f,
                                const double z_0,const double z_f,const unsigned N_x,const unsigned N_y,const unsigned N_z);
        ~Integral3DGaussLegendre();

        void set_limits(const double x_0,const double x_f,const double y_0,const double y_f,const double z_0,const double z_f);
       
        template <class M> double integration(M&) const; //* In class M the operator() must be overloaded -> (x,y,z)
};


//* For class M (containing the integrand function)
//* The operator(x,y,z) must be overloaded to return a double
template <class M>
inline double Integral3DGaussLegendre::integration(M &funct)const{
    
    double sum=0;
    #pragma omp parallel for collapse (3) reduction (+:sum)
    for(unsigned i=0;i<N_x;i++)
        for(unsigned j=0;j<N_y;j++)
            for(unsigned k=0;k<N_z;k++)
                sum+=w_x[i]*w_y[j]*w_z[k]*funct(get_zero_x(i),get_zero_y(j),get_zero_z(k)); 
    
    return sum * integral_factor();
}



//**************************************************
//*           3D integration GaussHermite
//**************************************************
//* 3D integration!
//* (x,y,z)
//* \int_{-\inf}^{\inf}   f(x) exp(-2 alpha|r-r0|^2 )

class Integral3DGaussHermite{
    
    private:
        unsigned N_x,N_y,N_z;
        //* gaussian centroid and alpha(width)
        double x_0,y_0,z_0;
        double alpha;

        //* Gauss-Hermite weights and zeros for x, y and z
        double *w_x,*w_y,*w_z;
        double *zeros_x,*zeros_y,*zeros_z;


        double get_zero_x(const unsigned) const;
        double get_zero_y(const unsigned) const;
        double get_zero_z(const unsigned) const;

        double integral_factor() const;
    
    public:
        
        Integral3DGaussHermite(const double alpha,const double x_0,const double y_0,const double z_0,
                               const unsigned N_x,const unsigned N_y,const unsigned N_z);
        ~Integral3DGaussHermite();

        void set_centroid(const double x_0,const double y_0,const double z_0);
        void set_alpha(const double alpha);
       
        template <class M> double integration(M&) const; //* In class M the operator() must be overloaded -> (x,y,z)
};


//* For class M (containing the integrand function)
//* The operator(x,y,z) must be overloaded to return a double
template <class M>
inline double Integral3DGaussHermite::integration(M &funct)const{
    
    double sum=0;
    #pragma omp parallel for collapse (3) reduction (+:sum)
    for(unsigned i=0;i<N_x;i++)
        for(unsigned j=0;j<N_y;j++)
            for(unsigned k=0;k<N_z;k++)
                sum+=w_x[i]*w_y[j]*w_z[k]*funct(get_zero_x(i),get_zero_y(j),get_zero_z(k)); 
    
    return sum * integral_factor();
}



//**************************************************
//*                3D integration
//**************************************************


//* 3D integration!
//* x: Gauss-Laguerre ----> int_0^inf  x^alpha exp(-x) f(x) dx (no limits for the integral over x in the constructor)
//* y: Gauss-Legendre 
//* z: Gauss-Legendre

class Integral3DGaussLegendre_Laguerre{
    
    private:
        unsigned N_x,N_y,N_z;
        double *w_x,*w_y,*w_z;
        double *zeros_x,*zeros_y,*zeros_z;    
    
    public:
        
        Integral3DGaussLegendre_Laguerre(const double,const double,const double,const double, const int,const int,const int);
        ~Integral3DGaussLegendre_Laguerre();
        void set_alpha(const double);
        template <class M> double integration(M&) const; //* In class M the operator() must be overloaded -> (r,u,ph)
};

//* For class M (containing the integrand function)
//* The operator(x,y,z) must be overloaded to return a double
template <class M>
inline double Integral3DGaussLegendre_Laguerre::integration(M &funct)const{
    
    double sum=0;
    #pragma omp parallel for collapse (3) reduction (+:sum)
    for(unsigned i=0;i<N_x;i++)
        for(unsigned j=0;j<N_y;j++)
            for(unsigned k=0;k<N_z;k++)
                sum+=w_x[i]*w_y[j]*w_z[k]*funct(zeros_x[i],zeros_y[j],zeros_z[k]);
    
    return sum;
}




