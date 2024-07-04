#include <gauss_wgts.h>






//Auxiliar function
double gammln(const double xx) {
    int j;
    double x,tmp,y,ser;
    static const double cof[14]={57.1562356658629235,-59.5979603554754912,
                      14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
                      .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
                      -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
                      .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
    if (xx <= 0) throw("bad arg in gammln");
    y=x=xx;
    tmp = x+5.24218750000000000;
    tmp = (x+0.5)*log(tmp)-tmp;
    ser = 0.999999999999997092;
    for (j=0;j<14;j++) ser += cof[j]/++y;
    return tmp+log(2.5066282746310005*ser/x);
}

//*************************************************
//*                With Pointers!!
//*************************************************

//* Gauss-Legendre (gauleg) int_x1^x2
//* Numerical Recipies page 183
void gauleg(const double x1, const double x2, double *x, double *w,const int n){
    
    const double EPS=1.0e-10;
    double z1,z,xm,xl,pp,p3,p2,p1;
    int m=(n+1)/2;

    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    
    for (int i=0;i<m;i++) {
        z=cos(3.141592654*(i+0.75)/(n+0.5));
        do {
            p1=1.0;
            p2=0.0;
            
            for (int j=0;j<n;j++) {
                p3=p2;
                p2=p1;
                p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
            }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;

        } while (fabs(z-z1) > EPS);
        x[i]=xm-xl*z;
        x[n-1-i]=xm+xl*z;
        w[i]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n-1-i]=w[i];
    }
}

//* Gauss-Laguerre
//* improper integral_0^inf  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=x^alpha exp(-x)
//* Numerical Recipies page 184
void gaulag(double *x, double *w, const double alf, const int n){
    
    const int MAXIT=100;
    const double EPS=1.0e-10;
    int i,its,j;
    double ai,p1,p2,p3,pp,z,z1;
    
    for (i=0;i<n;i++) {
        if (i == 0) 
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        else if (i == 1) 
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        else {
            ai=i-1;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                 (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=0;its<MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=0;j<n;j++) {
                p3=p2;
                p2=p1;
                p1=((2.0*j+1.0+alf-z)*p2-(j+alf)*p3)/(double)(j+1.0);
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its >= MAXIT) throw("too many iterations in gaulag");
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
    
}

//* Gauss-Hermite
//* improper integral_-inf^inf  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=exp(-x^2)
//* Numerical Recipies page 185

void gauher(double *x, double *w, const int n){
    
    const double EPS=1.0e-10,PIM4=0.7511255444649425;
    const int MAXIT=10;
    int i,its,j,m;
    double p1,p2,p3,pp,z,z1;
    m=(n+1)/2;
    for (i=0;i<m;i++) {
        if (i == 0)
            z=sqrt(double(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
        else if (i == 1) z -= 1.14*pow((double)n,0.426)/z;
        else if (i == 2) z=1.86*z-0.86*x[0];
        else if (i == 3) z=1.91*z-0.91*x[1];
        else z=2.0*z-x[i-2];
        for (its=0;its<MAXIT;its++) {
            p1=PIM4;
            p2=0.0;
            for (j=0;j<n;j++) {
                p3=p2;
                p2=p1;
                p1=z*sqrt(2.0/(j+1))*p2-sqrt((double)j/(j+1))*p3;
            }
            pp=sqrt((double)(2*n))*p2;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }

        if (its >= MAXIT) throw("too many iterations in gauher");
        x[i]=z;
        x[n-1-i] = -z;
        w[i]=2.0/(pp*pp);
        w[n-1-i]=w[i];
    }
}



//* Gauss-Jacobi
//* improper integral_-1^1  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=(1-x)^alpha (1+x)^beta
//* Numerical Recipies page 186
void gaujac(double *x, double *w, const double alf, const double bet, const int n){
    
    const int MAXIT=10;
    const double EPS=1.0e-10;
    int i,its,j;
    double alfbet,an,bn,r1,r2,r3;
    double a,b,c,p1,p2,p3,pp,temp,z,z1;
    
    for (i=0;i<n;i++) {
        if (i == 0) {
            an=alf/n;
            bn=bet/n;
            r1=(1.0+alf)*(2.78/(4.0+n*n)+0.768*an/n);
            r2=1.0+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
            z=1.0-r1/r2;
        }else if (i == 1) {
            r1=(4.1+alf)/((1.0+alf)*(1.0+0.156*alf));
            r2=1.0+0.06*(n-8.0)*(1.0+0.12*alf)/n;
            r3=1.0+0.012*bet*(1.0+0.25*fabs(alf))/n;
            z -= (1.0-z)*r1*r2*r3;
        }else if (i == 2) {
            r1=(1.67+0.28*alf)/(1.0+0.37*alf);
            r2=1.0+0.22*(n-8.0)/n;
            r3=1.0+8.0*bet/((6.28+bet)*n*n);
            z -= (x[0]-z)*r1*r2*r3;
        }else if (i == n-2) {
            r1=(1.0+0.235*bet)/(0.766+0.119*bet);
            r2=1.0/(1.0+0.639*(n-4.0)/(1.0+0.71*(n-4.0)));
            r3=1.0/(1.0+20.0*alf/((7.5+alf)*n*n));
            z += (z-x[n-4])*r1*r2*r3;
        }else if (i == n-1) {
            r1=(1.0+0.37*bet)/(1.67+0.28*bet);
            r2=1.0/(1.0+0.22*(n-8.0)/n);
            r3=1.0/(1.0+8.0*alf/((6.28+alf)*n*n));
            z += (z-x[n-3])*r1*r2*r3;
        }else z=3.0*x[i-1]-3.0*x[i-2]+x[i-3];
        alfbet=alf+bet;
        for (its=1;its<=MAXIT;its++) {
            temp=2.0+alfbet;
            p1=(alf-bet+temp*z)/2.0;
            p2=1.0;
            for (j=2;j<=n;j++) {
                p3=p2;
                p2=p1;
                temp=2*j+alfbet;
                a=2*j*(j+alfbet)*(temp-2.0);
                b=(temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z);
                c=2.0*(j-1+alf)*(j-1+bet)*temp;
                p1=(b*p2-c*p3)/a;
            }
            pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z));
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) throw("too many iterations in gaujac");
        x[i]=z;
        w[i]=exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0)-gammln(n+alfbet+1.0))*temp*pow(2.0,alfbet)/(pp*p2);
    }
}


//*************************************************
//*            With vector<double>
//*************************************************

//* Gauss-Legendre (gauleg) int_x1^x2
//* Numerical Recipies page 183
void gauleg(const double x1, const double x2, std::vector<double> &x, std::vector<double> &w){
    
    const int n=x.size();
    const double EPS=1.0e-10;
    double z1,z,xm,xl,pp,p3,p2,p1;
    int m=(n+1)/2;

    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    
    for (int i=0;i<m;i++) {
        z=cos(3.141592654*(i+0.75)/(n+0.5));
        do {
            p1=1.0;
            p2=0.0;
            
            for (int j=0;j<n;j++) {
                p3=p2;
                p2=p1;
                p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
            }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;

        } while (fabs(z-z1) > EPS);
        x[i]=xm-xl*z;
        x[n-1-i]=xm+xl*z;
        w[i]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n-1-i]=w[i];
    }
}





//* Gauss-Laguerre
//* improper integral_0^inf  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=x^alpha exp(-x)
//* Numerical Recipies page 184
void gaulag(std::vector<double> &x, std::vector<double> &w, const double alf){
    
    const int n=x.size();
    const int MAXIT=100;
    const double EPS=1.0e-10;
    int i,its,j;
    double ai,p1,p2,p3,pp,z,z1;
    
    for (i=0;i<n;i++) {
        if (i == 0) 
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        else if (i == 1) 
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        else {
            ai=i-1;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                 (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=0;its<MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=0;j<n;j++) {
                p3=p2;
                p2=p1;
                p1=((2.0*j+1.0+alf-z)*p2-(j+alf)*p3)/(double)(j+1.0);
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its >= MAXIT) throw("too many iterations in gaulag");
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
    
}

//* Gauss-Hermite
//* improper integral_-inf^inf  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=exp(-x^2)
//* Numerical Recipies page 185

void gauher(std::vector<double> &x, std::vector<double> &w){
    
    const int n=x.size();
    const double EPS=1.0e-10,PIM4=0.7511255444649425;
    const int MAXIT=10;
    int i,its,j,m;
    double p1,p2,p3,pp,z,z1;
    m=(n+1)/2;
    for (i=0;i<m;i++) {
        if (i == 0)
            z=sqrt(double(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
        else if (i == 1) z -= 1.14*pow((double)n,0.426)/z;
        else if (i == 2) z=1.86*z-0.86*x[0];
        else if (i == 3) z=1.91*z-0.91*x[1];
        else z=2.0*z-x[i-2];
        for (its=0;its<MAXIT;its++) {
            p1=PIM4;
            p2=0.0;
            for (j=0;j<n;j++) {
                p3=p2;
                p2=p1;
                p1=z*sqrt(2.0/(j+1))*p2-sqrt((double)j/(j+1))*p3;
            }
            pp=sqrt((double)(2*n))*p2;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }

        if (its >= MAXIT) throw("too many iterations in gauher");
        x[i]=z;
        x[n-1-i] = -z;
        w[i]=2.0/(pp*pp);
        w[n-1-i]=w[i];
    }
}



//* Gauss-Jacobi
//* improper integral_-1^1  with integrand of the form w(x)f(x)
//* The weight function is of the form ----> w(x)=(1-x)^alpha (1+x)^beta
//* Numerical Recipies page 186
void gaujac(std::vector<double> &x,std::vector<double> &w, const double alf, const double bet){
    
    const int n=x.size();
    const int MAXIT=10;
    const double EPS=1.0e-10;
    int i,its,j;
    double alfbet,an,bn,r1,r2,r3;
    double a,b,c,p1,p2,p3,pp,temp,z,z1;
    
    for (i=0;i<n;i++) {
        if (i == 0) {
            an=alf/n;
            bn=bet/n;
            r1=(1.0+alf)*(2.78/(4.0+n*n)+0.768*an/n);
            r2=1.0+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
            z=1.0-r1/r2;
        }else if (i == 1) {
            r1=(4.1+alf)/((1.0+alf)*(1.0+0.156*alf));
            r2=1.0+0.06*(n-8.0)*(1.0+0.12*alf)/n;
            r3=1.0+0.012*bet*(1.0+0.25*fabs(alf))/n;
            z -= (1.0-z)*r1*r2*r3;
        }else if (i == 2) {
            r1=(1.67+0.28*alf)/(1.0+0.37*alf);
            r2=1.0+0.22*(n-8.0)/n;
            r3=1.0+8.0*bet/((6.28+bet)*n*n);
            z -= (x[0]-z)*r1*r2*r3;
        }else if (i == n-2) {
            r1=(1.0+0.235*bet)/(0.766+0.119*bet);
            r2=1.0/(1.0+0.639*(n-4.0)/(1.0+0.71*(n-4.0)));
            r3=1.0/(1.0+20.0*alf/((7.5+alf)*n*n));
            z += (z-x[n-4])*r1*r2*r3;
        }else if (i == n-1) {
            r1=(1.0+0.37*bet)/(1.67+0.28*bet);
            r2=1.0/(1.0+0.22*(n-8.0)/n);
            r3=1.0/(1.0+8.0*alf/((6.28+alf)*n*n));
            z += (z-x[n-3])*r1*r2*r3;
        }else z=3.0*x[i-1]-3.0*x[i-2]+x[i-3];
        alfbet=alf+bet;
        for (its=1;its<=MAXIT;its++) {
            temp=2.0+alfbet;
            p1=(alf-bet+temp*z)/2.0;
            p2=1.0;
            for (j=2;j<=n;j++) {
                p3=p2;
                p2=p1;
                temp=2*j+alfbet;
                a=2*j*(j+alfbet)*(temp-2.0);
                b=(temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z);
                c=2.0*(j-1+alf)*(j-1+bet)*temp;
                p1=(b*p2-c*p3)/a;
            }
            pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z));
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) throw("too many iterations in gaujac");
        x[i]=z;
        w[i]=exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0)-gammln(n+alfbet+1.0))*temp*pow(2.0,alfbet)/(pp*p2);
    }
}



//**************************************************
//*         3D integration Gauss-Legendre     
//**************************************************
//* 3D integration!
//* (x,y,z)


//* Constructor for Gauss-Legendre algorithm over x, y and z
Integral3DGaussLegendre::Integral3DGaussLegendre(const double x_0,const double x_f,const double y_0,const double y_f,
                                                 const double z_0,const double z_f,const unsigned N_x,const unsigned N_y,const unsigned N_z):  
                                                 N_x(N_x),N_y(N_y),N_z(N_z),x_0(x_0),x_f(x_f),y_0(y_0),y_f(y_f),z_0(z_0),z_f(z_f){
    
    try{
    
        //* x arrays are alwas allocated in memory
        w_x=new double[N_x];
        zeros_x=new double[N_x];

        //* Allways getting the zeros and weights for the x component
        gauleg(-1.0,1.0,zeros_x,w_x,N_x);//* x component


        
        //* y arrays are only allocated in memory when Nx != Ny
        if(N_x==N_y){
            w_y=w_x;
            zeros_y=zeros_x;
        }else{ 
            w_y=new double[N_y];
            zeros_y=new double[N_y];

            //* Getting the zeros and weights for the y component
            gauleg(-1.0,1.0,zeros_y,w_y,N_y);//* y component

        }

        //* z arrays are only allocated in memory when Nx != Nz and Ny!= Nz
        if(N_x==N_z){
            w_z=w_x;
            zeros_z=zeros_x;
        }else{

            if(N_y==N_z){
                w_z=w_y;
                zeros_z=zeros_y;
            }else{

                w_z=new double[N_z];
                zeros_z=new double[N_z];

                //* Getting the zeros and weights for the z component
                gauleg(-1.0,1.0,zeros_z,w_z,N_z);//* z component
            }
        }
            
    }catch(const char *err){
        std::cerr<<err<<std::endl;
    }                                                        
}

void Integral3DGaussLegendre::set_limits(const double x_0,const double x_f,const double y_0,const double y_f,const double z_0,const double z_f){

    this->x_0=x_0;
    this->x_f=x_f;
    this->y_0=y_0;
    this->y_f=y_f;
    this->z_0=z_0;
    this->z_f=z_f;
}



double new_zero(const double old_zero, const double inf_lim, const double sup_lim){

    return 0.5*(1.0-old_zero)*inf_lim+0.5*(1.0+old_zero)*sup_lim; 
}



double Integral3DGaussLegendre::get_zero_x(const unsigned i) const{

    return new_zero(zeros_x[i],x_0,x_f);
}


double Integral3DGaussLegendre::get_zero_y(const unsigned i) const{

    return new_zero(zeros_y[i],y_0,y_f);
}


double Integral3DGaussLegendre::get_zero_z(const unsigned i) const{

    return new_zero(zeros_z[i],z_0,z_f);
}


double Integral3DGaussLegendre::integral_factor() const{
    
    return (x_f-x_0)*(y_f-y_0)*(z_f-z_0)/8.0;
}





Integral3DGaussLegendre::~Integral3DGaussLegendre(){
    
    delete[] w_x;
    delete[] zeros_x;

    if(N_x!=N_y){
        delete[] w_y;
        delete[] zeros_y;
    }

    if(N_x!=N_z && N_y!=N_z){
        delete[] w_z;
        delete[] zeros_z;
    }
}


//**************************************************
//*           3D integration GaussHermite
//**************************************************
//* 3D integration!
//* (x,y,z)
//* \int_{-\inf}^{\inf}   f(x) exp(-2 alpha|r-r0|^2 )


//* Constructor for Gauss-Legendre algorithm over x, y and z
Integral3DGaussHermite::Integral3DGaussHermite(const double alpha,const double x_0,const double y_0,const double z_0,
                                               const unsigned N_x,const unsigned N_y,const unsigned N_z):  
                                               N_x(N_x),N_y(N_y),N_z(N_z),x_0(x_0),y_0(y_0),z_0(z_0),alpha(alpha){
    
    try{
    
        //* x arrays are alwas allocated in memory
        w_x=new double[N_x];
        zeros_x=new double[N_x];

        //* Allways getting the zeros and weights for the x component
        gauher(zeros_x,w_x,N_x); //* x component

        //* y arrays are only allocated in memory when Nx != Ny
        if(N_x==N_y){
            w_y=w_x;
            zeros_y=zeros_x;
        }else{ 
            w_y=new double[N_y];
            zeros_y=new double[N_y];

            //* Getting the zeros and weights for the y component
            gauher(zeros_y,w_y,N_y);//* y component

        }

        //* z arrays are only allocated in memory when Nx != Nz and Ny!= Nz
        if(N_x==N_z){
            w_z=w_x;
            zeros_z=zeros_x;
        }else{

            if(N_y==N_z){
                w_z=w_y;
                zeros_z=zeros_y;
            }else{

                w_z=new double[N_z];
                zeros_z=new double[N_z];

                //* Getting the zeros and weights for the z component
                gauher(zeros_z,w_z,N_z);//* z component
            }
        }
            
    }catch(const char *err){
        std::cerr<<err<<std::endl;
    }                                                        
}

void Integral3DGaussHermite::set_centroid(const double x_0,const double y_0,const double z_0){

    this->x_0=x_0;
    this->y_0=y_0;
    this->z_0=z_0;
}

void Integral3DGaussHermite::set_alpha(const double alpha){
    this->alpha=alpha;
}



double new_zero_hermite(const double old_zero, const double xyz_0, const double alpha){

    return old_zero/sqrt(2.0*alpha)+xyz_0; 
}



double Integral3DGaussHermite::get_zero_x(const unsigned i) const{

    return new_zero_hermite(zeros_x[i],x_0,alpha);
}


double Integral3DGaussHermite::get_zero_y(const unsigned i) const{

    return new_zero_hermite(zeros_y[i],y_0,alpha);
}


double Integral3DGaussHermite::get_zero_z(const unsigned i) const{

    return new_zero_hermite(zeros_z[i],z_0,alpha);
}


double Integral3DGaussHermite::integral_factor() const{
    
    return 1.0/pow(2.0*alpha,1.5);
}


Integral3DGaussHermite::~Integral3DGaussHermite(){
    
    delete[] w_x;
    delete[] zeros_x;

    if(N_x!=N_y){
        delete[] w_y;
        delete[] zeros_y;
    }

    if(N_x!=N_z && N_y!=N_z){
        delete[] w_z;
        delete[] zeros_z;
    }
}



//**************************************************
//*                3D integration
//**************************************************


//* 3D integration!
//* x: Gauss-Laguerre ----> int_0^inf  x^alpha exp(-x) f(x) dx (no limits for the integral over x in the constructor)
//* y: Gauss-Legendre 
//* z: Gauss-Legendre


//* Constructor for Gauss-Laguerre algorithm over x and Gauss-Legendre over y and z
//* Integration over x goes from 0 to inf
//* Integral limits must be given for y and z
Integral3DGaussLegendre_Laguerre::Integral3DGaussLegendre_Laguerre(const double y_0,const double y_f,const double z_0,const double z_f,
                                       const int N_x,const int N_y,const int N_z) : N_x(N_x), N_y(N_y), N_z(N_z){
    
    w_x=new double[N_x];
    zeros_x=new double[N_x];
 
    w_y=new double[N_y];
    zeros_y=new double[N_y];
    //---------------------
    w_z=new double[N_z];
    zeros_z=new double[N_z];
        
    try{
        
        //* Gauss-Legendre over y and z
        #pragma omp parallel sections
        {
            #pragma omp section
                gauleg(y_0,y_f,zeros_y,w_y,N_y);//* y component
            #pragma omp section               
                gauleg(z_0,z_f,zeros_z,w_z,N_z); //* z component
        }
        
    }catch(const char *err){
        std::cerr<<err<<std::endl;
    }
}


Integral3DGaussLegendre_Laguerre::~Integral3DGaussLegendre_Laguerre(){
    
    delete[] w_x;
    delete[] zeros_x;
    delete[] w_y;
    delete[] zeros_y;
    delete[] w_z;
    delete[] zeros_z;
}


//* The integrand must have a term of the form: exp(-alpha*x)
//! This method must be called before calling "integration" when the method Gauss-Laguerre is used 
//! For Gauss-Legendre this method is not used
void Integral3DGaussLegendre_Laguerre::set_alpha(const double alpha){
    
    //* Gauss-Laguerre over "x"
    try{
        gaulag(zeros_x,w_x,alpha,N_x);
    }catch(const char *err){
        std::cerr<<err<<std::endl;
    }
    
}









