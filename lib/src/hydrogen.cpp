#include <hydrogen.h>
#include <cmath>



double Hydrogen::wave_function_radial(const double r) const{

    return 2.0/(n*n)*exp(0.5*(gammal(static_cast<int>(n-l)) - gammal(static_cast<int>(n+l+1)) ) -r/n)*pow(2.0*r/n,l)* 
            std::assoc_laguerre(static_cast<int>(n-l-1),static_cast<int>(2*l+1), 2.0*r/n);

}

double Hydrogen::wave_function_theta(const double th) const{
    return sqrt((2.0*l+1.0)/(4.0*M_PI)) * exp(0.5*( gammal(static_cast<int>(l-abs(m)+1)) - gammal(static_cast<int>(l+abs(m)+1)) )) *
           std::assoc_legendre(l,abs(m),cos(th));
}



double Hydrogen::wave_function_phi(const double ph, const bool real_part) const{
    return ((real_part)? cos(m*ph) : sin(m*ph)); //* exp(i * m * phi)
}


double Hydrogen::wave_function_angular(const double th, const double ph, const bool real_part) const{

    return  wave_function_theta(th)*wave_function_phi(ph,real_part);
}


double Hydrogen::wave_function(const double r, const double th,const double ph,const bool real_part) const{

    return wave_function_radial(r)*wave_function_angular(th,ph,real_part);

}





