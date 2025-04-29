#pragma once
#include <complex>


//* Hydrogen atom (energies and normalized eigenstates)
class Hydrogen{

    private:
        unsigned n;
        unsigned l;
        int m;
        double rMax;

    public:

        Hydrogen(unsigned n, unsigned l, int m): n(n), l(l), m(m), rMax(3.0*n*n) {}

        unsigned get_n() const { return n; }
        unsigned get_l() const { return l; }
        unsigned get_m() const { return m; }
        double get_rMax() const { return rMax; }

        void set_n(const unsigned n) { this->n = n; rMax=3.0*n*n;}
        void set_l(const unsigned l) { this->l = l; }
        void set_m(const int m) { this->m = m; }


        //! All in atomic units!
        double energy() const { return -0.5/(n*n); } 
        
        //* Wave function (Radial component)
        double wave_function_radial(const double r) const; //* including radial normalization constant

        //* Wave function theta dependence only

        double wave_function_theta(const double th) const; //* including angular normalization constant


        //* Wave function phi dependence only
        double wave_function_phi(const double ph, const bool real_part) const; //* exp(i*m*phi)
        

        //* Angular component (both th and ph are included)
        double wave_function_angular(const double th, const double ph, const bool real_part) const; //* normalized real_part= true(real_part) or false(imag_part)

        //* Whole wave function!
        double wave_function(const double r, const double th,const double ph,const bool real_part) const; //* normalized wave function
        std::complex<double> wave_function(const double r, const double th,const double ph) const;

};

