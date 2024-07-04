#include <iostream>
#include <fstream>
#include <hydrogen.h>
#include <gauss_wgts.h>


class InnerProduct{

    private:
        Hydrogen *func1;
        Hydrogen *func2;
        bool real_part;

    public:
        InnerProduct() :  func1(nullptr), func2(nullptr) {}
        InnerProduct(Hydrogen &func1,Hydrogen &func2) : func1(&func1), func2(&func2) {}
        void set_func1(Hydrogen &f1) { func1 = &f1; }
        void set_func2(Hydrogen &f2) { func2 = &f2; }

        void set_real(const bool real_part){this->real_part = real_part;}

        double operator()(const double r, const double th, const double ph) const;
};


double InnerProduct::operator()(const double r, const double th, const double ph) const{

    //* Inner product Real part!
    return (func1->wave_function(r,th,ph,true)  * func2->wave_function(r,th,ph,true) + 
            func1->wave_function(r,th,ph,false) * func2->wave_function(r,th,ph,false) )* r*r*sin(th) ;
}




void file_wf_radial(Hydrogen &hydrogen, std::string &&filename){

    const unsigned N=1000;
    double dr = hydrogen.get_rMax()/N;
    std::ofstream file(filename.c_str());

    for(unsigned i=0;i<N;i++)
        file<<dr*i<<'\t'<<hydrogen.wave_function_radial(dr*i)<<'\n';
    
    file.close();

}


void file_wf_theta(Hydrogen &hydrogen, std::string &&filename){

    const unsigned N=1000;
    double dth = M_PI/N;
    std::ofstream file(filename.c_str());

    for(unsigned i=0;i<N;i++)
        file<<dth*i<<'\t'<<hydrogen.wave_function_theta(dth*i)<<'\n';
    
    file.close();

}



int main(){
    const int N_quadratures = 100;
    std::string output_file="output/hydrogen";
    
    Hydrogen wf(10,6,-1);
    InnerProduct norm(wf,wf);

    Integral3DGaussLegendre integral(0,wf.get_rMax(),0.0,M_PI,0.0,2.0*M_PI,N_quadratures,N_quadratures,N_quadratures);
    
    std::cout<<"Inner Product: "<<integral.integration(norm)<<'\n';

    std::cout<<"Angular l=5, m=-4: \t" <<wf.wave_function_theta(1.645)<<'\n';

    file_wf_radial(wf,output_file+".r_wf");
    file_wf_theta(wf,output_file+".th_wf");





    return 0;
}