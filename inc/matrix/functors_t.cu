#include <iostream>
#include "functors.h"

int main()
{
    std::cout << "Test correct evaluation of dividied differences formulas\n";

    std::vector<double> limits = {1,0.5,1./6.,1./24.};
    std::vector<double> evals_small = {5e-11,1.66666666e-11,4.16666666e-12,8.3333333e-13};
    std::vector<double> evals_half = {5.0000166667e-6,1.6666708333e-6,4.1666750000e-7,8.33334722224206e-8};
    std::vector<double> evals_third = {0.000500166783,0.0001667083416680,0.000041675001389,8.3347224206597e-6};
    std::vector<double> evals_one = {0.7182818284590452,0.2182818284590452,0.05161516179237857,0.009948495125711902};
    std::vector<double> evals_ten = {2201.546579480672,219.6546579480672,21.79879912814005,2.138213246147338};

    std::vector<std::function<double(double)>> funcs = {dg::mat::phi1<double>, dg::mat::phi2<double>, dg::mat::phi3<double>, dg::mat::phi4<double>};
    for( unsigned u=0; u<4; u++)
    {
        std::cout << "Function number "<<u+1<<std::endl;
        std::cout << "Ana "<<limits[u] << "\nNum "<<funcs[u](0)<<std::endl;
        std::cout << "Ana "<<evals_small[u] << "\nNum "<<funcs[u](1e-10)-funcs[u](0)<<std::endl;
        std::cout << "Error " <<evals_small[u] +funcs[u](0) -funcs[u](1e-10)<<std::endl;
        std::cout << "Ana "<<evals_half[u] << "\nNum "<<funcs[u](1e-5)-funcs[u](0)<<std::endl;
        std::cout << "Error " <<evals_half[u] +funcs[u](0) -funcs[u](1e-5)<<std::endl;
        std::cout << "Ana "<<evals_third[u] << "\nNum "<<funcs[u](1e-3)-funcs[u](0)<<std::endl;
        std::cout << "Error " <<evals_third[u] +funcs[u](0) -funcs[u](1e-3)<<std::endl;
        std::cout << "Ana "<<evals_one[u] << "\nNum "<<funcs[u](1.0)-funcs[u](0)<<std::endl;
        std::cout << "Error " <<evals_one[u] +funcs[u](0) -funcs[u](1.0)<<std::endl;
        std::cout << "Ana "<<evals_ten[u] << "\nNum "<<funcs[u](10)-funcs[u](0)<<std::endl;
        std::cout << "Error " <<evals_ten[u] +funcs[u](0) -funcs[u](10)<<std::endl;
    }

    return 0;
}
