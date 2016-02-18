#include <iostream>
#include <complex>

int main()
{
    // This array represents array of complex numbers
    // 1 + 2i, 3 + 4i, 5 + 6i
    double test[3][2] = {{1.0, 2.0}, {3, 4}, {5, 6}};
    // Cast to a pointer to a length 2 array
    double (*ccomp)[2] = static_cast<double(*)[2]>(test);

    for(int ii=0;ii<3;++ii)
        std::cout << ccomp[ii][0] << " + " << ccomp[ii][1] << "i, ";

    std::cout << std::endl;

    std::complex<double>* ccpp = reinterpret_cast<std::complex<double>*>(ccomp);

    for(int ii=0;ii<3;++ii)
        std::cout << ccpp[ii].real() << " + " << ccpp[ii].imag() << "i, ";

    std::cout << std::endl;

    double* ccpp2 = reinterpret_cast<double*>(ccomp);

    for(int ii=0;ii<3;++ii)
        std::cout << ccpp2[2*ii] << " + " << ccpp2[2*ii+1] << "i, ";

    std::cout << std::endl;


    // double (*ccomp2)[2] = reinterpret_cast<double(*)[2]>(ccpp);

    return 0;
}
