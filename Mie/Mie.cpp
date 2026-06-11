// Mie.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include"Mie.h"

#include <iostream>
#include<vector>
#include<cmath>
#include<numeric>



bool testLogarithmicDerivatives()
{
    std::complex<double> z(8.0798703484513279);
    size_t N(14);
    auto testForward = getLogarithmicDerivativesForward(z, N);
    auto testBackward = getLogarithmicDerivativesBackward(z, N);

    auto testForward2 = getLogarithmicDerivativesForward(50.0 * 4.0 / 3.0, 66);
    auto testBackward2 = getLogarithmicDerivativesBackward(50.0 * 4.0 / 3.0, 66);
    return true;
}

void testPrahl()
{
    double diameter(1.05);
    double wavelength(0.6328);
    double x(sci::m_pi * diameter / wavelength); //circumference/wavelength
    //MieFlt x(50); //circumference/wavelength
    //MieComplex refractiveIndex(1.55);
    double refractiveIndex(1.55);
    sci::GridData<double, 1> mus{ 1.0, -0.5, -0.5 };//cosines of scattering angle

    double extinctionEfficiency;
    double scatteringEfficiency;
    double asymmetryParameter;
    double backscatterEfficiency;
    sci::GridData<std::complex<double>, 1> s1;
    sci::GridData<std::complex<double>, 1> s2;

    mie(x, refractiveIndex, mus, extinctionEfficiency, scatteringEfficiency, backscatterEfficiency, asymmetryParameter, s1, s2);


    //compare to Scott Prahl code's values
    if (std::abs(extinctionEfficiency - 3.1054257433224577) > 0.00001)
    {
        std::cout << "Extinction efficiency failed\n";
    }

    if (std::abs(scatteringEfficiency - 3.1054257433224577) > 0.00001)
    {
        std::cout << "Scattering efficiency failed\n";
    }

    if (std::abs(asymmetryParameter - 0.63313677398198109) > 0.00001)
    {
        std::cout << "Asymmetry parameter failed\n";
    }

    if (std::abs(backscatterEfficiency - 2.9253412092248068) > 0.00001)
    {
        std::cout << "Asymmetry parameter failed\n";
    }


    std::array<std::complex<double>, 3> prahlS1{ std::complex<double>(21.096312269429383, -8.5770007912199411),
        std::complex<double>(-0.93011284575019837, 1.3792933605425168),
        std::complex<double>(-0.93011284575018260, 1.3792933605424973) };
    std::array<std::complex<double>, 3> prahlS2{ std::complex<double>(21.096312269429383, -8.5770007912199411),
        std::complex<double>(-1.9234842479740848, 0.44338173235166872),
        std::complex<double>(-1.9234842479740211, 0.44338173235171563) };

    if (s1.size() != prahlS1.size())
        std::cout << "S1 is the wrong size\n";
    if (s2.size() != prahlS2.size())
        std::cout << "S1 is the wrong size\n";

    for (size_t i = 0; i < std::min(prahlS1.size(), mus.size()); ++i)
    {
        if (std::abs(s1[i] - prahlS1[i]) > 0.00001)
            std::cout << "S1[" << i << "] failed\n";
    }
    for (size_t i = 0; i < std::min(prahlS2.size(), mus.size()); ++i)
    {
        if (std::abs(s2[i] - prahlS2[i]) > 0.00001)
            std::cout << "S2[" << i << "] failed\n";
    }

    std::cout << "Prahl test completed, if you saw no errors above then it passed\n";
}

int main()
{
    constexpr double testBessel = getBesselMinusHalfOverPlusHalfRatio<double>(1.0, 9, 0.0); //test, check it matches with the Lentz paper
    if (std::abs(testBessel - double(18.95228198)) > 0.00000001)
    {
        std::cout << "Bessel test failed. Result was " << testBessel << " when it should have been 18.95228198" << std::endl;
        return 1;
    }

    getBesselMinusHalfOverPlusHalfRatio<std::complex<double>>(std::complex(1.0, 0.01), 9, 0.0); //just checking it compiles really

    testPrahl();

    testLogarithmicDerivatives();

    
}
