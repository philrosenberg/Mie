// Mie.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include<vector>
#include<cmath>
#include<numeric>
#include<complex>

//https://github.com/philrosenberg/sci
//elementwise maths on arrays
#include<scieng/grid.h>
#include<scieng/gridtransformview.h>

//maths constants
#include<scieng/math.h>

#include"bh.h"

using MieFlt = double;
using MieComplex = std::complex<MieFlt>;
constexpr MieComplex im(MieFlt(0), MieFlt(1));

inline constexpr MieFlt innerProduct(MieComplex v1, MieComplex v2)
{
    return v1.real() * v2.real() + v1.imag() * v2.imag();
}

inline constexpr MieFlt innerProduct(MieFlt v1, MieFlt v2)
{
    return v1 * v2;
}

/*
template<class X>
inline constexpr MieFlt legendreDash0(X x)
{
    return MieFlt(0.0);
}

template<int n, class X>
auto constexpr eigenPi(X x)
{
    auto s = x * eigenPi<n-1>(x);
    auto t = s - eigenPi<n - 2>(x);
}

template<int n, class X>
auto constexpr eigenTau(X x)
{
    auto s = x * eigenPi<n>(x);
    auto t = s - eigenPi<n - 1>(x);
    return MieFlt(n) * t - eigenPi<n - 1>(x);
}

template<class X>
double constexpr eigenPi<0, X>(X x)
{
    return legendreDash0(x);
}
*/

constexpr size_t nTerms(MieFlt x)
{
    if (x <= MieFlt(0.02))
        return size_t(0);
    if (x <= MieFlt(8.0))
        return size_t(x + MieFlt(4.0) * std::pow(x, MieFlt(1.0) / MieFlt(3.0))) + size_t(1);
    if( x < MieFlt(4200))
        return size_t(x + MieFlt(4.05) * std::pow(x, MieFlt(1.0) / MieFlt(3.0))) + size_t(2);
    if (x <= MieFlt(20000.0))
        return size_t(x + MieFlt(4.0) * std::pow(x, MieFlt(1.0) / MieFlt(3.0))) + size_t(2);
    return size_t(0);
}

template<class FLT>
struct Complex
{
    inline static constexpr FLT realUnity = 1.0;

};

template<class FLT>
struct Complex<std::complex<FLT>>
{
    inline static constexpr FLT realUnity = 1.0;
};

template<class FLT>
constexpr inline bool continueTest(FLT ratio, FLT accuracySquared)
{
    FLT diffFromUnity = ratio - FLT(1);
    return diffFromUnity * diffFromUnity > accuracySquared;
}


template<class FLT>
constexpr inline bool continueTest(std::complex<FLT> ratio, FLT accuracySquared)
{
    std::complex<FLT> diffFromUnitySquared = std::norm(ratio) - FLT(1);
    return std::abs(diffFromUnitySquared) > accuracySquared;
}

template<class FLT, class ACC>
constexpr FLT getBesselMinusHalfOverPlusHalfRatio(FLT xInverse, size_t order, ACC accuracy)
{
    //From Lentz(1976) Applied Optics Vol 15, Issue 3,pp. 668-671 Generating Bessel functions in Mie scattering calculations using continued fractions
    //
    //calculates J[n-1/2](x)/J[n+1/2](x) as per the paper above
    //
    //xInverse is 1/x, it may be complex
    //order is n
    //accuracy is the fractional accuracy required to stop the iteration. Note that the
    // norm of this accuracy is tested against the norm of the result. This means that
    // this accuracy could be complex or real, independant of whether x is complex or real
    //
    // If xInverse is complex, then the result will be complex, if xInverse is real the result will be real

    FLT nu = FLT(order) + FLT(0.5);
    auto accuracySquared = std::norm(accuracy);

    FLT top = FLT(2) * nu * xInverse;
    FLT result = top;

    FLT bottom = FLT(2) * (nu + FLT(1)) * xInverse;
    FLT add = FLT(1);
    FLT a = FLT(2) * (nu + add) * xInverse;
    top = a - Complex<FLT>::realUnity / top;
    FLT ratio = top / bottom;
    result *= ratio;

    FLT normRatio = std::norm(ratio);
    while (continueTest(ratio, accuracySquared))
    {
        add += FLT(1);
        a = FLT(2) * (nu + add) * xInverse;
        bottom = a - Complex<FLT>::realUnity / bottom;
        top = a - Complex<FLT>::realUnity / top;
        ratio = top / bottom;
        FLT normRatio = std::norm(ratio);
        result *= ratio;
    }

    return result;
}

template<class FLT1, class FLT2>
constexpr FLT1 getLogarithmicDerivative(FLT1 x, size_t n, FLT2 accuracy)
{
    return -FLT1(n) / x + getBesselMinusHalfOverPlusHalfRatio(MieFlt(1)/x, n, accuracy);
}

//z is the size parameter (x) multplied by the refractive index. It can be real or complex
template<class FLT>
std::vector<MieComplex> getLogarithmicDerivativesForward(FLT z, size_t N)
{
    std::vector<MieComplex>A(N);

    MieComplex exp = std::exp(-MieFlt(2) * im * z);
    A[1] = -MieFlt(1) / z + (MieFlt(1) - exp) / ((MieFlt(1) - exp) / z - im * (MieFlt(1) + exp));
    for (size_t i = 2; i < A.size(); ++i)
        //A[i] = MieFlt(1) / (MieFlt(i + 1) / z - A[i - 1]) - MieFlt(i + 1) / z;
        A[i] = MieFlt(1) / (MieFlt(i) / z - A[i - 1]) - MieFlt(i) / z;

    return A;
}

//z is the size parameter (x) multplied by the refractive index. It can be real or complex
template<class FLT>
std::vector<MieComplex> getLogarithmicDerivativesBackward(FLT z, size_t N)
{
    std::vector<MieComplex>A(N);

    A.back() = getLogarithmicDerivative(z, N, 0.0000001);
    for (size_t i = N - 2; i != size_t(-1); --i)
    {
        auto nOverZ = MieFlt(i + 1) / z;
        A[i] = nOverZ - MieFlt(1) / (A[i + 1] + nOverZ);
    }

    return A;
}


template<class FLT>
std::vector<MieComplex> getLogarithmicDerivativesFastestStable(FLT refractiveIndex, MieFlt x, size_t N)
{
    //calculate A. In some situations this can be done with forward recursion which
    //is faster, If not, then we calculate the last A using the Lentz method and then
    //do a backwards recursion which is always stable.
    //wiscombe found a limiting relationship from m.real from 1.05 to 9.25 and x from 1 to 10,000
    //I'm not sure if this relationship holds outside these bounds, but let's be conservative and
    //use the more stable method.
    auto z = refractiveIndex * x;
    if (refractiveIndex.real() < MieFlt(1.05) || refractiveIndex.real() > MieFlt(9.25) ||
        x < MieFlt(1) || x > MieFlt(10000) ||
        std::abs(refractiveIndex.imag()) * x > (MieFlt(13.78) * refractiveIndex.real() - MieFlt(10.8)) * refractiveIndex.real() + MieFlt(3.9))
        //Use the stable but slower backwards recursion
        return getLogarithmicDerivativesBackward(z, N);
    else
    
        //use the faster forwards recursion
        return getLogarithmicDerivativesForward(z, N);
}

bool testLogarithmicDerivatives()
{
    MieComplex z(8.0798703484513279);
    size_t N(14);
    auto testForward = getLogarithmicDerivativesForward(z, N);
    auto testBackward = getLogarithmicDerivativesBackward(z, N);

    auto testForward2 = getLogarithmicDerivativesForward(50.0 * 4.0 / 3.0, 66);
    auto testBackward2 = getLogarithmicDerivativesBackward(50.0 * 4.0 / 3.0, 66);
    return true;
}

int main()
{
    constexpr MieFlt testBessel = getBesselMinusHalfOverPlusHalfRatio<double>(1.0, 9, 0.0); //test, check it matches with the Lentz paper
    if (std::abs(testBessel - MieFlt(18.95228198)) > 0.00000001)
    {
        std::cout << "Bessel test failed. Result was " << testBessel << " when it should have been 18.95228198" << std::endl;
        return 1;
    }

    getBesselMinusHalfOverPlusHalfRatio<std::complex<double>>(std::complex(1.0, 0.01), 9, 0.0); //just checking it compiles really

    callBH();

    testLogarithmicDerivatives();

    MieFlt diameter(1.05);
    MieFlt wavelength(0.6328);
    MieFlt x(sci::m_pi * diameter / wavelength); //circumference/wavelength
    //MieFlt x(50); //circumference/wavelength
    MieComplex refractiveIndex(1.55);

    using COMPLEX = decltype(refractiveIndex);

    auto z = x * refractiveIndex;
    size_t N = nTerms(x);
    sci::GridData<MieFlt, 1> mus{ 1.0, -0.5, -0.5 };//cosines of scattering angle
    if (N == 0)
        return 1;

    //logarithmic derivative, called A in Wiscombe or D in Bohren anf Huffman
    std::vector<MieComplex>A(getLogarithmicDerivativesFastestStable(refractiveIndex, x, N));


    sci::GridData<MieFlt, 1> eigenPi(mus.size());
    sci::GridData<MieFlt, 1> eigenTau(mus.size());
    sci::GridData<MieFlt, 1> s(mus.size());
    sci::GridData<MieFlt, 1> t(mus.size());
    MieFlt psi;
    MieComplex xi;
    MieComplex a;
    MieComplex b;


    sci::GridData<MieFlt, 1> eigenPiPrev(mus.size());
    sci::GridData<MieFlt, 1> eigenTauPrev(mus.size());
    sci::GridData<MieFlt, 1> sPrev(mus.size());
    sci::GridData<MieFlt, 1> tPrev(mus.size());
    sci::GridData<MieFlt, 1> eigenPiNext(mus.size());
    MieFlt psiPrev;
    MieComplex xiPrev;
    MieFlt psiNext;
    MieComplex xiNext;







    //N will alsways be at least 2, so we can safely directly define the first 2 eigenvalues
    //from Wiscombe equations 3 and 4
    eigenPiPrev = MieFlt(0.0);
    eigenPi = MieFlt(1.0);
    eigenTauPrev = MieFlt(0.0);
    s = mus * eigenPi;
    t = s - eigenPiPrev;
    eigenTau = t - eigenPiPrev;
    eigenPiNext = s + MieFlt(2) * t;
    psiPrev = std::sin(x);
    psi = std::sin(x) / x - std::cos(x);
    xiPrev = psiPrev + im * std::cos(x);
    xi = psi + im * (std::cos(x) / x + std::sin(x));

    //calculate the first a and b terms from the values above
    //Note that the summing over a and b starts with index 1
    //and the prev suffixed variables above have index 0
    COMPLEX aFirstTerm = A[1] / refractiveIndex + MieFlt(1) / x;
    a = (aFirstTerm * psi - psiPrev) / (aFirstTerm * xi - xiPrev);
    COMPLEX bFirstTerm = A[1] * refractiveIndex + MieFlt(1) / x;
    b = (bFirstTerm * psi - psiPrev) / (bFirstTerm * xi - xiPrev);

    //start the summations
    MieFlt commonFactor(3);
    MieFlt extinctionEfficiency(commonFactor * (a.real() + b.real()));
    auto norma = std::norm(a);
    auto normb = std::norm(b);
    MieFlt scatteringEfficiency(commonFactor * (std::norm(a) + std::norm(b)));
    MieFlt asymmetryParameter = commonFactor / MieFlt(2) * innerProduct(a, b); //note first term is 0 for n=1
    MieComplex backscatterTemp = -commonFactor * (a - b);
    MieFlt sign(-1);

    sci::GridData<MieComplex, 1> sPlus = MieFlt(1.5) * (a + b) * (eigenPi + eigenTau);
    sci::GridData<MieComplex, 1> sMinus = MieFlt(1.5) * (a - b) * (eigenPi - eigenTau);

    for (size_t i = 2; i < N; ++i)
    {
        //MieFlt n = MieFlt(i + 1);
        psiNext = commonFactor * psi / x - psiPrev;
        //psiNext = xiNext.real();
        xiNext = commonFactor * xi / x - xiPrev;
        //shuffle the next to current and current to next. This is more efficient if we use swap as it avoids assignments
        std::swap(psiPrev, psi);
        std::swap(psi, psiNext);
        std::swap(xiPrev, xi);
        std::swap(xi, xiNext);

        std::swap(eigenPiPrev, eigenPi);
        std::swap(eigenPi, eigenPiNext);
        s = mus * eigenPi;
        t = s - eigenPiPrev;
        eigenPiNext = s + t *MieFlt(i + 1) / MieFlt(i);
        std::swap(eigenTauPrev, eigenTau);
        eigenTau = MieFlt(i) * t - eigenPiPrev;

        //calculate new a and b
        MieComplex aPrev = a;
        MieComplex bPrev = b;
        sign = sign * MieFlt(-1);
        auto aFirstTerm = A[i] / refractiveIndex + MieFlt(i) / x;
        a = (aFirstTerm * psi - psiPrev) / (aFirstTerm * xi - xiPrev);
        auto bFirstTerm = A[i] * refractiveIndex + MieFlt(i) / x;
        b = (bFirstTerm * psi - psiPrev) / (bFirstTerm * xi - xiPrev);

        commonFactor = (MieFlt(2) * i + MieFlt(1));

        extinctionEfficiency += commonFactor * (a.real() + b.real());
        scatteringEfficiency += commonFactor * (std::norm(a) + std::norm(b));
        asymmetryParameter += (MieFlt(i * i) - MieFlt(1.0)) / MieFlt(i) * (innerProduct(aPrev, a) + innerProduct( bPrev, b))
            + commonFactor / (MieFlt(i) * MieFlt(i + 1)) * innerProduct(a, b);
        backscatterTemp += sign * commonFactor * (a - b);

        sPlus += MieFlt(2 * i + 1) / MieFlt(i * (i + 1)) * (a + b) * (eigenPi + eigenTau);
        sMinus += MieFlt(2 * i + 1) / MieFlt(i * (i + 1)) * (a - b) * (eigenPi - eigenTau);
        //The variables below were just for debugging
        //sci::GridData<MieComplex, 1> s1 = (sPlus + sMinus) / MieFlt(2);
        //sci::GridData<MieComplex, 1> s2 = (sPlus - sMinus) / MieFlt(2);
    }

    MieFlt invXSquared = MieFlt(1) / (x * x);
    extinctionEfficiency *= MieFlt(2) * invXSquared;
    scatteringEfficiency *= MieFlt(2) * invXSquared;
    asymmetryParameter *= MieFlt(4) / scatteringEfficiency * invXSquared;
    MieFlt backscatterEfficiency = std::norm(backscatterTemp) * invXSquared;
    sci::GridData<MieComplex, 1> s1 = (sPlus + sMinus) / MieFlt(2);
    sci::GridData<MieComplex, 1> s2 = (sPlus - sMinus) / MieFlt(2);

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

    std::array<MieComplex, 3> prahlS1{ MieComplex(21.096312269429383, -8.5770007912199411),
        MieComplex(-0.93011284575019837, 1.3792933605425168),
        MieComplex(-0.93011284575018260, 1.3792933605424973) };
    std::array<MieComplex, 3> prahlS2{ MieComplex(21.096312269429383, -8.5770007912199411),
        MieComplex(-1.9234842479740848, 0.44338173235166872),
        MieComplex(-1.9234842479740211, 0.44338173235171563) };
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
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
