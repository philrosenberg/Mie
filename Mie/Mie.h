#pragma once
#include<complex>

//https://github.com/philrosenberg/sci
//elementwise maths on arrays
#include<scieng/grid.h>
#include<scieng/gridtransformview.h>

//maths constants
#include<scieng/math.h>


template<class T>
concept IsComplex = requires(T t)
{
    t.real();
    t.imag();
};

template<IsComplex COMPLEX>
inline constexpr auto innerProduct(COMPLEX v1, COMPLEX v2)
{
    return v1.real() * v2.real() + v1.imag() * v2.imag();
}

template<class T>
inline constexpr auto innerProduct(T v1, T v2)
{
    return v1 * v2;
}

template<class FLT>
constexpr size_t nTerms(FLT x)
{
    //to match prahl
    return size_t(x + FLT(4.05) * std::pow(x, FLT(1.0) / FLT(3.0))) + size_t(2);

    if (x <= FLT(0.02))
        return size_t(0);
    if (x <= FLT(8.0))
        return size_t(x + FLT(4.0) * std::pow(x, FLT(1.0) / FLT(3.0))) + size_t(1);
    if (x < FLT(4200))
        return size_t(x + FLT(4.05) * std::pow(x, FLT(1.0) / FLT(3.0))) + size_t(2);
    if (x <= FLT(20000.0))
        return size_t(x + FLT(4.0) * std::pow(x, FLT(1.0) / FLT(3.0))) + size_t(2);
    return size_t(0);
}

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

template<class MAYBECOMPLEX, class ACC>
constexpr MAYBECOMPLEX getBesselMinusHalfOverPlusHalfRatio(MAYBECOMPLEX xInverse, size_t order, ACC accuracy)
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

    using FLT = decltype(std::real(xInverse));

    MAYBECOMPLEX nu = static_cast<FLT>(order) + FLT(0.5);
    auto accuracySquared = std::norm(accuracy);

    MAYBECOMPLEX top = FLT(2) * nu * xInverse;
    MAYBECOMPLEX result = top;

    MAYBECOMPLEX bottom = FLT(2) * (nu + FLT(1)) * xInverse;
    FLT add = FLT(1);
    MAYBECOMPLEX a = FLT(2) * (nu + add) * xInverse;
    top = a - FLT(1) / top;
    MAYBECOMPLEX ratio = top / bottom;
    result *= ratio;

    FLT normRatio = std::norm(ratio);
    while (continueTest(ratio, accuracySquared))
    {
        add += FLT(1);
        a = FLT(2) * (nu + add) * xInverse;
        bottom = a - FLT(1) / bottom;
        top = a - FLT(1) / top;
        ratio = top / bottom;
        FLT normRatio = std::norm(ratio);
        result *= ratio;
    }

    return result;
}

template<class MAYBECOMPLEX, class FLT>
constexpr MAYBECOMPLEX getLogarithmicDerivative(MAYBECOMPLEX x, size_t n, FLT accuracy)
{
    return -MAYBECOMPLEX(static_cast<FLT>(n)) / x + getBesselMinusHalfOverPlusHalfRatio(FLT(1) / x, n, accuracy);
}

//z is the size parameter (x) multplied by the refractive index. It can be real or complex
template<class COMPLEX>
sci::GridData<COMPLEX, 1> getLogarithmicDerivativesForward(COMPLEX z, size_t N)
    requires (IsComplex<COMPLEX>)
{
    sci::GridData<COMPLEX, 1>A(N);
    using FLT = decltype(z.real());
    const COMPLEX im(0, 1);

    COMPLEX exp = std::exp(-FLT(2) * im * z);
    A[1] = -FLT(1) / z + (FLT(1) - exp) / ((FLT(1) - exp) / z - im * (FLT(1) + exp));
    for (size_t i = 2; i < A.size(); ++i)
        A[i] = FLT(1) / (FLT(i) / z - A[i - 1]) - FLT(i) / z;

    return A;
}

template<class FLT>
sci::GridData<FLT, 1> getLogarithmicDerivativesForward(FLT z, size_t N)
    requires (!IsComplex<FLT>)
{
    sci::GridData<FLT, 1>A(N);

    FLT tanz = std::tan(z);
    A[1] = -FLT(1) / z + (z * tanz) / (tanz - z);
    for (size_t i = 2; i < A.size(); ++i)
        A[i] = FLT(1) / (FLT(i) / z - A[i - 1]) - FLT(i) / z;

    return A;
}

//z is the size parameter (x) multplied by the refractive index. It can be real or complex
template<class MAYBECOMPLEX>
sci::GridData<MAYBECOMPLEX, 1> getLogarithmicDerivativesBackward(MAYBECOMPLEX z, size_t N)
{
    sci::GridData<MAYBECOMPLEX, 1>A(N);
    using FLT = decltype(std::real(z));

    A.back() = getLogarithmicDerivative(z, N, 0.0000001);
    for (size_t i = N - 2; i != size_t(-1); --i)
    {
        auto nOverZ = FLT(i + 1) / z;
        A[i] = nOverZ - FLT(1) / (A[i + 1] + nOverZ);
    }

    return A;
}


template<class MAYBECOMPLEX, class FLT>
sci::GridData<MAYBECOMPLEX, 1> getLogarithmicDerivativesFastestStable(MAYBECOMPLEX refractiveIndex, FLT x, size_t N)
{
    //calculate A. In some situations this can be done with forward recursion which
    //is faster, If not, then we calculate the last A using the Lentz method and then
    //do a backwards recursion which is always stable.
    //wiscombe found a limiting relationship from m.real from 1.05 to 9.25 and x from 1 to 10,000
    //I'm not sure if this relationship holds outside these bounds, but let's be conservative and
    //use the more stable method.

    if (std::real(refractiveIndex) == std::numeric_limits<double>::infinity())
        return sci::GridData<MAYBECOMPLEX, 1>(N, MAYBECOMPLEX(1.0)); //avoid nans when refractiveIndex.real() is infinity. A cancels out anyway in this case

    auto z = refractiveIndex * x;
    if (std::real(refractiveIndex) < FLT(1.05) || std::real(refractiveIndex) > FLT(9.25) ||
        x < FLT(1) || x > FLT(10000) ||
        std::abs(std::imag(refractiveIndex)) * x > (FLT(13.78) * std::real(refractiveIndex) - FLT(10.8)) * std::real(refractiveIndex) + FLT(3.9))
        //Use the stable but slower backwards recursion
        return getLogarithmicDerivativesBackward(z, N);
    else
        //use the faster forwards recursion
        return getLogarithmicDerivativesForward(z, N);
}

template<class FLT, class COMPLEX, class RI>
void rayleigh(FLT x, RI refractiveIndex, const sci::GridData<FLT, 1>& mus,
    sci::GridData<COMPLEX, 1>& s1, sci::GridData<COMPLEX, 1>& s2,
    FLT& extinctionEfficiency, FLT& scatteringEfficiency,
    FLT& backscatterEfficiency, FLT& asymmetryParameter)
{
    using MAYBECOMPLEX = decltype(refractiveIndex);
    auto refractiveIndexSquared = refractiveIndex * refractiveIndex;
    FLT xSquared = x * x;

    FLT aHat2Top = FLT(1) - (FLT(1) / FLT(14)) * xSquared;
    MAYBECOMPLEX aHat2Bottom = FLT(2) * refractiveIndexSquared + FLT(3) - (FLT(2) * refractiveIndexSquared - FLT(7)) / FLT(14) * xSquared;
    COMPLEX aHat2 = COMPLEX(0.0, 1.0) * (xSquared * (refractiveIndexSquared - FLT(1.0)) / FLT(15.0) * aHat2Top / aHat2Bottom);

    MAYBECOMPLEX bHat1Top = FLT(1) + (FLT(2) * refractiveIndexSquared - FLT(5)) / FLT(70) * xSquared;
    MAYBECOMPLEX bHat1Bottom = FLT(1) - (FLT(2) * refractiveIndexSquared - FLT(5)) / FLT(30) * xSquared;
    COMPLEX bHat1 = COMPLEX(0.0, 1.0) * (xSquared * (refractiveIndexSquared - FLT(1)) / FLT(45) * bHat1Top / bHat1Bottom);

    COMPLEX D = refractiveIndexSquared + FLT(2) + (FLT(1) - FLT(7) / FLT(10) * refractiveIndexSquared) * xSquared
        - (FLT(8) * refractiveIndexSquared * refractiveIndexSquared - FLT(385) * refractiveIndexSquared + FLT(350)) / FLT(1400) * xSquared * xSquared
        + COMPLEX(0.0, 1.0) * (FLT(2) * (refractiveIndexSquared - FLT(1)) / FLT(3) * x * xSquared * (FLT(1) - xSquared / FLT(10)));
    
    MAYBECOMPLEX aHat1Top = FLT(1) - xSquared / FLT(10) + (FLT(4) * refractiveIndexSquared + FLT(5)) / FLT(1400) * xSquared * xSquared;
    COMPLEX aHat1 = COMPLEX(0.0, 1.0) * (FLT(2) * (refractiveIndexSquared - FLT(1)) * aHat1Top / (FLT(3) * D));

    FLT T = std::norm(aHat1) + std::norm(bHat1) + (FLT(5) / FLT(3)) * std::norm(aHat2);

    s1 = FLT(1.5) * x * xSquared * (aHat1 + (bHat1 + (FLT(5) / FLT(3)) * aHat2) * mus);
    s2 = FLT(1.5) * x * xSquared * (bHat1 + (bHat1 + aHat1 * mus + (FLT(5) / FLT(3)) * aHat2) * (FLT(2) * mus - FLT(1)));

    asymmetryParameter = std::real(aHat1 * std::conj(aHat2 + bHat1)) / T;

    scatteringEfficiency = FLT(6) * xSquared * xSquared * T;
    extinctionEfficiency = FLT(6) * x * std::real(aHat1 + bHat1 + (FLT(5) / FLT(3)) * aHat2);

    COMPLEX ss1 = 1.5 * xSquared * (aHat1 - bHat1 - (FLT(5.0) / FLT(3.0)) * aHat2);
    backscatterEfficiency = FLT(4) * std::norm(ss1);
}

//Set the template argument LARGE_RI to false for intermediate values of RI and true
//for large values of RI. When set to true the calculation uses maths that is stable
//as RI -> infinity (conducting sphere). When set to false the calculation uses maths
//that is faster
template<class FLT, class COMPLEX, class RI>
void mie(FLT x, RI refractiveIndex, const sci::GridData<FLT, 1>& mus,
    sci::GridData<COMPLEX, 1>& s1, sci::GridData<COMPLEX, 1>& s2,
    FLT& extinctionEfficiency, FLT& scatteringEfficiency,
    FLT& backscatterEfficiency, FLT& asymmetryParameter)
{
    using MAYBECOMPLEX = decltype(refractiveIndex);

    if (std::norm(refractiveIndex) * x * x < 0.01)
    {
        rayleigh(x, refractiveIndex, mus, s1, s2, extinctionEfficiency, scatteringEfficiency, backscatterEfficiency, asymmetryParameter);
        return;
    }

    auto z = x * refractiveIndex;
    size_t N = nTerms(x);
    //sci::GridData<MieFlt, 1> mus{ 1.0, -0.5, -0.5 };//cosines of scattering angle
    if (N < 2)
        throw("Failure in Mie code - the number of terms needed was less than 2. Did you pass in negative sizes or something?");

    //logarithmic derivative, called A in Wiscombe or D in Bohren anf Huffman
    //In some circumstances it is stable to calculate this forward, but some 
    //circumstances need us to calculate the last element, then work backwards.
    //So we pre-calculate this to satisfy both scenarios
    sci::GridData<MAYBECOMPLEX, 1>A(getLogarithmicDerivativesFastestStable(refractiveIndex, x, N));




    //**************************************************
    //Declare all the variables needed for the sums
    //and assign them using the values from the n=1 term
    //**************************************************

    //These variables are needed for working out the overall scatering parameters
    FLT psiPrev = std::sin(x);
    FLT psi = std::sin(x) / x - std::cos(x);
    FLT psiNext(0); //this isn't used until it is defined in the loop
    COMPLEX xiPrev = psiPrev + COMPLEX(FLT(0), std::cos(x));
    COMPLEX xi = psi + COMPLEX(FLT(0), (std::cos(x) / x + std::sin(x)));
    COMPLEX xiNext(FLT(0)); //this isn't used until it is defined in the loop
    FLT sign(-1);


    //These variables are needed for the angular resolved scallering
    //
    // eigenPi = 1 for all angles and eigenPiPrev = 0 for all angles
    // hence
    //t = s = mus
    // hence 
    // eigenTau = mus
    // eigenPiNext = 3*mus
    //eigenTauPrev = MieFlt(0.0);
    sci::GridData<FLT, 1>  eigenPiPrev(mus.size(), FLT(0.0));
    sci::GridData<FLT, 1>  eigenPi(mus.size(), FLT(1.0));
    sci::GridData<FLT, 1>  eigenTau(mus);
    sci::GridData<FLT, 1>  eigenPiNext = FLT(3) * mus;

    auto reciprocalRefractiveIndex = FLT(1.0) / refractiveIndex;
    //a and b are used for both angular resolved and overall scattering properties
    // A will alsways have a size of at least 2, so this is safe
    //calculate the first a and b terms from the values above
    //Note that the summing over a and b starts with index 1
    //and the prev suffixed variables above have index 0
    MAYBECOMPLEX aFirstTerm = A[1] * reciprocalRefractiveIndex + FLT(1) / x;
    COMPLEX a = (aFirstTerm * psi - psiPrev) / (aFirstTerm * xi - xiPrev);
    //this formulation is more stable when using large refractive indices
    MAYBECOMPLEX bFirstTerm = A[1] + FLT(1) * reciprocalRefractiveIndex / x;
    COMPLEX b = (bFirstTerm * psi - psiPrev * reciprocalRefractiveIndex) / (bFirstTerm * xi - xiPrev * reciprocalRefractiveIndex);
    


    //These are the values we will be summing
    FLT commonFactor(3);
    extinctionEfficiency = (commonFactor * (a.real() + b.real()));
    scatteringEfficiency = (commonFactor * (std::norm(a) + std::norm(b)));
    asymmetryParameter = commonFactor / FLT(2) * innerProduct(a, b); //note first term is 0 for n=1
    COMPLEX backscatterTemp = -commonFactor * (a - b);

    sci::GridData<COMPLEX, 1> sPlus = FLT(1.5) * (a + b) * (eigenPi + eigenTau);
    sci::GridData<COMPLEX, 1> sMinus = FLT(1.5) * (a - b) * (eigenPi - eigenTau);

    //s and t are intermediate calculations for calculating eigenPi and eigenTau
    sci::GridData<FLT, 1> s(mus.size());
    sci::GridData<FLT, 1> t(mus.size());

    //loop through each of the remaining sum terms
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
        eigenPiNext = s + t * FLT(i + 1) / FLT(i);
        eigenTau = FLT(i) * t - eigenPiPrev;

        //calculate new a and b
        COMPLEX aPrev = a;
        COMPLEX bPrev = b;
        sign = sign * FLT(-1);
        aFirstTerm = A[i] * reciprocalRefractiveIndex + FLT(i) / x;
        a = (aFirstTerm * psi - psiPrev) / (aFirstTerm * xi - xiPrev);
        
        //this formulation is more stable when using large refractive indices
        bFirstTerm = A[i] + FLT(i) * reciprocalRefractiveIndex / x;
        b = (bFirstTerm * psi - psiPrev * reciprocalRefractiveIndex) / (bFirstTerm * xi - xiPrev * reciprocalRefractiveIndex);
        

        commonFactor = (FLT(2) * i + FLT(1));

        extinctionEfficiency += commonFactor * (a.real() + b.real());
        scatteringEfficiency += commonFactor * (std::norm(a) + std::norm(b));
        asymmetryParameter += (FLT(i * i) - FLT(1.0)) / FLT(i) * (innerProduct(aPrev, a) + innerProduct(bPrev, b))
            + commonFactor / (FLT(i) * FLT(i + 1)) * innerProduct(a, b);
        backscatterTemp += sign * commonFactor * (a - b);

        sPlus += FLT(2 * i + 1) / FLT(i * (i + 1)) * (a + b) * (eigenPi + eigenTau);
        sMinus += FLT(2 * i + 1) / FLT(i * (i + 1)) * (a - b) * (eigenPi - eigenTau);
        //The variables below were just for debugging
        //sci::GridData<MieComplex, 1> s1 = (sPlus + sMinus) / MieFlt(2);
        //sci::GridData<MieComplex, 1> s2 = (sPlus - sMinus) / MieFlt(2);
    }

    FLT invXSquared = FLT(1) / (x * x);
    extinctionEfficiency *= FLT(2) * invXSquared;
    scatteringEfficiency *= FLT(2) * invXSquared;
    asymmetryParameter *= FLT(4) / scatteringEfficiency * invXSquared;
    backscatterEfficiency = std::norm(backscatterTemp) * invXSquared;
    s1 = (sPlus + sMinus) / FLT(2);
    s2 = (sPlus - sMinus) / FLT(2);
}