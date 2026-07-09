#pragma once

//Copyright 2026 Philip Rosenberg
//p.d.rosenberg@gmail.com


#include<complex>

//https://github.com/philrosenberg/sci
//elementwise maths on arrays
#include<scieng/grid.h>
#include<scieng/gridview.h>
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
inline constexpr size_t nTerms(FLT x)
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
inline constexpr bool continueTest(FLT ratio, FLT accuracySquared)
{
    FLT diffFromUnity = ratio - FLT(1);
    return diffFromUnity * diffFromUnity > accuracySquared;
}


template<class FLT>
inline constexpr bool continueTest(std::complex<FLT> ratio, FLT accuracySquared)
{
    std::complex<FLT> diffFromUnitySquared = std::norm(ratio) - FLT(1);
    return std::abs(diffFromUnitySquared) > accuracySquared;
}

template<class MAYBECOMPLEX, class ACC>
inline constexpr MAYBECOMPLEX getBesselMinusHalfOverPlusHalfRatio(MAYBECOMPLEX xInverse, size_t order, ACC accuracy)
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
inline constexpr MAYBECOMPLEX getLogarithmicDerivative(MAYBECOMPLEX x, size_t n, FLT accuracy)
{
    return -MAYBECOMPLEX(static_cast<FLT>(n)) / x + getBesselMinusHalfOverPlusHalfRatio(FLT(1) / x, n, accuracy);
}

//z is the size parameter (x) multplied by the refractive index. It can be real or complex
template<class COMPLEX, class A_TYPE>
inline void getLogarithmicDerivativesForward(COMPLEX z, size_t N, A_TYPE& A)
    requires (IsComplex<COMPLEX>)
{
    using FLT = decltype(z.real());
    const COMPLEX im(0, 1);

    COMPLEX exp = std::exp(-FLT(2) * im * z);
    A[1] = -FLT(1) / z + (FLT(1) - exp) / ((FLT(1) - exp) / z - im * (FLT(1) + exp));
    for (size_t i = 2; i < A.size(); ++i)
        A[i] = FLT(1) / (FLT(i) / z - A[i - 1]) - FLT(i) / z;

}

template<class FLT, class A_TYPE>
inline void getLogarithmicDerivativesForward(FLT z, size_t N, A_TYPE& A)
    requires (!IsComplex<FLT>)
{
    FLT tanz = std::tan(z);
    A[1] = -FLT(1) / z + (z * tanz) / (tanz - z);
    for (size_t i = 2; i < A.size(); ++i)
        A[i] = FLT(1) / (FLT(i) / z - A[i - 1]) - FLT(i) / z;
}

//z is the size parameter (x) multplied by the refractive index. It can be real or complex
template<class MAYBECOMPLEX, class A_TYPE>
inline void getLogarithmicDerivativesBackward(MAYBECOMPLEX z, size_t N, A_TYPE& A)
{
    using FLT = decltype(std::real(z));

    A[N-1] = getLogarithmicDerivative(z, N, 0.0000001);
    for (size_t i = N - 2; i != size_t(-1); --i)
    {
        auto nOverZ = FLT(i + 1) / z;
        A[i] = nOverZ - FLT(1) / (A[i + 1] + nOverZ);
    }
}


template<class MAYBECOMPLEX, class FLT, class A_TYPE>
inline void getLogarithmicDerivativesFastestStable(MAYBECOMPLEX refractiveIndex, FLT x, size_t N, A_TYPE &A)
{
    //calculate A. In some situations this can be done with forward recursion which
    //is faster, If not, then we calculate the last A using the Lentz method and then
    //do a backwards recursion which is always stable.
    //wiscombe found a limiting relationship from m.real from 1.05 to 9.25 and x from 1 to 10,000
    //I'm not sure if this relationship holds outside these bounds, but let's be conservative and
    //use the more stable method.

    if (std::real(refractiveIndex) == std::numeric_limits<double>::infinity())
        for (auto& a : A)
            a = MAYBECOMPLEX(1.0); //avoid nans when refractiveIndex.real() is infinity. A cancels out anyway in this case

    auto z = refractiveIndex * x;
    if (std::real(refractiveIndex) < FLT(1.05) || std::real(refractiveIndex) > FLT(9.25) ||
        x < FLT(1) || x > FLT(10000) ||
        std::abs(std::imag(refractiveIndex)) * x > (FLT(13.78) * std::real(refractiveIndex) - FLT(10.8)) * std::real(refractiveIndex) + FLT(3.9))
        //Use the stable but slower backwards recursion
        getLogarithmicDerivativesBackward(z, N, A);
    else
        //use the faster forwards recursion
        getLogarithmicDerivativesForward(z, N, A);
}

template<class FLT, class COMPLEX, class RI>
inline void rayleigh(FLT x, RI refractiveIndex, const sci::GridData<FLT, 1>& mus,
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

class PreAllocator
{
public:
    PreAllocator(size_t nBytes)
        :m_currentByte(0), m_data(nBytes)
    {
    }
    template<class T>
    std::span<T> getSpan(size_t nElements)
    {
        T* start = (T*)&m_data[m_currentByte];
        //get a length in Bytes that is an integer number of the size of a pointer
        size_t lengthBytes = (nElements * sizeof(T) / sizeof(void*) + 1) *sizeof(void*);
        if (m_currentByte + lengthBytes > m_data.size())
            throw(std::bad_alloc());
        m_currentByte += lengthBytes;
        return std::span<T>(start, nElements);
    }
    template<class T>
    std::span<T> getSpan(size_t nElements, const T& value)
    {
        std::span<T> result = getSpan<T>(nElements);
        for (T& r : result)
            r = value;
        return result;
    }
    template<class T, class RANGE>
    std::span<T> getSpan(const RANGE& range)
    {
        size_t nElements = std::size(range);
        std::span<T> result = getSpan<T>(nElements);
        auto dest = result.begin();
        auto src = std::begin(range);
        for (; dest != result.end(); ++dest, ++src)
            *dest = *src;
        return result;
    }
    void clear()
    {
        m_currentByte = 0;
    }
    bool isClear() const
    {
        return m_currentByte == 0;
    }
private:
    std::vector<char> m_data;
    size_t m_currentByte;
};

template<class FLT, class COMPLEX>
class MiePreAllocator : public PreAllocator
{
public:
    MiePreAllocator(size_t nAngles, FLT largestX)
        :PreAllocator((nAngles + sizeof(void*))* (6 * sizeof(FLT) + 2 * sizeof(COMPLEX)) + (nTerms(largestX)+sizeof(void*))*sizeof(COMPLEX))
    {
    }
};

class PreAllocatorClearer
{
public:
    PreAllocatorClearer(PreAllocator *preAllocator)
        :m_preAllocator(preAllocator)
    {
        if(!m_preAllocator->isClear())
            throw(std::bad_alloc());
    }
    ~PreAllocatorClearer()
    {
        m_preAllocator->clear();
    }
private:
    PreAllocator* m_preAllocator;
};

template<class FLT, class COMPLEX, class MAYBECOMPLEX>
inline void mie(FLT x, MAYBECOMPLEX refractiveIndex, const sci::GridData<FLT, 1>& mus,
    sci::GridData<COMPLEX, 1>& s1, sci::GridData<COMPLEX, 1>& s2,
    FLT& extinctionEfficiency, FLT& scatteringEfficiency,
    FLT& backscatterEfficiency, FLT& asymmetryParameter)
{
    MiePreAllocator<FLT, COMPLEX> preAllocator(mus.size(), x);
    mie(x, refractiveIndex, mus, s1, s2, extinctionEfficiency, scatteringEfficiency, backscatterEfficiency, asymmetryParameter, preAllocator);
}


template<class FLT, class COMPLEX, class MAYBECOMPLEX>
inline void mie(FLT x, MAYBECOMPLEX refractiveIndex, const sci::GridData<FLT, 1>& mus,
    sci::GridData<COMPLEX, 1>& s1, sci::GridData<COMPLEX, 1>& s2,
    FLT& extinctionEfficiency, FLT& scatteringEfficiency,
    FLT& backscatterEfficiency, FLT& asymmetryParameter,
    PreAllocator &preAllocator)
{
    //creating this object checks the preallocated memory is empty and
    // when it goes out of scope, it clears the memory ready for reuse
    PreAllocatorClearer memoryClearer(&preAllocator);

    if (std::norm(refractiveIndex) * x * x < 0.01)
    {
        rayleigh(x, refractiveIndex, mus, s1, s2, extinctionEfficiency, scatteringEfficiency, backscatterEfficiency, asymmetryParameter);
        return;
    }

    auto z = x * refractiveIndex;
    size_t N = nTerms(x);
    if (N < 2)
        throw("Failure in Mie code - the number of terms needed was less than 2. Did you pass in negative sizes or something?");

    //logarithmic derivative, called A in Wiscombe or D in Bohren anf Huffman
    //In some circumstances it is stable to calculate this forward, but some 
    //circumstances need us to calculate the last element, then work backwards.
    //So we pre-calculate A to satisfy both scenarios
    auto  A = sci::views::make_grid_view_1d(preAllocator.getSpan<MAYBECOMPLEX>(N));
    getLogarithmicDerivativesFastestStable(refractiveIndex, x, N, A);
    //sci::GridData<MAYBECOMPLEX, 1>A(getLogarithmicDerivativesFastestStable(refractiveIndex, x, N));




    //**************************************************
    //Declare all the variables needed for the sums
    //and assign them using the values from the n=1 term
    //**************************************************
    auto reciprocalRefractiveIndex = FLT(1.0) / refractiveIndex;
    FLT reciprocalX = FLT(1.0) / x;
    FLT sinX = std::sin(x);
    FLT cosX = std::cos(x);
    const long nAngles = (long)mus.size();

    //These variables are needed for working out the overall scatering parameters
    FLT psiPrev = sinX;
    FLT psi = sinX * reciprocalX - cosX;
    FLT psiNext(0); //this isn't used until it is defined in the loop
    COMPLEX xiPrev = psiPrev + COMPLEX(FLT(0), cosX);
    COMPLEX xi = psi + COMPLEX(FLT(0), (cosX * reciprocalX + sinX));
    COMPLEX xiNext(FLT(0)); //this isn't used until it is defined in the loop
    FLT sign(-1.0);

    
    //These variables are needed for the angular resolved scattering
    //
    // initially, eigenPi = 1 for all angles and eigenPiPrev = 0 for all angles
    // hence
    //t = s = mus
    // hence 
    // eigenTau = mus
    // eigenPiNext = 3*mus
    //eigenTauPrev = MieFlt(0.0);
    auto  eigenPiPrev = sci::views::make_grid_view_1d(preAllocator.getSpan<FLT>(mus.size(), FLT(0.0)));
    auto  eigenPi = sci::views::make_grid_view_1d(preAllocator.getSpan<FLT>(mus.size(), FLT(1.0)));
    auto  eigenTau = sci::views::make_grid_view_1d(preAllocator.getSpan<FLT>(mus));
    auto  eigenPiNext = sci::views::make_grid_view_1d(preAllocator.getSpan<FLT>(FLT(3) * mus));
    //sci::GridData<FLT, 1>  eigenPiPrev(mus.size(), FLT(0.0));
    //sci::GridData<FLT, 1>  eigenPi(mus.size(), FLT(1.0));
    //sci::GridData<FLT, 1>  eigenTau(mus);
    //sci::GridData<FLT, 1>  eigenPiNext = FLT(3) * mus;
    static_assert (sci::IsGrid<decltype(eigenPi)>);
    
    auto  sPlus = sci::views::make_grid_view_1d(preAllocator.getSpan<COMPLEX>(mus.size()));
    auto  sMinus = sci::views::make_grid_view_1d(preAllocator.getSpan<COMPLEX>(mus.size()));
    //sci::GridData<COMPLEX, 1> sPlus(mus.size());
    //sci::GridData<COMPLEX, 1> sMinus(mus.size());
    
    
    //a and b are used for both angular resolved and overall scattering properties
    // A will alsways have a size of at least 2, so this is safe
    //calculate the first a and b terms from the values above
    //Note that the summing over a and b starts with index 1
    //and the prev suffixed variables above have index 0
    MAYBECOMPLEX aFirstTerm = A[1] * reciprocalRefractiveIndex + FLT(1) * reciprocalX;
    COMPLEX a = (aFirstTerm * psi - psiPrev) / (aFirstTerm * xi - xiPrev);
    //this formulation is more stable when using large refractive indices
    MAYBECOMPLEX bFirstTerm = A[1] + FLT(1) * reciprocalRefractiveIndex * reciprocalX;
    COMPLEX b = (bFirstTerm * psi - psiPrev * reciprocalRefractiveIndex) / (bFirstTerm * xi - xiPrev * reciprocalRefractiveIndex);
    

    
    
    //These are the values we will be summing
    FLT commonFactor(3);
    extinctionEfficiency = (commonFactor * (a.real() + b.real()));
    scatteringEfficiency = (commonFactor * (std::norm(a) + std::norm(b)));
    asymmetryParameter = commonFactor / FLT(2) * innerProduct(a, b); //note first term is 0 for n=1
    COMPLEX backscatterTemp = -commonFactor * (a - b);
    for (size_t i = 0; i < mus.size(); ++i)
    {
        sPlus[i] = FLT(1.5) * (a + b) * (eigenPi[i] + eigenTau[i]);
        sMinus[i] = FLT(1.5) * (a - b) * (eigenPi[i] - eigenTau[i]);
    }

    //loop through each of the remaining sum terms
    for (size_t i = 2; i < N; ++i)
    {
        psiNext = commonFactor * psi * reciprocalX - psiPrev;
        //psiNext = xiNext.real();
        xiNext = commonFactor * xi * reciprocalX - xiPrev;
        //shuffle the next to current and current to next. This is more efficient if we use swap as it avoids assignments
        
        psiPrev = psi;
        psi = psiNext;
        xiPrev = xi;
        xi = xiNext;
        decltype(eigenPiPrev) eigenPiTemp(eigenPiPrev);
        eigenPiPrev.retarget(eigenPi);
        eigenPi.retarget(eigenPiNext);
        eigenPiNext.retarget(eigenPiTemp);
        for (long j = 0; j < nAngles; ++j)
        {
            auto s = mus[j] * eigenPi[j];
            auto t = s - eigenPiPrev[j];
            eigenPiNext[j] = s + t * FLT(i + 1) / FLT(i);
            eigenTau[j] = FLT(i) * t - eigenPiPrev[j];
        }
        
        //calculate new a and b
        COMPLEX aPrev = a;
        COMPLEX bPrev = b;
        sign = sign * FLT(-1);
        aFirstTerm = A[i] * reciprocalRefractiveIndex + FLT(i) * reciprocalX;
        a = (aFirstTerm * psi - psiPrev) / (aFirstTerm * xi - xiPrev);
        
        //this formulation is more stable when using large refractive indices
        bFirstTerm = A[i] + FLT(i) * reciprocalRefractiveIndex * reciprocalX;
        b = (bFirstTerm * psi - psiPrev * reciprocalRefractiveIndex) / (bFirstTerm * xi - xiPrev * reciprocalRefractiveIndex);
        
        
        commonFactor = (FLT(2) * i + FLT(1));

        extinctionEfficiency += commonFactor * (a.real() + b.real());
        scatteringEfficiency += commonFactor * (std::norm(a) + std::norm(b));
        asymmetryParameter += (FLT(i * i) - FLT(1.0)) / FLT(i) * (innerProduct(aPrev, a) + innerProduct(bPrev, b))
            + commonFactor / (FLT(i) * FLT(i + 1)) * innerProduct(a, b);
        backscatterTemp += sign * commonFactor * (a - b);
        
        for (long j = 0; j < nAngles; ++j)
        {
            sPlus[j] += FLT(2 * i + 1) / FLT(i * (i + 1)) * (a + b) * (eigenPi[j] + eigenTau[j]);
            sMinus[j] += FLT(2 * i + 1) / FLT(i * (i + 1)) * (a - b) * (eigenPi[j] - eigenTau[j]);
        }
        //The variables below were just for debugging
        //sci::GridData<MieComplex, 1> s1 = (sPlus + sMinus) / MieFlt(2);
        //sci::GridData<MieComplex, 1> s2 = (sPlus - sMinus) / MieFlt(2);
    }
    
    FLT invXSquared = FLT(1) / (x * x);
    extinctionEfficiency *= FLT(2) * invXSquared;
    scatteringEfficiency *= FLT(2) * invXSquared;
    asymmetryParameter *= FLT(4) / scatteringEfficiency * invXSquared;
    backscatterEfficiency = std::norm(backscatterTemp) * invXSquared;
    s1.resize(mus.size());
    s2.resize(mus.size());
    for (long i = 0; i < nAngles; ++i)
    {
        s1[i] = (sPlus[i] + sMinus[i]) / FLT(2);
        s2[i] = (sPlus[i] - sMinus[i]) / FLT(2);
    }
    
    preAllocator.clear();
}