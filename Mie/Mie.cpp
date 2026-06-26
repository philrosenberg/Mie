// mie.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include"mie.h"

#include <iostream>
#include<vector>
#include<cmath>
#include<numeric>
#include<random>



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

    mie(x, refractiveIndex, mus, s1, s2, extinctionEfficiency, scatteringEfficiency, backscatterEfficiency, asymmetryParameter);


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

std::complex<double> Lentz_Dn(std::complex<double> z, long n)
{
	return getLogarithmicDerivative(z, n, 0.0000001);
}

void Dn_up(std::complex<double> z, long n, sci::GridData< std::complex<double>, 1>& D)
{
	D = getLogarithmicDerivativesForward(z, n);
}

void Dn_down(std::complex<double> z, long n, sci::GridData< std::complex<double>, 1>& D)
{
	D = getLogarithmicDerivativesBackward(z, n);
}

double
LegendrePn(long n, double x)
{
	double		pk, pkp1, pkm1;
	long		k;

	if (n <= 0)
		return 1.0;
	if (n == 1)
		return x;

	if (x >= 1.0)
		return 1.0;
	if (x <= -1.0)
		return (n % 2) ? -1.0 : 1.0;

	pk = x;
	pkm1 = 1.0;
	for (k = 1; k < n; k++) {
		pkp1 = ((2 * k + 1) * x * pk - k * pkm1) / (k + 1);
		pkm1 = pk;
		pk = pkp1;
	}

	return pk;
}




void
LegendrePn_and_Pnm1(long n, double x, double* Pnm1Val, double* PnVal)
{
	long		k;
	double		Pk, Pkp1;
	double		Pkm1 = 1.0;

	*Pnm1Val = 1.0;
	*PnVal = 1.0;
	if (x >= 1.0)
		x = 1.0;
	if (x <= -1.0)
		x = -1.0;

	Pk = x;

	for (k = 1; k < n; k++) {
		Pkp1 = ((2 * k + 1) * x * Pk - k * Pkm1) / (k + 1);
		Pkm1 = Pk;
		Pk = Pkp1;
	}

	*Pnm1Val = Pkm1;
	*PnVal = Pk;
}





double
LegendrePnd(long n, double x)
{
	double		p, pminus, pplus;
	long		i;

	if (n <= 0)
		return 0;
	if (n == 1)
		return 1;

	if (x > 1.0)
		x = 1.0;
	if (x < -1.0)
		x = -1.0;

	pminus = 0;
	p = 1;

	for (i = 1; i < n; i++) {
		pplus = ((2 * i + 1) * x * p - (i + 1) * pminus) / i;
		pminus = p;
		p = pplus;
	}
	return p;
}

double
LegendrePndd(long n, double x)
{
	double		p, pminus, pplus;
	long		m;

	if (n <= 1)
		return 0;
	if (n == 2)
		return 3;

	if (x > 1.0)
		x = 1.0;
	if (x < -1.0)
		x = -1.0;

	pminus = 0;
	p = 3;

	for (m = 2; m < n; m++) {
		pplus = ((2 * m + 1) * x * p - (m + 2) * pminus) / (m - 1);
		pminus = p;
		p = pplus;
	}
	return p;

}

void		    bracketroot(double (*fx) (double), double x1, double x2, long n,
	double	    xb1[], double xb2[], long* nrequested) {
	long		nfound, i;
	double		x, fp, fc, dx;

	if ((n <= 0) | (*nrequested <= 0))
		return;

	nfound = 0;
	dx = (x2 - x1) / n;
	x = x1;
	fp = (*fx) (x);
	for (i = 0; i < n; i++) {
		x += dx;
		fc = (*fx) (x);
		if (((fc < 0.0) && (fp > 0.0)) || ((fp < 0.0) && (fc > 0.0))) {
			nfound++;
			xb1[nfound - 1] = x - dx;
			xb2[nfound - 1] = x;
			if (*nrequested == nfound)
				return;
		}
		fp = fc;
	}
	*nrequested = nfound;
}




double		    saferoot(void (*funcd) (double, double*, double*), double x1, double x2, double xacc) {
	double		df, dx, dxold, f, fh, fl;
	double		temp, xh, xl, rts;
	double		temp1, temp2;
	long		j, MAXIT = 100;

	(*funcd) (x1, &fl, &df);
	if (fl == 0.0)
		return x1;

	(*funcd) (x2, &fh, &df);
	if (fh == 0.0)
		return x2;

	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		printf("saferoot -- Root must be bracketed.\n");
		printf("saferoot -- x1  = %10.5f; f(x1) = %10.5f \n", x1, fl);
		printf("saferoot -- x2  = %10.5f; f(x2) = %10.5f \n", x2, fh);
		exit(1);
	}
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	}
	else {
		xh = x1;
		xl = x2;
	}

	rts = 0.5 * (x1 + x2);
	dxold = fabs(x2 - x1);
	dx = dxold;
	(*funcd) (rts, &f, &df);

	for (j = 1; j <= MAXIT; j++) {

		temp1 = (rts - xh) * df - f;
		temp2 = (rts - xl) * df - f;

		if ((temp1 * temp2 >= 0.0) || (fabs(2.0 * f) > fabs(dxold * df))) {
			dxold = dx;
			dx = 0.5 * (xh - xl);
			rts = xl + dx;
			if (xl == rts)
				return rts;
		}
		else {
			dxold = dx;
			dx = f / df;
			temp = rts;
			rts -= dx;
			if (temp == rts)
				return rts;
		}

		if (fabs(dx) < xacc)
			return rts;

		(*funcd) (rts, &f, &df);

		if (f < 0.0)
			xl = rts;
		else
			xh = rts;
	}

	printf("saferoot -- Root cannot be found.\n");
	printf("saferoot -- Executed %ld iterations. \n", j);
	exit(1);
	return 0.0;
}

#define NSLICES 1000
#define EPS 1e-16 \

static long	    Lobatto_n_minus_1;



static void
Lobatto_error(const char* s)
{
	printf("%s\n", s);
	exit(1);
}




static double
Lobatto_fn1(double x)
{
	return LegendrePnd(Lobatto_n_minus_1, x);
}




static void
Lobatto_fn2(double x, double* f, double* df)
{
	*f = LegendrePnd(Lobatto_n_minus_1, x);
	*df = LegendrePndd(Lobatto_n_minus_1, x);
}

void
Lobatto(double a, double b, double* x, double* w, long n)
{
	long		nby2, n_odd, i;
	double		xm, xl, pnval;

	if (n < 3)
		Lobatto_error("Number of Lobatto quadrature points less than 3");

	if (x == NULL)
		Lobatto_error("NULL value passed for x array to Lobatto");

	if (w == NULL)
		Lobatto_error("NULL value passed for w array to Lobatto");

	x[n - 1] = 1.0;
	w[n - 1] = 2.0 / n / (n - 1);
	nby2 = n / 2 - 1;
	n_odd = n % 2;

	switch (n) {
	case 4:
		x[2] = 0.4472135954999579;

		w[2] = 0.8333333333333333;
		break;


	case 8:
		x[6] = 0.8717401485096066;
		x[5] = 0.5917001814331423;
		x[4] = 0.2092992179024789;

		w[6] = 0.2107042271435061;
		w[5] = 0.3411226924835043;
		w[4] = 0.4124587946587038;
		break;


	case 16:
		x[14] = 0.9695680462702180;
		x[13] = 0.8992005330934720;
		x[12] = 0.7920082918618151;
		x[11] = 0.6523887028824931;
		x[10] = 0.4860594218871376;
		x[9] = 0.2998304689007632;
		x[8] = 0.1013262735219495;

		w[14] = 0.0508503610059200;
		w[13] = 0.0893936973259308;
		w[12] = 0.1242553821325141;
		w[11] = 0.1540269808071643;
		w[10] = 0.1774919133917041;
		w[9] = 0.1936900238252036;
		w[8] = 0.2019583081782299;
		break;


	default:
	{
		long		nb, ndiv, size;
		double		z;

		Lobatto_n_minus_1 = n - 1;
		size = NSLICES;
		sci::GridData<double, 1> xb1(size);
		sci::GridData<double, 1> xb2(size);


		ndiv = nby2;
		do {
			ndiv *= 2;
			if (ndiv >= NSLICES)
				ndiv = NSLICES - 1;
			nb = nby2;
			bracketroot(Lobatto_fn1, 0.0, 1.0, ndiv, &xb1[0], &xb2[0], &nb);
		} while (nb < nby2 && ndiv < NSLICES - 1);

		if (nb < nby2)
			Lobatto_error("Cannot find enough roots for Lobatto quadrature");



		for (i = 0; i < nby2; i++) {
			z = saferoot(Lobatto_fn2, xb1[i], xb2[i], EPS);
			x[n - nby2 + i - 1] = z;
			pnval = LegendrePn(n - 1, z);
			w[n - nby2 + i - 1] = w[n - 1] / (pnval * pnval);
		}
		break;
	}


	}

	if (n_odd) {
		i = nby2 + 1;
		x[i] = 0.0;
		pnval = LegendrePn(n - 1, 0.0);
		w[i] = 2 / (n * (n - 1) * pnval * pnval);
	}
	for (i = 0; i <= nby2; i++) {
		w[i] = w[n - i - 1];
		x[i] = -x[n - i - 1];
	}



	if ((a != -1.0) | (b != 1.0)) {
		xm = (b + a) / 2.0;
		xl = (b - a) / 2.0;

		for (i = 0; i < n; i++) {
			x[i] = xm - xl * x[i];
			w[i] = xl * w[i];
		}
	}
}

void TestsFromPrahl()
{
	{
		long		    nstop;

		printf("\n***********************************************\n");
		printf("Zeroth test for logarithmic derivative\n");
		printf("   The result should for D_9(1.0) = 9.95228198\n");
		std::complex<double> z = std::complex<double>(1.0, 0.0);
		std::complex<double> y = Lentz_Dn(z, 9L);
		printf("   The actual value               = %11.8f +i%12.8f\n\n", y.real(), y.imag());

		z = std::complex<double>(62 * 1.28, -62 * 1.37);
		nstop = 50;

		sci::GridData< std::complex<double>, 1> D(nstop);
		printf("   For n = %ld \n", nstop);
		printf("   For j = %ld \n", 10L);
		printf("   For z = %10.5f +i %10.5f\n", z.real(), z.imag());

		printf("   Mathematica                 gives %10.6f +i %10.6f\n", 0.004087, 1.0002620);
		y = Lentz_Dn(z, 10L);
		printf("   Dn[10] continued fraction   gives %10.6f +i %10.6f\n", y.real(), y.imag());

		Dn_up(z, nstop, D);
		printf("   Dn[10] upwards recurrence   gives %10.6f +i %10.6f\n", D[10].real(), D[10].imag());

		Dn_down(z, nstop, D);
		printf("   Dn[10] downwards recurrence gives %10.6f +i %10.6f\n", D[10].real(), D[10].imag());

	}


	
	{
		long		    nangles, i;
		double x, rho, qext, qsca, qback, g;

		nangles = 10;
		std::complex<double> m(1.55,0.0);
		x = 5.213;
		rho = 2 * x * (m.real() - 1);

		printf("\n***********************************************\n");
		printf("First mie Test -- cf. Bohren and Huffman pg 482\n");
		printf("    index of medium      %7.4f\n", 1.0);
		printf("    real index of sphere %7.4f\n", m.real());
		printf("    imag index of sphere %7.4f\n", m.imag());
		printf("    quadrature angles    %ld\n\n", nangles);

		sci::GridData< std::complex<double>, 1> s1(nangles);
		sci::GridData< std::complex<double>, 1> s2(nangles);
		sci::GridData< double, 1> mu(nangles);
		sci::GridData< double, 1> w(nangles);

		Lobatto(0.0, 3.1415926535, &mu[0], &w[0], nangles);
		for (i = 0; i < nangles; ++i)
			mu[i] = cos(mu[i]);

		mie(x, m, mu, s1, s2, qext, qsca, qback, g);

		printf("          x        Qsca       Qext      Qback    g\n");
		printf("BH    %7.3f %10.6f %10.6f %10.6f \n", 5.213, 3.105430, 3.10543, 2.92534);
		printf("Prahl %7.3f %10.6f %10.6f %10.6f %10.6f \n", 5.213, 3.104996, 3.104996, 2.92534, 0.633104);
		printf("This  %7.3f %10.6f %10.6f %10.6f %10.6f \n", x, qsca, qext, qback, g);
	}

	

	{
		double		    pi = 3.14159265358979;
		long		    nangles = 0;
		sci::GridData<std::complex<double>, 1> s1;
		sci::GridData<std::complex<double>, 1> s2;
		sci::GridData<double, 1> mu;
		std::complex<double>    m;
		double		    x = 50.0 * pi;
		double		    qext, qsca, qback, g;

		printf("\n***********************************************\n");
		printf("Second mie Test -- Dave Table 2\n");
		printf("          n                 Qa           Dave         Prahl\n");

		m = std::complex<double>(1.342, 0.0);
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%10.5g%-+7.4fi    %10.5f    %10.5f    %10.5f\n", m.real(), m.imag(), qext - qsca, 0.0, 0.0);

		m = std::complex<double>(1.342, -0.0001);
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%10.5g%-+7.4fi    %10.5f    %10.5f    %10.5f\n", m.real(), m.imag(), qext - qsca, 0.0535, 0.05355);

		m = std::complex<double>(1.342, -0.01);
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%10.5g%-+7.4fi    %10.5f    %10.5f    %10.5f\n", m.real(), m.imag(), qext - qsca, 0.9649, 0.96495);

		m = std::complex<double>(1.342, -0.2);
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%10.5g%-+7.4fi    %10.5f    %10.5f    %10.5f\n", m.real(), m.imag(), qext - qsca, 0.9542, 0.95419);

		m = std::complex<double>(1.342, -0.4);
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%10.5g%-+7.4fi    %10.5f    %10.5f    %10.5f\n", m.real(), m.imag(), qext - qsca, 0.9221, 0.92111);

		m = std::complex<double>(1.342, -0.6);
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%10.5g%-+7.4fi    %10.5f    %10.5f    %10.5f\n", m.real(), m.imag(), qext - qsca, 0.8808, 0.88081);

		m = std::complex<double>(1.342, -0.8);
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%10.5g%-+7.4fi    %10.5f    %10.5f    %10.5f\n", m.real(), m.imag(), qext - qsca, 0.8369, 0.83686);

		m = std::complex<double>(1.342, -1.0);
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%10.5g%-+7.4fi    %10.5f    %10.5f    %10.5f\n", m.real(), m.imag(), qext - qsca, 0.7910, 0.79097);

		printf("\n");
	}


	

	{
		long		    nangles = 0;
		sci::GridData<std::complex<double>, 1> s1;
		sci::GridData<std::complex<double>, 1> s2;
		sci::GridData<double, 1> mu;
		std::complex<double>    m;
		double		    x = 20;
		double		    qext, qsca, qback, g;

		printf("\n***********************************************\n");
		printf("Third mie Test -- van de Hulst page 161\n");
		printf(" x          Qs   Prahl    vdH       Qs*g   Prahl   vdH\n");

		m = std::complex<double>(std::numeric_limits<double>::infinity(), 0.0);
		x = 0.3;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.1f   %7.3f %7.3f %7.3f   %7.3f %7.3f %7.3f\n", x, qsca, 0.028, 0.028, qsca * g, -0.010, -0.011);
		//ez_Mie(x, 0.0, qsca, g);
		//printf("%4.1f   %7.3f %7.3f %7.3f   %7.3f %7.3f %7.3f [ez mie]\n", x, qsca, 0.028, 0.028, qsca * g, -0.10, -0.011);

		x = 1.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.1f   %7.3f %7.3f %7.3f   %7.3f %7.3f %7.3f\n", x, qsca, 2.036, 2.036, qsca * g, -0.384, -0.385);

		x = 1.5;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.1f   %7.3f %7.3f %7.3f   %7.3f %7.3f %7.3f\n", x, qsca, 2.154, 2.155, qsca * g, 0.156, 0.156);

		x = 5.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.1f   %7.3f %7.3f %7.3f   %7.3f %7.3f %7.3f\n", x, qsca, 2.116, 2.116, qsca * g, 0.965, 0.965);

		printf("\n");
	}


	{
		long		    nangles = 0;
		sci::GridData<std::complex<double>, 1> s1;
		sci::GridData<std::complex<double>, 1> s2;
		sci::GridData<double, 1> mu;
		std::complex<double>    m;
		double		    x = 20;
		double		    qext, qsca, qback, g;

		printf("\n***********************************************\n");
		printf("Fourth mie Test -- van de Hulst page 277\n");
		printf(" x          Qs   Prahl    vdH       Qs*g   Prahl   vdH\n");

		m = std::complex<double>(3.41, -1.94);
		x = 1.3;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.1f   %7.3f %7.3f %7.3f   %7.2f %7.2f %7.2f\n", x, qsca, 1.670, 1.669, qsca * g, 0.38, 0.30);

		m = std::complex<double>(7.20, -2.65);
		x = 1.3;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.1f   %7.3f %7.3f %7.3f   %7.2f %7.2f %7.2f\n", x, qsca, 1.861, 1.860, qsca * g, 0.20, 0.31);

		m = std::complex<double>(std::numeric_limits<double>::infinity(), 0.0);
		x = 1.3;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.1f   %7.3f %7.3f %7.3f   %7.2f %7.2f %7.2f\n", x, qsca, 2.265, 2.266, qsca * g, -0.11, -0.05);

		printf("\n");
	}

	{
		long		    nangles = 0;
		sci::GridData<std::complex<double>, 1> s1;
		sci::GridData<std::complex<double>, 1> s2;
		sci::GridData<double, 1> mu;
		std::complex<double>    m;
		double		    x;
		double		    qext, qsca, qback, g;

		printf("\n***********************************************\n");
		printf("Fifth mie Test -- Wiscombe\n");


		printf("\nNon-Absorbing Spheres m=(0.75+0.0i)\n");
		printf("               Calc.     Wiscombe     Calc     Wiscombe\n");
		printf("   x            Qs          Qs          g          g\n");
		m = std::complex<double>(0.75, 0.0);
		x = 0.099;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 0.000007, g, 0.001448);
		//ez_Mie(x, 0.75, qsca, g);
		//printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f [ez mie]\n", x, qsca, 0.000007, g, 0.001448);

		m = std::complex<double>(0.75, 0.0);
		x = 0.101;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 0.000008, g, 0.001507);

		m = std::complex<double>(0.75, 0.0);
		x = 10.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 2.232265, g, 0.896473);

		m = std::complex<double>(0.75, 0.0);
		x = 1000.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 1.997908, g, 0.844944);
		//ez_Mie(x, 0.75, qsca, g);
		//printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f [ez mie]\n", x, qsca, 1.997908, g, 0.844944);



		printf("\nAbsorbing Water Spheres m=(1.33-0.00001i)\n");
		printf("               Calc.     Wiscombe     Calc     Wiscombe\n");
		printf("   x            Qs          Qs          g          g\n");
		m = std::complex<double>(1.33, -0.00001);
		x = 1.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 0.093923, g, 0.184517);

		x = 100.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 2.096594, g, 0.868959);

		x = 10000.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 1.723857, g, 0.907840);



		printf("\nAbsorbing Spheres m=(1.5-i)\n");
		printf("               Calc.     Wiscombe     Calc     Wiscombe\n");
		printf("   x            Qs          Qs          g          g\n");
		m = std::complex<double>(1.5, -1.00);

		x = 0.055;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 0.000011, g, 0.000491);

		x = 0.056;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 0.000012, g, 0.000509);

		x = 1.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 0.6634538, g, 0.192136);

		x = 100.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 1.283697, g, 0.850252);

		x = 10000.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 1.236574, g, 0.846310);



		printf("\n Yet More Absorbing Spheres m=(10-10i)\n");
		printf("               Calc.     Wiscombe     Calc     Wiscombe\n");
		printf("   x            Qs          Qs          g          g\n");
		m = std::complex<double>(10.0, -10.00);

		x = 1.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 2.049405, g, -0.110664);

		x = 100.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 1.836785, g, 0.556215);

		x = 10000.0;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 1.795393, g, 0.548194);



		printf("\nPerfectly Conducting Spheres\n");
		printf("               Calc.     Wiscombe     Calc     Wiscombe\n");
		printf("   x            Qs          Qs          g          g\n");
		m = std::complex<double>(std::numeric_limits<double>::infinity(), 0.0);
		x = 0.099;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 0.000321, g, -0.397357);

		x = 0.101;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 0.000348, g, -0.397262);

		x = 100;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 2.008102, g, 0.500926);
		//ez_Mie(x, 0.0, qsca, g);
		//printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f [ez mie]\n", x, qsca, 2.008102, g, 0.500926);

		x = 10000;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n", x, qsca, 2.000289, g, 0.500070);



		printf("\n");
	}

	

	{
		double		    pi = 3.14159265358979;
		long		    nangles = 7;
		sci::GridData<std::complex<double>, 1> s1;
		sci::GridData<std::complex<double>, 1> s2;
		sci::GridData<double, 1> mu(nangles);
		std::complex<double>    m;
		double		    x;
		double		    qext, qsca, qback, g;
		char form[] = "%7.4f %8.5f%+-8.5fi    %8.5f%+-8.5fi  Calc\n";
		char form2[] = "%7.4f %8.5f%+-8.5fi    %8.5f%+-8.5fi   Wiscombe\n\n";
		long		    i;

		s1 = sci::GridData< std::complex<double>, 1>(nangles);
		s2 = sci::GridData< std::complex<double>, 1>(nangles);

		for (i = 0; i < nangles; i++)
			mu[i] = cos(pi * i / 6.0);

		printf("\n***********************************************\n");
		printf("Sixth mie Test -- Wiscombe\n");
		printf("   angle       S1                    S2         \n");
		x = 1.0;
		m = std::complex<double>(1.5, -1.0);
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);

		printf(form, mu[0], s1[0].real(), s1[0].imag(), s2[0].real(), s2[0].imag());
		printf(form2, mu[0], 5.84080E-01, 1.90515E-01, 5.84080E-01, 1.90515E-01);
		printf(form, mu[1], s1[1].real(), s1[1].imag(), s2[1].real(), s2[1].imag());
		printf(form2, mu[1], 5.65702E-01, 1.87200E-01, 5.00161E-01, 1.45611E-01);
		printf(form, mu[2], s1[2].real(), s1[2].imag(), s2[2].real(), s2[2].imag());
		printf(form2, mu[2], 5.17525E-01, 1.78443E-01, 2.87964E-01, 4.10540E-02);
		printf(form, mu[3], s1[3].real(), s1[3].imag(), s2[3].real(), s2[3].imag());
		printf(form2, mu[3], 4.56340E-01, 1.67167E-01, 3.62285E-02, -6.18265E-02);
		printf(form, mu[4], s1[4].real(), s1[4].imag(), s2[4].real(), s2[4].imag());
		printf(form2, mu[4], 4.00212E-01, 1.56643E-01, -1.74875E-01, -1.22959E-01);
		printf(form, mu[5], s1[5].real(), s1[5].imag(), s2[5].real(), s2[5].imag());
		printf(form2, mu[5], 3.62157E-01, 1.49391E-01, -3.05682E-01, -1.43846E-01);
		printf(form, mu[6], s1[6].real(), s1[6].imag(), s2[6].real(), s2[6].imag());
		printf(form2, mu[6], 3.48844E-01, 1.46829E-01, -3.48844E-01, -1.46829E-01);

		printf("\n");
	}

	

	{
		long		    nangles = 0;
		sci::GridData<std::complex<double>, 1> s1;
		sci::GridData<std::complex<double>, 1> s2;
		sci::GridData<double, 1> mu;
		std::complex<double>    m;
		double		    x = 20;
		double		    qext, qsca, qback, g;

		printf("\n***********************************************\n");
		printf("Small mie Test\n");
		printf("           calc    Wiscombe       calc   Wiscombe\n");
		printf(" X         Qsca    Qsca            g         g\n");

		m = std::complex<double>(0.75, 0.0);
		x = 0.099;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n", x, qext, 0.000007, g, 0.001448);
		x = 0.101;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n", x, qext, 0.000008, g, 0.001507);

		m = std::complex<double>(1.5, -1.0);
		x = 0.055;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("the following line does not agree with the Prahl\nor Wiscombe values for g.\nThe difference is because for this example I use\nthe more stable down recursion (due to stricter\ntests) and Prahl uses less stable up recursion.\nForcing Prahl to use down recursion gives\nmatching answers. Hence I am potentially correct.\n");
		printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n", x, qext, 0.101491, g, 0.000491);
		x = 0.056;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n", x, qext, 0.103347, g, 0.000509);

		m = std::complex<double>(1e-10, -1e10);
		x = 0.099;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n", x, qext, 0.000321, g, -0.397357);
		x = 0.101;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n", x, qext, 0.000348, g, -0.397262);

		m = std::complex<double>(0.0, -1e10);
		x = 0.099;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n", x, qext, 0.000321, g, -0.397357);
		x = 0.101;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n", x, qext, 0.000348, g, -0.397262);

		printf("\n");
	}

	

	{
		long		    nangles = 0;
		sci::GridData<std::complex<double>, 1> s1;
		sci::GridData<std::complex<double>, 1> s2;
		sci::GridData<double, 1> mu;
		std::complex<double>    m;
		double		    x, qext, qsca, qback, g, ref;

		printf("\n***********************************************\n");
		printf("Backscattering Efficiency\n");
		printf("                                   Calc          reference\n");
		printf("    X         m.real()     m.imag()        Qsca          qsca          ratio\n");

		m = std::complex<double>(1.55, 0.0);
		x = 2 * 3.1415926535 * 0.525 / 0.6328;
		ref = 2.92534;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8g\n", x, m.real(), m.imag(), qback, ref, qback / ref);

		m = std::complex<double>(0.0, -1000.0);
		x = 0.099;
		ref = (4.77373E-07 * 4.77373E-07 + 1.45416E-03 * 1.45416E-03) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.2f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);
		x = 0.101;
		ref = (5.37209E-07 * 5.37209E-07 + 1.54399E-03 * 1.54399E-03) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.2f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);
		x = 100;
		ref = (4.35251E+01 * 4.35251E+01 + 2.45587E+01 * 2.45587E+01) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.2f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);

		m = std::complex<double>(0.75, 0.0);
		x = 0.099;
		ref = (1.81756E-08 * 1.81756E-08 + 1.64810E-04 * 1.64810E-04) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);
		x = 0.101;
		ref = (2.04875E-08 * 2.04875E-08 + 1.74965E-04 * 1.74965E-04) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);
		x = 10.0;
		ref = (1.07857E+00 * 1.07857E+00 + 3.60881E-02 * 3.60881E-02) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);
		x = 1000.0;
		ref = (1.70578E+01 * 1.70578E+01 + 4.84251E+02 * 4.84251E+02) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);

		m = std::complex<double>(1.33, -0.00001);
		x = 1.0;
		ref = (2.24362E-02 * 2.24362E-02 + 1.43711E-01 * 1.43711E-01) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);
		x = 100.0;
		ref = (5.65921E+01 * 5.65921E+01 + 4.65097E+01 * 4.65097E+01) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);

		m = std::complex<double>(1.5, -1.0);
		x = 0.055;
		ref = (7.66140E-05 * 7.66140E-05 + 8.33814E-05 * 8.33814E-05) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);
		x = 0.056;
		ref = (8.08721E-05 * 8.08721E-05 + 8.80098E-05 * 8.80098E-05) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);
		x = 1.0;
		ref = (3.48844E-01 * 3.48844E-01 + 1.46829E-01 * 1.46829E-01) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);
		x = 100.0;
		ref = (2.02936E+01 * 2.02936E+01 + 4.38444E+00 * 4.38444E+00) / x / x * 4;
		mie(x, m, mu, s1, s2, qext, qsca, qback, g);
		printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n", x, m.real(), m.imag(), qback, ref, qback / ref);


		printf("\n");
	}
	

}

void speedTest()
{
	srand((unsigned int)time(NULL));

	sci::GridData<std::complex<double>, 1> s1(100);
	sci::GridData<std::complex<double>, 1> s2(100);
	sci::GridData<double, 1> mus (100);
	for (size_t i = 0; i < mus.size(); ++i)
		mus[i] = double(i) / double(mus.size() - 1) * 2.0;
	MiePreAllocator<double, std::complex<double>> preAllocator(mus.size());
	std::complex<double>    refractiveIndex(1.5, 0.00001);
	double x;
	double extinctionEfficiency;
	double scatteringEfficiency;
	double backscatterEfficiency;
	double asymmetryParameter;

	double extinctionEfficiencyMean = 0.0;
	double scatteringEfficiencyMean = 0.0;
	double backscatterEfficiencyMean = 0.0;
	double asymmetryParameterMean = 0.0;

	size_t n = 100000;
	auto start = clock();
	for (size_t i = 0; i < n; ++i)
	{
		x = double(rand()) / double(RAND_MAX) + 1.0;
		mie(x, refractiveIndex, mus, s1, s2, extinctionEfficiency, scatteringEfficiency, backscatterEfficiency, asymmetryParameter, preAllocator);
		extinctionEfficiencyMean += extinctionEfficiency;
		scatteringEfficiencyMean += scatteringEfficiency;
		backscatterEfficiencyMean += backscatterEfficiency;
		asymmetryParameterMean += asymmetryParameter;
	}
	auto end = clock();
	std::cout << "Doing " << n << " Mie calculations with x between 1 and 2 took " << double(end - start) / double(CLOCKS_PER_SEC) << " s\n";

	extinctionEfficiencyMean = 0.0;
	scatteringEfficiencyMean = 0.0;
	backscatterEfficiencyMean = 0.0;
	asymmetryParameterMean = 0.0;

	start = clock();
	for (size_t i = 0; i < n; ++i)
	{
		x = 9.0 * double(rand()) / double(RAND_MAX) + 1.0;
		mie(x, refractiveIndex, mus, s1, s2, extinctionEfficiency, scatteringEfficiency, backscatterEfficiency, asymmetryParameter, preAllocator);
		extinctionEfficiencyMean += extinctionEfficiency;
		scatteringEfficiencyMean += scatteringEfficiency;
		backscatterEfficiencyMean += backscatterEfficiency;
		asymmetryParameterMean += asymmetryParameter;
	}
	end = clock();
	std::cout << "Doing " << n << " Mie calculations with x between 1 and 10 took " << double(end - start) / double(CLOCKS_PER_SEC) << " s\n";

	extinctionEfficiencyMean = 0.0;
	scatteringEfficiencyMean = 0.0;
	backscatterEfficiencyMean = 0.0;
	asymmetryParameterMean = 0.0;

	start = clock();
	for (size_t i = 0; i < n; ++i)
	{
		x = 90.0 * double(rand()) / double(RAND_MAX) + 1.0;
		mie(x, refractiveIndex, mus, s1, s2, extinctionEfficiency, scatteringEfficiency, backscatterEfficiency, asymmetryParameter, preAllocator);
		extinctionEfficiencyMean += extinctionEfficiency;
		scatteringEfficiencyMean += scatteringEfficiency;
		backscatterEfficiencyMean += backscatterEfficiency;
		asymmetryParameterMean += asymmetryParameter;
	}
	end = clock();
	std::cout << "Doing " << n << " Mie calculations with x between 1 and 100 took " << double(end - start) / double(CLOCKS_PER_SEC) << " s\n";

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

	TestsFromPrahl();

    testPrahl();

    testLogarithmicDerivatives();
	speedTest();
}
