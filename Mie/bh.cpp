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

double cabs(std::complex<double> x)
{
	return std::sqrt(std::norm(x));
}

void bhMie(double x, std::complex<double> refRel, size_t nAng, std::vector<std::complex<double>>& s1, std::vector<std::complex<double>>& s2, double& qExt, double& qSca, double& qBack)
{
	std::vector<double> amu(100);
	std::vector<double> theta(100);
	std::vector<double> pi(100);
	std::vector<double> tau(100);
	std::vector<double> pi0(100);
	std::vector<double> pi1(100);

	std::vector<std::complex<double>>D(3000);
	std::complex<double> y;
	std::complex<double> xi;
	std::complex<double> xi0;
	std::complex<double> xi1;
	std::complex<double> an;
	std::complex<double> bn;
	s1.resize(200);
	s2.resize(200);

	double psi0;
	double psi1;
	double psi;
	double dn;
	double dx;

	dx = x;
	y = x * refRel;
	size_t xStop = size_t(std::round(x + 4 * std::pow(x, 0.3333) + 2.0));
	size_t nStop = xStop;
	double yMod = cabs(y);
	size_t nMx = std::max(xStop, size_t(std::round(yMod))) + 15;
	double dAng = sci::m_pi_2 / double(nAng - 1);
	for (size_t j = 0; j < nAng; ++j)
	{
		theta[j] = double(j - 1) * dAng;
		amu[j] = std::cos(theta[j]);
	}

	D[nMx-1] = 0.0; //line 100
	size_t nn = nMx - 1;
	for (size_t j = 0; j < nn; ++j)
	{
		size_t n = j + 1;
		double rn = double(nMx - n + 1);
		D[nMx - n - 1] = rn / y - (1.0 / (D[nMx - n] + rn / y));
	}

	for (size_t j = 0; j < nAng; ++j)
	{
		pi0[j] = 0.0;
		pi1[j] = 1.0;
	}

	nn = 2 * nAng - 1;
	for (size_t j = 0; j < nn; ++j)
	{
		s1[j] = 0.0;
		s2[j] = 0.0;
	}

	//line 116
	psi0 = std::cos(dx);
	psi1 = std::sin(dx);
	double chi0 = -std::sin(x);
	double chi1 = cos(x);
	double apsi0 = psi0;
	double apsi1 = psi1;
	xi0 = std::complex(apsi0, -chi0);
	xi1 = std::complex(apsi1, -chi1);
	qSca = 0.0;
	size_t n = 1;
	do //200, line126
	{
		dn = double(n);
		double rn = double(n);
		double fn = (2.0 * rn + 1.0) / (rn * (rn + 1.0));
		psi = (2.0 * dn - 1.0) * psi1 / dx - psi0;
		double apsi = psi;
		double chi = (2.0 * rn - 1.0) * chi1 / x - chi0;
		xi = std::complex(apsi, -chi);
		an = (D[n - 1] / refRel + rn / x) * apsi - apsi1;
		an = an / ((D[n - 1] / refRel + rn / x) * xi - xi1);
		bn = (refRel * D[n - 1] + rn / x) * apsi - apsi1;
		bn = bn / ((refRel * D[n - 1] + rn / x) * xi - xi1);
		qSca = qSca + (2.0 * rn + 1.0) * (std::norm(an) + std::norm(bn));
		for (size_t j = 0; j < nAng; ++j) //line138
		{
			size_t jj = 2 * nAng - j;
			pi[j] = pi1[j];
			tau[j] = rn * amu[j] * pi[j] - (rn + 1.0) * pi0[j];
			double p = std::pow(-1.0, n - 1);
			s1[j] = s1[j] + fn * (an * pi[j] + bn * tau[j]);
			double t = std::pow(-1.0, n);
			s2[j] = s2[j] + fn * (an * tau[j] + bn * pi[j]);
			//146
			if (j != jj)
			{
				s1[jj] = s1[jj] + fn * (an * pi[j] * p + bn * tau[j] * t);
				s2[jj] = s2[jj] + fn * (an * tau[j] * t + bn * pi[j] * p);
			}
		} //789, line 149
		psi0 = psi1;
		psi1 = psi;
		apsi1 = psi1;
		chi0 = chi1;
		chi1 = chi;
		xi1 = std::complex(apsi1, -chi1);
		n = n + 1;
		rn = double(n);
		for (size_t j = 0; j < nAng; ++j)
		{
			pi1[j] = ((2.0 * rn - 1.0) / (rn - 1.0)) * amu[j] * pi[j];
			pi1[j] = pi1[j] - rn * pi0[j] / (rn - 1.0);
			pi0[j] = pi1[j];
		}
	} while (n - 1 - nStop < 0);
	qSca = (2.0 / (x * x)) * qSca;
	qExt = (4.0 / (x * x)) * s1[1].real();
	qBack = (4.0 / (x * x)) * std::norm(s1[2 * nAng - 1]);
}

void callBH()
{
	double refMed = 1.0;
	double refRe = 1.55;
	double refIm = 0.0;
	std::complex<double> refRel = std::complex<double>(refRe, refIm) / refMed;
	double rad = 0.525;
	double wavel = 0.6328;
	double x = sci::m_pi * 2.0 * rad * refMed / wavel;

	size_t nAng = 11;
	double dAng = sci::m_pi_2 / double(nAng - 1);

	std::vector<std::complex<double>> s1;
	std::vector<std::complex<double>> s2;
	double qExt;
	double qSc;
	double qBack;
	bhMie(x, refRel, nAng, s1, s2, qExt, qSc, qBack);

	std::cout << "qSca = " << qSc << "\nqExt = " << qExt << "\nqBack = " << qBack << "\n";
	std::cout << "\nangle\ts11\tpol\ts33\ts34\n\n";

	double s11Nor = 0.5 * std::norm(s2[0]) + std::norm(s1[0]);
	size_t nAn = 2 * nAng - 1;
	for (size_t j = 0; j < nAn; ++j)
	{
		double aj = double(j);
		double s11 = 0.5 * std::norm(s2[j]);
		s11 = s11 + 0.5 * std::norm(s1[j]);
		double s12 = 0.5 * std::norm(s2[j]);
		s12 = s12 - 0.5 * std::norm(s1[j]);
		double pol = -s12 / s11;
		double s33 = (s2[j] * std::conj(s1[j])).real();
		s33 = s33 / s11;
		double s34 = (s2[j] * std::conj(s1[j])).imag();
		s34 = s34 / s11;
		s11 = s11 / s11Nor;
		double ang = dAng * aj * 57.2958;
		std::cout << ang << "\t" << s11 << "\t" << pol << "\t" << s33 << "\t" << s34 << std::endl;
	}
}