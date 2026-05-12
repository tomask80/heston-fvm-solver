#pragma once

#include <complex>

using Complex = std::complex<double>;

class HestonAnalytical
{
public:
	double kappa, theta, sigma, rho, r;
	HestonAnalytical(double k, double th, double sig, double rh, double r_val)
		: kappa(k), theta(th), sigma(sig), rho(rh), r(r_val) {
	}
	double getCallPrice(double S, double v, double tau, double E) {
		double P1 = calculateP(1, S, v, tau, E);
		double P2 = calculateP(2, S, v, tau, E);

		return S * P1 - E * std::exp(-r * tau) * P2;
	}

private:
	// charaktreristicke funkice
	Complex characteristicFunction(int j, double x, double v, double tau, double phi);
	double integrand(int j, double S, double v, double tau, double E, double phi);
	double calculateP(int j, double S, double v, double tau, double E);


};

