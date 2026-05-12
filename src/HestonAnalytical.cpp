#include "HestonAnalytical.h"


Complex HestonAnalytical::characteristicFunction(int j, double x, double v, double tau, double phi)
{
    Complex I(0.0, 1.0); 
    double a = kappa * theta;

    Complex b_j;
    double w_j;

    if (j == 1) {
        b_j = kappa - rho * sigma;
        w_j = 0.5;
    }
    else { // j == 2
        b_j = kappa;
        w_j = -0.5;
    }

    
    Complex d_inner = std::pow(I * rho * sigma * phi - b_j, 2.0) -
        std::pow(sigma, 2.0) * (2.0 * w_j * I * phi - phi * phi);
    Complex d_j = std::sqrt(d_inner);

    Complex c_j = (b_j - rho * sigma * I * phi - d_j) /
        (b_j - rho * sigma * I * phi + d_j);


    Complex C_j = r * I * phi * tau +
        (a / std::pow(sigma, 2.0)) * ((b_j - rho * sigma * I * phi - d_j) * tau -
            2.0 * std::log((1.0 - c_j * std::exp(-d_j * tau)) / (1.0 - c_j)));

    Complex D_j = ((b_j - rho * sigma * I * phi - d_j) / std::pow(sigma, 2.0)) *
        ((1.0 - std::exp(-d_j * tau)) / (1.0 - c_j * std::exp(-d_j * tau)));

    return std::exp(C_j + D_j * v + I * phi * x);
}


double HestonAnalytical::integrand(int j, double S, double v, double tau, double E, double phi) {
    Complex I(0.0, 1.0); // Imaginarna jednotka
    double x = std::log(S);

    Complex f_j = characteristicFunction(j, x, v, tau, phi);
    Complex numerator = std::exp(-I * phi * std::log(E)) * f_j;
    Complex denominator = I * phi;

    return std::real(numerator / denominator);
}

//P_j (Gauss-Legendre)

double HestonAnalytical::calculateP(int j, double S, double v, double tau, double E) {
    // 20-bodové Gauss-Legendre uzly (na intervale -1 až 1)
    const double nodes[20] = {
        -0.9931285991850949, -0.9639719272779138, -0.9122344282513259,
        -0.8391169718222188, -0.7463319064601508, -0.6360536807265150,
        -0.5108670019508271, -0.3737060887154196, -0.2277858511416451,
        -0.07652652113349733, 0.07652652113349733, 0.2277858511416451,
        0.3737060887154196, 0.5108670019508271, 0.6360536807265150,
        0.7463319064601508, 0.8391169718222188, 0.9122344282513259,
        0.9639719272779138, 0.9931285991850949
    };

    // 20-bodové Gauss-Legendre váhy
    const double weights[20] = {
        0.017614007139152, 0.040601429800387, 0.062672048334109, 0.083276741576705, 0.101930119817240,
        0.118194531961518, 0.131688638449177, 0.142096109318382, 0.149172986472604, 0.152753387130726,
        0.152753387130726, 0.149172986472604, 0.142096109318382, 0.131688638449177, 0.118194531961518,
        0.101930119817240, 0.083276741576705, 0.062672048334109, 0.040601429800387, 0.017614007139152
    };

    double phi_max = 200.0;
    int num_subintervals = 100;
    double sub_length = phi_max / num_subintervals; 

    double integral_sum = 0.0;

    for (int step = 0; step < num_subintervals; ++step) {
        double a = step * sub_length;
        double b = (step + 1) * sub_length;

        double half_width = (b - a) / 2.0;
        double center = (a + b) / 2.0;

        // Vnútorný cyklus: 20-bodová Gaussova kvadratúra na aktuálnom podintervale
        for (int k = 0; k < 20; ++k) {
            double phi_k = half_width * nodes[k] + center;
            integral_sum += (half_width * weights[k]) * integrand(j, S, v, tau, E, phi_k);
        }
    }

    double pi = 3.14159265358979323846;
    return 0.5 + (1.0 / pi) * integral_sum;
}