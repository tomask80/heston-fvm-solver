#pragma once

struct HestonParams
{
    double kappa   = 0.0; // mean reversion speed
    double theta   = 0.0; // long-term variance
    double sigma   = 0.0; // volatility of volatility
    double rho     = 0.0; // correlation between asset and variance
    double r       = 0.0; // risk-free rate
    double T       = 0.0; // time to maturity
    double X_l     = 0.0, X_r = 0.0; // spatial domain limits
    double y_maxl  = 0.0; // maximum variance level
};

