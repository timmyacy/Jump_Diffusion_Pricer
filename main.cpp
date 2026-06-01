#include <cmath>
#include <complex>
#include <iostream>

#include "black_scholes.h"
#include "carr_madan.h"

// Standard normal CDF using erfc
double N(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// Black-Scholes closed form for a European call
double bsClosedForm(double S0, double K, double r, double sigma, double T) {
  double d1 = (std::log(S0 / K) + (r + 0.5 * sigma * sigma) * T) /
              (sigma * std::sqrt(T));
  double d2 = d1 - sigma * std::sqrt(T);
  return S0 * N(d1) - K * std::exp(-r * T) * N(d2);
}

int main() {
  // Model parameters
  double r = 0.05;
  double sigma = 0.20;
  double S0 = 100.0;
  double K = 100.0;
  double T = 1.0;

  // Build model and pricer
  BlackScholes bs(r, sigma);
  CarrMadan pricer(bs, r, T, 4096, 0.25, 1.5);

  // Price both ways
  double cm = pricer.price(S0, K);
  double bf = bsClosedForm(S0, K, r, sigma, T);

  std::cout << "Carr-Madan : " << cm << "\n";
  std::cout << "Closed form: " << bf << "\n";
  std::cout << "Difference : " << std::abs(cm - bf) << "\n";

  return 0;
}
