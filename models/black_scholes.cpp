#include "black_scholes.h"
#include <cmath>
#include <random>
using namespace std::complex_literals;
BlackScholes::BlackScholes(double r, double sigma) : r_(r), sigma_(sigma) {}

std::complex<double>
BlackScholes::characteristicFunction(std::complex<double> u, double T) const {

  std::complex<double> i = 1i;
  return std::exp(i * u * (r_ - 0.5 * sigma_ * sigma_) * T -
                  0.5 * sigma_ * sigma_ * u * u * T);
}

double BlackScholes::simulate(double S0, double T, std::mt19937 &rng) const {
  std::normal_distribution<double> norm(0.0, 1.0);
  double z = norm(rng);
  return S0 *
         std::exp((r_ - 0.5 * sigma_ * sigma_) * T + sigma_ * std::sqrt(T) * z);
}
