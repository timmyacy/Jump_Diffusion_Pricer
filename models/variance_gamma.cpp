#include "variance_gamma.h"
#include <cmath>
#include <random>
#include <stdexcept>

using cd = std::complex<double>;
VarianceGamma::VarianceGamma(double sigma, double nu, double theta, double r,
                             double q)
    : sigma_(sigma), nu_(nu), theta_(theta), r_(r), q_(q) {

  double arg = 1.0 - theta_ * nu_ - 0.5 * sigma_ * sigma_ * nu_;

  omega_ = std::log(arg) / nu_;
}

cd VarianceGamma::characteristicFunction(cd u, double T) const {

  cd drift = cd(0.0, 1.0) * u * (r_ - q_ + omega_) * T;

  cd inner = 1.0 - cd(0.0, 1.0) * u * theta_ * nu_ +
             0.5 * u * u * sigma_ * sigma_ * nu_;
  cd vg_part = std::pow(inner, -T / nu_);

  return std::exp(drift) * vg_part;
}

double VarianceGamma::simulate(double S0, double T, std::mt19937 &rng) const {
  std::normal_distribution<double> normal(0.0, 1.0);
  std::gamma_distribution<double> gamma_dist(T / nu_, nu_);
  double dG = gamma_dist(rng);
  double dW = normal(rng);
  return S0 * std::exp((r_ - q_ + omega_) * T + theta_ * dG +
                       sigma_ * std::sqrt(dG) * dW);
}
