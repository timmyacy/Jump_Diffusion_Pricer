#include "merton.h"
#include <cmath>
#include <complex>
#include <random>

Merton::Merton(double S0, double sigma, double r, double lambda, double muJ,
               double sigmaJ)
    : S0_(S0), sigma_(sigma), r_(r), lambda_(lambda), muJ_(muJ),
      sigmaJ_(sigmaJ) {
  kBar_ = std::exp(muJ_ + 0.5 * sigmaJ_ * sigmaJ_) - 1.0;
}

std::complex<double> Merton::characteristicFunction(std::complex<double> u,
                                                    double T) const {
  using namespace std::complex_literals;
  const std::complex<double> i(0.0, 1.0);

  std::complex<double> diffusion =
      i * u * (r_ - 0.5 * sigma_ * sigma_ - lambda_ * kBar_) -
      0.5 * u * u * sigma_ * sigma_;

  std::complex<double> jumpCF =
      std::exp(i * u * muJ_ - 0.5 * u * u * sigmaJ_ * sigmaJ_) - 1.0;

  std::complex<double> exponent = T * (diffusion + lambda_ * jumpCF);

  return std::exp(i * u * std::log(S0_) + exponent);
}

double Merton::simulate(double S0, double T, std::mt19937 &rng) const {
  std::normal_distribution<double> stdNorm(0.0, 1.0);
  std::poisson_distribution<int> poisson(lambda_ * T);

  int N = poisson(rng);

  // Sum of N normal jump sizes
  double jumpSum = 0.0;
  std::normal_distribution<double> jumpDist(muJ_, sigmaJ_);
  for (int j = 0; j < N; ++j)
    jumpSum += jumpDist(rng);

  // Diffusion component
  double Z = stdNorm(rng);
  double logS = std::log(S0) +
                (r_ - 0.5 * sigma_ * sigma_ - lambda_ * kBar_) * T +
                sigma_ * std::sqrt(T) * Z + jumpSum;

  return std::exp(logS);
}

double Merton::bsCall(double S, double K, double r, double sigma,
                      double T) const {
  if (sigma <= 0.0 || T <= 0.0)
    return std::max(S - K * std::exp(-r * T), 0.0);

  double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) /
              (sigma * std::sqrt(T));
  double d2 = d1 - sigma * std::sqrt(T);

  auto N = [](double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); };

  return S * N(d1) - K * std::exp(-r * T) * N(d2);
}

double Merton::semiAnalyticPrice(double S0, double K, double T,
                                 int nTerms) const {
  double lambdaPrime = lambda_ * (1.0 + kBar_);
  double price = 0.0;
  double logFactorial = 0.0;
  for (int n = 0; n < nTerms; ++n) {
    if (n > 0)
      logFactorial += std::log(static_cast<double>(n));

    double weight = std::exp(-lambdaPrime * T + n * std::log(lambdaPrime * T) -
                             logFactorial);

    double sigmaN = std::sqrt(sigma_ * sigma_ + n * sigmaJ_ * sigmaJ_ / T);
    double rN = r_ - lambda_ * kBar_ + n * (muJ_ + 0.5 * sigmaJ_ * sigmaJ_) / T;

    price += weight * bsCall(S0, K, rN, sigmaN, T);
  }

  return price;
}
