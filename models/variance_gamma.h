#pragma once
#include "model.h"
#include <complex>

class VarianceGamma : public Model {
public:
  VarianceGamma(double sigma, double nu, double theta, double r,
                double q = 0.0);

  std::complex<double> characteristicFunction(std::complex<double> u,
                                              double T) const override;

  double simulate(double S0, double T, std::mt19937 &rng) const override;

  // Getters
  double sigma() const { return sigma_; }
  double nu() const { return nu_; }
  double theta() const { return theta_; }
  double omega() const { return omega_; }

private:
  double sigma_;
  double nu_;
  double theta_;
  double r_;
  double q_;
  double omega_;
};
