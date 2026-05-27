#pragma once
#include "model.h"
#include <complex>

class BlackScholes : public Model {
public:
  BlackScholes(double r, double sigma);

  std::complex<double> characteristicFunction(std::complex<double> u,
                                              double T) const override;

  double simulate(double S0, double T, std::mt19937 &rng) const override;

private:
  double r_;
  double sigma_;
};
