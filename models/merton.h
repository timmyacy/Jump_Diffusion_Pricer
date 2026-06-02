#pragma once
#include "model.h"
#include <cmath>

class Merton : public Model {
public:
  Merton(double S0, double sigma, double r, double lambda, double muJ,
         double sigmaJ);

  std::complex<double> characteristicFunction(std::complex<double> u,
                                              double T) const override;

  double simulate(double S0, double T, std::mt19937 &rng) const override;

  double semiAnalyticPrice(double S0, double K, double T,
                           int nTerms = 20) const;

private:
  double S0_, sigma_, r_, lambda_, muJ_, sigmaJ_, kBar_;

  double bsCall(double S, double K, double r, double sigma, double T) const;
};
