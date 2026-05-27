#pragma once
#include <complex>
#include <random>

class Model {
public:
  virtual std::complex<double> characteristicFunction(std::complex<double> u,
                                                      double T) const = 0;

  virtual double simulate(double S0, double T, std::mt19937 &rng) const = 0;

  virtual ~Model() = default;
};
