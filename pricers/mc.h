#pragma once
#include "model.h"
#include <random>

class MonteCarlo {
public:
  MonteCarlo(const Model &model, double r, double T, int numPaths,
             int seed = 42);

  double price(double S0, double K) const;

private:
  const Model &model_;
  double r_, T_;
  int numPaths_;
  int seed_;
};
