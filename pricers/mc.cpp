#include "mc.h"
#include <cmath>

MonteCarlo::MonteCarlo(const Model &model, double r, double T, int numPaths,
                       int seed)
    : model_(model), r_(r), T_(T), numPaths_(numPaths), seed_(seed) {}

double MonteCarlo::price(double S0, double K) const {
  std::mt19937 rng(seed_);
  double sumPayoff = 0.0;

  for (int i = 0; i < numPaths_; ++i) {
    double ST = model_.simulate(S0, T_, rng);
    double payoff = std::max(ST - K, 0.0);
    sumPayoff += payoff;
  }

  return std::exp(-r_ * T_) * sumPayoff / numPaths_;
}
