#pragma once
#include "../models/model.h"
#include <vector>

class CarrMadan {
public:
  CarrMadan(const Model &model, double r, double T, int N, double eta,
            double alpha);

  double price(double S0, double K) const;

private:
  const Model &model_;
  double r_, T_, eta_, alpha_;
  int N_;
};
