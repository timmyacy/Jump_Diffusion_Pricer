#include "carr_madan.h"
#include "fft.h"
CarrMadan::CarrMadan(const Model &model, double r, double T, int N, double eta,
                     double alpha)
    : model_(model), r_(r), T_(T), N_(N), eta_(eta), alpha_(alpha) {};

double CarrMadan::price(double S0, double K) const {

  double b = 2.0 * M_PI / eta_;
  std::vector<std::complex<double>> input(N_);
  for (int i = 0; i < N_; i++) {

    double v = i * eta_;

    // shift v down the imaginary axis by (alpha+1)
    std::complex<double> u(v, -(alpha_ + 1.0));

    // Evaluate the CF at this argument
    std::complex<double> phi = model_.characteristicFunction(u, T_);

    std::complex<double> denom =
        std::complex<double>(alpha_, v) * std::complex<double>(alpha_ + 1.0, v);

    // Psi(v): the Carr-Madan integrand
    std::complex<double> psi =
        std::exp(std::complex<double>(-r_ * T_, 0.0)) * phi / denom;

    // Simpson weights: 1/3 on endpoints, alternating 4/3 and 2/3 inside
    double w;
    if (i == 0 || i == N_ - 1)
      w = 1.0 / 3.0;
    else if (i % 2 == 1)
      w = 4.0 / 3.0;
    else
      w = 2.0 / 3.0;

    // Phase shift: centres the log-strike grid around zero
    std::complex<double> phase =
        std::exp(std::complex<double>(0.0, v * b / 2.0));

    input[i] = phase * psi * w * eta_;
  }

  auto output = fft::fft(input);

  double lambda = (2.0 * M_PI) / (N_ * eta_);
  double k = std::log(K / S0);

  int idx = static_cast<int>(std::round((k + b / 2.0) / lambda));

  if (idx < 0 || idx >= N_)
    return 0.0;

  double call = S0 * (std::exp(-alpha_ * k) / M_PI) * output[idx].real();
  return call;
};
