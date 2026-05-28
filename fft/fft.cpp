#include "fft.h"
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>
namespace fft {
std::vector<std::complex<double>> fft(std::vector<std::complex<double>> x) {
  std::vector<std::complex<double>> even;
  std::vector<std::complex<double>> odd;
  if (x.size() == 1) {
    return x;
  }
  for (size_t i = 0; i < x.size(); i++) {
    if (i % 2 == 0)
      even.push_back(x[i]);
    else
      odd.push_back(x[i]);
  };
  auto f_even = fft(even);
  auto f_odd = fft(odd);

  std::vector<std::complex<double>> output(x.size());
  for (size_t k = 0; k < x.size() / 2; k++) {

    std::complex<double> W =
        std::exp(std::complex<double>(0, -2.0 * M_PI * k / x.size()));

    output[k] = f_even[k] + W * f_odd[k];

    output[k + x.size() / 2] = f_even[k] - W * f_odd[k];
  }

  return output;
}
} // namespace fft
