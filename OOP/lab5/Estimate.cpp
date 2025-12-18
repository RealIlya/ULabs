#include "Estimate.hpp"

#include <cmath>
#include <stdexcept>

namespace {

constexpr int kMaxIters = 200;
constexpr ld kEps = 1e-12L;
constexpr ld kInf = 1e300L;

inline uint32_t min_u32(uint32_t a, uint32_t b) { return (a < b) ? a : b; }

inline void sort_nstu_vector(nstu::vector& v) {
  for (uint32_t i = 1; i < v.size(); ++i) {
    ld key = v[i];
    uint32_t j = i;
    while (j > 0 && v[j - 1] > key) {
      v[j] = v[j - 1];
      --j;
    }
    v[j] = key;
  }
}

}  // namespace

Estimate::Estimate(Empiric& data, ld sigma, ld c, int starts)
    : m_data(data), m_sigma(sigma), m_c(c), m_starts(starts), m_mu(0) {
  if (m_data.empty()) throw std::invalid_argument("Пустая выборка");
  if (!(m_sigma > 0)) throw std::invalid_argument("sigma должен быть > 0");
  if (!(m_c > 0)) throw std::invalid_argument("c должен быть > 0");
  if (m_starts < 1) throw std::invalid_argument("starts должен быть >= 1");

  estimate();
  m_data.attach(this);
}

Estimate::~Estimate() { m_data.detach(this); }

void Estimate::update() { estimate(); }

ld Estimate::mu() const noexcept { return m_mu; }
ld Estimate::sigma() const noexcept { return m_sigma; }
ld Estimate::c() const noexcept { return m_c; }
int Estimate::starts() const noexcept { return m_starts; }

void Estimate::setSigma(ld sigma) {
  if (!(sigma > 0)) throw std::invalid_argument("sigma должен быть > 0");
  m_sigma = sigma;
  estimate();
}

void Estimate::setC(ld c) {
  if (!(c > 0)) throw std::invalid_argument("c должен быть > 0");
  m_c = c;
  estimate();
}

void Estimate::setStarts(int starts) {
  if (starts < 1) throw std::invalid_argument("starts должен быть >= 1");
  m_starts = starts;
  estimate();
}

ld Estimate::weightZ(ld z) const noexcept {
  const ld c2 = m_c * m_c;
  return c2 / (c2 + z * z);
}

ld Estimate::rhoZ(ld z) const {
  const ld c2 = m_c * m_c;
  return std::log(1 + (z * z) / c2);
}

void Estimate::estimate() {
  const auto& s = m_data.sample();
  if (s.empty()) throw std::invalid_argument("Пустая выборка");

  nstu::vector xs = s;
  sort_nstu_vector(xs);

  const int total_starts =
      (m_starts < static_cast<int>(xs.size())) ? m_starts : xs.size();
  ld best_mu = xs[xs.size() / 2];
  ld best_obj = kInf;

  for (int k = 0; k < total_starts; ++k) {
    const uint32_t idx = static_cast<uint32_t>(
        ((static_cast<long double>(k) + 0.5L) * xs.size()) / total_starts);
    ld mu = xs[min_u32(idx, xs.size() - 1)];

    for (int iter = 0; iter < kMaxIters; ++iter) {
      ld w_sum = 0;
      ld xw_sum = 0;
      for (const ld x : s) {
        const ld z = (x - mu) / m_sigma;
        const ld w = weightZ(z);
        w_sum += w;
        xw_sum += x * w;
      }
      if (!(w_sum > 0)) break;
      const ld next_mu = xw_sum / w_sum;
      const ld diff = std::fabs(next_mu - mu);
      mu = next_mu;
      if (diff <= kEps * (1 + std::fabs(mu))) break;
    }

    ld obj = 0;
    for (const ld x : s) {
      const ld z = (x - mu) / m_sigma;
      obj += rhoZ(z);
    }

    if (obj < best_obj) {
      best_obj = obj;
      best_mu = mu;
    }
  }

  m_mu = best_mu;
}

ld Estimate::weightForObservation(ld x) const {
  const ld z = (x - m_mu) / m_sigma;
  return weightZ(z);
}

nstu::vector Estimate::weights() const {
  const auto& s = m_data.sample();
  nstu::vector w(s.size());
  for (uint32_t i = 0; i < s.size(); ++i) {
    w[i] = weightForObservation(s[i]);
  }
  return w;
}
