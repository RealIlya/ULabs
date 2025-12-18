#pragma once

#include <cstdint>

#include "../lab2/header.h"
#include "Empiric.hpp"
#include "IObserver.hpp"

class Estimate final : public IObserver {
 private:
  Empiric& m_data;
  ld m_sigma;
  ld m_c;
  int m_starts;
  ld m_mu;

  void estimate();
  ld weightZ(ld z) const noexcept;
  ld rhoZ(ld z) const;

 public:
  Estimate(Empiric& data, ld sigma, ld c, int starts = 100);
  ~Estimate() override;

  Estimate(const Estimate&) = delete;
  Estimate& operator=(const Estimate&) = delete;
  Estimate(Estimate&&) = delete;
  Estimate& operator=(Estimate&&) = delete;

  void update() override;

  ld mu() const noexcept;
  ld sigma() const noexcept;
  ld c() const noexcept;
  int starts() const noexcept;

  void setSigma(ld sigma);
  void setC(ld c);
  void setStarts(int starts);

  ld weightForObservation(ld x) const;
  nstu::vector weights() const;
};

