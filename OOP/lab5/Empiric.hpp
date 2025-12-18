#pragma once

#include <cstdint>
#include <stdexcept>

#include "../lab2/header.h"
#include "IObserver.hpp"

class Empiric {
 private:
  nstu::vector m_sample;
  nstu::vector_t<IObserver*> m_observers;

 public:
  Empiric() = default;
  explicit Empiric(const nstu::vector& sample);

  Empiric(const Empiric&) = delete;
  Empiric& operator=(const Empiric&) = delete;
  Empiric(Empiric&&) = delete;
  Empiric& operator=(Empiric&&) = delete;
  ~Empiric() = default;

  const nstu::vector& sample() const noexcept;
  uint32_t size() const noexcept;
  bool empty() const noexcept;

  ld at(uint32_t index) const;
  void changeData(uint32_t index, ld value);
  void setSample(const nstu::vector& sample);

  ld mean() const;

  void attach(IObserver* observer);
  void detach(IObserver* observer);
  void notify();
};
