#include "Empiric.hpp"

Empiric::Empiric(const nstu::vector& sample) : m_sample(sample) {
  if (m_sample.empty()) throw std::invalid_argument("Пустая выборка");
}

const nstu::vector& Empiric::sample() const noexcept { return m_sample; }

uint32_t Empiric::size() const noexcept { return m_sample.size(); }

bool Empiric::empty() const noexcept { return m_sample.empty(); }

ld Empiric::at(uint32_t index) const {
  if (index >= m_sample.size()) throw std::out_of_range("Индекс вне выборки");
  return m_sample[index];
}

void Empiric::changeData(uint32_t index, ld value) {
  if (index >= m_sample.size()) throw std::out_of_range("Индекс вне выборки");
  m_sample[index] = value;
}

void Empiric::setSample(const nstu::vector& sample) {
  if (sample.empty()) throw std::invalid_argument("Пустая выборка");
  m_sample = sample;
}

ld Empiric::mean() const {
  if (m_sample.empty()) throw std::invalid_argument("Пустая выборка");
  ld sum = 0;
  for (const ld x : m_sample) sum += x;
  return sum / m_sample.size();
}

void Empiric::attach(IObserver* observer) {
  if (!observer) throw std::invalid_argument("observer == nullptr");
  for (uint32_t i = 0; i < m_observers.size(); ++i) {
    if (m_observers[i] == observer) return;
  }
  m_observers.push_back(observer);
}

void Empiric::detach(IObserver* observer) {
  if (!observer) return;
  for (uint32_t i = 0; i < m_observers.size(); ++i) {
    if (m_observers[i] == observer) {
      for (uint32_t j = i + 1; j < m_observers.size(); ++j) {
        m_observers[j - 1] = m_observers[j];
      }
      m_observers.pop_back();
      return;
    }
  }
}

void Empiric::notify() {
  for (uint32_t i = 0; i < m_observers.size(); ++i) {
    if (m_observers[i]) {
      m_observers[i]->update();
    }
  }
}
