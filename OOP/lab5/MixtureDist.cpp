#include <cmath>
#include <fstream>
#include <memory>

#include "MixtureDist.hpp"

template <typename dist1_t, typename dist2_t>
MixtureDist<dist1_t, dist2_t>::MixtureDist() : m_p(0.5), m_dist1(), m_dist2() {}

template <typename dist1_t, typename dist2_t>
MixtureDist<dist1_t, dist2_t>::MixtureDist(ld p, const dist1_t& dist1,
                                           const dist2_t& dist2)
    : m_p(p), m_dist1(dist1), m_dist2(dist2) {
  if (p < 0 || p > 1)
    throw std::invalid_argument(
        "Параметр смеси должен быть в диапазоне [0, 1]");
}

template <typename dist1_t, typename dist2_t>
void MixtureDist<dist1_t, dist2_t>::save(const std::string& file_name) const {
  std::ofstream file(file_name);
  if (!file.is_open()) throw std::runtime_error("Не удалось открыть файл");
  save(file);
  file.close();
}

template <typename dist1_t, typename dist2_t>
void MixtureDist<dist1_t, dist2_t>::save(std::ofstream& file) const noexcept {
  file << "MixtureDist ";
  file << m_p << " ";
  m_dist1.save(file);
  m_dist2.save(file);
}

template <typename dist1_t, typename dist2_t>
void MixtureDist<dist1_t, dist2_t>::load(const std::string& file_name) {
  std::ifstream file(file_name);
  if (!file.is_open()) throw std::runtime_error("Не удалось открыть файл");
  load(file);
  file.close();
}

template <typename dist1_t, typename dist2_t>
void MixtureDist<dist1_t, dist2_t>::load(std::ifstream& file) {
  std::string d_name;
  file >> d_name;
  if (d_name == "MixtureDist") {
    file >> m_p;
    if (m_p < 0 || m_p > 1)
      throw std::invalid_argument(
          "Параметр смеси должен быть в диапазоне [0, 1]");
    m_dist1.load(file);
    m_dist2.load(file);
  } else
    throw std::runtime_error("Неверный тип распределения");
}

template <typename dist1_t, typename dist2_t>
MixtureDist<dist1_t, dist2_t>::MixtureDist(const MixtureDist& other)
    : m_p(other.m_p), m_dist1(other.m_dist1), m_dist2(other.m_dist2) {}

template <typename dist1_t, typename dist2_t>
MixtureDist<dist1_t, dist2_t>& MixtureDist<dist1_t, dist2_t>::operator=(
    const MixtureDist& other) {
  if (this != &other) {
    m_p = other.m_p;
    m_dist1 = other.m_dist1;
    m_dist2 = other.m_dist2;
  }
  return *this;
}

template <typename dist1_t, typename dist2_t>
MixtureDist<dist1_t, dist2_t>::MixtureDist(MixtureDist&& other) noexcept
    : m_p(other.m_p),
      m_dist1(std::move(other.m_dist1)),
      m_dist2(std::move(other.m_dist2)) {
  other.m_p = 0;
}
template <typename dist1_t, typename dist2_t>
MixtureDist<dist1_t, dist2_t>& MixtureDist<dist1_t, dist2_t>::operator=(
    MixtureDist&& other) noexcept {
  if (this != &other) {
    m_p = other.m_p;
    m_dist1 = std::move(other.m_dist1);
    m_dist2 = std::move(other.m_dist2);
    other.m_p = 0;
  }
  return *this;
}

template <typename dist1_t, typename dist2_t>
dist1_t& MixtureDist<dist1_t, dist2_t>::component1() noexcept {
  return m_dist1;
}
template <typename dist1_t, typename dist2_t>
const dist1_t& MixtureDist<dist1_t, dist2_t>::component1() const noexcept {
  return m_dist1;
}
template <typename dist1_t, typename dist2_t>
dist2_t& MixtureDist<dist1_t, dist2_t>::component2() noexcept {
  return m_dist2;
}
template <typename dist1_t, typename dist2_t>
const dist2_t& MixtureDist<dist1_t, dist2_t>::component2() const noexcept {
  return m_dist2;
}

template <typename dist1_t, typename dist2_t>
ld MixtureDist<dist1_t, dist2_t>::getP() const noexcept {
  return m_p;
}

template <typename dist1_t, typename dist2_t>
void MixtureDist<dist1_t, dist2_t>::setP(ld p) {
  if (p < 0 || p > 1)
    throw std::invalid_argument(
        "Параметр смеси должен быть в диапазоне [0, 1]");
  else
    m_p = p;
}

template <typename dist1_t, typename dist2_t>
ld MixtureDist<dist1_t, dist2_t>::density(ld x) const noexcept {
  return (1 - m_p) * m_dist1.density(x) + m_p * m_dist2.density(x);
}

template <typename dist1_t, typename dist2_t>
nstu::vector MixtureDist<dist1_t, dist2_t>::density(
    const nstu::vector& x) const {
  nstu::vector result(x);
  for (uint32_t i = 0; i < result.size(); ++i) {
    result[i] = (1 - m_p) * m_dist1.density(x[i]) + m_p * m_dist2.density(x[i]);
  }
  return result;
}

template <typename dist1_t, typename dist2_t>
ld MixtureDist<dist1_t, dist2_t>::M() const noexcept {
  return (1 - m_p) * m_dist1.M() + m_p * m_dist2.M();
}

template <typename dist1_t, typename dist2_t>
ld MixtureDist<dist1_t, dist2_t>::D() const noexcept {
  ld m1 = m_dist1.M(), m2 = m_dist2.M();
  return (1 - m_p) * (std::pow(m1, 2) + m_dist1.D()) +
         m_p * (std::pow(m2, 2) + m_dist2.D()) - std::pow(M(), 2);
}

template <typename dist1_t, typename dist2_t>
ld MixtureDist<dist1_t, dist2_t>::G1() const noexcept {
  ld m1 = m_dist1.M(), m2 = m_dist2.M(), d1 = m_dist1.D(), d2 = m_dist2.D();
  ld m = M();
  return ((1 - m_p) * (std::pow(m1 - m, 3) + 3 * (m1 - m) * d1 +
                       std::pow(d1, 1.5) * m_dist1.G1()) +
          m_p * (std::pow(m2 - m, 3) + 3 * (m2 - m) * d2 +
                 std::pow(d2, 1.5) * m_dist2.G1())) /
         std::pow(D(), 1.5);
}

template <typename dist1_t, typename dist2_t>
ld MixtureDist<dist1_t, dist2_t>::G2() const noexcept {
  ld m1 = m_dist1.M(), m2 = m_dist2.M(), d1 = m_dist1.D(), d2 = m_dist2.D();
  ld m = M();
  return ((1 - m_p) * (std::pow(m1 - m, 4) + 6 * std::pow(m1 - m, 2) * d1 +
                       4 * (m1 - m) * std::pow(d1, 1.5) * m_dist1.G1() +
                       std::pow(d1, 2) * (m_dist1.G2() + 3)) +
          m_p * (std::pow(m2 - m, 4) + 6 * std::pow(m2 - m, 2) * d2 +
                 4 * (m2 - m) * std::pow(d2, 1.5) * m_dist2.G1() +
                 std::pow(d2, 2) * (m_dist2.G2() + 3))) /
             std::pow(D(), 2) -
         3;
}

template <typename dist1_t, typename dist2_t>
ld MixtureDist<dist1_t, dist2_t>::Xi() const noexcept {
  return (nstu::randNum() < (1 - m_p)) ? m_dist1.Xi() : m_dist2.Xi();
}

template <typename dist1_t, typename dist2_t>
nstu::vector MixtureDist<dist1_t, dist2_t>::Xi(uint32_t size) const {
  nstu::vector result(size);
  for (ld& i : result) i = Xi();
  return result;
}