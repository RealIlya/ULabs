#include "MainDist.hpp"

#include <cmath>
#include <fstream>
#include <stdexcept>

MainDist::MainDist() : m_nu(1.0L), m_mu(0.0L), m_lambda(1.0L) {}

MainDist::MainDist(ld nu, ld mu, ld lambda)
    : m_nu(nu), m_mu(mu), m_lambda(lambda) {
  if (m_lambda == 0) throw std::invalid_argument("Лямбда не может быть 0");
  if (m_nu == 0) throw std::invalid_argument("Нулевое значение nu");
}

MainDist::MainDist(std::initializer_list<ld> list) {
  if (list.size() != 3)
    throw std::invalid_argument("Неверное количество аргументов");
  auto it = list.begin();
  m_nu = *(it++);
  m_mu = *(it++);
  m_lambda = *it;
  if (m_lambda == 0) throw std::invalid_argument("Лямбда не может быть 0");
  if (m_nu == 0) throw std::invalid_argument("Нулевое значение nu");
}

ld MainDist::getNu() const noexcept { return m_nu; }
ld MainDist::getMu() const noexcept { return m_mu; }
ld MainDist::getLambda() const noexcept { return m_lambda; }

void MainDist::setNu(ld nu) {
  if (nu == 0) throw std::invalid_argument("Нулевое значение nu");
  m_nu = nu;
}

void MainDist::setMu(ld mu) noexcept { m_mu = mu; }

void MainDist::setLambda(ld lambda) {
  if (lambda == 0) throw std::invalid_argument("Лямбда не может быть 0");
  m_lambda = lambda;
}

ld MainDist::density(ld x) const noexcept {
  const ld v = 1.0L;
  const ld k1 = std::cyl_bessel_k(static_cast<double>(v),
                                 static_cast<double>(m_nu));
  const ld coeff = 2 * m_lambda * std::sqrt(m_nu) * k1;
  const ld inner = 1 + std::pow((x - m_mu) / m_lambda, 2) / m_nu;
  const ld exponent = std::exp(-m_nu * std::sqrt(inner));
  return exponent / coeff;
}

ld MainDist::M() const noexcept { return m_mu; }

ld MainDist::D() const noexcept {
  const ld k1 = std::cyl_bessel_k(1.0, static_cast<double>(m_nu));
  const ld k2 = std::cyl_bessel_k(2.0, static_cast<double>(m_nu));
  return std::pow(m_lambda, 2) * k2 / k1;
}

ld MainDist::G1() const noexcept { return 0; }

ld MainDist::G2() const noexcept {
  const ld k1 = std::cyl_bessel_k(1.0, static_cast<double>(m_nu));
  const ld k2 = std::cyl_bessel_k(2.0, static_cast<double>(m_nu));
  const ld k3 = std::cyl_bessel_k(3.0, static_cast<double>(m_nu));
  return 3 * k3 * k1 / (k2 * k2) - 3;
}

ld MainDist::Xi() const noexcept {
  ld r1, r2, delta, t;
  do {
    r1 = nstu::randNum();
    r2 = nstu::randNum();
    delta = 2 * (std::sqrt(1 + std::pow(m_nu, 2)) - 1) / m_nu;
    t = -2 * std::log(r1) / delta;
  } while (-std::log(r2) <= (m_nu - delta) * t / 2 + m_nu / (2 * t) -
                               std::sqrt(m_nu * (m_nu - delta)));

  const ld r3 = nstu::randNum();
  const ld r4 = nstu::randNum();
  const ld z = std::sqrt(-2 * std::log(r3)) * std::cos(2 * nstu::pi * r4);
  return z * std::sqrt(t) * m_lambda + m_mu;
}

nstu::vector MainDist::Xi(uint32_t size) const {
  nstu::vector result(size);
  for (auto& x : result) x = Xi();
  return result;
}

void MainDist::save(const std::string& file_name) const {
  std::ofstream file(file_name);
  if (!file.is_open()) throw std::runtime_error("Не удалось открыть файл");
  file << "MainDist " << m_nu << " " << m_mu << " " << m_lambda << " ";
}

void MainDist::load(const std::string& file_name) {
  std::ifstream file(file_name);
  if (!file.is_open()) throw std::runtime_error("Не удалось открыть файл");
  std::string tag;
  file >> tag;
  if (tag != "MainDist") throw std::runtime_error("Неверный тип распределения");
  file >> m_nu >> m_mu >> m_lambda;
  if (m_lambda == 0) throw std::invalid_argument("Лямбда не может быть 0");
  if (m_nu == 0) throw std::invalid_argument("Нулевое значение nu");
}
