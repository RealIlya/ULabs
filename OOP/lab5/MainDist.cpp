#include "MainDist.hpp"

#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <fstream>

MainDist::MainDist() : m_nu(1.0), m_mu(0.0), m_lambda(1.0) {}

MainDist::MainDist(ld nu, ld mu, ld lambda)
    : m_nu(nu), m_mu(mu), m_lambda(lambda) {
  if (lambda == 0) throw std::invalid_argument("Лямбда не может быть 0");
  if (nu == 0) throw std::invalid_argument("Нулевое значение nu");
}

void MainDist::load(const std::string& file_name) {
  std::ifstream file(file_name);
  if (!file.is_open()) throw std::runtime_error("Не удалось открыть файл");
  load(file);
  file.close();
}

void MainDist::load(std::ifstream& file) {
  std::string d_name;
  file >> d_name;
  if (d_name != "MainDist")
    throw std::runtime_error("Неверный тип распределения");
  file >> m_nu >> m_mu >> m_lambda;
  if (m_lambda == 0) throw std::invalid_argument("Лямбда не может быть 0");
  if (m_nu == 0) throw std::invalid_argument("Нулевое значение nu");
}

void MainDist::save(const std::string& file_name) const {
  std::ofstream file(file_name);
  if (!file.is_open()) throw std::runtime_error("Не удалось открыть файл");
  save(file);
  file.close();
}

void MainDist::save(std::ofstream& file) const noexcept {
  file << "MainDist ";
  file << m_nu << " " << m_mu << " " << m_lambda << " ";
}

MainDist::MainDist(std::initializer_list<ld> list) {
  if (list.size() != 3)
    throw std::invalid_argument("Неверное количество аргументов");
  auto i = list.begin();
  m_nu = *(i++);
  m_mu = *(i++);
  m_lambda = *i;
  if (m_lambda == 0) throw std::invalid_argument("Лямбда не может быть 0");
  if (m_nu == 0) throw std::invalid_argument("Нулевое значение nu");
}

MainDist::MainDist(const MainDist& other)
    : m_nu(other.m_nu), m_mu(other.m_mu), m_lambda(other.m_lambda) {}

MainDist& MainDist::operator=(const MainDist& other) {
  if (this != &other) {
    m_nu = other.m_nu;
    m_mu = other.m_mu;
    m_lambda = other.m_lambda;
  }
  return *this;
}

MainDist::MainDist(MainDist&& other) noexcept
    : m_nu(other.m_nu), m_mu(other.m_mu), m_lambda(other.m_lambda) {
  other.m_nu = 0;
  other.m_mu = 0;
  other.m_lambda = 0;
}

MainDist& MainDist::operator=(MainDist&& other) noexcept {
  if (this != &other) {
    m_nu = other.m_nu;
    m_mu = other.m_mu;
    m_lambda = other.m_lambda;
    other.m_nu = 0;
    other.m_mu = 0;
    other.m_lambda = 0;
  }
  return *this;
}

ld MainDist::getNu() const noexcept { return m_nu; }
ld MainDist::getMu() const noexcept { return m_mu; }
ld MainDist::getLambda() const noexcept { return m_lambda; }

void MainDist::setNu(ld nu) {
  if (nu == 0)
    throw std::invalid_argument("Нулевое значение nu");
  else
    m_nu = nu;
}
void MainDist::setMu(ld mu) noexcept { m_mu = mu; }
void MainDist::setLambda(ld lambda) {
  if (lambda == 0)
    throw std::invalid_argument("Лямбда не может быть 0");
  else
    m_lambda = lambda;
}

ld MainDist::density(ld x) const noexcept {
  ld coeff = 2 * m_lambda * sqrt(m_nu) * boost::math::cyl_bessel_k(1, m_nu);
  ld exponent = exp(-m_nu * sqrt(1 + pow((x - m_mu) / m_lambda, 2) / m_nu));
  return exponent / coeff;
}

nstu::vector MainDist::density(const nstu::vector& x) const {
  nstu::vector result(x.size());
  for (uint32_t i = 0; i < x.size(); ++i) {
    result[i] = density(x[i]);
  }
  return result;
}

ld MainDist::M() const noexcept { return m_mu; }

ld MainDist::D() const noexcept {
  return pow(m_lambda, 2) * boost::math::cyl_bessel_k(2, m_nu) /
         boost::math::cyl_bessel_k(1, m_nu);
}

ld MainDist::G1() const noexcept { return 0; }

ld MainDist::G2() const noexcept {
  return 3 * boost::math::cyl_bessel_k(3, m_nu) *
             boost::math::cyl_bessel_k(1, m_nu) /
             pow(boost::math::cyl_bessel_k(2, m_nu), 2) -
         3;
}

ld MainDist::Xi() const noexcept {
  ld r1, r2, delta, t;
  do {
    r1 = nstu::randNum();
    r2 = nstu::randNum();
    delta = 2 * (sqrt(1 + pow(m_nu, 2)) - 1) / m_nu;
    t = -2 * log(r1) / delta;
  } while (-log(r2) <= (m_nu - delta) * t / 2 + m_nu / (2 * t) -
                           sqrt(m_nu * (m_nu - delta)));
  ld r3 = nstu::randNum(), r4 = nstu::randNum();
  ld z = sqrt(-2 * log(r3)) * cos(2 * nstu::pi * r4);
  // z = sqrt(-2 * log(r3)) * sin(2 * pi * r4);
  return z * sqrt(t) * m_lambda + m_mu;
}

nstu::vector MainDist::Xi(uint32_t size) const {
  nstu::vector result(size);
  for (auto& i : result) i = Xi();
  return result;
}