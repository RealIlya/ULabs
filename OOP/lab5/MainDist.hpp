#pragma once

#include <cstdint>
#include <initializer_list>
#include <string>

#include "../lab2/header.h"

class MainDist {
 private:
  ld m_nu;
  ld m_mu;
  ld m_lambda;

 public:
  MainDist();
  explicit MainDist(ld nu, ld mu, ld lambda);
  explicit MainDist(std::initializer_list<ld> list);

  ld getNu() const noexcept;
  ld getMu() const noexcept;
  ld getLambda() const noexcept;

  void setNu(ld nu);
  void setMu(ld mu) noexcept;
  void setLambda(ld lambda);

  ld density(ld x) const noexcept;
  ld M() const noexcept;
  ld D() const noexcept;
  ld G1() const noexcept;
  ld G2() const noexcept;

  ld Xi() const noexcept;
  nstu::vector Xi(uint32_t size) const;

  void save(const std::string& file_name) const;
  void load(const std::string& file_name);
};

