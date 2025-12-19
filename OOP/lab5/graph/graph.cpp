#include <iostream>

#include "../../lab2/header.h"
#include "../EmpiricDist.hpp"
#include "../MainDist.hpp"
#include "../MixtureDist.hpp"

using namespace std;

void testDistribution_3_1_1() {
  ld mu = 0, lambda = 1, nu = 1;
  ld x = 0;
  MainDist main_dist(nu, mu, lambda);
  cout << "-----Тест 3.1.1-----" << endl;
  cout << "Тест основного распределения:" << endl;
  cout << "nu = " << nu << ", mu = " << mu << ", lambda = " << lambda
       << ", x = " << x << endl;
  cout << "Плотность вероятности f(x): " << main_dist.density(x) << endl;
  cout << "M = " << main_dist.M() << "; D = " << main_dist.D()
       << "; G1 = " << main_dist.G1() << "; G2 = " << main_dist.G2() << endl;
  auto fout = ofstream("test311.csv");
  fout << "x,f(x),xi,f(xi)" << endl;
  for (ld i = -15; i <= 15; i += 0.5) {
    ld xi = main_dist.Xi();
    ld density = main_dist.density(i);
    ld density_xi = main_dist.density(xi);
    fout << i << "," << density << "," << xi << "," << density_xi << endl;
  }
}

void testDistribution_3_1_2() {
  ld mu = 0, lambda = 2, nu = 1;
  ld x = 0;
  MainDist main_dist(nu, mu, lambda);
  cout << "\n-----Тест 3.1.2-----" << endl;
  cout << "Тест основного распределения:" << endl;
  cout << "nu = " << nu << ", mu = " << mu << ", lambda = " << lambda
       << ", x = " << x << endl;
  cout << "Плотность вероятности f(x): " << main_dist.density(x) << endl;
  cout << "M = " << main_dist.M() << "; D = " << main_dist.D()
       << "; G1 = " << main_dist.G1() << "; G2 = " << main_dist.G2() << endl;
  auto fout = ofstream("test312.csv");
  fout << "x,f(x),xi,f(xi)" << endl;
  for (ld i = -15; i <= 15; i += 0.5) {
    ld xi = main_dist.Xi();
    ld density = main_dist.density(i);
    ld density_xi = main_dist.density(xi);
    fout << i << "," << density << "," << xi << "," << density_xi << endl;
  }
}

void testDistribution_3_1_3() {
  ld mu = 10, lambda = 2, nu = 1;
  ld x = 0;
  MainDist main_dist(nu, mu, lambda);
  cout << "\n-----Тест 3.1.3-----" << endl;
  cout << "Тест основного распределения:" << endl;
  cout << "nu = " << nu << ", mu = " << mu << ", lambda = " << lambda
       << ", x = " << x << endl;
  cout << "Плотность вероятности f(x): " << main_dist.density(x) << endl;
  cout << "M = " << main_dist.M() << "; D = " << main_dist.D()
       << "; G1 = " << main_dist.G1() << "; G2 = " << main_dist.G2() << endl;
  auto fout = ofstream("test313.csv");
  fout << "x,f(x),xi,f(xi)" << endl;
  for (ld i = -15; i <= 15; i += 0.5) {
    ld xi = main_dist.Xi();
    ld density = main_dist.density(i);
    ld density_xi = main_dist.density(xi);
    fout << i << "," << density << "," << xi << "," << density_xi << endl;
  }
}

void testDistribution_3_2_1() {
  ld mu1 = 10, mu2 = mu1, lambda1 = 2, lambda2 = lambda1, nu1 = 1, nu2 = nu1;
  ld p = 0.5;
  ld x = 0;
  MixtureDist<MainDist, MainDist> mixt_dist(p, MainDist{nu1, mu1, lambda1},
                                            MainDist{nu2, mu2, lambda2});
  cout << "\n-----Тест 3.2.1-----" << endl;
  cout << "Тест распределения смеси:" << endl;
  cout << "nu1 = " << nu1 << "mu1 = " << mu1 << ", lambda1 = " << lambda1
       << '\n'
       << "nu2 = " << nu2 << "mu2 = " << mu2 << ", lambda2 = " << lambda2
       << '\n'
       << "p = " << p << ", x = " << x << endl;
  cout << "Плотность вероятности f(x): " << mixt_dist.density(x) << endl;
  cout << "M = " << mixt_dist.M() << "; D = " << mixt_dist.D()
       << "; G1 = " << mixt_dist.G1() << "; G2 = " << mixt_dist.G2() << endl;
  auto fout = ofstream("test321.csv");
  fout << "x,f(x),xi,f(xi)" << endl;
  for (ld i = -5; i <= 25; i += 0.5) {
    auto xi = mixt_dist.Xi();
    fout << i << "," << mixt_dist.density(i) << "," << xi << ","
         << mixt_dist.density(xi) << endl;
  }
}
void testDistribution_3_2_2() {
  ld mu1 = 0, mu2 = 2, lambda1 = 1, lambda2 = lambda1, nu1 = 1, nu2 = nu1;
  ld p = 0.75;
  ld x = 0;
  MixtureDist<MainDist, MainDist> mixt_dist(p, MainDist{nu1, mu1, lambda1},
                                            MainDist{nu2, mu2, lambda2});
  cout << "\n-----Тест 3.2.2-----" << endl;
  cout << "Тест распределения смеси:" << endl;
  cout << "nu1 = " << nu1 << "mu1 = " << mu1 << ", lambda1 = " << lambda1
       << '\n'
       << "nu2 = " << nu2 << "mu2 = " << mu2 << ", lambda2 = " << lambda2
       << '\n'
       << "p = " << p << ", x = " << x << endl;
  cout << "Плотность вероятности f(x): " << mixt_dist.density(x) << endl;
  cout << "M = " << mixt_dist.M() << "; D = " << mixt_dist.D()
       << "; G1 = " << mixt_dist.G1() << "; G2 = " << mixt_dist.G2() << endl;
  auto fout = ofstream("test322.csv");
  fout << "x,f(x),xi,f(xi)" << endl;
  for (ld i = -15; i <= 15; i += 0.5) {
    auto xi = mixt_dist.Xi();
    fout << i << "," << mixt_dist.density(i) << "," << xi << ","
         << mixt_dist.density(xi) << endl;
  }
}

void testDistribution_3_2_3() {
  ld mu1 = 0, mu2 = mu1, lambda1 = 1, lambda2 = 3, nu1 = 1, nu2 = nu1;
  ld p = 0.5;
  ld x = 0;
  MixtureDist mixt_dist(p, MainDist{nu1, mu1, lambda1},
                        MainDist{nu2, mu2, lambda2});
  cout << "\n-----Тест 3.2.3-----" << endl;
  cout << "Тест распределения смеси:" << endl;
  cout << "nu1 = " << nu1 << "mu1 = " << mu1 << ", lambda1 = " << lambda1
       << '\n'
       << "nu2 = " << nu2 << "mu2 = " << mu2 << ", lambda2 = " << lambda2
       << '\n'
       << "p = " << p << ", x = " << x << endl;
  cout << "Плотность вероятности f(x): " << mixt_dist.density(x) << endl;
  cout << "M = " << mixt_dist.M() << "; D = " << mixt_dist.D()
       << "; G1 = " << mixt_dist.G1() << "; G2 = " << mixt_dist.G2() << endl;
  auto fout = ofstream("test323.csv");
  fout << "x,f(x),xi,f(xi)" << endl;
  for (ld i = -15; i <= 15; i += 0.5) {
    auto xi = mixt_dist.Xi();
    fout << i << "," << mixt_dist.density(i) << "," << xi << ","
         << mixt_dist.density(xi) << endl;
  }
}

void testDistribution_3_2_4() {
  ld mu1 = 0, mu2 = mu1, lambda1 = 1, lambda2 = lambda1, nu1 = 0.1, nu2 = 30;
  ld p = 0.5;
  ld x = 0;
  MixtureDist mixt_dist(p, MainDist{nu1, mu1, lambda1},
                        MainDist{nu2, mu2, lambda2});
  cout << "\n-----Тест 3.2.4-----" << endl;
  cout << "Тест распределения смеси:" << endl;
  cout << "nu1 = " << nu1 << "mu1 = " << mu1 << ", lambda1 = " << lambda1
       << '\n'
       << "nu2 = " << nu2 << "mu2 = " << mu2 << ", lambda2 = " << lambda2
       << '\n'
       << "p = " << p << ", x = " << x << endl;
  cout << "Плотность вероятности f(x): " << mixt_dist.density(x) << endl;
  cout << "M = " << mixt_dist.M() << "; D = " << mixt_dist.D()
       << "; G1 = " << mixt_dist.G1() << "; G2 = " << mixt_dist.G2() << endl;
  auto fout = ofstream("test324.csv");
  fout << "x,f(x),xi,f(xi)" << endl;
  for (ld i = -15; i <= 15; i += 0.5) {
    auto xi = mixt_dist.Xi();
    fout << i << "," << mixt_dist.density(i) << "," << xi << ","
         << mixt_dist.density(xi) << endl;
  }
}

void testDistribution_3_3_1_1() {
  ld mu = 0, lambda = 1, nu = 1;
  ld x = 0;
  uint32_t n = 10000;
  MainDist dist(nu, mu, lambda);
  EmpiricDist e_dist(n, dist);
  cout << "\n-----Тест 3.3.1.1-----" << endl;
  cout << "Тест эмпирического распределения:" << endl;
  cout << "Распределение построено на основе выборки " << n
       << " элементов из основного распределения со следующими параметрами:"
       << endl;
  cout << "nu = " << nu << ", mu = " << mu << ", lambda = " << lambda << endl;
  cout << "Выборка:\nf(0) = " << dist.density(x);
  cout << ", M = " << dist.M() << ", D = " << dist.D() << ", G1 = " << dist.G1()
       << ", G2 = " << dist.G2() << endl;
  cout << "Эмпирическое распределение:\nf(0) = " << e_dist.density(x);
  cout << ", M = " << e_dist.M() << "; D = " << e_dist.D()
       << "; G1 = " << e_dist.G1() << "; G2 = " << e_dist.G2() << endl;
  auto fout = ofstream("test3311.csv");
  fout << "x,f(x),fe(x),xi,fe(xi)" << endl;
  for (ld i = -15; i <= 15; i += 0.5) {
    ld xi = e_dist.Xi();
    fout << i << "," << dist.density(i) << "," << e_dist.density(i) << "," << xi
         << "," << e_dist.density(xi) << endl;
  }
}

void testDistribution_3_3_1_2() {
  ld mu = 10, lambda = 2, nu = 3;
  ld x = 0;
  uint32_t n = 10000;
  MainDist dist(nu, mu, lambda);
  EmpiricDist e_dist(n, dist);
  cout << "\n-----Тест 3.3.1.2-----" << endl;
  cout << "Тест эмпирического распределения:" << endl;
  cout << "Распределение построено на основе выборки " << n
       << " элементов из основного распределения со следующими параметрами:"
       << endl;
  cout << "nu = " << nu << ", mu = " << mu << ", lambda = " << lambda << endl;
  cout << "Выборка:\nf(0) = " << dist.density(x);
  cout << ", M = " << dist.M() << ", D = " << dist.D() << ", G1 = " << dist.G1()
       << ", G2 = " << dist.G2() << endl;
  cout << "Эмпирическое распределение:\nf(0) = " << e_dist.density(x);
  cout << ", M = " << e_dist.M() << "; D = " << e_dist.D()
       << "; G1 = " << e_dist.G1() << "; G2 = " << e_dist.G2() << endl;
  auto fout = ofstream("test3312.csv");
  fout << "x,f(x),fe(x),xi,fe(xi)" << endl;
  for (ld i = -5; i <= 25; i += 0.5) {
    ld xi = e_dist.Xi();
    fout << i << "," << dist.density(i) << "," << e_dist.density(i) << "," << xi
         << "," << e_dist.density(xi) << endl;
  }
}

void testDistribution_3_3_1_3() {
  ld mu1 = 0, mu2 = mu1, lambda1 = 1, lambda2 = 2, nu1 = 1, nu2 = nu1;
  ld x = 0, p = 0.5;
  uint32_t n = 10000;
  MixtureDist<MainDist, MainDist> mixture(p, MainDist{nu1, mu1, lambda1},
                                          MainDist{nu2, mu2, lambda2});
  EmpiricDist e_dist(n, mixture);
  cout << "\n-----Тест 3.3.1.3-----" << endl;
  cout << "Тест эмпирического распределения:" << endl;
  cout << "Распределение построено на основе выборки " << n
       << " элементов из распределения смеси со следующими параметрами:"
       << endl;
  cout << "nu1 = " << nu1 << "mu1 = " << mu1 << ", lambda1 = " << lambda1
       << '\n'
       << "nu2 = " << nu2 << "mu2 = " << mu2 << ", lambda2 = " << lambda2
       << '\n'
       << "p = " << p << endl;
  cout << "Выборка:\nf(0) = " << mixture.density(x);
  cout << ", M = " << mixture.M() << ", D = " << mixture.D()
       << ", G1 = " << mixture.G1() << ", G2 = " << mixture.G2() << endl;
  cout << "Эмпирическое распределение:\nf(0) = " << e_dist.density(x);
  cout << ", M = " << e_dist.M() << "; D = " << e_dist.D()
       << "; G1 = " << e_dist.G1() << "; G2 = " << e_dist.G2() << endl;
  auto fout = ofstream("test3313.csv");
  fout << "x,f(x),fe(x),xi,fe(xi)" << endl;
  for (ld i = -15; i <= 15; i += 0.5) {
    ld xi = e_dist.Xi();
    fout << i << "," << mixture.density(i) << "," << e_dist.density(i) << ","
         << xi << "," << e_dist.density(xi) << endl;
  }
}
void testDistribution_3_3_1_4() {
  ld mu1 = 0, mu2 = 3, lambda1 = 1, lambda2 = 2, nu1 = 1, nu2 = 3;
  ld x = 0, p = 0.4;
  uint32_t n = 10000;
  MixtureDist<MainDist, MainDist> mixture(p, MainDist{nu1, mu1, lambda1},
                                          MainDist{nu2, mu2, lambda2});
  EmpiricDist e_dist(n, mixture);
  cout << "\n-----Тест 3.3.1.4-----" << endl;
  cout << "Тест эмпирического распределения:" << endl;
  cout << "Распределение построено на основе выборки " << n
       << " элементов из распределения смеси со следующими параметрами:"
       << endl;
  cout << "nu1 = " << nu1 << "mu1 = " << mu1 << ", lambda1 = " << lambda1
       << '\n'
       << "nu2 = " << nu2 << "mu2 = " << mu2 << ", lambda2 = " << lambda2
       << '\n'
       << "p = " << p << endl;
  cout << "Выборка:\nf(0) = " << mixture.density(x);
  cout << ", M = " << mixture.M() << ", D = " << mixture.D()
       << ", G1 = " << mixture.G1() << ", G2 = " << mixture.G2() << endl;
  cout << "Эмпирическое распределение:\nf(0) = " << e_dist.density(x);
  cout << ", M = " << e_dist.M() << "; D = " << e_dist.D()
       << "; G1 = " << e_dist.G1() << "; G2 = " << e_dist.G2() << endl;
  auto fout = ofstream("test3314.csv");
  fout << "x,f(x),fe(x),xi,fe(xi)" << endl;
  for (ld i = -15; i <= 15; i += 0.5) {
    ld xi = e_dist.Xi();
    fout << i << "," << mixture.density(i) << "," << e_dist.density(i) << ","
         << xi << "," << e_dist.density(xi) << endl;
  }
}

void testDistribution_3_3_2() {
  ld mu = 0, lambda = 1, nu = 1;
  ld x = 0;
  uint32_t n = 10000;
  MainDist dist(nu, mu, lambda);
  EmpiricDist e_dist(n, dist);
  cout << "\n-----Тест 3.3.2-----" << endl;
  cout << "Тест эмпирического распределения:" << endl;
  cout << "Распределение построено на основе выборки " << n
       << " элементов из основного распределения со следующими параметрами:"
       << endl;
  cout << "nu = " << nu << ", mu = " << mu << ", lambda = " << lambda << endl;
  cout << "Выборка:\nf(0) = " << dist.density(x);
  cout << ", M = " << dist.M() << ", D = " << dist.D() << ", G1 = " << dist.G1()
       << ", G2 = " << dist.G2() << endl;
  cout << "Эмпирическое распределение:\nf(0) = " << e_dist.density(x);
  cout << ", M = " << e_dist.M() << "; D = " << e_dist.D()
       << "; G1 = " << e_dist.G1() << "; G2 = " << e_dist.G2() << endl;
  EmpiricDist e2_dist(n, e_dist);
  cout << "Эмпирическое распределение2:\nf(0) = " << e2_dist.density(x);
  cout << ", M = " << e2_dist.M() << ", D = " << e2_dist.D()
       << ", G1 = " << e2_dist.G1() << ", G2 = " << e2_dist.G2() << endl;
  auto fout = ofstream("test332.csv");
  fout << "x,f(x),fe(x),xi,fe(xi),fee(x),xie,fee(xie)" << endl;
  for (ld i = -15; i <= 15; i += 0.5) {
    ld xi = e_dist.Xi();
    ld xie = e2_dist.Xi();
    fout << i << "," << dist.density(i) << "," << e_dist.density(i) << "," << xi
         << "," << e_dist.density(xi) << "," << e2_dist.density(i) << "," << xie
         << "," << e2_dist.density(xie) << endl;
  }
}

int main() {
  testDistribution_3_1_1();
  testDistribution_3_1_2();
  testDistribution_3_1_3();
  testDistribution_3_2_1();
  testDistribution_3_2_2();
  testDistribution_3_2_3();
  testDistribution_3_2_4();
  testDistribution_3_3_1_1();
  testDistribution_3_3_1_2();
  testDistribution_3_3_1_3();
  testDistribution_3_3_1_4();
  testDistribution_3_3_2();
  return 0;
}