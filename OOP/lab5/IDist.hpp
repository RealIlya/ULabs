#pragma once

#include "../lab2/header.h"

class IDist {
 public:
  /**
   * Создание копии объекта и возвращение указателя на него
   */
  // virtual std::unique_ptr<IDist> copy() const noexcept = 0;

  /**
   * Вычисление плотности вероятности в точке x
   * @param x Точка, в которой вычисляется плотность вероятности
   * @return Значение плотности вероятности в точке x
   */
  virtual ld density(ld x) const noexcept = 0;
  /**
   * Вычисление плотности вероятности в каждой точке векторa x
   * @param x Вектор точек, в которых вычисляется плотность вероятности
   * @return Вектор значений плотности вероятности в каждой точке из вектора x
   */
  virtual nstu::vector density(const nstu::vector& x) const = 0;
  /**
   * Вычисление математического ожидания
   * @return Значение математического ожидания
   */
  virtual ld M() const noexcept = 0;
  /**
   * Вычисление дисперсии
   * @return Значение дисперсии
   */
  virtual ld D() const noexcept = 0;
  /**
   * Вычисление коэффициента асимметрии
   * @return Значение коэффициента асимметрии
   */
  virtual ld G1() const noexcept = 0;
  /**
   * Вычисление коэффициента эксцесса
   * @return Значение коэффициента эксцесса
   */
  virtual ld G2() const noexcept = 0;
  /**
   * Моделирование случайной величины
   * @return Случайная величина
   */
  virtual ld Xi() const noexcept = 0;
  /**
   * Моделирование выборки случайных величин
   * @param size Размер выборки случайных величин
   * @return Выборка случайных величин
   */
  virtual nstu::vector Xi(uint32_t size) const = 0;
};