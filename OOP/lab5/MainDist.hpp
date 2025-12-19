#pragma once

#include "../lab2/header.h"
#include "IDist.hpp"
#include "IPersistent.hpp"

class MainDist : public IDist, private IPersistent {
 private:
  ld m_nu;      // параметр формы
  ld m_mu;      // сдвиг
  ld m_lambda;  // масштаб

 public:
  // Конструкторы

  // Стандартный конструктор без параметров
  MainDist();
  /**
   * Конструктор с параметрами распределения
   * @param nu Параметр формы распределения (не может быть равен 0).
   * @param mu Параметр сдвига распределения.
   * @param lambda Параметр масштаба распределения (не может быть равен 0).
   */
  explicit MainDist(ld nu, ld mu, ld lambda);
  /**
   * Конструктор с параметрами распределения из списка инициализации
   * @param list Список инициализации из трех элементов: nu, mu, lambda - форма,
   * сдвиг и масштаб соответственно.
   */
  explicit MainDist(std::initializer_list<ld> list);
  // Конструктор копирования
  MainDist(const MainDist& other);
  // Оператор присваивания копированием
  MainDist& operator=(const MainDist& other);
  // Конструктор перемещения
  MainDist(MainDist&& other) noexcept;
  // Оператор присваивания перемещением
  MainDist& operator=(MainDist&& other) noexcept;
  // Деструктор
  ~MainDist() = default;

  // Персистентность

  void save(const std::string& file_name) const override;
  void save(std::ofstream& file) const noexcept override;
  void load(const std::string& file_name) override;
  void load(std::ifstream& file) override;

  // Геттеры

  /**
   * Получение параметра формы nu
   * @return Значение параметра nu
   */
  ld getNu() const noexcept;
  /**
   * Получение параметра сдвига mu
   * @return Значение параметра mu
   */
  ld getMu() const noexcept;
  /**
   * Получение параметра масштаба lambda
   * @return Значение параметра lambda
   */
  ld getLambda() const noexcept;

  // Сеттеры

  /**
   * Изменение параметра формы nu (не может быть 0)
   * @param nu Новое значение параметра nu
   */
  void setNu(ld nu);
  /**
   * Изменение параметра сдвига mu
   * @param mu Новое значение параметра mu
   */
  void setMu(ld mu) noexcept;
  /**
   * Изменение параметра масштаба lambda (не может быть 0)
   * @param lambda Новое значение параметра lambda
   */
  void setLambda(ld lambda);

  // Функции интерфейса

  ld density(ld x) const noexcept override;
  nstu::vector density(const nstu::vector& x) const override;
  ld M() const noexcept override;
  ld D() const noexcept override;
  ld G1() const noexcept override;
  ld G2() const noexcept override;
  ld Xi() const noexcept override;
  nstu::vector Xi(uint32_t size) const override;
};