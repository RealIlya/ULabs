#pragma once

// подключаем заголовочный файл с nstu::vector из предыдущей лабы
#include "../lab2/header.h"
// подключаем заголовочный файл с MainDist из предыдущей лабы
#include "IDist.hpp"
#include "IPersistent.hpp"
#include "IObserver.hpp"

class EmpiricDist : public IDist, public IPersistent {
 private:
  nstu::vector m_init_sample;    // Начальная выборка
  uint32_t m_k;                  // Количество интервалов
  nstu::vector m_interval_freq;  // Частота на интервале
  ld m_min;                      // Минимальное значение выборки
  ld m_max;                      // Максимальное значение выборки
  ld m_delta;                    // Ширина интервала

  // Наблюдатель
  IObserver** m_observers;
  uint32_t m_observers_count;
  uint32_t m_observers_capacity;

  void notify() const;
  void resizeObservers();

 public:
  // Конструкторы

  // Стандартный конструктор без параметров
  EmpiricDist();
  /**
   * Конструктор с параметром начальной выборки
   * @param init_sample Начальная выборка, которая используется для вычисления
   * интервалов и вероятностей.
   * @param k Количество интервалов. По умолчанию рассчитывается по формуле
   * Стёрджеса.
   */
  explicit EmpiricDist(const nstu::vector& init_sample, uint32_t k = 0);
  /**
   * Конструктор с параметром основного распределения
   * @param n Размер выборки, которая будет сгенерирована из распределения.
   * @param dist Основное распределение, из которого генерируется выборка.
   * @param k Количество интервалов. По умолчанию рассчитывается по формуле
   * Стёрджеса.
   */
  explicit EmpiricDist(uint32_t n, IDist& dist, uint32_t k = 0);

  // Персистентность

  void save(const std::string& file_name) const override;
  void save(std::ofstream& file) const noexcept override;
  void load(std::ifstream& file) override;
  void load(const std::string& file_name) override;

  // Конструктор копирования
  EmpiricDist(const EmpiricDist& other);
  // Оператор присваивания копированием
  EmpiricDist& operator=(const EmpiricDist& other);
  // Конструктор перемещения
  EmpiricDist(EmpiricDist&& other) noexcept;
  // Оператор присваивания перемещением
  EmpiricDist& operator=(EmpiricDist&& other) noexcept;
  // Деструктор
  ~EmpiricDist() = default;

  // Методы Subject
  void attach(IObserver* obs);
  void detach(IObserver* obs);

  // Метод для "засорения" данных (вызывает notify)
  void setSample(const nstu::vector& new_sample);

  // Геттеры

  // Полученить копию начальной выборки
  const nstu::vector& getInitSample() const;
  // Получить количество интервалов
  ld getK() const noexcept;
  // Получить минимальное значение выборки
  ld getMin() const noexcept;
  // Получить максимальное значение выборки
  ld getMax() const noexcept;

  // Сеттеры

  // Установить количество интервалов
  void setK(uint32_t k);
  // Установить начальную выборку
  void setInitSample(const nstu::vector& sample);

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