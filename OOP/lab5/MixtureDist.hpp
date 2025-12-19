#pragma once

// подключаем заголовочный файл с nstu::vector из предыдущей лабы
#include "../lab2/header.h"
#include "IDist.hpp"
#include "IPersistent.hpp"

template <typename dist1_t, typename dist2_t>
class MixtureDist : public IDist, public IPersistent {
  static_assert(std::is_base_of<IDist, dist1_t>::value,
                "dist1_t должен наследоваться от IDist");
  static_assert(std::is_base_of<IDist, dist2_t>::value,
                "dist2_t должен наследоваться от IDist");

 private:
  ld m_p;           // Параметр смеси
  dist1_t m_dist1;  // Первый компонент
  dist2_t m_dist2;  // Второй компонент

 public:
  // Конструкторы

  // Конструктор по умолчанию
  MixtureDist();
  /**
   * Конструктор с парметрами
   * @param p Параметр смеси (в диапазоне [0, 1])
   * @param dist1 Первый компонент
   * @param dist2 Второй компонент
   */
  explicit MixtureDist(ld p, const dist1_t& dist1, const dist2_t& dist2);
  // Конструктор копирования
  MixtureDist(const MixtureDist& other);
  // Оператор присваивания копированием
  MixtureDist& operator=(const MixtureDist& other);
  // Конструктор перемещения
  MixtureDist(MixtureDist&& other) noexcept;
  // Оператор присваивания перемещением
  MixtureDist& operator=(MixtureDist&& other) noexcept;
  // Деструктор
  ~MixtureDist() = default;

  // Персистентность

  void save(const std::string& file_name) const override;
  void save(std::ofstream& file) const noexcept override;
  void load(const std::string& file_name) override;
  void load(std::ifstream& file) override;

  // Геттеры

  // Получить первый компонент
  dist1_t& component1() noexcept;
  const dist1_t& component1() const noexcept;
  // Получить второй компонент
  dist2_t& component2() noexcept;
  const dist2_t& component2() const noexcept;
  // Получить параметр смеси p
  ld getP() const noexcept;

  // Сеттеры

  // Задать параметр смеси p (в диапазоне [0, 1])
  void setP(ld p);

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
