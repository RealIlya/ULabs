#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>  // для std::accumulate
#include <stdexcept>

#include "../EmpiricDist.hpp"
#include "../MainDist.hpp"
#include "../MixtureDist.hpp"

// Вспомогательная функция для сравнения чисел с плавающей точкой
bool areFloatsEqual(ld a, ld b, ld epsilon = 1e-9) {
  return std::fabs(a - b) < epsilon;
}

// Вспомогательная функция для сравнения векторов
bool areVectorsEqual(const nstu::vector& v1, const nstu::vector& v2) {
  if (v1.size() != v2.size()) return false;
  for (size_t i = 0; i < v1.size(); ++i) {
    if (!areFloatsEqual(v1[i], v2[i])) return false;
  }
  return true;
}

// Вспомогательная функция для печати состояния объекта EmpiricDist
void printEmpiricState(const EmpiricDist& dist, const std::string& name) {
  std::cout << "   - Состояние " << name << ":\n";
  std::cout << "     - k=" << dist.getK() << ", min=" << dist.getMin()
            << ", max=" << dist.getMax()
            << ", sample_size=" << dist.getInitSample().size() << "\n";
}

/**
 * @brief Основная функция, тестирующая все методы класса EmpiricDist.
 */
void testAllEmpiricDistMethods() {
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "--- НАЧАЛО ПОЛНОГО ТЕСТИРОВАНИЯ КЛАССА EMPIRICDIST ---\n\n";

  const unsigned int seed = 12345;
  srand(seed);

  // ======================================================================
  // 1. ТЕСТИРОВАНИЕ КОНСТРУКТОРОВ
  // ======================================================================
  std::cout << "1. Тестирование конструкторов\n";
  std::cout << "--------------------------------\n";

  std::cout << "1.1. Конструктор из готовой выборки nstu::vector\n";
  nstu::vector sample = {10, 20, 21, 30, 32,
                         35, 40, 41, 48, 50};  // 10 элементов
  // При k=0, k должно стать floor(log2(10)) + 1 = 3 + 1 = 4
  EmpiricDist ed_from_vec(sample);
  printEmpiricState(ed_from_vec, "ed_from_vec");
  if (ed_from_vec.getK() == 4 && areFloatsEqual(ed_from_vec.getMin(), 10) &&
      areFloatsEqual(ed_from_vec.getMax(), 50)) {
    std::cout << "   >>> ПОДТВЕРЖДЕНО: Параметры рассчитаны корректно.\n";
  } else {
    std::cout << "   >>> ОШИБКА: Некорректный расчет параметров.\n";
  }

  std::cout << "\n1.2. Конструктор из другого распределения (MainDist)\n";
  MainDist main_dist(5.0, 0.0, 2.0);
  EmpiricDist ed_from_main(100, main_dist);
  printEmpiricState(ed_from_main, "ed_from_main");

  std::cout << "\n1.3. Конструктор из другого распределения (MixtureDist)\n";
  MixtureDist<MainDist, MainDist> mix_dist(0.4, MainDist(5, -5, 1),
                                           MainDist(5, 5, 1));
  EmpiricDist ed_from_mix(100, mix_dist);
  printEmpiricState(ed_from_mix, "ed_from_mix");
  std::cout << "\n";

  // ======================================================================
  // 2. ТЕСТИРОВАНИЕ СОХРАНЕНИЯ И ЗАГРУЗКИ (SAVE/LOAD)
  // ======================================================================
  std::cout << "2. Тестирование сохранения и загрузки (save/load)\n";
  std::cout << "------------------------------------------------\n";
  const std::string save_filename = "empiric_dist.txt";
  ed_from_vec.save(save_filename);
  EmpiricDist ed_loaded;
  ed_loaded.load(save_filename);
  printEmpiricState(ed_loaded, "ed_loaded");
  if (areVectorsEqual(ed_from_vec.getInitSample(), ed_loaded.getInitSample())) {
    std::cout << "   >>> ПОДТВЕРЖДЕНО: Сохраненная и загруженная выборки "
                 "совпадают!\n";
  } else {
    std::cout
        << "   >>> ОШИБКА: Выборки не совпадают после сохранения/загрузки!\n";
  }
  std::cout << "\n";

  // ======================================================================
  // 3. ТЕСТИРОВАНИЕ СЕТТЕРОВ
  // ======================================================================
  std::cout << "3. Тестирование сеттеров\n";
  std::cout << "------------------------\n";
  EmpiricDist ed_set(sample);
  std::cout << "   - Начальное k=" << ed_set.getK() << "\n";
  ed_set.setK(5);  // Пересчет на 5 интервалов
  std::cout << "   - Новое k после setK(5): " << ed_set.getK() << "\n";

  nstu::vector new_sample = {1, 2, 3, 4, 5};
  ed_set.setInitSample(new_sample);
  std::cout << "   - Установлена новая выборка размером " << new_sample.size()
            << "\n";
  printEmpiricState(ed_set, "ed_set (после setInitSample)");
  if (areFloatsEqual(ed_set.getMin(), 1) &&
      areFloatsEqual(ed_set.getMax(), 5)) {
    std::cout << "   >>> ПОДТВЕРЖДЕНО: Объект корректно пересчитан.\n";
  } else {
    std::cout << "   >>> ОШИБКА: Объект не был корректно пересчитан.\n";
  }
  std::cout << "\n";

  // ======================================================================
  // 4. ТЕСТИРОВАНИЕ АНАЛИТИЧЕСКИХ МЕТОДОВ
  // ======================================================================
  std::cout << "4. Тестирование аналитических методов\n";
  std::cout << "-------------------------------------\n";
  nstu::vector simple_sample = {1, 2, 3, 4, 5};
  EmpiricDist ed_analysis(simple_sample);
  printEmpiricState(ed_analysis, "ed_analysis");
  // M = (1+2+3+4+5)/5 = 3
  // D = ((1-3)^2 + (2-3)^2 + (3-3)^2 + (4-3)^2 + (5-3)^2)/5 = (4+1+0+1+4)/5 = 2
  std::cout << "   - Математическое ожидание (M): " << ed_analysis.M()
            << " (ожидалось: 3.0)\n";
  std::cout << "   - Дисперсия (D): " << ed_analysis.D()
            << " (ожидалось: 2.0)\n";
  // Для симметричного распределения G1 = 0
  std::cout << "   - Асимметрия (G1): " << ed_analysis.G1()
            << " (ожидалось: 0.0)\n";
  ld x = 2.5;  // Эта точка попадает во второй интервал
  std::cout << "   - Плотность в точке x=" << x << ": "
            << ed_analysis.density(x) << "\n";
  std::cout << "\n";

  // ======================================================================
  // 5. ТЕСТИРОВАНИЕ ГЕНЕРАЦИИ СЛУЧАЙНЫХ ЧИСЕЛ (Xi)
  // ======================================================================
  std::cout << "5. Тестирование генерации случайных чисел (Xi)\n";
  std::cout << "----------------------------------------------\n";
  srand(seed);
  nstu::vector xi_sample1 = ed_from_vec.Xi(5);
  std::cout << "   - Запуск 1: ";
  for (const auto& val : xi_sample1) std::cout << val << " ";
  std::cout << "\n";

  srand(seed);  // Сброс
  nstu::vector xi_sample2 = ed_from_vec.Xi(5);
  std::cout << "   - Запуск 2: ";
  for (const auto& val : xi_sample2) std::cout << val << " ";
  std::cout << "\n";

  if (areVectorsEqual(xi_sample1, xi_sample2)) {
    std::cout << "   >>> ПОДТВЕРЖДЕНО: Воспроизводимость генерации для "
                 "EmpiricDist работает!\n";
  } else {
    std::cout << "   >>> ОШИБКА: Воспроизводимость генерации НЕ работает!\n";
  }
  std::cout << "\n";

  // ======================================================================
  // 6. ТЕСТИРОВАНИЕ ОБРАБОТКИ ИСКЛЮЧЕНИЙ
  // ======================================================================
  std::cout << "6. Тестирование обработки исключений\n";
  std::cout << "------------------------------------\n";
  try {
    nstu::vector empty_sample;
    EmpiricDist tmp(empty_sample);
  } catch (const std::invalid_argument& e) {
    std::cout << "   - Успех (пустая выборка в конструкторе): " << e.what()
              << std::endl;
  }
  try {
    EmpiricDist tmp;
    tmp.load("non_existent_file.txt");
  } catch (const std::runtime_error& e) {
    std::cout << "   - Успех (несуществующий файл): " << e.what() << std::endl;
  }

  std::cout << "\n--- ПОЛНОЕ ТЕСТИРОВАНИЕ ЗАВЕРШЕНО ---\n";
}

int main() {
  try {
    testAllEmpiricDistMethods();
  } catch (const std::exception& e) {
    std::cerr << "Критическая ошибка во время тестирования: " << e.what()
              << std::endl;
    return 1;
  }
  return 0;
}