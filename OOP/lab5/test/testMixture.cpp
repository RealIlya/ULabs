#include <cmath>    // для fabs
#include <fstream>  // для работы с файлами
#include <iomanip>  // для std::fixed, std::setprecision
#include <iostream>
#include <stdexcept>  // для std::invalid_argument

#include "../MainDist.hpp"
#include "../MixtureDist.hpp"

// Вспомогательная функция для сравнения двух векторов
bool areVectorsEqual(const nstu::vector& v1, const nstu::vector& v2) {
  if (v1.size() != v2.size()) {
    return false;
  }
  for (size_t i = 0; i < v1.size(); ++i) {
    if (std::fabs(v1[i] - v2[i]) > 1e-9) {
      return false;
    }
  }
  return true;
}

bool compareMixtures(const MixtureDist<MainDist, MainDist>& mix1,
                     const MixtureDist<MainDist, MainDist>& mix2) {
  return areVectorsEqual(
      {mix1.getP(), mix1.component1().getNu(), mix1.component1().getMu(),
       mix1.component1().getLambda()},
      {mix2.component2().getNu(), mix2.component2().getMu(),
       mix2.component2().getLambda()});
}

void printDistState(const MainDist& dist, const std::string& name) {
  std::cout << dist.getNu() << ", " << dist.getMu() << ", " << dist.getLambda()
            << "\n";
}

// Вспомогательная функция для печати состояния объекта MixtureDist
void printMixState(const MixtureDist<MainDist, MainDist>& mix,
                   const std::string& name) {
  std::cout << "   - Состояние " << name << ": p=" << mix.getP() << "\n";
  std::cout << "     - Comp1 (nu, mu, lambda): ";
  printDistState(mix.component1(), "comp1");
  std::cout << "     - Comp2 (nu, mu, lambda): " << mix.component2().getNu()
            << ", " << mix.component2().getMu() << ", "
            << mix.component2().getLambda() << "\n";
}

/**
 * @brief Основная функция, тестирующая все методы класса MixtureDist.
 */
void testAllMixtureDistMethods() {
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "--- НАЧАЛО ПОЛНОГО ТЕСТИРОВАНИЯ КЛАССА MIXTUREDIST ---\n\n";

  // Подготовим два базовых распределения для использования в тестах
  MainDist d1(5.0, -2.0, 1.0);
  MainDist d2(3.0, 2.0, 1.5);

  // ======================================================================
  // 1. ТЕСТИРОВАНИЕ КОНСТРУКТОРОВ
  // ======================================================================
  std::cout << "1. Тестирование конструкторов\n";
  std::cout << "--------------------------------\n";

  std::cout << "1.1. Конструктор по умолчанию\n";
  MixtureDist<MainDist, MainDist> mix_default;
  printMixState(mix_default, "mix_default");

  std::cout << "\n1.2. Конструктор с параметрами\n";
  MixtureDist mix_params(0.3, d1, d2);
  printMixState(mix_params, "mix_params");

  std::cout << "\n1.3. Конструктор копирования\n";
  MixtureDist mix_copy = mix_params;
  printMixState(mix_copy, "mix_copy");

  std::cout << "\n1.4. Конструктор перемещения\n";
  MixtureDist mix_move = MixtureDist(0.8, d1, d2);
  printMixState(mix_move, "mix_move");
  std::cout << "\n";

  // ======================================================================
  // 2. ТЕСТИРОВАНИЕ СОХРАНЕНИЯ И ЗАГРУЗКИ (SAVE/LOAD)
  // ======================================================================
  std::cout << "2. Тестирование сохранения и загрузки (save/load)\n";
  std::cout << "------------------------------------------------\n";
  const std::string save_filename = "mixture_params.txt";
  std::cout << "   - Сохраняем mix_params в файл '" << save_filename
            << "'...\n";
  mix_params.save(save_filename);

  std::cout << "   - Загружаем из файла в новый объект...\n";
  MixtureDist<MainDist, MainDist> mix_loaded;
  mix_loaded.load(save_filename);
  printMixState(mix_loaded, "mix_loaded");

  if (compareMixtures(mix_params, mix_loaded)) {
    std::cout
        << "   >>> ПОДТВЕРЖДЕНО: Сохраненные и загруженные данные совпадают!\n";
  } else {
    std::cout
        << "   >>> ОШИБКА: Данные после сохранения и загрузки не совпадают!\n";
  }
  std::cout << "\n";

  // ======================================================================
  // 3. ТЕСТИРОВАНИЕ ГЕТТЕРОВ И СЕТТЕРОВ
  // ======================================================================
  std::cout << "3. Тестирование геттеров и сеттеров\n";
  std::cout << "------------------------------------\n";
  MixtureDist mix_set(0.5, d1, d2);
  std::cout << "   - Начальное значение p: " << mix_set.getP() << "\n";
  mix_set.setP(0.9);
  std::cout << "   - Новое значение p после setP(0.9): " << mix_set.getP()
            << "\n";

  std::cout << "   - Проверка неконстантного геттера component1():\n";
  mix_set.component1().setMu(99.0);
  std::cout << "   - Значение mu у component1 после изменения через ссылку: "
            << mix_set.component1().getMu() << "\n";
  if (mix_set.component1().getMu() == 99.0) {
    std::cout
        << "   >>> ПОДТВЕРЖДЕНО: Неконстантный геттер работает корректно.\n";
  } else {
    std::cout << "   >>> ОШИБКА: Неконстантный геттер не работает.\n";
  }
  std::cout << "\n";

  // ======================================================================
  // 4. ТЕСТИРОВАНИЕ АНАЛИТИЧЕСКИХ МЕТОДОВ
  // ======================================================================
  std::cout << "4. Тестирование аналитических методов\n";
  std::cout << "-------------------------------------\n";
  printMixState(mix_params, "mix_params");
  std::cout << "   - Математическое ожидание (M): " << mix_params.M()
            << std::endl;
  std::cout << "   - Дисперсия (D): " << mix_params.D() << std::endl;
  std::cout << "   - Асимметрия (G1): " << mix_params.G1() << std::endl;
  std::cout << "   - Эксцесс (G2): " << mix_params.G2() << std::endl;
  ld x = 1.0;
  ld expected_density = (1 - mix_params.getP()) * d1.density(x) +
                        mix_params.getP() * d2.density(x);
  std::cout << "   - Плотность в точке x=" << x << ": " << mix_params.density(x)
            << " (ожидалось: " << expected_density << ")\n";
  std::cout << "\n";

  // ======================================================================
  // 5. ТЕСТИРОВАНИЕ ГЕНЕРАЦИИ СЛУЧАЙНЫХ ЧИСЕЛ (Xi)
  // ======================================================================
  std::cout << "5. Тестирование генерации случайных чисел (Xi)\n";
  std::cout << "----------------------------------------------\n";
  const unsigned int seed = 12345;
  const uint32_t sample_size = 5;

  std::cout << "   - Запуск 1 с seed = " << seed << "\n";
  srand(seed);
  nstu::vector sample1 = mix_params.Xi(sample_size);
  std::cout << "     Сгенерированные числа: ";
  for (const auto& val : sample1) std::cout << val << " ";
  std::cout << "\n";

  std::cout << "   - Запуск 2 с тем же seed = " << seed << "\n";
  srand(seed);  // Сброс генератора для воспроизводимости
  nstu::vector sample2 = mix_params.Xi(sample_size);
  std::cout << "     Сгенерированные числа: ";
  for (const auto& val : sample2) std::cout << val << " ";
  std::cout << "\n";

  if (areVectorsEqual(sample1, sample2)) {
    std::cout << "   >>> ПОДТВЕРЖДЕНО: Воспроизводимость генерации для "
                 "MixtureDist работает!\n";
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
    MixtureDist(-0.1, d1, d2);
  } catch (const std::invalid_argument& e) {
    std::cout << "   - Успех (p < 0 в конструкторе): " << e.what() << std::endl;
  }

  try {
    MixtureDist(1.1, d1, d2);
  } catch (const std::invalid_argument& e) {
    std::cout << "   - Успех (p > 1 в конструкторе): " << e.what() << std::endl;
  }

  MixtureDist<MainDist, MainDist> mix_exc;
  try {
    mix_exc.setP(-0.5);
  } catch (const std::invalid_argument& e) {
    std::cout << "   - Успех (p < 0 в сеттере): " << e.what() << std::endl;
  }

  try {
    mix_exc.setP(1.5);
  } catch (const std::invalid_argument& e) {
    std::cout << "   - Успех (p > 1 в сеттере): " << e.what() << std::endl;
  }

  // Тестируем, что исключение из MainDist пробрасывается через конструктор
  // MixtureDist
  std::ofstream bad_file("bad_mixture.txt");
  bad_file << "0.5\n5.0 -2.0 1.0\n3.0 2.0 0.0";  // lambda2 = 0
  bad_file.close();
  try {
    MixtureDist<MainDist, MainDist> mix;
    mix.load("bad_mixture.txt");
  } catch (const std::runtime_error& e) {
    std::cout << "   - Успех (невалидные параметры для MainDist в файле): "
              << e.what() << std::endl;
  }
  bad_file = std::ofstream("bad_mixture.txt");
  bad_file
      << "MixtureDist 0.5 MainDist 5.0 -2.0 1.0 MainDist 3.0 2.0 0.0";  // lambda2
                                                                        // = 0
  bad_file.close();
  try {
    MixtureDist<MainDist, MainDist> mix;
    mix.load("bad_mixture.txt");
  } catch (const std::invalid_argument& e) {
    std::cout << "   - Успех (невалидные параметры для MainDist в файле): "
              << e.what() << std::endl;
  }

  std::cout << "\n--- ПОЛНОЕ ТЕСТИРОВАНИЕ ЗАВЕРШЕНО ---\n";
}

int main() {
  try {
    testAllMixtureDistMethods();
  } catch (const std::exception& e) {
    std::cerr << "Критическая ошибка во время тестирования: " << e.what()
              << std::endl;
    return 1;
  }
  return 0;
}