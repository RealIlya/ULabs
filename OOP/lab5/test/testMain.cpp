#include <cmath>    // для fabs
#include <fstream>  // для работы с файлами
#include <iomanip>  // для std::fixed, std::setprecision
#include <iostream>
#include <stdexcept>  // для std::invalid_argument

#include "../MainDist.hpp"

// Вспомогательная функция для печати состояния объекта MainDist
void printDistState(const MainDist& dist, const std::string& name) {
  std::cout << "   - Состояние " << name
            << " (nu, mu, lambda): " << dist.getNu() << ", " << dist.getMu()
            << ", " << dist.getLambda() << std::endl;
}

// Вспомогательная функция для сравнения двух векторов
bool areVectorsEqual(const nstu::vector& v1, const nstu::vector& v2) {
  if (v1.size() != v2.size()) {
    return false;
  }
  for (size_t i = 0; i < v1.size(); ++i) {
    // Сравниваем числа с плавающей точкой с небольшой погрешностью
    if (std::fabs(v1[i] - v2[i]) > 1e-9) {
      return false;
    }
  }
  return true;
}

void testAllMainDistMethods() {
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "--- НАЧАЛО ПОЛНОГО ТЕСТИРОВАНИЯ КЛАССА MAINDIST ---\n\n";

  // ======================================================================
  // 1. ТЕСТИРОВАНИЕ КОНСТРУКТОРОВ И ГЕТТЕРОВ
  // ======================================================================
  std::cout << "1. Тестирование конструкторов и геттеров\n";
  std::cout << "-------------------------------------------\n";

  // 1.1. Конструктор с параметрами
  std::cout << "1.1. Конструктор с параметрами: MainDist(5.0, 1.0, 2.0)\n";
  MainDist dist_params(5.0, 1.0, 2.0);
  printDistState(dist_params, "dist_params");

  // 1.2. Конструктор из файла
  std::cout << "\n1.2. Конструктор из файла (test_params.txt)\n";
  MainDist dist_file_params(5.0, 1.0, 2.0);
  dist_file_params.save("test_params.txt");
  MainDist dist_file;
  dist_file.load("test_params.txt");
  printDistState(dist_file_params, "dist_file_params");
  printDistState(dist_file, "dist_file");

  // 1.3. Конструктор с initializer_list
  std::cout << "\n1.3. Конструктор с initializer_list: {5.0, 1.0, 2.0}\n";
  MainDist dist_list{5.0, 1.0, 2.0};
  printDistState(dist_list, "dist_list");

  // 1.4. Конструктор копирования
  std::cout
      << "\n1.4. Конструктор копирования: MainDist dist_copy = dist_params\n";
  MainDist dist_copy = dist_params;
  printDistState(dist_copy, "dist_copy");
  // Убедимся, что это глубокая копия (хотя для ld это очевидно)
  dist_params.setNu(99.0);
  std::cout << "   Изменили nu у оригинала до 99.0\n";
  printDistState(dist_params, "dist_params (измененный)");
  printDistState(dist_copy, "dist_copy (неизмененный)");
  dist_params.setNu(5.0);  // Вернем обратно

  // 1.5. Конструктор перемещения
  std::cout
      << "\n1.5. Конструктор перемещения: MainDist dist_move = MainDist(...)\n";
  MainDist dist_move = MainDist(10.0, -5.0, 1.0);
  printDistState(dist_move, "dist_move");
  // Состояние перемещенного объекта не определено стандартом, но в вашей
  // реализации оно обнуляется

  std::cout << "\n";

  // ======================================================================
  // 2. ТЕСТИРОВАНИЕ ОПЕРАТОРОВ ПРИСВАИВАНИЯ
  // ======================================================================
  std::cout << "2. Тестирование операторов присваивания\n";
  std::cout << "--------------------------------------\n";

  // 2.1. Оператор присваивания копированием
  std::cout << "2.1. Присваивание копированием\n";
  MainDist assign_target(1, 1, 1);
  MainDist assign_source(8, 8, 8);
  printDistState(assign_target, "assign_target (до)");
  printDistState(assign_source, "assign_source");
  assign_target = assign_source;
  printDistState(assign_target, "assign_target (после)");

  // 2.2. Оператор присваивания перемещением
  std::cout << "\n2.2. Присваивание перемещением\n";
  MainDist move_target(1, 1, 1);
  MainDist move_source(9, 9, 9);
  printDistState(move_target, "move_target (до)");
  printDistState(move_source, "move_source");
  move_target = std::move(move_source);
  printDistState(move_target, "move_target (после)");
  printDistState(move_source, "move_source (после перемещения, обнулен)");

  std::cout << "\n";

  // ======================================================================
  // 3. ТЕСТИРОВАНИЕ СЕТТЕРОВ
  // ======================================================================
  std::cout << "3. Тестирование сеттеров\n";
  std::cout << "------------------------\n";
  MainDist dist_set(1, 1, 1);
  printDistState(dist_set, "dist_set (начальное)");
  dist_set.setNu(10.0);
  dist_set.setMu(20.0);
  dist_set.setLambda(30.0);
  printDistState(dist_set, "dist_set (конечное)");

  std::cout << "\n";

  // ======================================================================
  // 4. ТЕСТИРОВАНИЕ АНАЛИТИЧЕСКИХ МЕТОДОВ
  // ======================================================================
  std::cout << "4. Тестирование аналитических методов\n";
  std::cout << "-------------------------------------\n";
  MainDist dist_analysis(5.0, 1.0, 2.0);
  printDistState(dist_analysis, "dist_analysis");
  std::cout << "   - Математическое ожидание (M): " << dist_analysis.M()
            << std::endl;
  std::cout << "   - Дисперсия (D): " << dist_analysis.D() << std::endl;
  std::cout << "   - Асимметрия (G1): " << dist_analysis.G1() << std::endl;
  std::cout << "   - Эксцесс (G2): " << dist_analysis.G2() << std::endl;
  ld x_point = 1.5;
  std::cout << "   - Плотность в точке x = " << x_point << ": "
            << dist_analysis.density(x_point) << std::endl;

  nstu::vector x_vec = {0.5, 1.0, 1.5, 2.0};
  nstu::vector densities = dist_analysis.density(x_vec);
  std::cout << "   - Плотности для вектора {0.5, 1.0, 1.5, 2.0}:\n     { ";
  for (const auto& d : densities) {
    std::cout << d << " ";
  }
  std::cout << "}\n\n";

  // ======================================================================
  // 5. ТЕСТИРОВАНИЕ ГЕНЕРАЦИИ СЛУЧАЙНЫХ ЧИСЕЛ (Xi)
  // ======================================================================
  std::cout << "5. Тестирование генерации случайных чисел (Xi)\n";
  std::cout << "----------------------------------------------\n";
  const unsigned int seed = 12345;
  const uint32_t sample_size = 5;

  srand(seed);
  ld single_xi = dist_analysis.Xi();
  std::cout << "   - Одиночное случайное число Xi(): " << single_xi << "\n";

  srand(seed);  // Сброс для воспроизводимости
  nstu::vector sample1 = dist_analysis.Xi(sample_size);
  std::cout << "   - Вектор Xi(" << sample_size << ") [Запуск 1]: { ";
  for (const auto& val : sample1) std::cout << val << " ";
  std::cout << "}\n";

  srand(seed);  // Снова сброс
  nstu::vector sample2 = dist_analysis.Xi(sample_size);
  std::cout << "   - Вектор Xi(" << sample_size << ") [Запуск 2]: { ";
  for (const auto& val : sample2) std::cout << val << " ";
  std::cout << "}\n";

  if (areVectorsEqual(sample1, sample2)) {
    std::cout << "   >>> ПОДТВЕРЖДЕНО: Воспроизводимость генерации работает!\n";
  } else {
    std::cout << "   >>> ОШИБКА: Воспроизводимость генерации НЕ работает!\n";
  }
  std::cout << "\n";

  // ======================================================================
  // 6. ТЕСТИРОВАНИЕ ОБРАБОТКИ ИСКЛЮЧЕНИЙ
  // ======================================================================
  std::cout << "6. Тестирование обработки исключений\n";
  std::cout << "------------------------------------\n";

  // 6.1. lambda = 0 в конструкторе
  try {
    MainDist(1, 1, 0);
  } catch (const std::invalid_argument& e) {
    std::cout << "   - Успех: перехвачено исключение для lambda=0: " << e.what()
              << std::endl;
  }

  // 6.2. Неверное число аргументов в initializer_list
  try {
    MainDist{1.0, 2.0};
  } catch (const std::invalid_argument& e) {
    std::cout
        << "   - Успех: перехвачено исключение для неверного числа аргументов: "
        << e.what() << std::endl;
  }

  // 6.3. Не удалось открыть файл
  try {
    MainDist tmp;
    tmp.load("non_existent_file.txt");
  } catch (const std::runtime_error& e) {
    std::cout
        << "   - Успех: перехвачено исключение для несуществующего файла: "
        << e.what() << std::endl;
  }

  // 6.4. lambda = 0 в файле
  std::ofstream bad_file("bad_params.txt");
  bad_file << "5.0 1.0 0.0";
  bad_file.close();
  try {
    MainDist tmp;
    tmp.load("bad_params.txt");
  } catch (const std::runtime_error& e) {
    std::cout
        << "   - Успех: перехвачено исключение для неверного типа в файле: "
        << e.what() << std::endl;
  }
  // 6.4. lambda = 0 в файле
  bad_file = std::ofstream("bad_params.txt");
  bad_file << "MainDist 5.0 1.0 0.0";
  bad_file.close();
  try {
    MainDist tmp;
    tmp.load("bad_params.txt");
  } catch (const std::invalid_argument& e) {
    std::cout << "   - Успех: перехвачено исключение для lambda=0 в файле: "
              << e.what() << std::endl;
  }
  // 6.5. Установка lambda = 0 через сеттер
  MainDist dist_exc(1, 1, 1);
  try {
    dist_exc.setLambda(0);
  } catch (const std::invalid_argument& e) {
    std::cout << "   - Успех (setLambda(0)): " << e.what() << std::endl;
  }

  // 6.6. Установка nu = 0 через сеттер
  try {
    dist_exc.setNu(0);
  } catch (const std::invalid_argument& e) {
    std::cout << "   - Успех (setNu(0)): " << e.what() << std::endl;
  }
  std::cout << "\n--- ПОЛНОЕ ТЕСТИРОВАНИЕ ЗАВЕРШЕНО ---\n";
}

int main() {
  try {
    testAllMainDistMethods();
  } catch (const std::exception& e) {
    std::cerr << "Критическая ошибка во время тестирования: " << e.what()
              << std::endl;
    return 1;
  }
  return 0;
}