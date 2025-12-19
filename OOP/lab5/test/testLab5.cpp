#include <iostream>
#include <cmath>
#include <iomanip>
#include "../MainDist.hpp"
#include "../EmpiricDist.hpp"
#include "../CauchyEstimate.hpp"

// Функция для генерации засоренных данных
nstu::vector contaminate(const nstu::vector& clean, const IDist& poison_dist, double percent) {
  nstu::vector contaminated = clean;
  size_t n = clean.size();
  size_t n_poison = static_cast<size_t>(n * percent);
  
  // Заменяем последние n_poison элементов
  // (в реальной задаче можно выбирать случайно, но для теста сойдет конец)
  for (size_t i = n - n_poison; i < n; ++i) {
    contaminated[i] = poison_dist.Xi();
  }
  return contaminated;
}

int main() {
  std::cout << "--- Observer And M-Estimation ---" << std::endl;
  
  // 1. Подготовка распределений
  // Идеальное: N(0, 1) -> Nu=infinity (approx 100), Mu=0, Lambda=1
  MainDist ideal_dist(100.0, 0.0, 1.0); 
  
  // Засоряющее: Сдвинутое N(5, 1)
  MainDist poison_dist(100.0, 5.0, 1.0);

  // Размер выборки
  uint32_t N = 500;
  
  // Генерируем чистую выборку
  nstu::vector clean_sample = ideal_dist.Xi(N);

  // 2. Создаем Субъект (EmpiricDist)
  EmpiricDist empiric(clean_sample);
  std::cout << "Чистая выборка создана. Истинное среднее ~ 0." << std::endl;
  std::cout << "Среднее арифметическое (Empiric): " << empiric.M() << std::endl;

  // 3. Создаем Наблюдателя (CauchyEstimate)
  // Параметры c для теста: 0.75, 1.0, 1.725, 2.385
  ld test_c = 1.0;
  CauchyEstimate m_estimator(empiric, test_c);
  
  std::cout << "M-оценка Коши (c=" << test_c << "): " << m_estimator.getEstimate() << std::endl;
  std::cout << "------------------------------------------" << std::endl;

  // 4. Тест паттерна Observer (Засорение данных)
  std::cout << ">>> Применяем засорение 25% данными из N(5, 1)..." << std::endl;
  
  nstu::vector dirty_sample = contaminate(clean_sample, poison_dist, 0.25);
  
  // Здесь срабатывает Observer
  // Мы меняем данные в empiric, а m_estimator должен сам пересчитать значение.
  empiric.setSample(dirty_sample);

  std::cout << "Данные обновлены в EmpiricDist." << std::endl;
  std::cout << "Новое среднее арифметическое (чувствительно к выбросам): " << empiric.M() << std::endl;
  std::cout << "Новая M-оценка (автоматически пересчитана): " << m_estimator.getEstimate() << std::endl;
  
  if (std::abs(m_estimator.getEstimate()) < std::abs(empiric.M())) {
    std::cout << "РЕЗУЛЬТАТ: M-оценка оказалась устойчивее!" << std::endl;
  } else {
    std::cout << "РЕЗУЛЬТАТ: Что-то пошло не так." << std::endl;
  }
  std::cout << "------------------------------------------" << std::endl;

  // 5. Вычислительный эксперимент (подбор C)
  std::cout << ">>> Подбор оптимального параметра C (MAE)" << std::endl;
  
  ld c_values[] = {0.75, 1.0, 1.725, 2.385};
  ld true_mu = 0.0; // Мы знаем, что идеал это 0

  for (ld c : c_values) {
    m_estimator.setC(c);  // Меняем параметр, пересчет происходит внутри
    ld est = m_estimator.getEstimate();
    ld error = std::abs(est - true_mu);
    std::cout << "c = " << std::setw(5) << c 
              << " | Estimate = " << std::setw(8) << est 
              << " | Error = " << error << std::endl;
  }

  // Создаем объект с оптимальным параметром (из теста c=0.75)
  empiric.setSample(dirty_sample);
  ld best_c = 0.75;
  CauchyEstimate estimator(empiric, best_c);

  std::cout << "Расчет с оптимальным C=" << best_c << ". Оценка: " << estimator.getEstimate() << std::endl;

  estimator.saveDebugData("debug_data.txt");

  return 0;
}