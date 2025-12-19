#include "CauchyEstimate.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

CauchyEstimate::CauchyEstimate(EmpiricDist& data, ld c) 
  : m_data(data), m_c(c), m_mu(0), m_scale(1) {
  m_data.attach(this); // Подписываемся на уведомления
  estimate();          // Первичный расчет
}

CauchyEstimate::~CauchyEstimate() {
  m_data.detach(this); // Отписываемся при уничтожении
}

void CauchyEstimate::update() {
  // Этот метод вызывается автоматически, когда EmpiricDist вызывает notify()
  // std::cout << "Наблюдатель уведомлен! Пересчет оценки...\n"; 
  estimate();
}

void CauchyEstimate::bubbleSort(nstu::vector& v) {
  size_t n = v.size();
  for (size_t i = 0; i < n - 1; ++i) {
    for (size_t j = 0; j < n - i - 1; ++j) {
      if (v[j] > v[j + 1]) {
        ld temp = v[j];
        v[j] = v[j + 1];
        v[j + 1] = temp;
      }
    }
  }
}

ld CauchyEstimate::calculateMedian(const nstu::vector& sample) {
  if (sample.empty()) return 0;
  // Делаем копию для сортировки, чтобы не менять исходные данные
  nstu::vector sorted_sample = sample; 
  bubbleSort(sorted_sample);

  size_t n = sorted_sample.size();
  if (n % 2 == 0) {
    return (sorted_sample[n / 2 - 1] + sorted_sample[n / 2]) / 2.0;
  } else {
    return sorted_sample[n / 2];
  }
}

ld CauchyEstimate::calculateMAD(const nstu::vector& sample, ld median) {
  nstu::vector abs_diffs(sample.size());
  for (size_t i=0; i < sample.size(); ++i) {
    abs_diffs[i] = std::abs(sample[i] - median);
  }
  // MAD = median(|x_i - median(X)|)
  return calculateMedian(abs_diffs);
}

ld CauchyEstimate::cauchyWeight(ld z) const {
  // Формула: w(z) = c^2 / (c^2 + z^2)
  // Это то же самое, что 1 / (1 + (z/c)^2)
  
  if (m_c == 0) return 0;
  
  // z - это уже нормированный остаток (x - mu)/s
  // Нам нужно применить формулу Коши
  
  ld c2 = m_c * m_c;
  ld z2 = z * z;
  
  return c2 / (c2 + z2);
}

void CauchyEstimate::estimate() {
  const nstu::vector& X = m_data.getInitSample();
  if (X.empty()) return;

  // Начальное приближение - Медиана
  ld current_mu = calculateMedian(X);
  
  // Оценка масштаба S = 1.4826 * MAD
  ld mad = calculateMAD(X, current_mu);
  m_scale = 1.4826 * mad;
  
  if (m_scale < 1e-9) {
    // Если разброс нулевой, оценка равна медиане
    m_mu = current_mu;
    return;
  }

  // Итерационный процесс
  int max_iter = 50;
  const ld EPS = 1e-6;

  for (int k = 0; k < max_iter; ++k) {
    ld sum_w_xi = 0;
    ld sum_w = 0;

    for (size_t i = 0; i < X.size(); ++i) {
      ld z = (X[i] - current_mu) / m_scale;
      ld w = cauchyWeight(z);
      sum_w_xi += w * X[i];
      sum_w += w;
    }

    if (sum_w == 0) break; 

    ld next_mu = sum_w_xi / sum_w;

    if (std::abs(next_mu - current_mu) < EPS) {
      current_mu = next_mu;
      break;
    }
    current_mu = next_mu;
  }

  m_mu = current_mu;
}

ld CauchyEstimate::getEstimate() const {
  return m_mu;
}

void CauchyEstimate::setC(ld c) {
  m_c = c;
  estimate(); // Пересчитываем с новым параметром
}

void CauchyEstimate::saveDebugData(const std::string& filename) const {
  const nstu::vector& X = m_data.getInitSample();
  
  std::ofstream file(filename);
  if (!file.is_open()) return;

  // Первая строка: метаданные (mu, scale, c)
  file << m_mu << " " << m_scale << " " << m_c << "\n";

  // Остальные строки: index value weight
  for (size_t i = 0; i < X.size(); ++i) {
    ld z = (X[i] - m_mu) / m_scale;
    ld w = cauchyWeight(z);
    file << i << " " << X[i] << " " << w << "\n";
  }
  
  file.close();
  std::cout << "Данные для графиков сохранены в " << filename << std::endl;
}