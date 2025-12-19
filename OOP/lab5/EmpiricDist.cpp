#include "EmpiricDist.hpp"

#include <cmath>

EmpiricDist::EmpiricDist()
    : m_init_sample(), m_k(1), m_interval_freq(), m_min(0), m_max(0), m_delta(1.0),
      m_observers(nullptr), m_observers_count(0), m_observers_capacity(0) {}

EmpiricDist::EmpiricDist(const nstu::vector& init_sample, uint32_t k)
    : m_init_sample(init_sample),
      m_k((k == 0) ? floor(log2(init_sample.size())) + 1 : 1),
      m_interval_freq(m_k),
      m_observers(nullptr), m_observers_count(0), m_observers_capacity(0) {
  if (m_init_sample.empty())
    throw std::invalid_argument("Начальная выборка не может быть пустой");
  m_min = m_max = m_init_sample[0];
  for (const ld i : m_init_sample) {
    if (i < m_min) m_min = i;
    if (i > m_max) m_max = i;
  }
  m_delta = (m_max - m_min) / m_interval_freq.size();
  static int n = 0;
  if (n == 0) {
    std::ofstream file("logini.txt");
    file << m_init_sample.size() << "\n";
    file << m_k << "\n";
    file << m_min << "\n";
    file << m_max << "\n";
    file << m_delta << "\n";
    for (const ld i : m_interval_freq) file << i << " ";
    n++;
  }
  for (const ld x : m_init_sample) {
    uint32_t index = std::min(static_cast<uint32_t>((x - m_min) / m_delta),
                              m_interval_freq.size() - 1);
    m_interval_freq[index]++;
  }
  for (uint32_t i = 0; i < m_interval_freq.size(); i++) {
    m_interval_freq[i] /= m_init_sample.size();
  }

  if (m_init_sample.empty()) throw std::invalid_argument("Empty sample");
    m_min = m_max = m_init_sample[0];
    for (const ld i : m_init_sample) {
      if (i < m_min) m_min = i;
      if (i > m_max) m_max = i;
    }
    m_delta = (m_max - m_min) / m_interval_freq.size();
    if (m_delta == 0) m_delta = 1.0; // Защита от деления на 0

    for (const ld x : m_init_sample) {
      uint32_t index = std::min(static_cast<uint32_t>((x - m_min) / m_delta),
                                (uint32_t)m_interval_freq.size() - 1);
      m_interval_freq[index]++;
    }
    for (uint32_t i = 0; i < m_interval_freq.size(); i++) {
      m_interval_freq[i] /= m_init_sample.size();
    }
}

EmpiricDist::EmpiricDist(uint32_t n, IDist& dist, uint32_t k)
    : EmpiricDist(dist.Xi(n), k) {}

void EmpiricDist::load(const std::string& file_name) {
  std::ifstream file(file_name);
  if (!file.is_open())
    throw std::runtime_error("Не удалось открыть файл " + file_name);
  load(file);
  file.close();
}

// Реализация Subject
void EmpiricDist::resizeObservers() {
  uint32_t new_cap = (m_observers_capacity == 0) ? 2 : m_observers_capacity * 2;
  IObserver** new_arr = new IObserver*[new_cap];
  for(uint32_t i = 0; i < m_observers_count; ++i) {
    new_arr[i] = m_observers[i];
  }
  delete[] m_observers;
  m_observers = new_arr;
  m_observers_capacity = new_cap;
}

void EmpiricDist::attach(IObserver* obs) {
  if (m_observers_count == m_observers_capacity) {
    resizeObservers();
  }
  m_observers[m_observers_count++] = obs;
}

void EmpiricDist::detach(IObserver* obs) {
  for (uint32_t i = 0; i < m_observers_count; ++i) {
    if (m_observers[i] == obs) {
      // Сдвигаем элементы
      for (uint32_t j = i; j < m_observers_count - 1; ++j) {
        m_observers[j] = m_observers[j + 1];
      }
      m_observers_count--;
      return;
    }
  }
}

void EmpiricDist::notify() const {
  for (uint32_t i = 0; i < m_observers_count; ++i) {
    if (m_observers[i])
      m_observers[i]->update();
  }
}

// Изменение данных (Засорение)
void EmpiricDist::setSample(const nstu::vector& new_sample) {
  m_init_sample = new_sample;
  // Пересчет параметров гистограммы
  *this = EmpiricDist(m_init_sample, m_k);
  
  // Уведомляем наблюдателей
  notify();
}

void EmpiricDist::load(std::ifstream& file) {
  std::string name;
  file >> name;
  if (name != "EmpiricDist")
    throw std::runtime_error("Неверный тип распределения");
  file >> m_k;
  uint32_t sample_size;
  file >> sample_size;
  if (sample_size < 1)
    throw std::invalid_argument(
        "Размер начальной выборки должен быть больше 0");
  m_init_sample = nstu::vector(sample_size);
  for (ld& i : m_init_sample) {
    file >> i;
  }
  *this = EmpiricDist(m_init_sample, m_k);
}

void EmpiricDist::save(const std::string& file_name) const {
  std::ofstream file(file_name);
  if (!file.is_open())
    throw std::runtime_error("Не удалось открыть файл " + file_name);
  save(file);
  file.close();
}

void EmpiricDist::save(std::ofstream& file) const noexcept {
  file << "EmpiricDist ";
  file << m_k << " ";
  file << m_init_sample.size() << " ";
  for (const ld i : m_init_sample) {
    file << i << " ";
  }
}

EmpiricDist::EmpiricDist(const EmpiricDist& other) 
    : m_init_sample(other.m_init_sample), m_k(other.m_k),
      m_interval_freq(other.m_interval_freq), m_min(other.m_min),
      m_max(other.m_max), m_delta(other.m_delta),
      m_observers(nullptr), m_observers_count(0), m_observers_capacity(0) {
    // Наблюдателей не копируем, новая копия дистрибутива не имеет подписчиков
}

EmpiricDist& EmpiricDist::operator=(const EmpiricDist& other) {
  if (this != &other) {
    m_init_sample = other.m_init_sample;
    m_k = other.m_k;
    m_interval_freq = other.m_interval_freq;
    m_min = other.m_min;
    m_max = other.m_max;
    m_delta = other.m_delta;
  }
  return *this;
}

EmpiricDist::EmpiricDist(EmpiricDist&& other) noexcept
    : m_init_sample(std::move(other.m_init_sample)),
      m_k(other.m_k),
      m_interval_freq(std::move(other.m_interval_freq)),
      m_min(other.m_min),
      m_max(other.m_max),
      m_delta(other.m_delta),
      m_observers(other.m_observers),
      m_observers_count(other.m_observers_count),
      m_observers_capacity(other.m_observers_capacity)  {
  other.m_k = 0;
  other.m_min = 0;
  other.m_max = 0;
  other.m_delta = 0;
  other.m_observers = nullptr;
  other.m_observers_count = 0;
  other.m_observers_capacity = 0;
}

EmpiricDist& EmpiricDist::operator=(EmpiricDist&& other) noexcept {
  if (this != &other) {
    delete[] m_observers;

    m_init_sample = std::move(other.m_init_sample);
    m_k = other.m_k;
    m_interval_freq = std::move(other.m_interval_freq);
    m_min = other.m_min;
    m_max = other.m_max;
    m_delta = other.m_delta;

    m_observers = other.m_observers;
    m_observers_count = other.m_observers_count;
    m_observers_capacity = other.m_observers_capacity;

    other.m_k = 0;
    other.m_min = 0;
    other.m_max = 0;
    other.m_delta = 0;
    other.m_observers = nullptr;
    other.m_observers_count = 0;
    other.m_observers_capacity = 0;
  }
  return *this;
}

const nstu::vector& EmpiricDist::getInitSample() const { return m_init_sample; }

ld EmpiricDist::getK() const noexcept { return m_k; }

ld EmpiricDist::getMin() const noexcept { return m_min; }

ld EmpiricDist::getMax() const noexcept { return m_max; }

void EmpiricDist::setK(uint32_t k = 0) {
  m_k = (k == 0) ? floor(log2(m_init_sample.size())) + 1 : k;
  *this = EmpiricDist(m_init_sample, m_k);
}

void EmpiricDist::setInitSample(const nstu::vector& sample) {
  m_init_sample = sample;
  *this = EmpiricDist(m_init_sample, m_k);
}

ld EmpiricDist::density(ld x) const noexcept {
  if (x < m_min || x > m_max) return 0;
  uint32_t index =
      std::min(static_cast<uint32_t>((x - m_min) / m_delta), m_k - 1);
  static int n = 0;
  if (n == 0) {
    std::ofstream file("log" + std::to_string(n++) + ".txt");
    for (uint32_t i = 0; i < m_k; ++i) {
      file << m_interval_freq[i] << "\n";
    }
  }
  return m_interval_freq[index] / m_delta;
}

nstu::vector EmpiricDist::density(const nstu::vector& x) const {
  nstu::vector result(x.size());
  for (uint32_t i = 0; i < result.size(); ++i) {
    result[i] = density(x[i]);
  }
  return result;
}

ld EmpiricDist::M() const noexcept {
  ld sum = 0;
  for (const ld x : m_init_sample) sum += x;
  return sum / m_init_sample.size();
}

ld EmpiricDist::D() const noexcept {
  ld m = M();
  ld sum = 0;
  for (const ld x : m_init_sample) sum += std::pow(x - m, 2);
  return sum / m_init_sample.size();
}

ld EmpiricDist::G1() const noexcept {
  ld m = M();
  ld d = D();
  ld sum = 0;
  for (const ld x : m_init_sample) sum += std::pow(x - m, 3);
  return (sum / m_init_sample.size()) / std::pow(d, 1.5);
}

ld EmpiricDist::G2() const noexcept {
  ld m = M();
  ld d = D();
  ld sum = 0;
  for (const ld x : m_init_sample) sum += std::pow(x - m, 4);
  return (sum / m_init_sample.size()) / (d * d) - 3;
}

ld EmpiricDist::Xi() const noexcept {
  ld r = nstu::randNum();
  ld cumulative = 0;
  for (uint32_t i = 0; i < m_k; ++i) {
    cumulative += m_interval_freq[i];
    if (r < cumulative) {
      return m_min + i * m_delta + m_delta * nstu::randNum();
    }
  }
  return m_max;
}

nstu::vector EmpiricDist::Xi(uint32_t size) const {
  nstu::vector result(size);
  for (ld& i : result) i = Xi();
  return result;
}
