#pragma once
#include "EmpiricDist.hpp"
#include "IObserver.hpp"
#include <cmath>
#include <algorithm>

class CauchyEstimate : public IObserver {
 private:
  EmpiricDist& m_data;   // Агрегация по ссылке
  ld m_c;                // Параметр весовой функции
  ld m_mu;               // Текущая оценка (сдвиг)
  ld m_scale;            // Оценка масштаба (MAD)

  void bubbleSort(nstu::vector& v);
  
  // Медиана
  ld calculateMedian(const nstu::vector& sample);
  
  // MAD (Median Absolute Deviation)
  ld calculateMAD(const nstu::vector& sample, ld median);

  // Весовая функция Коши
  ld cauchyWeight(ld z) const;

  // Основной метод пересчета
  void estimate();
    
 public:
  CauchyEstimate(EmpiricDist& data, ld c);
  ~CauchyEstimate();

  // Реализация интерфейса IObserver
  void update() override;

  // Геттер оценки
  ld getEstimate() const;
  
  // Изменить параметр c
  void setC(ld c);

  void saveDebugData(const std::string& filename) const;
};