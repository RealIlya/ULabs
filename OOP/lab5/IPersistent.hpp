#pragma once

#include <fstream>

#include "../lab2/header.h"

class IPersistent {
 public:
  /**
   * Сохранение в файл
   * @param file_name Имя файла
   */
  virtual void save(const std::string& file_name) const = 0;
  /**
   * Чтение из файла
   * @param file_name Имя файла
   */
  virtual void load(const std::string& file_name) = 0;
  /**
   * Сохранение в файл
   * @param file Поток открытого файла
   */
  virtual void save(std::ofstream& file) const noexcept = 0;
  /**
   * Чтение из файла
   * @param file Поток открытого файла
   */
  virtual void load(std::ifstream& file) = 0;
};