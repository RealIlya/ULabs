import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

def plot_lab5_results():
    # Чтение данных из файла
    try:
        with open("debug_data.txt", "r") as f:
            # Читаем первую строку с параметрами
            header = f.readline().split()
            mu_est = float(header[0])
            scale_est = float(header[1])
            c_param = float(header[2])
            
            # Читаем данные: index, value, weight
            indices = []
            values = []
            weights = []
            for line in f:
                parts = line.split()
                indices.append(int(parts[0]))
                values.append(float(parts[1]))
                weights.append(float(parts[2]))
    except FileNotFoundError:
        print("Ошибка: Файл debug_data.txt не найден. Запустите сначала C++ программу.")
        return

    # --- ГРАФИК 1: Плотности и Весовая функция ---
    plt.figure(figsize=(10, 6))
    
    # Диапазон x для построения линий
    x = np.linspace(-4, 10, 1000)
    
    # Параметры распределений (из задания)
    p = 0.25 # уровень засорения
    
    # Идеальное N(0, 1)
    pdf_ideal = norm.pdf(x, 0, 1)
    plt.plot(x, pdf_ideal, 'r-', label='Идеальное (N(0,1))', linewidth=1.5)
    
    # Засоряющее N(5, 1)
    pdf_poison = norm.pdf(x, 5, 1)
    plt.plot(x, pdf_poison, 'b-', label='Засоряющее (N(5,1))', linewidth=1.5)
    
    # Смесь (Засоренное)
    pdf_mixture = (1 - p) * pdf_ideal + p * pdf_poison
    plt.plot(x, pdf_mixture, 'g-', label='Смесь (Засоренная)', linewidth=2)
    
    # Весовая функция Коши
    # w(z) = 1 / (1 + (z/c)^2), где z = (x - mu)/s
    z = (x - mu_est) / scale_est
    w_values = 1 / (1 + (z / c_param)**2)
    
    # Масштабируем весовую функцию для наглядности (чтобы была высотой с графики плотности)
    # Обычно умножают на макс. высоту основного графика (около 0.4)
    scaling_factor = 0.4 
    plt.plot(x, w_values * scaling_factor, 'm--', label=f'Весовая функция (c={c_param}) [scaled]', linewidth=2)

    plt.title('Плотности распределений и весовая функция Коши')
    plt.xlabel('x')
    plt.ylabel('Плотность / Вес (усл.)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('plot_densities.png') # Сохранить картинку
    plt.show()

    # --- ГРАФИК 2: Веса наблюдений (точки) ---
    plt.figure(figsize=(10, 6))
    
    # Рисуем точки
    plt.scatter(indices, weights, s=10, c='blue', alpha=0.6, label='Веса наблюдений')
    
    plt.title(f'Веса наблюдений (c={c_param})')
    plt.xlabel('Индекс наблюдения')
    plt.ylabel('Вес w_i')
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    
    # Добавим линию отсечки (опционально)
    plt.axhline(y=0, color='black', linewidth=0.5)
    
    plt.legend()
    plt.savefig('plot_weights.png') # Сохранить картинку
    plt.show()

if __name__ == "__main__":
    plot_lab5_results()