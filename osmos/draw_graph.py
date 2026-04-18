import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# Получаем имя эксперимента: либо из аргументов командной строки, либо спрашиваем
if len(sys.argv) == 2:
    exp_name = sys.argv[1]
else:
    exp_name = input("Введите название эксперимента: ")

# Путь к файлу с данными
csv_file = os.path.join("cmake-build-debug", "graphics", exp_name, "pressure_diff.csv")

# Проверяем, существует ли файл
if not os.path.exists(csv_file):
    print(f"Ошибка: файл {csv_file} не найден.")
    print("Убедитесь, что симуляция завершена и имя введено правильно.")
    sys.exit(1)

# Читаем данные
data = pd.read_csv(csv_file)

# Создаём график
fig, ax1 = plt.subplots(figsize=(10, 6))

# Разность давлений (левая ось)
ax1.plot(data['time_ps'], data['dP'], color='red', linewidth=1.5, label='dP (bar)')
ax1.set_xlabel('Время (пс)')
ax1.set_ylabel('dP (bar)', color='red')
ax1.tick_params(axis='y', labelcolor='red')
ax1.grid(True, linestyle='--', alpha=0.5)

# Температура (правая ось)
ax2 = ax1.twinx()
ax2.plot(data['time_ps'], data['T'], color='blue', linewidth=1.5, linestyle='--', label='T (K)')
ax2.set_ylabel('T (K)', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')

# Заголовок с именем эксперимента
plt.title(f'Осмотическое давление и температура — {exp_name}')

# Легенда
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

# Сохраняем картинку в ту же папку
output_path = os.path.join("cmake-build-debug", "graphics", exp_name, "pressure_plot.png")
plt.savefig(output_path, dpi=150, bbox_inches='tight')
print(f"График сохранён: {output_path}")

# Показываем на экране
plt.show()
