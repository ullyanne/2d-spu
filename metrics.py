import os
import re
from decimal import Decimal, getcontext
import sys

getcontext().prec = 50

def extract_values(content, pattern):
    values = [Decimal(match.group(1)) for match in re.finditer(pattern, '\n'.join(content))]
    return values[-20:]

def calculate_average(values):
    return sum(values) / Decimal(len(values)) if values else None

def process_folder(directory):
    min_heights, avg_heights, max_heights = [], [], []
    avg_times = []

    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".log"):
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.readlines()

                heights = extract_values(content, r"Avaliação = (\d+\.?\d*)")
                times = extract_values(content, r"Tempo de execução: (\d+\.?\d*)")

                if heights:
                    min_heights.append(min(heights))
                    avg_heights.append(calculate_average(heights))
                    max_heights.append(max(heights))
                if times:
                    sum_times = sum(times)
                    avg_times.append(sum_times / Decimal(len(times)))

    return (
        min_heights, avg_heights, max_heights, avg_times
    )

def process_all_folders(base_dir):
    folders = [folder for folder in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, folder))]

    all_min_heights, all_avg_heights, all_max_heights = [], [], []
    all_avg_times = []

    for folder in sorted(folders):
        path = os.path.join(base_dir, folder)
        print(f"\n📂 Processing folder: {folder}")

        min_heights, avg_heights, max_heights, avg_times = [], [], [], []

        for file in sorted(os.listdir(path)):
            if file.endswith(".log"):
                file_path = os.path.join(path, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.readlines()

                heights = extract_values(content, r"Avaliação = (\d+\.?\d*)")
                times = extract_values(content, r"Tempo de execução: (\d+\.?\d*)")

                if heights:
                    min_heights.append(min(heights))
                    avg_heights.append(calculate_average(heights))
                    max_heights.append(max(heights))
                if times:
                    sum_times = sum(times)
                    avg_times.append(sum_times / Decimal(len(times)))

        if min_heights:
            avg_min = calculate_average(min_heights)
            avg_avg = calculate_average(avg_heights)
            avg_max = calculate_average(max_heights)
            avg_time = calculate_average(avg_times)

            print(f"Height ({folder}) - Average of minimums: {avg_min:.20f}, "
                  f"Average of averages: {avg_avg:.20f}, Average of maximums: {avg_max:.20f}")
            print(f"Execution time ({folder}) - Average of times: {avg_time:.20f}")

            all_min_heights.extend(min_heights)
            all_avg_heights.extend(avg_heights)
            all_max_heights.extend(max_heights)
            all_avg_times.extend(avg_times)

    print("\n📊 Global statistics:")
    if all_min_heights:
        print(f"Height (Global) - Average of minimums: {calculate_average(all_min_heights):.20f}, "
              f"Average of averages: {calculate_average(all_avg_heights):.20f}, "
              f"Average of maximums: {calculate_average(all_max_heights):.20f}")
    if all_avg_times:
        print(f"Execution time (Global) - Average of times: {calculate_average(all_avg_times):.20f}")

if __name__ == "__main__":
    package_name = sys.argv[1]
    if package_name == "2lcvrp":
        package_name = ""

    base_directory = "./logs/" + package_name
    process_all_folders(base_directory)
