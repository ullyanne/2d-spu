import os
import re

def process_logs(directory):
    melhor_alturas = []

    # Percorrer recursivamente os arquivos no diretório
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".log"):
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.readlines()

                # Procurar pela última ocorrência de "Melhor altura = {número}"
                for line in reversed(content):
                    match = re.search(r"Avaliação = (\d+\.?\d*)", line)
                    if match:
                        melhor_alturas.append(float(match.group(1)))
                        break

    # Calcular a média, o menor e o maior valor
    if melhor_alturas:
        media = sum(melhor_alturas) / len(melhor_alturas)
        menor = min(melhor_alturas)
        maior = max(melhor_alturas)
        print(f"Menor: {menor:.2f}, Média: {media:.2f}, Maior: {maior:.2f}")
    else:
        print("Nenhuma altura encontrada nos arquivos .log")

if __name__ == "__main__":
    diretorio_base = "/home/ullyanne/Documents/2spp/logs/2lcvrp"
    process_logs(diretorio_base)
