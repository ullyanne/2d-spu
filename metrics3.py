import os
import re

def extrair_valores(content, padrao):
    valores = [float(match.group(1)) for match in re.finditer(padrao, '\n'.join(content))]
    return valores[-20:]  # Pegar as últimas 20 ocorrências

def calcular_estatisticas(valores):
    if valores:
        return min(valores), sum(valores) / len(valores), max(valores)
    return None, None, None

def process_logs(directory):
    alturas = []
    tempos = []

    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".log"):
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.readlines()

                # Extrair as últimas 20 alturas e tempos de execução
                alturas.extend(extrair_valores(content, r"Avaliação = (\d+\.?\d*)"))
                tempos.extend(extrair_valores(content, r"Tempo de execução: (\d+\.?\d*)"))

    # Calcular e exibir estatísticas
    for categoria, valores in [("Altura", alturas), ("Tempo de execução", tempos)]:
        menor, media, maior = calcular_estatisticas(valores)
        if menor is not None:
            print(f"{categoria} - Menor: {menor:.2f}, Média: {media:.2f}, Maior: {maior:.2f}")
        else:
            print(f"Nenhum dado encontrado para {categoria}")

if __name__ == "__main__":
    diretorio_base = "/home/ullyanne/Documents/2spp/logs/bke/T20"
    process_logs(diretorio_base)
