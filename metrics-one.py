import os
import re
from decimal import Decimal, getcontext

# Define a precisão (20 casas decimais)
getcontext().prec = 20

def process_logs(directory):
    avaliacoes = []
    tempos_execucao = []

    # Percorrer recursivamente os arquivos no diretório
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".log"):
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.readlines()

                # Procurar pela última ocorrência de "Avaliação = {número}"
                for line in reversed(content):
                    match_avaliacao = re.search(r"Avaliação = (\d+\.?\d*)", line)
                    if match_avaliacao:
                        avaliacoes.append(Decimal(match_avaliacao.group(1)))
                        break

                # Procurar pela última ocorrência de "Tempo de execução: {número}"
                for line in reversed(content):
                    match_tempo = re.search(r"Tempo de execução: (\d+\.?\d*)", line)
                    if match_tempo:
                        tempos_execucao.append(Decimal(match_tempo.group(1)))
                        break

    # Calcular e exibir as médias
    if avaliacoes:
        media_avaliacao = sum(avaliacoes) / Decimal(len(avaliacoes))
        print(f"Média da Avaliação: {media_avaliacao:.10f}")
    else:
        print("Nenhuma avaliação encontrada nos arquivos .log")

    if tempos_execucao:
        media_tempo = sum(tempos_execucao) / Decimal(len(tempos_execucao))
        print(f"Média do Tempo de Execução: {media_tempo:.10f} segundos")
    else:
        print("Nenhum tempo de execução encontrado nos arquivos .log")

if __name__ == "__main__":
    diretorio_base = "/home/ullyanne/Documents/2spp/logs/bke/T20/"
    process_logs(diretorio_base)
