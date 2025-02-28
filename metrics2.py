import os
import re

def process_logs(directory):
    melhor_alturas = []
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
                    match = re.search(r"Avaliação = (\d+\.?\d*)", line)
                    if match:
                        melhor_alturas.append(float(match.group(1)))
                        break

                # Procurar pela última ocorrência de "Tempo de execução: {número}"
                for line in reversed(content):
                    match = re.search(r"Tempo de execução: (\d+\.?\d*)", line)
                    if match:
                        tempos_execucao.append(float(match.group(1)))
                        break

    # Calcular a média, o menor e o maior valor
    if melhor_alturas:
        media = sum(melhor_alturas) / len(melhor_alturas)
        menor = min(melhor_alturas)
        maior = max(melhor_alturas)
        media_tempo = sum(tempos_execucao) / len(tempos_execucao) if tempos_execucao else 0
        return media, menor, maior, media_tempo
    else:
        return None

def process_all_directories(base_directory):
    todas_medias = []
    todas_menores = []
    todas_maiores = []
    todos_tempos = []

    for subdir in sorted(os.listdir(base_directory)):
        subdir_path = os.path.join(base_directory, subdir)
        if os.path.isdir(subdir_path):
            resultado = process_logs(subdir_path)
            if resultado:
                media, menor, maior, media_tempo = resultado
                todas_medias.append(media)
                todas_menores.append(menor)
                todas_maiores.append(maior)
                todos_tempos.append(media_tempo)
                print(f"{subdir} -> Menor: {menor:.2f}, Média: {media:.2f},  Maior: {maior:.2f}, Média Tempo: {media_tempo:.2f}s")

    if todas_medias:
        media_geral = sum(todas_medias) / len(todas_medias)
        menor_geral = sum(todas_menores) / len(todas_menores)
        maior_geral = sum(todas_maiores) / len(todas_maiores)
        media_tempo_geral = sum(todos_tempos) / len(todos_tempos) if todos_tempos else 0
        print("\nMédias Gerais:")
        print(f"Média do Menor: {menor_geral:.2f}")
        print(f"Média da Média: {media_geral:.2f}")
        print(f"Média do Maior: {maior_geral:.2f}")
        print(f"Média do Tempo de Execução: {media_tempo_geral:.2f}s")
    else:
        print("Nenhum dado encontrado em nenhum diretório.")

if __name__ == "__main__":
    diretorio_base = "/home/ullyanne/Documents/2spp/logs/bke"
    process_all_directories(diretorio_base)
