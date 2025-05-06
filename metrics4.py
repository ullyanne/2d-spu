import os
import re
from decimal import Decimal, getcontext

# Define a precisão que você quiser. Aqui está altíssima:
getcontext().prec = 50

def extrair_valores(content, padrao):
    """Extrai os últimos 20 valores que correspondem ao padrão."""
    valores = [Decimal(match.group(1)) for match in re.finditer(padrao, '\n'.join(content))]
    return valores[-20:]  # Pegamos apenas os últimos 20 valores

def calcular_media(valores):
    """Calcula a média de uma lista de valores."""
    return sum(valores) / Decimal(len(valores)) if valores else None

def processar_pasta(directory):
    """Processa todos os logs dentro de uma pasta e retorna as estatísticas da pasta."""
    minimos_alturas, medias_alturas, maximos_alturas = [], [], []
    medias_tempos = []

    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".log"):
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.readlines()

                alturas = extrair_valores(content, r"Avaliação = (\d+\.?\d*)")
                tempos = extrair_valores(content, r"Tempo de execução: (\d+\.?\d*)")

                if alturas:
                    minimos_alturas.append(min(alturas))
                    medias_alturas.append(calcular_media(alturas))
                    maximos_alturas.append(max(alturas))
                if tempos:
                    soma_tempos = sum(tempos)
                    medias_tempos.append(soma_tempos / Decimal(len(tempos)))  # média dos 20 tempos

    return (
        minimos_alturas, medias_alturas, maximos_alturas, medias_tempos
    )

def processar_todas_as_pastas(base_dir):
    """Processa todas as pastas dentro do diretório base e calcula estatísticas globais."""
    pastas = [pasta for pasta in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, pasta))]

    todos_minimos_alturas, todas_medias_alturas, todos_maximos_alturas = [], [], []
    todas_medias_tempos = []

    for pasta in sorted(pastas):
        caminho = os.path.join(base_dir, pasta)
        print(f"\n📂 Processando pasta: {pasta}")

        minimos_alturas, medias_alturas, maximos_alturas, medias_tempos = [], [], [], []

        for file in sorted(os.listdir(caminho)):
            if file.endswith(".log"):
                file_path = os.path.join(caminho, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.readlines()

                alturas = extrair_valores(content, r"Avaliação = (\d+\.?\d*)")
                tempos = extrair_valores(content, r"Tempo de execução: (\d+\.?\d*)")

                if alturas:
                    minimos_alturas.append(min(alturas))
                    medias_alturas.append(calcular_media(alturas))
                    maximos_alturas.append(max(alturas))
                    media_altura = calcular_media(alturas)
                if tempos:
                    soma_tempos = sum(tempos)
                    medias_tempos.append(soma_tempos / Decimal(len(tempos)))

        if minimos_alturas:
            media_minimos = calcular_media(minimos_alturas)
            media_medias = calcular_media(medias_alturas)
            media_maximos = calcular_media(maximos_alturas)
            media_tempos = calcular_media(medias_tempos)

            print(f"Altura ({pasta}) - Média dos menores: {media_minimos:.20f}, "
                  f"Média das médias: {media_medias:.20f}, Média dos maiores: {media_maximos:.20f}")
            print(f"Tempo de execução ({pasta}) - Média dos tempos: {media_tempos:.20f}")

            todos_minimos_alturas.extend(minimos_alturas)
            todas_medias_alturas.extend(medias_alturas)
            todos_maximos_alturas.extend(maximos_alturas)
            todas_medias_tempos.extend(medias_tempos)

    # Cálculo das estatísticas globais
    print("\n📊 Estatísticas globais:")
    if todos_minimos_alturas:
        print(f"Altura (Global) - Média dos menores: {calcular_media(todos_minimos_alturas):.20f}, "
              f"Média das médias: {calcular_media(todas_medias_alturas):.20f}, "
              f"Média dos maiores: {calcular_media(todos_maximos_alturas):.20f}")
    if todas_medias_tempos:
        print(f"Tempo de execução (Global) - Média dos tempos: {calcular_media(todas_medias_tempos):.20f}")

if __name__ == "__main__":
    diretorio_base = "/home/ullyanne/Documents/2spp/logs/"
    processar_todas_as_pastas(diretorio_base)
