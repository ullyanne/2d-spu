import os
import re

def extrair_valores(content, padrao):
    """Extrai os últimos 20 valores que correspondem ao padrão."""
    valores = [float(match.group(1)) for match in re.finditer(padrao, '\n'.join(content))]
    return valores[-20:]  # Pegamos apenas os últimos 20 valores

def calcular_media(valores):
    """Calcula a média de uma lista de valores."""
    return sum(valores) / len(valores) if valores else None

def processar_pasta(directory):
    """Processa todos os logs dentro de uma pasta e retorna as estatísticas da pasta."""
    medias_alturas = []
    medias_tempos = []

    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".log"):
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.readlines()

                # Extrair os últimos 20 valores e calcular a média para cada log
                alturas = extrair_valores(content, r"Avaliação = (\d+\.?\d*)")
                tempos = extrair_valores(content, r"Tempo de execução: (\d+\.?\d*)")

                if alturas:
                    medias_alturas.append(calcular_media(alturas))
                if tempos:
                    medias_tempos.append(calcular_media(tempos))

    # Calcular estatísticas para a pasta
    def calcular_estatisticas(medias):
        return (min(medias), calcular_media(medias), max(medias)) if medias else (None, None, None)

    return calcular_estatisticas(medias_alturas), calcular_estatisticas(medias_tempos)

def processar_todas_as_pastas(base_dir):
    """Processa todas as pastas dentro do diretório base e calcula estatísticas globais."""
    pastas = [pasta for pasta in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, pasta))]

    todos_minimos_alturas, todas_medias_alturas, todos_maximos_alturas = [], [], []
    todos_minimos_tempos, todas_medias_tempos, todos_maximos_tempos = [], [], []

    for pasta in sorted(pastas):
        caminho = os.path.join(base_dir, pasta)
        print(f"\n📂 Processando pasta: {pasta}")

        estatisticas_alturas, estatisticas_tempos = processar_pasta(caminho)

        # Exibir estatísticas da pasta
        for categoria, (menor, media, maior) in [("Altura", estatisticas_alturas), ("Tempo de execução", estatisticas_tempos)]:
            if menor is not None:
                print(f"{categoria} ({pasta}) - Menor: {menor:.2f}, Média: {media:.2f}, Maior: {maior:.2f}")

        # Armazena os valores para o cálculo global
        if estatisticas_alturas[0] is not None:
            todos_minimos_alturas.append(estatisticas_alturas[0])
            todas_medias_alturas.append(estatisticas_alturas[1])
            todos_maximos_alturas.append(estatisticas_alturas[2])

        if estatisticas_tempos[0] is not None:
            todos_minimos_tempos.append(estatisticas_tempos[0])
            todas_medias_tempos.append(estatisticas_tempos[1])
            todos_maximos_tempos.append(estatisticas_tempos[2])

    # Cálculo das estatísticas globais
    print("\n📊 Estatísticas globais:")
    for categoria, minimos, medias, maximos in [
        ("Altura", todos_minimos_alturas, todas_medias_alturas, todos_maximos_alturas),
        ("Tempo de execução", todos_minimos_tempos, todas_medias_tempos, todos_maximos_tempos),
    ]:
        if minimos:
            print(f"{categoria} (Global) - Menor dos menores: {min(minimos):.2f}, Média geral: {calcular_media(medias):.2f}, Maior dos maiores: {max(maximos):.2f}")

if __name__ == "__main__":
    diretorio_base = "/home/ullyanne/Documents/2spp/logs/bke"
    processar_todas_as_pastas(diretorio_base)
