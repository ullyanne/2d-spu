import os

def modify_txt_file(file_path):
    # Leitura do conteúdo do arquivo
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Deleta a primeira linha
    del lines[0]

    # Remove o primeiro número da terceira linha e concatena
    lines[2] = lines[2].split()[1] + "\n"
    # Escreve o conteúdo modificado de volta no arquivo
    with open(file_path, 'w') as file:
        file.writelines(lines)

def modify_all_txt_files_in_directory(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            if file.lower().endswith('.txt'):
                file_path = os.path.join(root, file)
                modify_txt_file(file_path)
                print(f'Arquivo {file_path} modificado com sucesso!')

# Exemplo de uso:
if __name__ == "__main__":
    directory_path = 'instances'  # Substitua pelo caminho da sua pasta principal
    modify_all_txt_files_in_directory(directory_path)
    print(f'Todos os arquivos .txt em {directory_path} foram modificados com sucesso!')
