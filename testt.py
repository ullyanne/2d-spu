import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import sys

class DraggableManager:
    def __init__(self, ax, rectangles):
        self.ax = ax
        self.rectangles = rectangles
        self.selected_rect = None
        self.offset = (0, 0)
        self.legend_visible = False
        self.legend = None

        # Conectar eventos
        self.cid_press = ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_motion = ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.cid_release = ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_key = ax.figure.canvas.mpl_connect('key_press_event', self.toggle_legend)
        self.cid_scroll = ax.figure.canvas.mpl_connect('scroll_event', self.on_scroll)  # Evento de zoom

    def on_press(self, event):
        if event.inaxes != self.ax:
            return
        for rect, text in self.rectangles:
            contains, _ = rect.contains(event)
            if contains:
                self.selected_rect = (rect, text)
                x0, y0 = rect.xy
                self.offset = (event.xdata - x0, event.ydata - y0)
                break

    def on_motion(self, event):
        if self.selected_rect and event.inaxes == self.ax:
            rect, text = self.selected_rect
            new_x = event.xdata - self.offset[0]
            new_y = event.ydata - self.offset[1]
            rect.set_xy((new_x, new_y))
            text.set_position((new_x + rect.get_width() / 2, new_y + rect.get_height() / 2))
            self.ax.figure.canvas.draw_idle()

    def on_release(self, event):
        self.selected_rect = None

    def toggle_legend(self, event):
        if event.key == 'l':  # Pressionar "l" para alternar legenda
            if self.legend_visible:
                self.legend.remove()
                self.legend_visible = False
            else:
                self.legend = self.ax.legend(handles=legend_patches, title="Clientes", loc='upper right')
                self.legend_visible = True
            self.ax.figure.canvas.draw_idle()

    def on_scroll(self, event):
        """Aplica zoom quando o mouse é rolado."""
        if event.inaxes != self.ax:
            return

        # Fator de zoom
        zoom_factor = 1.1 if event.button == 'down' else 0.9

        # Obtém os limites atuais dos eixos
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        # Calcula novos limites com base no fator de zoom
        new_xlim = [
            event.xdata - (event.xdata - xlim[0]) * zoom_factor,
            event.xdata + (xlim[1] - event.xdata) * zoom_factor
        ]
        new_ylim = [
            event.ydata - (event.ydata - ylim[0]) * zoom_factor,
            event.ydata + (ylim[1] - event.ydata) * zoom_factor
        ]

        # Aplica os novos limites
        self.ax.set_xlim(new_xlim)
        self.ax.set_ylim(new_ylim)

        # Redesenha o canvas
        self.ax.figure.canvas.draw_idle()

# Função para ler os dados do arquivo txt
def read_data_from_txt(file_path):
    with open(file_path, 'r') as file:
        lines = [line.strip() for line in file if line.strip()]
    fator_x = int(lines[0])
    data = np.array(lines[1:], dtype=str).reshape(-1, 4)

    rectangles_data = []
    for row in data:
        bottom_point = tuple(map(float, row[0].split()))
        top_point = tuple(map(float, row[1].split()))
        item_index = int(row[2])
        item_client = int(row[3])
        rectangles_data.append({"bottom_point": bottom_point, "top_point": top_point, "item_index": item_index, "item_client": item_client})
    return rectangles_data, fator_x

# Gerar cores únicas por cliente
# def generate_colors(clients):
#     cmap = plt.get_cmap("tab20")
#     num_colors = 20  # Usaremos todas as 20 cores da paleta
#     colors = {}

#     sorted_clients = sorted(clients)  # Ordenamos os clientes para manter consistência
#     used_colors = [-1]  # Começamos com um valor inválido para evitar problemas de índice

#     for i, client in enumerate(sorted_clients):
#         # Lista de cores disponíveis, evitando repetir a última usada
#         available_colors = [c for c in range(num_colors) if c != used_colors[-1]]

#         # Escolhemos a cor na ordem, garantindo que não repetimos a última usada
#         chosen_color = available_colors[i % len(available_colors)]
#         used_colors.append(chosen_color)

#         colors[client] = cmap(chosen_color / num_colors)

#     return colors

def generate_colors(clients):
    cmap = plt.get_cmap("viridis")
    return {client: cmap(i / len(clients)) for i, client in enumerate(clients)}


solfile = sys.argv[1]
file_path = f'/home/ullyanne/Documents/2spp/sol/bke/T60/N{solfile}BurkeSol.txt'

print(file_path)
rectangles_data, fator_x = read_data_from_txt(file_path)
clients = {rect["item_client"] for rect in rectangles_data}
colors = generate_colors(clients)
legend_patches = [patches.Patch(color=colors[client], label=f"Cliente {client}") for client in clients]

fig, ax = plt.subplots()
ax.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
rectangles = []

for rect_data in rectangles_data:
    color = colors[rect_data["item_client"]]
    bottom_point, top_point = rect_data["bottom_point"], rect_data["top_point"]
    width, height = top_point[0] - bottom_point[0], top_point[1] - bottom_point[1]
    rect = patches.Rectangle(bottom_point, width, height, linewidth=1, edgecolor='r', facecolor=color, alpha=0.5)
    ax.add_patch(rect)
    text = ax.text((bottom_point[0] + top_point[0]) / 2, (bottom_point[1] + top_point[1]) / 2, str(rect_data["item_index"]), ha='center', va='center', fontsize=10, color='black')
    rectangles.append((rect, text))

# Criar gerenciador de eventos
draggable_manager = DraggableManager(ax, rectangles)

# Definir limites do gráfico
all_x = np.array([p for rect in rectangles_data for p in (rect["bottom_point"][0], rect["top_point"][0])])
all_y = np.array([p for rect in rectangles_data for p in (rect["bottom_point"][1], rect["top_point"][1])])
margin_x, margin_y = (all_x.max() - all_x.min()) * 0.2, (all_y.max() - all_y.min()) * 0.2
ax.set_xlim(all_x.min() - margin_x, all_x.max() + margin_x + 200)
ax.set_ylim(all_y.min() - margin_y, all_y.max() + margin_y)
ax.set_aspect('equal', adjustable='box')
ax.plot(fator_x, 0, 'ro', markersize=8, label="Fator-X")

plt.show()