import tkinter as tk
from tkinter import filedialog, Button, Entry, Listbox, MULTIPLE, ttk
import comtrade
import json
import os
from processamentodados import processamento

class JanelaSelecaoArquivos:
    def __init__(self, root, parametros):
        self.root = root
        self.parametros = parametros
        self.arquivo1 = tk.StringVar(value=parametros.get('arquivo1', ''))
        self.arquivo2 = tk.StringVar(value=parametros.get('arquivo2', ''))
        self.freq_amostragem = tk.DoubleVar(value=parametros.get('freq_amostragem', 0.0))
        self.freq_corte_max = tk.DoubleVar(value=parametros.get('freq_corte_max', 0.0))
        self.R1 = tk.DoubleVar(value=parametros.get('dadoslinha', {}).get('R1', 0.0))
        self.X1 = tk.DoubleVar(value=parametros.get('dadoslinha', {}).get('X1', 0.0))
        self.R0 = tk.DoubleVar(value=parametros.get('dadoslinha', {}).get('R0', 0.0))
        self.X0 = tk.DoubleVar(value=parametros.get('dadoslinha', {}).get('X0', 0.0))
        self.L = tk.DoubleVar(value=parametros.get('dadoslinha', {}).get('L', 0.0))

        self.entries_caracteristicas = {}  # Para armazenar as caixas de entrada de características

        self.colunas = []  # Para armazenar os nomes das colunas

        self.criar_interface()

    def criar_interface(self):
        self.root.title("Selecionar Arquivos COMTRADE")
        self.root.state('zoomed')  # Maximiza a janela

        # Configura o restante da interface conforme já feito
        label_instrucao = tk.Label(self.root, text="Selecione ou insira manualmente os caminhos dos arquivos COMTRADE")
        label_instrucao.pack(pady=10)

        frame_arquivo1 = tk.Frame(self.root)
        frame_arquivo1.pack(pady=10)
        label_arquivo1 = tk.Label(frame_arquivo1, text="Arquivo .CFG:")
        label_arquivo1.grid(row=0, column=0, padx=5, sticky="w")
        self.entry_arquivo1 = tk.Entry(frame_arquivo1, textvariable=self.arquivo1, width=50)
        self.entry_arquivo1.grid(row=0, column=1, padx=5, sticky="w")
        button_selecionar1 = tk.Button(frame_arquivo1, text="Selecionar", command=self.selecionar_arquivo1)
        button_selecionar1.grid(row=0, column=2, padx=5)

        frame_arquivo2 = tk.Frame(self.root)
        frame_arquivo2.pack(pady=10)
        label_arquivo2 = tk.Label(frame_arquivo2, text="Arquivo .DAT:")
        label_arquivo2.grid(row=0, column=0, padx=5, sticky="w")
        self.entry_arquivo2 = tk.Entry(frame_arquivo2, textvariable=self.arquivo2, width=50)
        self.entry_arquivo2.grid(row=0, column=1, padx=5, sticky="w")
        button_selecionar2 = tk.Button(frame_arquivo2, text="Selecionar", command=self.selecionar_arquivo2)
        button_selecionar2.grid(row=0, column=2, padx=5)

        # Botão de Carregar Arquivo logo abaixo das entradas
        button_carregar = tk.Button(self.root, text="Carregar Arquivo", command=self.carregar_colunas)
        button_carregar.pack(pady=10)

        # Adiciona caixas de texto e seus rótulos com as características da linha
        frame_caracteristicas = tk.Frame(self.root)
        frame_caracteristicas.pack(pady=10)
        tk.Label(frame_caracteristicas, text="Características da Linha").grid(row=0, column=0, columnspan=3, pady=5)

        for i, (label, unit) in enumerate([("R1", "Ω/km"), ("X1", "Ω/km"), ("R0", "Ω/km"), ("X0", "Ω/km"), ("L", "km")]):
            tk.Label(frame_caracteristicas, text=f"{label}").grid(row=i+1, column=0, padx=5, sticky="w")
            entry = tk.Entry(frame_caracteristicas, width=15)
            entry.grid(row=i+1, column=1, padx=5, sticky="w")
            tk.Label(frame_caracteristicas, text=unit).grid(row=i+1, column=2, padx=5, sticky="w")
            self.entries_caracteristicas[label] = entry  # Salva a referência à caixa de entrada

        # Caixas de entrada de frequências
        frame_frequencias = tk.Frame(self.root)
        frame_frequencias.pack(pady=10)
        label_freq_amostragem = tk.Label(frame_frequencias, text="Frequência de Amostragem (Hz):")
        label_freq_amostragem.grid(row=0, column=0, padx=5, sticky="w")
        entry_freq_amostragem = tk.Entry(frame_frequencias, textvariable=self.freq_amostragem, width=15)
        entry_freq_amostragem.grid(row=0, column=1, padx=5, sticky="w")

        label_freq_max = tk.Label(frame_frequencias, text="Frequência de Corte Máxima (Hz):")
        label_freq_max.grid(row=1, column=0, padx=5, sticky="w")
        entry_freq_max = tk.Entry(frame_frequencias, textvariable=self.freq_corte_max, width=15)
        entry_freq_max.grid(row=1, column=1, padx=5, sticky="w")

        # Lista suspensa para número de terminais
        self.opcao_terminal = tk.StringVar(value="")
        label_opcao_terminal = tk.Label(self.root, text="Selecione o número de terminais:")
        label_opcao_terminal.pack(pady=10)
        self.combo_terminal = ttk.Combobox(self.root, textvariable=self.opcao_terminal, values=["1 Terminal", "2 Terminais"], state="readonly")
        self.combo_terminal.pack(pady=10)

        # Adicionando as listas suspensas (inicia com valor vazio)
        self.frames_listas = []
        self.listas_tensao = []
        self.listas_corrente = []

        for i in range(6):
            frame_lista = tk.Frame(self.root)
            self.frames_listas.append(frame_lista)
            lista_tensao = tk.StringVar()
            lista_corrente = tk.StringVar()
            self.listas_tensao.append(lista_tensao)
            self.listas_corrente.append(lista_corrente)
            rotulo_tensao = ["øA", "øB", "øC", "øA", "øB", "øC"][i]  # Nomes alternativos
            rotulo_corrente = ["øA", "øB", "øC", "øA", "øB", "øC"][i]  # Nomes alternativos
            tk.Label(frame_lista, text=f"Tensão {rotulo_tensao}").grid(row=0, column=0)
            opcao_tensao = ttk.Combobox(frame_lista, textvariable=lista_tensao, state="readonly")
            opcao_tensao.grid(row=0, column=1, padx=5)
            tk.Label(frame_lista, text=f"Corrente {rotulo_corrente}").grid(row=0, column=2)
            opcao_corrente = ttk.Combobox(frame_lista, textvariable=lista_corrente, state="readonly")
            opcao_corrente.grid(row=0, column=3, padx=5)
            frame_lista.opcao_tensao = opcao_tensao
            frame_lista.opcao_corrente = opcao_corrente

        # Adicionando Checkboxes
        self.checkbox_frame = tk.Frame(self.root)
        self.checkbox_frame.pack_forget()  # Esconde o frame inicialmente

        self.label_checkboxes = tk.Label(self.checkbox_frame, text="Seleção de gráficos:")
        self.label_checkboxes.grid(row=0, column=0, columnspan=3)

        # Criação dos checkboxes
        self.plotar_filtro = tk.BooleanVar(value=False)
        self.plotar_rms = tk.BooleanVar(value=False)
        self.plotar_fasores = tk.BooleanVar(value=False)

        self.checkbox_filtro = tk.Checkbutton(self.checkbox_frame, text="Sinal pós filtro", variable=self.plotar_filtro)
        self.checkbox_rms = tk.Checkbutton(self.checkbox_frame, text="RMS", variable=self.plotar_rms)
        self.checkbox_fasores = tk.Checkbutton(self.checkbox_frame, text="Fasores", variable=self.plotar_fasores)

        # Posiciona os checkboxes
        self.checkbox_filtro.grid(row=1, column=0, padx=5, pady=5)
        self.checkbox_rms.grid(row=1, column=1, padx=5, pady=5)
        self.checkbox_fasores.grid(row=1, column=2, padx=5, pady=5)

        # Botão Confirmar, inicialmente escondido
        self.botao_confirmar = tk.Button(self.root, text="Confirmar", command=self.confirmar)
        self.botao_confirmar.pack(pady=20)
        self.botao_confirmar.pack_forget()  # Esconde o botão

        # Configurando para chamar 'atualizar_interface' quando o valor da Combobox mudar
        self.combo_terminal.bind("<<ComboboxSelected>>", self.atualizar_interface)

    def atualizar_interface(self, *args):
        opcao = self.opcao_terminal.get()

        # Limpa todos os elementos antes de reconfigurar
        for frame in self.frames_listas:
            frame.pack_forget()
        if hasattr(self, 'label_terminal1'):
            self.label_terminal1.pack_forget()
        if hasattr(self, 'label_terminal2'):
            self.label_terminal2.pack_forget()
        self.botao_confirmar.pack_forget()

        if opcao == "1 Terminal":
            self.exibir_listas_terminal(3, titulo_terminal1=True)
        elif opcao == "2 Terminais":
            self.exibir_listas_terminal(6, titulo_terminal1=True, titulo_terminal2=True)

        # Exibe o botão Confirmar e o frame com os checkboxes
        self.botao_confirmar.pack(pady=20)
        self.checkbox_frame.pack(pady=10)

    def exibir_listas_terminal(self, num_listas, titulo_terminal1=False, titulo_terminal2=False):
        if titulo_terminal1:
            if not hasattr(self, 'label_terminal1'):
                self.label_terminal1 = tk.Label(self.root, text="Terminal 1", font=("Helvetica", 12, "bold"))
            self.label_terminal1.pack(pady=10)

        for i in range(min(num_listas, 3)):
            self.frames_listas[i].pack(pady=5)

        if titulo_terminal2:
            if not hasattr(self, 'label_terminal2'):
                self.label_terminal2 = tk.Label(self.root, text="Terminal 2", font=("Helvetica", 12, "bold"))
            self.label_terminal2.pack(pady=10)
            for i in range(3, num_listas):
                self.frames_listas[i].pack(pady=5)

    def carregar_colunas(self):
        arquivo1 = self.arquivo1.get()
        arquivo2 = self.arquivo2.get()
        freq_amostragem = float(self.freq_amostragem.get())  # Frequência de amostragem fornecida na interface
        freq_corte_max = float(self.freq_corte_max.get())    # Frequência de corte máxima fornecida na interface

        rec = comtrade.Comtrade()
        try:
            rec.load(arquivo1, arquivo2)
            print("001_____ARQUIVOS IMPORTADOS COM SUCESSO")
            self.colunas = rec.analog_channel_ids  # Obtendo os nomes das colunas
            self.atualizar_listas_suspensas()  # Atualiza as listas suspensas com os nomes das colunas
        except TypeError:
            print("Caminhos não foram definidos corretamente.")
            return

    def atualizar_listas_suspensas(self):
        for i in range(6):
            self.frames_listas[i].opcao_tensao['values'] = self.colunas
            self.frames_listas[i].opcao_corrente['values'] = self.colunas

    def selecionar_arquivo1(self):
        arquivo = filedialog.askopenfilename(title="Selecione o primeiro arquivo", filetypes=[("Arquivos CFG", "*.cfg")])
        if arquivo:
            self.arquivo1.set(arquivo)

    def selecionar_arquivo2(self):
        arquivo = filedialog.askopenfilename(title="Selecione o segundo arquivo", filetypes=[("Arquivos DAT", "*.dat")])
        if arquivo:
            self.arquivo2.set(arquivo)
    
    def confirmar(self):
        # Cria um dicionário para armazenar os valores de 'dadoslinha'
        dados_linha = {
            'R1': self.R1.get(),
            'X1': self.X1.get(),
            'R0': self.R0.get(),
            'X0': self.X0.get(),
            'L': self.L.get()
        }

        # Adiciona os dados ao dicionário de parâmetros
        self.parametros['dadoslinha'] = dados_linha

        # Atualiza os outros parâmetros
        self.parametros['arquivo1'] = self.arquivo1.get()
        self.parametros['arquivo2'] = self.arquivo2.get()
        self.parametros['freq_amostragem'] = self.freq_amostragem.get()
        self.parametros['freq_corte_max'] = self.freq_corte_max.get()

        # Coleta as colunas selecionadas das listas suspensas
        colunas_selecionadas = []
        for i in range(6):
            coluna_tensao = self.frames_listas[i].opcao_tensao.get()
            coluna_corrente = self.frames_listas[i].opcao_corrente.get()
            if coluna_tensao:
                colunas_selecionadas.append(coluna_tensao)
            if coluna_corrente:
                colunas_selecionadas.append(coluna_corrente)

        self.parametros['colunas'] = colunas_selecionadas

        # Adiciona os estados dos checkboxes aos parâmetros
        self.parametros['plotar_filtro'] = self.plotar_filtro.get()
        self.parametros['plotar_rms'] = self.plotar_rms.get()
        self.parametros['plotar_fasores'] = self.plotar_fasores.get()

        # Chama o método para salvar os parâmetros
        processamento.salvar_parametros(self.parametros)
        self.initial()


    def initial(self):
        processamento.process(self, self.parametros, self.colunas)
