import tkinter as tk
from tkinter import filedialog, Button, Entry, Listbox, MULTIPLE, Toplevel, OptionMenu
import comtrade
import json
import os
from plot import plotar_sinais
from processamentodados import processamento

class JanelaSelecaoArquivos:
    def __init__(self, root, parametros):
        self.root = root
        self.arquivo1 = tk.StringVar(value=parametros['arquivo1'])
        self.arquivo2 = tk.StringVar(value=parametros['arquivo2'])
        self.freq_amostragem = tk.DoubleVar(value=parametros['freq_amostragem'])
        self.freq_corte_max = tk.DoubleVar(value=parametros['freq_corte_max'])
        self.parametros = parametros
        self.colunas = []  # Para armazenar os nomes das colunas

        self.criar_interface()

    def criar_interface(self):
        self.root.title("Selecionar Arquivos COMTRADE")

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

        # Adicionando as checkboxes
        self.opcao = tk.IntVar(value=0)
        self.frame_checkboxes = tk.Frame(self.root)
        self.frame_checkboxes.pack(padx=10, pady=10, expand=True)
        self.checkbox1 = tk.Checkbutton(self.root, text="1 Terminal", variable=self.opcao, onvalue=1, offvalue=0, command=self.atualizar_interface)
        self.checkbox1.pack(padx=10, pady=10, anchor = "center")
        self.checkbox2 = tk.Checkbutton(self.root, text="2 Terminais", variable=self.opcao, onvalue=2, offvalue=0, command=self.atualizar_interface)
        self.checkbox2.pack(padx=10, pady=10, anchor = "center")
        self.frame_checkboxes.pack_propagate(False)
        

        # Adicionando as listas suspensas (desabilitadas inicialmente)
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
            tk.Label(frame_lista, text=f"Tensão {i+1}").grid(row=0, column=0)
            opcao_tensao = tk.OptionMenu(frame_lista, lista_tensao, "")
            opcao_tensao.grid(row=0, column=1, padx=5)
            tk.Label(frame_lista, text=f"Corrente {i+1}").grid(row=0, column=2)
            opcao_corrente = tk.OptionMenu(frame_lista, lista_corrente, "")
            opcao_corrente.grid(row=0, column=3, padx=5)
            frame_lista.opcao_tensao = opcao_tensao
            frame_lista.opcao_corrente = opcao_corrente

        # Botão de Carregar
        button_carregar = tk.Button(self.root, text="Carregar Arquivo", command=self.carregar_colunas)
        button_carregar.pack(pady=10)

        # self.label_tipo_falta = tk.Label(self.root, text="Tipo de Falta Identificada:")
        # self.label_tipo_falta.pack(pady=10)

        # self.label_porcentagem_falta = tk.Label(self.root, text="Porcentagem da Linha de Transmissão:")
        # self.label_porcentagem_falta.pack(pady=10)

        botao_confirmar = tk.Button(self.root, text="Confirmar", command=self.confirmar)
        botao_confirmar.pack(pady=20)

    def atualizar_interface(self):
        opcao = self.opcao.get()
        
        for i, frame in enumerate(self.frames_listas):
            if opcao == 1 and i < 3:
                frame.pack(pady=5)
            elif opcao == 2:
                frame.pack(pady=5)
            else:
                frame.pack_forget()

    def carregar_colunas(self):
        # Carrega os arquivos e extrai as colunas conforme o código fornecido
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
            # Atualizando os menus de Tensão e Corrente
            menu_tensao = self.frames_listas[i].opcao_tensao["menu"]
            menu_corrente = self.frames_listas[i].opcao_corrente["menu"]
            
            menu_tensao.delete(0, 'end')
            menu_corrente.delete(0, 'end')
            
            for coluna in self.colunas:
                menu_tensao.add_command(label=coluna, command=tk._setit(self.listas_tensao[i], coluna))
                menu_corrente.add_command(label=coluna, command=tk._setit(self.listas_corrente[i], coluna))

    def selecionar_arquivo1(self):
        arquivo = filedialog.askopenfilename(title="Selecione o primeiro arquivo", filetypes=[("Arquivos CFG", "*.cfg")])
        if arquivo:
            self.arquivo1.set(arquivo)

    def selecionar_arquivo2(self):
        arquivo = filedialog.askopenfilename(title="Selecione o segundo arquivo", filetypes=[("Arquivos DAT", "*.dat")])
        if arquivo:
            self.arquivo2.set(arquivo)

    def confirmar(self):
        self.parametros['arquivo1'] = self.arquivo1.get()
        self.parametros['arquivo2'] = self.arquivo2.get()
        self.parametros['freq_amostragem'] = self.freq_amostragem.get()
        self.parametros['freq_corte_max'] = self.freq_corte_max.get()
        processamento.salvar_parametros(self.parametros)
        self.abrir_janela_selecao_colunas()

    def abrir_janela_selecao_colunas(self):
        self.janela_colunas = Toplevel(self.root)
        self.janela_colunas.title("Selecionar Colunas para Plotagem")
        label_instrucao = tk.Label(self.janela_colunas, text="Selecione as colunas que deseja plotar:")
        label_instrucao.pack(pady=10)

        self.listbox_colunas = Listbox(self.janela_colunas, selectmode=MULTIPLE, width=50)
        self.listbox_colunas.pack(pady=10)

        for i, canal in enumerate(self.colunas):
            self.listbox_colunas.insert(tk.END, f"{i} - {canal}")

        botao_confirmar_colunas = tk.Button(self.janela_colunas, text="Confirmar", command=self.confirmar_colunas)
        botao_confirmar_colunas.pack(pady=20)

    def confirmar_colunas(self):
        colunas_selecionadas = [int(i) for i in self.listbox_colunas.curselection()]
        max_index = len(self.colunas) - 1
        colunas_validas = [i for i in colunas_selecionadas if 0 <= i <= max_index]
        if len(colunas_validas) != len(colunas_selecionadas):
            print("Alguns índices selecionados são inválidos e foram removidos.")
        self.parametros['colunas'] = colunas_validas
        processamento.salvar_parametros(self.parametros)
        self.janela_colunas.destroy()
        self.initial()

    def initial(self):
        plotar_sinais(self.parametros, self.colunas, self.label_tipo_falta, self.label_porcentagem_falta)

if __name__ == "__main__":
    parametros = {
        'arquivo1': '',
        'arquivo2': '',
        'freq_amostragem': 0.0,
        'freq_corte_max': 0.0
    }

    root = tk.Tk()
    app = JanelaSelecaoArquivos(root, parametros)
    root.mainloop()
