import tkinter as tk
from tkinter import filedialog, Button, Entry, Listbox, MULTIPLE, Toplevel
import comtrade
import json
import os
from plot import plotar_sinais
from processamentodados import salvar_parametros, process_signal, calculate_vrms, detectar_tipo_falta

class JanelaSelecaoArquivos:
    def __init__(self, root, parametros):
        self.root = root
        self.arquivo1 = tk.StringVar(value=parametros['arquivo1'])
        self.arquivo2 = tk.StringVar(value=parametros['arquivo2'])
        self.freq_amostragem = tk.DoubleVar(value=parametros['freq_amostragem'])
        self.freq_corte_max = tk.DoubleVar(value=parametros['freq_corte_max'])
        self.parametros = parametros

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

        self.label_tipo_falta = tk.Label(self.root, text="Tipo de Falta Identificada:")
        self.label_tipo_falta.pack(pady=10)

        self.label_porcentagem_falta = tk.Label(self.root, text="Porcentagem da Linha de Transmissão:")
        self.label_porcentagem_falta.pack(pady=10)

        botao_confirmar = tk.Button(self.root, text="Confirmar", command=self.confirmar)
        botao_confirmar.pack(pady=20)

    def selecionar_arquivo1(self):
        arquivo = filedialog.askopenfilename(title="Selecione o primeiro arquivo")
        if arquivo:
            self.arquivo1.set(arquivo)

    def selecionar_arquivo2(self):
        arquivo = filedialog.askopenfilename(title="Selecione o segundo arquivo")
        if arquivo:
            self.arquivo2.set(arquivo)

    def confirmar(self):
        self.parametros['arquivo1'] = self.arquivo1.get()
        self.parametros['arquivo2'] = self.arquivo2.get()
        self.parametros['freq_amostragem'] = self.freq_amostragem.get()
        self.parametros['freq_corte_max'] = self.freq_corte_max.get()
        salvar_parametros(self.parametros)
        self.abrir_janela_selecao_colunas()

    def abrir_janela_selecao_colunas(self):
        self.janela_colunas = Toplevel(self.root)
        self.janela_colunas.title("Selecionar Colunas para Plotagem")
        label_instrucao = tk.Label(self.janela_colunas, text="Selecione as colunas que deseja plotar:")
        label_instrucao.pack(pady=10)

        self.listbox_colunas = Listbox(self.janela_colunas, selectmode=MULTIPLE, width=50)
        self.listbox_colunas.pack(pady=10)

        rec = comtrade.Comtrade()
        rec.load(self.parametros['arquivo1'], self.parametros['arquivo2'])
        self.canais = rec.analog_channel_ids
        for i, canal in enumerate(self.canais):
            self.listbox_colunas.insert(tk.END, f"{i} - {canal}")

        botao_confirmar_colunas = tk.Button(self.janela_colunas, text="Confirmar", command=self.confirmar_colunas)
        botao_confirmar_colunas.pack(pady=20)

    def confirmar_colunas(self):
        colunas_selecionadas = [int(i) for i in self.listbox_colunas.curselection()]
        max_index = len(self.canais) - 1
        colunas_validas = [i for i in colunas_selecionadas if 0 <= i <= max_index]
        if len(colunas_validas) != len(colunas_selecionadas):
            print("Alguns índices selecionados são inválidos e foram removidos.")
        self.parametros['colunas'] = colunas_validas
        salvar_parametros(self.parametros)
        self.janela_colunas.destroy()
        self.plotar_sinais()

    def plotar_sinais(self):
        plotar_sinais(self.parametros, self.canais, self.label_tipo_falta, self.label_porcentagem_falta)
