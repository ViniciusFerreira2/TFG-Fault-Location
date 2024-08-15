import numpy as np
import cmath
from scipy.fft import fft, ifft, fftfreq
from scipy.signal import butter, filtfilt, resample
import matplotlib.pyplot as plt
import comtrade
import tkinter as tk
from tkinter import filedialog, Button, Entry, Listbox, MULTIPLE, Toplevel
import json
import os

# Nome do arquivo de configuração
CONFIG_FILE = 'config.json'

def carregar_parametros():
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as f:
            return json.load(f)
    return {
        'arquivo1': '',
        'arquivo2': '',
        'freq_corte_min': 45.0,
        'freq_corte_max': 75.0,
        'colunas': []
    }

def salvar_parametros(parametros):
    with open(CONFIG_FILE, 'w') as f:
        json.dump(parametros, f)

class JanelaSelecaoArquivos:
    def __init__(self, root, parametros):
        self.root = root
        self.arquivo1 = tk.StringVar(value=parametros['arquivo1'])
        self.arquivo2 = tk.StringVar(value=parametros['arquivo2'])
        self.freq_corte_min = tk.DoubleVar(value=parametros['freq_corte_min'])
        self.freq_corte_max = tk.DoubleVar(value=parametros['freq_corte_max'])
        self.parametros = parametros

        self.criar_interface()

    def criar_interface(self):
        self.root.title("Selecionar Arquivos COMTRADE")

        label_instrucao = tk.Label(self.root, text="Selecione ou insira manualmente os caminhos dos arquivos COMTRADE")
        label_instrucao.pack(pady=10)

        frame_arquivo1 = tk.Frame(self.root)
        frame_arquivo1.pack(pady=10)
        label_arquivo1 = tk.Label(frame_arquivo1, text="Primeiro arquivo:")
        label_arquivo1.grid(row=0, column=0, padx=5, sticky="w")
        self.entry_arquivo1 = tk.Entry(frame_arquivo1, textvariable=self.arquivo1, width=50)
        self.entry_arquivo1.grid(row=0, column=1, padx=5, sticky="w")
        button_selecionar1 = tk.Button(frame_arquivo1, text="Selecionar", command=self.selecionar_arquivo1)
        button_selecionar1.grid(row=0, column=2, padx=5)

        frame_arquivo2 = tk.Frame(self.root)
        frame_arquivo2.pack(pady=10)
        label_arquivo2 = tk.Label(frame_arquivo2, text="Segundo arquivo:")
        label_arquivo2.grid(row=0, column=0, padx=5, sticky="w")
        self.entry_arquivo2 = tk.Entry(frame_arquivo2, textvariable=self.arquivo2, width=50)
        self.entry_arquivo2.grid(row=0, column=1, padx=5, sticky="w")
        button_selecionar2 = tk.Button(frame_arquivo2, text="Selecionar", command=self.selecionar_arquivo2)
        button_selecionar2.grid(row=0, column=2, padx=5)

        frame_frequencias = tk.Frame(self.root)
        frame_frequencias.pack(pady=10)
        label_freq_min = tk.Label(frame_frequencias, text="Frequência de Corte Mínima (Hz):")
        label_freq_min.grid(row=0, column=0, padx=5, sticky="w")
        entry_freq_min = tk.Entry(frame_frequencias, textvariable=self.freq_corte_min, width=15)
        entry_freq_min.grid(row=0, column=1, padx=5, sticky="w")

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
        self.parametros['freq_corte_min'] = self.freq_corte_min.get()
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

        # Preencher Listbox com os nomes das colunas dos dados COMTRADE
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

    def downsample_signal(self, signal, original_rate, target_rate):
        num_samples = int(len(signal) * target_rate / original_rate)
        downsampled_signal = resample(signal, num_samples)
        return downsampled_signal

    def butter_lowpass_filter(self, data, cutoff, fs, order=2):
        nyquist = 0.5 * fs
        normal_cutoff = cutoff / nyquist
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        y = filtfilt(b, a, data)
        return y

    def process_signal(self, signal, original_rate=1000000, target_rate=1920, cutoff_freq=100):
        # Reamostragem do sinal
        downsampled_signal = self.downsample_signal(signal, original_rate, target_rate)
        
        # Aplicação do filtro passa-baixa de Butterworth de 2ª ordem
        filtered_signal = self.butter_lowpass_filter(downsampled_signal, cutoff_freq, target_rate)
        
        return filtered_signal

    # Função para calcular Vrms
    # Função para calcular Vrms
    def calculate_vrms(self, signal, sampling_rate, offset, scale_factor):
        vrms_values = []
        time_new = []
        squared_values = []

        # Certifica-se de que há pelo menos duas amostras
        if len(signal) < 2:
            raise ValueError("O sinal deve ter pelo menos duas amostras para calcular o RMS entre intervalos.")

        # Itera sobre os pares de amostras consecutivas
        for i in range(len(signal) - 1):
            # Extrai os dois valores do sinal
            x1 = signal[i]
            x2 = signal[i + 1]
            # Aplicar fator de multiplicação e offset
            #x1 = (x1 * scale_factor) + offset
            
            #x2 = (x2 * scale_factor) + offset
            # Calcula o quadrado dos valores
            squared_x1 = x1 ** 2
            squared_x2 = x2 ** 2
            
            # Calcula o RMS para o intervalo entre x1 e x2
            area = (squared_x1 + squared_x2)*10**(-6)/ 2
            squared_values.append(area)
        # Limite para a lógica de somatório
        limite = 16666
        
        # Parte 1: até a linha 16666
        for i in range(limite):
            if i == 0:
                somatorio = 0
            else:
                somatorio += squared_values[i]
            vrms = np.sqrt(somatorio*60)/1000
            vrms_values.append(vrms)
            time_new.append(i / sampling_rate)
        
        # Parte 2: após a linha 16666
        for i in range(limite, len(squared_values)):
            somatorio -= squared_values[i - limite]
            somatorio += squared_values[i]
            vrms = np.sqrt(somatorio*60)/1000
            vrms_values.append(vrms)
            time_new.append(i / sampling_rate)

        return vrms_values, time_new
    
    def detectar_tipo_falta(self, vrms_values):
        if len(vrms_values) < 12:  # Ajuste no número de valores necessários
            print("Não há dados suficientes para detectar a falta.")
            return None

        # Converte os valores de tensão e corrente para arrays NumPy
        vrms_A = np.array(vrms_values[0])
        vrms_B = np.array(vrms_values[1])
        vrms_C = np.array(vrms_values[2])
        irms_A = np.array(vrms_values[3])
        irms_B = np.array(vrms_values[4])
        irms_C = np.array(vrms_values[5])

        # Impedâncias da linha
        Z_pos = 3.70678  # Impedância da sequência positiva
        Z_neg = 3.70678  # Impedância da sequência negativa
        Z_zero = 98.74326  # Impedância da sequência zero

        # Matriz de transformação e sua inversa
        a = cmath.exp(2j * cmath.pi / 3)

        A = np.array([
            [1, 1, 1],
            [1, a, a**2],
            [1, a**2, a]
        ])
        A_inv = (1/3) * np.array([
            [1, 1, 1],
            [1, a**2, a],
            [1, a, a**2]
        ])

        # Vetor de impedâncias de sequência
        Z_seq = np.array([Z_pos, Z_neg, Z_zero])

        # Converter impedâncias de sequência para impedâncias de fase
        Z_fase = np.dot(A, Z_seq)
        print("Impedâncias de fase calculadas:", Z_fase)

        # Calcula a impedância na falta (usando operações element-wise)
        z_falta_fase_A = vrms_A / irms_A
        z_falta_fase_B = vrms_B / irms_B
        z_falta_fase_C = vrms_C / irms_C

        # Calcular correntes de falta usando a impedância de fase
        z_fase_A = z_falta_fase_A / Z_fase[0]
        z_fase_B = z_falta_fase_B / Z_fase[1]
        z_fase_C = z_falta_fase_C / Z_fase[2]

        # Corrente de sequência zero
        I_fase_zero = vrms_A / Z_zero

        # Verificar se a falta é trifásica ou trifásica para terra
        if (np.allclose(vrms_A, vrms_B, rtol=0.1) and np.allclose(vrms_A, vrms_C, rtol=0.1)):
            tipo_falta = "Falta Trifásica"
            porcentagem_falta = 100.0 * z_fase_A
        else:
            tipo_falta = "Falta Trifásica para Terra"
            porcentagem_falta = 100.0 * z_fase_A

        # Atualiza os labels da interface
        self.label_tipo_falta.config(text=f"Tipo de Falta Identificada: {tipo_falta}")
        self.label_porcentagem_falta.config(text=f"Porcentagem da Linha de Transmissão: {porcentagem_falta:.2f}%")
        
        return tipo_falta, porcentagem_falta


    def plotar_sinais(self):
        vrms_final = []
        arquivo1 = self.parametros['arquivo1']
        arquivo2 = self.parametros['arquivo2']
        freq_corte_min = float(self.parametros['freq_corte_min'])
        freq_corte_max = float(self.parametros['freq_corte_max'])
        colunas_selecionadas = self.parametros['colunas']

        rec = comtrade.Comtrade()
        try:
            rec.load(arquivo1, arquivo2)
        except TypeError:
            print("Caminhos não foram definidos corretamente.")
            return
        signals = np.array(rec.analog)
        original_rate = 1000000
        target_rate = 1920
        tipodefalta = 0
        porcentagemfalta = 0
        cutoff_freq = freq_corte_max  # Usando freq_corte_max para o filtro

        for coluna in colunas_selecionadas:
            sinal = signals[coluna]
            channel_index = coluna  # Assumindo que a coluna selecionada é a primeira na lista
            # Ajuste os atributos corretos aqui
            scale_factor = rec.cfg.analog_channels[channel_index].a
            offset = rec.cfg.analog_channels[channel_index].b  
            
            vrms_values, time_new = self.calculate_vrms(sinal, target_rate, offset, scale_factor)
            vrms_final.append(vrms_values)
            sinal_processado = self.process_signal(vrms_values, original_rate, target_rate, cutoff_freq)
            # Criar um vetor de tempo para sinal_processado com o mesmo comprimento
            num_samples_processado = len(sinal_processado)
            time_processado = np.linspace(0, num_samples_processado / target_rate, num_samples_processado)

            
            # Criar uma figura e eixos
            plt.figure(figsize=(12, 6))
            
            # Plotar o Vrms
            plt.subplot(2, 1, 1)
            plt.plot(time_new, vrms_values)
            plt.title(f'Vrms - Canal {self.canais[coluna]}')
            plt.xlabel('Tempo (s)')
            plt.ylabel('Vrms')
            
            # Plotar o sinal processado
            plt.subplot(2, 1, 2)
            plt.plot(time_processado, sinal_processado)
            plt.title(f'Sinal Processado - Canal {self.canais[coluna]}')
            plt.xlabel('Tempo (s)')
            plt.ylabel('Sinal Processado')
            
            # Ajustar layout e mostrar a figura
            plt.tight_layout()
            plt.show()

        #print(vrms_final[1])
        #tipodefalta, porcentagemfalta = self.detectar_tipo_falta(vrms_final)




def main():
    parametros = carregar_parametros()
    root = tk.Tk()
    app = JanelaSelecaoArquivos(root, parametros)
    root.mainloop()

if __name__ == "__main__":
    main()
