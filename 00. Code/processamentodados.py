import json
import os
import numpy as np
from scipy.signal import butter, filtfilt, resample
import cmath
import math
import matplotlib.pyplot as plt

# Nome do arquivo de configuração
CONFIG_FILE = 'config.json'

class processamento():
    def carregar_parametros():
        parametros_padroes = {
            'arquivo1': '',
            'arquivo2': '',
            'freq_amostragem': 1000.0,
            'freq_corte_max': 75.0,
            'colunas': []
        }

        if os.path.exists(CONFIG_FILE):
            with open(CONFIG_FILE, 'r') as f:
                parametros = json.load(f)
                # Atualiza os parâmetros com os valores padrão para chaves ausentes
                for chave, valor_padrao in parametros_padroes.items():
                    parametros.setdefault(chave, valor_padrao)
                return parametros
        else:
            return parametros_padroes

    def salvar_parametros(parametros):
        with open(CONFIG_FILE, 'w') as f:
            json.dump(parametros, f)

    def filter_signal(signal, original_rate, target_rate, cutoff_freq):
        # Reamostragem do sinal
        num_samples = int(len(signal) * target_rate / original_rate)
        downsampled_signal = resample(signal, num_samples)
        
        # Aplicação do filtro passa-baixa de Butterworth de 2ª ordem
        nyquist = 0.5 * target_rate
        normal_cutoff = cutoff_freq / nyquist
        b, a = butter(2, normal_cutoff, btype='low', analog=False)
        filtered_signal = filtfilt(b, a, downsampled_signal)
        num_samples_processado = len(filtered_signal)
        time_processado = np.linspace(0, num_samples_processado / target_rate, num_samples_processado)
        sample_rate = num_samples_processado

        return filtered_signal, time_processado, sample_rate
        
    def vrm_per_phase(signal, time):
        vrms_values = []
        time_new = []
        squared_values = []

        if len(signal) < 2:
            raise ValueError("O sinal deve ter pelo menos duas amostras para calcular o RMS entre intervalos.")

        for i in range(len(signal) - 1):
            x1 = signal[i]
            x2 = signal[i + 1]
            squared_x1 = x1 ** 2
            squared_x2 = x2 ** 2
            area = (squared_x1 + squared_x2) * (time[i+1]-time[i]) / 2
            squared_values.append(area)
            
        limite = (1/60)/(time[-1]/len(time))
 
        limite = math.floor(limite)
        #print(limite)

        for i in range(limite):
            if i == 0:
                somatorio = 0
            else:
                somatorio += squared_values[i]
            vrms = np.sqrt(somatorio * 60)
            vrms_values.append(vrms)
            time_new.append(time[i])

        for i in range(limite, len(squared_values)):
            somatorio -= squared_values[i - limite]
            somatorio += squared_values[i]
            vrms = np.sqrt(somatorio * 60)
            vrms_values.append(vrms)
            time_new.append(time[i])

        return vrms_values, time_new
    
    def fasor(signal, time):
        # Frequência angular
        omega = 2 * np.pi * 60  # 60 Hz
        sample_rate = math.floor((1/60)/(time[-1]/len(time)))

        yk = [] 
        xt = []
        mod_values = []
        ang_values = []
        R_values = []
        I_values = []
        ang_ant = 0 
        vrms_ant = 0

        # Pré-preenchimento para as primeiras amostras
        for i in range(sample_rate):
            yk.append(signal[i])
            xt.append([1, math.sin(omega*time[i]), math.cos(omega*time[i]), time[i]])
            mod_values.append(0)
            ang_values.append(0)

        # Processamento das amostras restantes
        for j in range(sample_rate, len(signal)):
            xt_array = np.array(xt)
            yk_array = np.array(yk)

            xt_transpoto = xt_array.T
            yk_transpoto = yk_array.T

            produt1 = np.dot(xt_transpoto, xt_array)
            if np.linalg.det(produt1) == 0:
                ang = ang_ant
                vrms = vrms_ant
            else:
                inversa1 = np.linalg.inv(produt1)
                produt2 = np.dot(inversa1, xt_transpoto)
                deltas = np.dot(produt2, yk_transpoto)

                value = complex(deltas[1], deltas[2])

                vrms = abs(value)
                ang = math.degrees(cmath.phase(value))

            vrms_ant = vrms
            ang_ant = ang

            mod_values.append(vrms)
            ang_values.append(ang)
            
            # Atualiza as listas de yk e xt
            del yk[0]  # Remove a primeira amostra
            del xt[0]  # Remove a primeira amostra

            yk.append(signal[j])  # Adiciona a nova amostra
            xt.append([1, math.sin(omega*time[j]), math.cos(omega*time[j]), time[j]])

        # Conversão para numpy array
        signal_modulo = np.array(mod_values)
        signal_ang = np.array(ang_values)
        # Parte real e imaginária
        real_values = signal_modulo * np.cos(signal_ang)
        imag_values = signal_modulo * np.sin(signal_ang)

        # Sinal complexo (opcional, se quiser juntar as partes real e imaginária em um array complexo)
        sinal_complexo = real_values + 1j * imag_values

        #signal_ang_mtz = signal_ang.reshape(960,1)
        return signal_modulo, signal_ang, sinal_complexo
        
    

    # def plot_XR(sinal_complexo):
    #     # Separar a parte real e imaginária do sinal complexo
    #     parte_real = np.real(sinal_complexo)
    #     parte_imaginaria = np.imag(sinal_complexo)

    #     # Plotar o gráfico
    #     plt.figure(figsize=(8, 6))
    #     plt.plot(parte_real, parte_imaginaria, 'bo-', label='Sinal Complexo')
        
    #     # Títulos e labels
    #     plt.title('Gráfico da Parte Real vs Parte Imaginária do Sinal')
    #     plt.xlabel('Parte Real')
    #     plt.ylabel('Parte Imaginária')
        
    #     # Adicionar grade e legenda
    #     plt.grid(True)
    #     plt.legend()

    #     # Mostrar o gráfico
    #     plt.show()

    def detectar_tipo_falta(vrms_values):
        if len(vrms_values) < 12:
            print("Não há dados suficientes para detectar a falta.")
            return None

        vrms_A = np.array(vrms_values[0])
        vrms_B = np.array(vrms_values[1])
        vrms_C = np.array(vrms_values[2])
        irms_A = np.array(vrms_values[3])
        irms_B = np.array(vrms_values[4])
        irms_C = np.array(vrms_values[5])

        Z_pos = 3.70678
        Z_neg = 3.70678
        Z_zero = 98.74326

        a = cmath.exp(2j * cmath.pi / 3)

        A = np.array([
            [1, 1, 1],
            [1, a, a**2],
            [1, a**2, a]
        ])
        A_inv = (1 / 3) * np.array([
            [1, 1, 1],
            [1, a**2, a],
            [1, a, a**2]
        ])

        Z_seq = np.array([Z_pos, Z_neg, Z_zero])
        Z_fase = np.dot(A, Z_seq)

        z_falta_fase_A = vrms_A / irms_A
        z_falta_fase_B = vrms_B / irms_B
        z_falta_fase_C = vrms_C / irms_C

        z_fase_A = z_falta_fase_A / Z_fase[0]
        z_fase_B = z_falta_fase_B / Z_fase[1]
        z_fase_C = z_falta_fase_C / Z_fase[2]

        I_fase_zero = vrms_A / Z_zero

        if (np.allclose(vrms_A, vrms_B, rtol=0.1) and np.allclose(vrms_A, vrms_C, rtol=0.1)):
            tipo_falta = "Falta Trifásica"
            porcentagem_falta = 100.0 * z_fase_A
        else:
            tipo_falta = "Falta Trifásica para Terra"
            porcentagem_falta = 100.0 * z_fase_A

        return tipo_falta, porcentagem_falta
