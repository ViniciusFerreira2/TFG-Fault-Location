import json
import os
import numpy as np
from scipy.signal import butter, filtfilt, resample

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

def process_signal(signal, original_rate=1000000, target_rate=1920, cutoff_freq=100):
    # Reamostragem do sinal
    num_samples = int(len(signal) * target_rate / original_rate)
    downsampled_signal = resample(signal, num_samples)
    
    # Aplicação do filtro passa-baixa de Butterworth de 2ª ordem
    nyquist = 0.5 * target_rate
    normal_cutoff = cutoff_freq / nyquist
    b, a = butter(2, normal_cutoff, btype='low', analog=False)
    filtered_signal = filtfilt(b, a, downsampled_signal)
    
    return filtered_signal

def calculate_vrms(signal, sampling_rate, offset, scale_factor):
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
        area = (squared_x1 + squared_x2) * 10**(-6) / 2
        squared_values.append(area)
    limite = 16666

    for i in range(limite):
        if i == 0:
            somatorio = 0
        else:
            somatorio += squared_values[i]
        vrms = np.sqrt(somatorio * 60) / 1000
        vrms_values.append(vrms)
        time_new.append(i / sampling_rate)

    for i in range(limite, len(squared_values)):
        somatorio -= squared_values[i - limite]
        somatorio += squared_values[i]
        vrms = np.sqrt(somatorio * 60) / 1000
        vrms_values.append(vrms)
        time_new.append(i / sampling_rate)

    return vrms_values, time_new

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
