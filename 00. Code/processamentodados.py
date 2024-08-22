import json
import os
import numpy as np
from scipy.signal import butter, filtfilt, resample

# Nome do arquivo de configuração
CONFIG_FILE = 'config.json'

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

def calculate_vrms(signal, original_rate):
    # Conversão da taxa de amostragem de µs para segundos
    #t = np.arange(0, len(signal)) / original_rate  # Tempo em segundos

    # Frequência angular
    omega = 2 * np.pi * 60  # 60 Hz

    # Inicializar listas para armazenar os valores de θ, módulo e ângulo
    theta_list = []
    theta_list2 = []
    modulus_list = []
    modulus_values = []
    angle_list = []
    vrms_values = []
    time_new = []

    # Inicializar variáveis para o cálculo de VRMS
    limite = 16666
    squared_values = []
    squared_values2 = []
    t=0.000001

    # Cálculo linha por linha
    for i in range(len(signal)-1):
        # Construir a matriz X para a linha i (deve ser uma matriz 2D)
        X_t_i = np.array([[1, np.sin(omega * t), np.cos(omega * t), t]])
        # O valor Y correspondente é o valor do sinal na linha i
        Y_i = np.array([signal[i]])
        # Calcular θ para a linha i usando pseudo-inversa
       # theta_i = np.linalg.pinv(X_t_i.T @ X_t_i) @ X_t_i.T @ Y_i
        theta_i =np.linalg.pinv(X_t_i) @ Y_i 
        # Adicionar o valor de θ à lista
        theta_list.append(theta_i.flatten())  # Transformar para 1D para armazenar

        # Extraindo Theta2 e Theta3
        Theta21 = theta_i[1]
        Theta31 = theta_i[2]
        # Calculando Theta2 + Theta3 * j
        complex_number = Theta21 + Theta31 * 1j
        squared1=complex_number**2
       # if i == 16667:
         #   print(complex_number)
        # Adicionar o módulo e o módulo ao quadrado à lista
        # modulus_list.append(modulus)
        squared_values.append(squared1)
        angle_list.append(np.angle(complex_number,True))
        t+=0.000001
        
    somatorio = 0
    for i in range(limite):
        if i == 0:
            somatorio = 0
        else:
            somatorio += squared_values[i]
        vrms = np.sqrt(somatorio*60)
        vrms_values.append(vrms)  # Calcular a raiz quadrada da média
        modulus_values.append(somatorio)
        #if i == 16665:
         #   print("modulo:", vrms)
        time_new.append(i / original_rate)

    for i in range(limite, len(squared_values)):
        somatorio -= squared_values[i - limite]
        somatorio += squared_values[i]
        vrms_values.append(np.sqrt(somatorio*60))  # Calcular a raiz quadrada da média
        modulus_values.append(somatorio)
        time_new.append(i / original_rate)

    # Converter as listas para arrays numpy
    vrms_array = np.array(vrms_values)
    modulus_array = np.array(modulus_values)
    angle_array = np.array(angle_list)
    time_new_array = np.array(time_new)

    return vrms_array, time_new_array, modulus_array, angle_array


    '''
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
'''
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
