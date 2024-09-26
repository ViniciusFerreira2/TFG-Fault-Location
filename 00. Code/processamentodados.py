import json
import os
import numpy as np
from scipy.signal import butter, filtfilt, resample
import cmath
import math
import comtrade
from plot import*
from datetime import datetime
import matplotlib.pyplot as plt

# Nome do arquivo de configuração
CONFIG_FILE = 'config.json'

class processamento():
    def process(self, parametros, canais):
        processamento.init_process(self, parametros, canais)
        print("003_____ REALIZANDO FILTRO")
        processamento.general_filter(self, parametros, canais)
        print("004_____ REALIZANDO RMS")
        processamento.general_rms(self, parametros, canais)
        print("004_____ REALIZANDO FASOR")
        signal = processamento.general_fasor(self, parametros, canais)
        print(f"XXX_____REALIZANDO COMPONENTES SIMETRICAS")
        processamento.general_symmetrical_components(self, parametros, signal)
        #processamento.symmetrical_componentes(parametros, signal)
        print("005_____ REALIZANDO IMPENDANCIA")
        #processamento.general_impendace(self, parametros, canais)
       
    def init_process(self, parametros, canais):
        #global modulo, angulo
        vrms_final = []
        
        arquivo1 = parametros['arquivo1']
        arquivo2 = parametros['arquivo2']
        freq_amostragem = float(parametros['freq_amostragem'])
        freq_corte_max = float(parametros['freq_corte_max'])
        colunas_selecionadas = parametros['colunas']
        print("COLUNAS", colunas_selecionadas)
        

        self.rec = comtrade.Comtrade()
        try:
            self.rec.load(arquivo1, arquivo2)
            print("001_____ARQUIVOS IMPORTADOS COM SUCESSO")
        except TypeError:
            print("Caminhos não foram definidos corretamente.")
            return

        signals = np.array(self.rec.analog)
        original_rate = 1000000
        self.original_rate = original_rate
        target_rate = freq_amostragem
        cutoff_freq = freq_corte_max

        now = datetime.now()
        self.timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")
        current_directory = os.path.dirname(os.path.abspath(__file__))
        self.pasta_nome = os.path.join(current_directory, f"Graficos_{self.timestamp}")
        os.makedirs(self.pasta_nome, exist_ok=True)
        print(f"002_____CRIAÇÃO DA PASTA: {self.pasta_nome}")

    def general_filter(self, parametros, canais):
        freq_amostragem = float(parametros['freq_amostragem'])
        freq_corte_max = float(parametros['freq_corte_max'])
        colunas_selecionadas = parametros['colunas']
        signals = np.array(self.rec.analog)
        original_rate = self.original_rate
        target_rate = freq_amostragem
        cutoff_freq = freq_corte_max
        plotar_filtro = parametros.get('plotar_filtro', False)
        
        self.sinal_filtrado = []
        self.time_filtrado = []
        self.sample_rate = []
        index = 0


        for coluna in colunas_selecionadas:
            
            coluna_index = self.rec.analog_channel_ids.index(coluna)
            sinal = signals[coluna_index]
            scale_factor = self.rec.cfg.analog_channels[coluna_index].a
            offset = self.rec.cfg.analog_channels[coluna_index].b 
            time = np.linspace(0, len(sinal) / original_rate, len(sinal))

            sinal_filtrado, time_filtrado, sample_rate = processamento.filter_signal(sinal, original_rate, target_rate, cutoff_freq)
            self.sinal_filtrado.append(sinal_filtrado)
            self.time_filtrado.append(time_filtrado)
            self.sample_rate.append(sample_rate)
            # Plotar o sinal filtrado somente se o checkbox correspondente estiver marcado
            if plotar_filtro:
                plot_filter(sinal, canais, coluna_index, self.time_filtrado[index], self.sinal_filtrado[index], original_rate, self.timestamp, self.pasta_nome)
            index = index + 1

    def general_rms(self, parametros, canais):
        freq_amostragem = float(parametros['freq_amostragem'])
        freq_corte_max = float(parametros['freq_corte_max'])
        colunas_selecionadas = parametros['colunas']
        signals = np.array(self.rec.analog)
        original_rate = self.original_rate
        plotar_rms = parametros.get('plotar_rms', False)

        self.sinal_vrms = []
        self.time_vrms = []
        index = 0

        for coluna in colunas_selecionadas:
            coluna_index = self.rec.analog_channel_ids.index(coluna)
            sinal = signals[coluna_index]
            scale_factor = self.rec.cfg.analog_channels[coluna_index].a
            offset = self.rec.cfg.analog_channels[coluna_index].b 
            time = np.linspace(0, len(sinal) / original_rate, len(sinal))

            # Plotar o RMS somente se o checkbox correspondente estiver marcado
            sinal_vrms, time_vrms = processamento.vrm_per_phase(self.sinal_filtrado[index], self.time_filtrado[index])

            self.sinal_vrms.append(sinal_vrms)
            self.time_vrms.append(time_vrms)
            if plotar_rms:
                plot_vrms(self.time_filtrado[index], self.sinal_filtrado[index], canais, coluna_index, self.time_vrms[index], self.sinal_vrms[index], original_rate, self.timestamp, self.pasta_nome)
            index = index + 1

    def general_fasor(self, parametros, canais):
        self.modulo = []
        self.angulo = []
        self.complexo = []
        freq_amostragem = float(parametros['freq_amostragem'])
        freq_corte_max = float(parametros['freq_corte_max'])
        colunas_selecionadas = parametros['colunas']
        signals = np.array(self.rec.analog)
        plotar_fasores = parametros.get('plotar_fasores', False)
        index = 0

        for coluna in colunas_selecionadas:
            coluna_index = self.rec.analog_channel_ids.index(coluna)
            sinal = signals[coluna_index]

            # Plotar os fasores somente se o checkbox correspondente estiver marcado
            signal_modulo, signal_ang, signal_complexo = processamento.fasor(self.sinal_filtrado[index], self.time_filtrado[index]) #RADIANOS
            self.modulo.append(signal_modulo)
            self.angulo.append(signal_ang)
            self.complexo.append(signal_complexo)
            index = index + 1

        print("004.1_____ CORRIGIR ANGULO")  
        self.angulo = processamento.ang_correction(self.angulo)
        if plotar_fasores:
            plot_fasor(self.modulo, self.angulo, self.pasta_nome)

        return self.complexo

    def general_symmetrical_components(self, parametros, signal):

        seq_V, seq_I = processamento.symmetrical_componentes(signal)
        processamento.impedance_symmetrical_components(seq_V, seq_I)
        
        # if plotar_fasores and self.modulo and self.angulo:
        #     plot_fasor([self.angulo[0], self.angulo[2], self.angulo[4]], [self.angulo[1], self.angulo[3], self.angulo[5]], self.pasta_nome) #[angulo[1], angulo[3], angulo[5]], pasta_nome)
        #     #plot_fasor(, angulo, pasta_nome)
        #     #plotar_fasores(modulo, angulo, 18)
            
    def carregar_parametros():

        parametros_padroes = {
            'arquivo1': '',
            'arquivo2': '',
            'freq_amostragem': 1000.0,
            'freq_corte_max': 75.0,
            'colunas': [],
            'dadoslinha': {
                'R1': 0.0,
                'X1': 0.0,
                'R0': 0.0,
                'X0': 0.0,
                'L': 0.0
            }
        }

        if os.path.exists(CONFIG_FILE):
            with open(CONFIG_FILE, 'r') as f:
                parametros = json.load(f)
                # Atualiza os parâmetros com os valores padrão para chaves ausentes
                for chave, valor_padrao in parametros_padroes.items():
                    if chave == 'dadoslinha':
                        # Para 'dadoslinha', atualiza as subchaves
                        parametros[chave] = {**parametros_padroes[chave], **parametros.get(chave, {})}
                    else:
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
            ang_values.append(ang) #GRAUS
            
            # Atualiza as listas de yk e xt
            del yk[0]  # Remove a primeira amostra
            del xt[0]  # Remove a primeira amostra

            yk.append(signal[j])  # Adiciona a nova amostra
            xt.append([1, math.sin(omega*time[j]), math.cos(omega*time[j]), time[j]])

        # Conversão para numpy array
        signal_modulo = np.array(mod_values)
        signal_ang = np.array(ang_values) #GRAUS
        # Parte real e imaginária
        real_values = signal_modulo * np.cos((signal_ang*math.pi/180)) #RADIANOS
        imag_values = signal_modulo * np.sin((signal_ang*math.pi/180))

        # Sinal complexo (opcional, se quiser juntar as partes real e imaginária em um array complexo)
        sinal_complexo = real_values + 1j * imag_values

        #signal_ang_mtz = signal_ang.reshape(960,1)
        return signal_modulo, signal_ang, sinal_complexo #SINAL COMPLEXO TA CERTO
    
    def ajustar_angulo(angulo):
        limiteang = 180
        if angulo > limiteang:
            angulo = -180 + (angulo - 180)
        elif angulo < -limiteang:
            angulo = 180 + (angulo + 180)
        return angulo

    def ang_correction(angulo):
        for i in range(0, len(angulo[0])):
            #angulo[0][i] = angulo[0][i] - angulo0[i]  
            # Ajuste das demais fases com relação à Fase A
            angulo[1][i] = processamento.ajustar_angulo(angulo[1][i] - angulo[0][i])
            angulo[2][i] = processamento.ajustar_angulo(angulo[2][i] - angulo[0][i]) 
            angulo[3][i] = processamento.ajustar_angulo(angulo[3][i] - angulo[0][i])
            angulo[4][i] = processamento.ajustar_angulo(angulo[4][i] - angulo[0][i])
            angulo[5][i] = processamento.ajustar_angulo(angulo[5][i] - angulo[0][i])
            angulo[0][i] = processamento.ajustar_angulo(angulo[0][i] - angulo[0][i])

        return angulo

    def calculo_impedancia(modulo, angulo):

        Zmod = [[0] * len(angulo[0]-120) for _ in range(int(round(len(angulo)/2)))]  # Inicializa uma lista 3x960
        Zang = [[0] * len(angulo[0]-120) for _ in range(int(round(len(angulo)/2)))]  # Inicializa uma lista 3x960
        Complexo = [[0] * len(angulo[0]-120) for _ in range(int(round(len(angulo)/2)))]  # Inicializa uma lista 3x960


        for i in range(120, len(angulo[0])):
            #Calculo dos modulos e angulos das impedancias
            Zmod[0][i] = modulo[0][i] / modulo[1][i]
            Zmod[1][i] = modulo[2][i] / modulo[3][i]
            Zmod[2][i] = modulo[4][i] / modulo[5][i]

            Zang[0][i] = processamento.ajustar_angulo(angulo[0][i] - angulo[1][i])
            Zang[1][i] = processamento.ajustar_angulo(angulo[2][i] - angulo[3][i])
            Zang[2][i] = processamento.ajustar_angulo(angulo[4][i] - angulo[5][i])
        
            Complexo[0][i] = (Zmod[0][i] * np.cos(Zang[0][i]*math.pi/180)) + 1j*(Zmod[0][i] * np.sin(Zang[0][i]*math.pi/180))
            Complexo[1][i] = (Zmod[1][i] * np.cos(Zang[1][i]*math.pi/180)) + 1j*(Zmod[1][i] * np.sin(Zang[1][i]*math.pi/180))
            Complexo[2][i] = (Zmod[2][i] * np.cos(Zang[2][i]*math.pi/180)) + 1j*(Zmod[1][i] * np.sin(Zang[2][i]*math.pi/180))


        return Complexo

    def line_data(parameters):
        # Extrai os valores de dadoslinha
        dadoslinha = parameters.get('dadoslinha', {})

        # Extrai e converte os valores para float, caso não sejam numéricos
        R1 = float(dadoslinha.get('R1', 0.0))
        X1 = float(dadoslinha.get('X1', 0.0))
        R0 = float(dadoslinha.get('R0', 0.0))
        X0 = float(dadoslinha.get('X0', 0.0))
        L = float(dadoslinha.get('L', 0.0))

        # Realiza os cálculos com os valores numéricos
        R1 = R1 * L
        X1 = X1 * L
        R0 = R0 * L
        X0 = X0 * L

    def symmetrical_componentes(complexo):

        mod_values = []
        ang_values = []
        # Inicialização de uma lista para modulo e angulo da impedancia
        seq_voltage = [[0] * len(complexo[0]) for _ in range(int(round(len(complexo))))]
        seq_current = [[0] * len(complexo[0]) for _ in range(int(round(len(complexo))))]   


        # Definição das matrizes de síntese e análise
        a = cmath.exp(2j * cmath.pi / 3)

        #Conversão da impedancia de entrada para forma polar
        for i in range(0, len(complexo[0])):

            A = np.array([
                [1, 1, 1],
                [1, a, a**2],
                [1, a**2, a]
            ])

            B_V = np.array([complexo[0][i], complexo[2][i], complexo[4][i]])

            B_I = np.array([complexo[1][i], complexo[3][i], complexo[5][i]])

            seq_V = (1/3)*(np.dot(A, B_V))
            seq_V = seq_V.T

            seq_I = (1/3)*(np.dot(A, B_I))
            seq_I = seq_I.T

            seq_voltage[0][i] = seq_V[0]
            seq_voltage[1][i] = seq_V[1]
            seq_voltage[2][i] = seq_V[2]

            seq_current[0][i] = seq_I[0]
            seq_current[1][i] = seq_I[1]
            seq_current[2][i] = seq_I[2]
        return seq_voltage, seq_current

    def impedance_symmetrical_components(seq_V, seq_I):
        
        seq_mod = [[0] * len(seq_V[0]) for _ in range(int(round(len(seq_V))))]
        seq_ang = [[0] * len(seq_V[0]) for _ in range(int(round(len(seq_V))))]
        Z_seq_mod = [[0] * len(seq_V[0]) for _ in range(int(round(len(seq_V))))]
        Z_seq_ang = [[0] * len(seq_V[0]) for _ in range(int(round(len(seq_V))))]
        Z_seq_real = [[0] * len(seq_V[0]) for _ in range(int(round(len(seq_V))))]
        Z_seq_imag = [[0] * len(seq_V[0]) for _ in range(int(round(len(seq_V))))]  

        for i in range(0, len(seq_V[0])):    
            seq_mod[0][i], seq_ang[0][i] = cmath.polar(seq_V[0][i])
            seq_mod[1][i], seq_ang[1][i] = cmath.polar(seq_V[1][i])
            seq_mod[2][i], seq_ang[2][i] = cmath.polar(seq_V[2][i])
            seq_mod[3][i], seq_ang[3][i] = cmath.polar(seq_I[0][i])
            seq_mod[4][i], seq_ang[4][i] = cmath.polar(seq_I[1][i])
            seq_mod[5][i], seq_ang[5][i] = cmath.polar(seq_I[2][i])   
            if i > 80:  # antes de 80 os valores estão como 0, estava dando erro por divisão por zero
                # Verificar se os divisores são iguais a zero e atribuir 1 se forem
                Z_seq_mod[0][i] = seq_mod[0][i] / (seq_mod[3][i])
                Z_seq_mod[1][i] = seq_mod[1][i] / (seq_mod[4][i])
                Z_seq_mod[2][i] = seq_mod[2][i] / (seq_mod[5][i])

                Z_seq_ang[0][i] = seq_ang[0][i] - seq_ang[3][i]
                Z_seq_ang[1][i] = seq_ang[1][i] - seq_ang[4][i]
                Z_seq_ang[2][i] = seq_ang[2][i] - seq_ang[5][i]
            else:
                Z_seq_mod[0][i] = 0
                Z_seq_mod[1][i] = 0
                Z_seq_mod[2][i] = 0

                Z_seq_ang[0][i] = 0
                Z_seq_ang[1][i] = 0
                Z_seq_ang[2][i] = 0

            # Parte real e imaginária
            Z_seq_real[0][i] = Z_seq_mod[0][i] * np.cos(Z_seq_ang[0][i])
            Z_seq_real[1][i] = Z_seq_mod[1][i] * np.cos(Z_seq_ang[1][i])
            Z_seq_real[2][i] = Z_seq_mod[2][i] * np.cos(Z_seq_ang[2][i])

            Z_seq_imag[0][i] = Z_seq_mod[0][i] * np.sin(Z_seq_ang[0][i])
            Z_seq_imag[1][i] = Z_seq_mod[1][i] * np.sin(Z_seq_ang[1][i])
            Z_seq_imag[2][i] = Z_seq_mod[2][i] * np.sin(Z_seq_ang[2][i])
        # Z_seq_r = np.array(Z_seq_real)
        # Z_seq_i = np.array(Z_seq_imag)
        # sinal_complex = Z_seq_r + 1j*Z_seq_i
        # print("SINAL100", sinal_complex[0][100])
       
        # plot_XR(sinal_complex,(0.0166*223.3),(0.2624*223.3),(0.4422*223.3),(223.3*1.3189))
        plot_Z_seq(Z_seq_mod, Z_seq_ang)

    
    def detectar_tipo_falta(rms_values):
        if len(rms_values) < 12:
            print("Não há dados suficientes para detectar a falta.")
            return None

        vrms_A = np.array(rms_values[0])
        vrms_B = np.array(rms_values[1])
        vrms_C = np.array(rms_values[2])
        irms_A = np.array(rms_values[3])
        irms_B = np.array(rms_values[4])
        irms_C = np.array(rms_values[5])

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
