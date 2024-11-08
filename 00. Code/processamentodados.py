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
import heapq
from sympy import Eq, solve, symbols

# Obtém o diretório onde o script está localizado
script_dir = os.path.dirname(os.path.abspath(__file__))

# Caminho completo para salvar o arquivo config.json na mesma pasta do script
CONFIG_FILE = os.path.join(script_dir, 'config.json')

class processamento():
    def process(self, parametros, canais):
        processamento.init_process(self, parametros, canais)
        print("002_____ REALIZANDO FILTRO")
        processamento.general_filter(self, parametros, canais)
        print("003_____ REALIZANDO RMS")
        rms=processamento.general_rms(self, parametros, canais)
        print("004_____ REALIZANDO FASOR")
        signal, modulorms = processamento.general_fasor(self, parametros, canais)
        print("005_____REALIZANDO COMPONENTES SIMETRICAS")
        seq_I,seq_I_2, Z_seq_r_1, Z_seq_i_1, Z_seq_r_2, Z_seq_i_2 = processamento.general_symmetrical_components(self, signal, modulorms, parametros)
        print("006_____REALIZANDO LOCALIZAÇÃO DA FALTA")
        processamento.general_fault_location(self, parametros, signal, rms, seq_I, seq_I_2, Z_seq_r_1, Z_seq_i_1, Z_seq_r_2, Z_seq_i_2)
        print("007_____ REALIZANDO IMPENDANCIA")
        #processamento.general_impendace(self, parametros, canais)
       
    def init_process(self, parametros, canais):
        #global modulo, angulo
        vrms_final = []
        
        arquivo1 = parametros['arquivo1']
        arquivo2 = parametros['arquivo2']
        freq_amostragem = float(parametros['freq_amostragem'])
        freq_corte_max = float(parametros['freq_corte_max'])
        colunas_selecionadas = parametros['colunas']  

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
        return self.sinal_vrms

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
        #self.angulo = processamento.ang_correction(self.angulo)
        if plotar_fasores:
            plot_fasor(self.modulo, self.angulo, self.pasta_nome)

        return self.complexo, self.modulo

    def general_symmetrical_components(self, signal, rms, parametros):

        seq_V, seq_I, seq_V_2, seq_I_2 = processamento.symmetrical_componentes(signal, rms)
        print("TENSAOVALORES", seq_I[1])
        Z_seq_r_1, Z_seq_i_1 = processamento.impedance_symmetrical_components(seq_V, seq_I, parametros)
        Z_seq_r_2, Z_seq_i_2 = processamento.impedance_symmetrical_components(seq_V_2, seq_I_2, parametros)
        #print("Z_SEQ_1", Z_seq_r_1[:300],Z_seq_i_1[:300])
        #print("Z_SEQ_2", Z_seq_r_2[:300],Z_seq_i_2[:300])
        return seq_I, seq_I_2, Z_seq_r_1, Z_seq_i_1, Z_seq_r_2, Z_seq_i_2

    def general_fault_location(self, parametros, signal, rms, seq_I, seq_I_2, Z_seq_r_1, Z_seq_i_1, Z_seq_r_2, Z_seq_i_2):
        # Obtém o valor da fase em falta do dicionário de parâmetros
        fases_afetadas = parametros.get('fase_em_falta', '')
        fault_time = processamento.fault_detection(rms)
        #processamento.SAHA_1_Terminal_Fault_Location(fases_afetadas)
        processamento.Takagi_single_terminal(parametros, fases_afetadas, signal, fault_time, seq_I)
        processamento.Tziouvaras_two_terminals(parametros, fases_afetadas, signal, fault_time, seq_I, seq_I_2, Z_seq_r_1, Z_seq_i_1, Z_seq_r_2, Z_seq_i_2)

    def salvar_parametros(parametros):
        # Salvar os parâmetros existentes no JSON
        with open('parametros.json', 'w') as f:
            json.dump(parametros, f, indent=4)

    def carregar_parametros():
        try:
            with open('parametros.json', 'r') as f:
                parametros = json.load(f)
                return parametros
        except FileNotFoundError:
            print("Arquivo de parâmetros não encontrado. Criando novo.")
            return {}
    
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
        #     mod_values.append(0)
        #     ang_values.append(0)

        # Processamento das amostras restantes
        for j in range(sample_rate, len(signal)):
        #for j in range(len(signal)):
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

    def symmetrical_componentes(complexo, rms):

        # Inicialização de uma lista para modulo e angulo da impedancia
        seq_voltage = [[0] * len(complexo[0]) for _ in range(int(round(len(complexo))))]
        seq_current = [[0] * len(complexo[0]) for _ in range(int(round(len(complexo))))]
        seq_voltage_2 = [[0] * len(complexo[0]) for _ in range(int(round(len(complexo))))]
        seq_current_2 = [[0] * len(complexo[0]) for _ in range(int(round(len(complexo))))]      
        
        #Definição da porcentagem do filtro
        filter_percentage = 0.10

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

            seq_voltage[0][i] = seq_V[0] if abs(seq_V[0]) > (rms[0][i]*filter_percentage) else 0
            seq_voltage[1][i] = seq_V[1] if abs(seq_V[1]) > (rms[2][i]*filter_percentage) else 0   
            seq_voltage[2][i] = seq_V[2] if abs(seq_V[2]) > (rms[4][i]*filter_percentage) else 0

            seq_current[0][i] = seq_I[0] if abs(seq_I[0]) > (rms[1][i]*filter_percentage) else 0
            seq_current[1][i] = seq_I[1] if abs(seq_I[1]) > (rms[3][i]*filter_percentage) else 0
            seq_current[2][i] = seq_I[2] if abs(seq_I[2]) > (rms[5][i]*filter_percentage) else 0

            if len(complexo) >=7:
                B_V_2 = np.array([complexo[6][i], complexo[8][i], complexo[10][i]])

                B_I_2 = np.array([complexo[7][i], complexo[9][i], complexo[11][i]])

                seq_V_2 = (1/3)*(np.dot(A, B_V_2))
                seq_V_2 = seq_V_2.T

                seq_I_2 = (1/3)*(np.dot(A, B_I_2))
                seq_I_2 = seq_I_2.T

                seq_voltage_2[0][i] = seq_V_2[0] if abs(seq_V_2[0]) > (rms[6][i]*filter_percentage) else 0
                seq_voltage_2[1][i] = seq_V_2[1] if abs(seq_V_2[1]) > (rms[8][i]*filter_percentage) else 0   
                seq_voltage_2[2][i] = seq_V_2[2] if abs(seq_V_2[2]) > (rms[10][i]*filter_percentage) else 0

                seq_current_2[0][i] = seq_I_2[0] if abs(seq_I_2[0]) > (rms[7][i]*filter_percentage) else 0
                seq_current_2[1][i] = seq_I_2[1] if abs(seq_I_2[1]) > (rms[9][i]*filter_percentage) else 0
                seq_current_2[2][i] = seq_I_2[2] if abs(seq_I_2[2]) > (rms[11][i]*filter_percentage) else 0

        return seq_voltage, seq_current, seq_voltage_2, seq_current_2

    def impedance_symmetrical_components(seq_V, seq_I, parametros):

        plotar_XR = parametros.get('plotar_XR', False)
        plotar_Z_seq = parametros.get('plotar_Z_seq', False)
        
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
                Z_seq_mod[0][i] = (seq_mod[0][i] / seq_mod[3][i]) if ((abs(seq_V[0][i]) != 0) and (abs(seq_I[0][i]) != 0)) else 0
                Z_seq_mod[1][i] = (seq_mod[1][i] / seq_mod[4][i]) if ((abs(seq_V[1][i]) != 0) and (abs(seq_I[1][i]) != 0)) else 0
                Z_seq_mod[2][i] = (seq_mod[2][i] / seq_mod[5][i]) if ((abs(seq_V[2][i]) != 0) and (abs(seq_I[2][i]) != 0)) else 0

                Z_seq_ang[0][i] = (seq_ang[0][i] - seq_ang[3][i]) if ((abs(seq_V[0][i]) != 0) and (abs(seq_I[0][i]) != 0)) else 0
                Z_seq_ang[1][i] = (seq_ang[1][i] - seq_ang[4][i]) if ((abs(seq_V[1][i]) != 0) and (abs(seq_I[1][i]) != 0)) else 0
                Z_seq_ang[2][i] = (seq_ang[2][i] - seq_ang[5][i]) if ((abs(seq_V[2][i]) != 0) and (abs(seq_I[2][i]) != 0)) else 0
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

        Z_seq_r = np.array(Z_seq_real)
        Z_seq_i = np.array(Z_seq_imag)
        sinal_complex = Z_seq_r + 1j*Z_seq_i

        #print(sinal_complex[2])
        #print("SINAL100", sinal_complex[0][100])
        if plotar_XR:
            plot_XR(sinal_complex, parametros)
        if plotar_Z_seq:
            plot_Z_seq(parametros, Z_seq_mod, Z_seq_ang)
        return Z_seq_r, Z_seq_i

    def fault_detection(rms):

        fault_filter = 1.15
        time = (80/4800)
        for i in range(len(rms[1])):
            time+=(1/4800)
            if (abs(rms[1][i+80])*fault_filter) < (abs(rms[1][i+88])):
                fault_time = i
                break
        print("TEMPO DA FALTA", time)
        return fault_time              

    def fault_type_detection(rms):
       
        current_A = heapq.nlargest(10, rms[1])
        current_B = heapq.nlargest(10, rms[3])
        current_C = heapq.nlargest(10, rms[5])

        current_A = sum(current_A) / len(current_A)
        current_B = sum(current_B) / len(current_B)
        current_C = sum(current_C) / len(current_C)
        current_N = current_A+current_B+current_C

        # Determinando o valor máximo das correntes
        max_current = max(current_A, current_B, current_C)

        # Porcentagens em relação à corrente máxima
        percentage_A = current_A / max_current
        percentage_B = current_B / max_current
        percentage_C = current_C / max_current

        # Definindo os limites
        limite_fase = 0.5  # Limite para identificar fases afetadas (50%)
        limite_terra_corrente = 0.9  # Limite para identificar aterramento (10%)

        # Verificando quais fases estão afetadas
        fases_afetadas = []
        if percentage_A > limite_fase:
            fases_afetadas.append('A')
        if percentage_B > limite_fase:
            fases_afetadas.append('B')
        if percentage_C > limite_fase:
            fases_afetadas.append('C')

        # Verificando se há falta para terra (corrente muito baixa em uma das fases)
        terra = True
        # Se qualquer fase tiver uma corrente extremamente baixa, consideramos que há falta para terra
        if current_A > (limite_terra_corrente*current_B) and current_A > (limite_terra_corrente*current_C):
            if current_N > (limite_terra_corrente*current_A):
                terra = False
        elif current_B > (limite_terra_corrente*current_A) and current_B > (limite_terra_corrente*current_C):
            if current_N > (limite_terra_corrente*current_B):
                terra = False
        elif current_C > (limite_terra_corrente*current_A) and current_C > (limite_terra_corrente*current_B):
            if current_N > (limite_terra_corrente*current_C):
                terra = False

        # Classificando a falta
        if len(fases_afetadas) == 3:
            if terra:
                tipo_falta = "FALTA TRIFÁSICA TERRA"
            else:
                tipo_falta = "FALTA TRIFÁSICA"
        elif len(fases_afetadas) == 2:
            if terra:
                tipo_falta = "FALTA BIFÁSICA TERRA"
            else:
                tipo_falta = "FALTA BIFÁSICA"
        elif len(fases_afetadas) == 1:
            tipo_falta = "FALTA MONOFÁSICA"
        else:
            tipo_falta = "NÃO HÁ FALTA"

        # Resultado
        print(tipo_falta)
        if fases_afetadas:
            print("Fases afetadas:", ', '.join(fases_afetadas))
        if terra:
            print("Falta para terra detectada.")

    def define_KF(fases_afetadas):
         
        if fases_afetadas == "A-T":
            # Matriz KF_AT
            KF = np.array([
                [1, 0, 0],
                [0, 0, 0],
                [0, 0, 0]
            ])
        elif fases_afetadas == "B-T":
            # Matriz KF_BT
            KF = np.array([
                [0, 0, 0],
                [0, 1, 0],
                [0, 0, 0]
            ])
        elif fases_afetadas == "C-T":
            # Matriz KF_CT
            KF = np.array([
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 1]
            ])
        elif fases_afetadas == "AB":
            # Matriz KF_AB
            KF = np.array([
                [1, -1, 0],
                [-1, 1, 0],
                [0, 0, 0]
            ])
        elif fases_afetadas == "BC":
            # Matriz KF_BC
            KF = np.array([
                [0, 0, 0],
                [0, 1, -1],
                [0, -1, 1]
            ])
        elif fases_afetadas == "CA":
            # Matriz KF_CA
            KF = np.array([
                [1, 0, -1],
                [0, 0, 0],
                [-1, 0, 1]
            ])
        elif fases_afetadas == "AB-T":
            # Matriz KF_ABT
            KF = np.array([
                [2, -1, 0],
                [-1, 2, 0],
                [0, 0, 0]
            ])

        elif fases_afetadas == "BC-T":
            # Matriz KF_BCT
            KF = np.array([
                [0, 0, 0],
                [0, 2, -1],
                [0, -1, 2]
            ])

        elif fases_afetadas == "CA-T":
            # Matriz KF_CAT
            KF = np.array([
                [2, 0, -1],
                [0, 0, 0],
                [-1, 0, 2]
            ])

        elif fases_afetadas == "ABC":
            # Matriz KF_ABC
            KF = np.array([
                [2, -1, -1],
                [-1, 2, -1],
                [-1, -1, 2]
            ])

        elif fases_afetadas == "ABC-T":
            # Matriz KF_ABCT
            KF = np.array([
                [3, -1, -1],
                [-1, 3, -1],
                [-1, -1, 3]
            ])

        else:
            raise ValueError(f"Fases afetadas '{fases_afetadas}' não reconhecido.")

        # Retorna a matriz KF correspondente
        return KF

    def SAHA_1_Terminal_Fault_Location(fases_afetadas, I_terminal):

        KF = processamento.define_KF(fases_afetadas)
        Z_L = 1
        I_G = I_terminal
        I_G_PRE = I_terminal
        Z_H = 1
        Z_G = 1
        V = 1
        V_G = 1
        IMP_GHL = Z_G + Z_H + Z_L
        CUR_IGPRE = I_G + I_G_PRE
       
        A = np.dot(Z_L, KF, Z_L, I_G)
        B = (np.dot(Z_L, KF, V,)+np.dot(Z_H, KF, Z_L, I_G))
        C = np.dot((Z_L+Z_H), KF, V_G)
        D = np.dot(IMP_GHL, CUR_IGPRE)

    def Takagi_single_terminal(parametros, fases_afetadas, RMS, fault_time, seq_I):
        local_falta = []
        plotar_LF = parametros.get('plotar_LF', False)
        freq_amostragem = float(parametros['freq_amostragem'])
        Z1_L = (parametros['dadoslinha']['R1']+1j* parametros['dadoslinha']['X1'])*parametros['dadoslinha']['L']
        Z0_L = (parametros['dadoslinha']['R0']+1j* parametros['dadoslinha']['X0'])*parametros['dadoslinha']['L']
        K= (Z0_L/Z1_L)-1

        cycle_ms = round(0.010/(1/freq_amostragem))
        prefault_ms = round(0.025/(1/freq_amostragem))
        fix_ms = round(0.025/(1/freq_amostragem))
        #print("inicio ig_preA", (cycle_ms)*0.000208)
        #print("fim ig_preA", (fault_time-45)*0.000208)
        # print("tamanho RMS 1:",len(RMS[1]))
        #print("faulttime-tenms", fault_time - ten_ms)
        if 'A' in fases_afetadas:
            Vg_A = RMS[0][fault_time+fix_ms:(len(RMS[1]) - fix_ms)]
            Ig_pre_A = RMS[1][cycle_ms:(fault_time-cycle_ms)]
            Ig_pre_avrg_A = sum(Ig_pre_A)/len(Ig_pre_A)
            Ig_A = RMS[1][fault_time+fix_ms:(len(RMS[1]) - fix_ms)]
            for i in range(len(Ig_A)):
                Ig_A[i] = Ig_A[i]-Ig_pre_avrg_A
        
        if 'B' in fases_afetadas:
            Vg_B = RMS[2][fault_time+fix_ms:]
            Ig_pre_B = RMS[3][130:((len(RMS[1]))-(fault_time-cycle_ms))]
            Ig_pre_avrg_B = sum(Ig_pre_B)/len(Ig_pre_B)
            Ig_B = RMS[3][fault_time+fix_ms:(len(RMS[1]) - fix_ms)]
            for i in range(len(Ig_B)):
                Ig_B[i] = Ig_B[i]-Ig_pre_avrg_B
        
        if 'C' in fases_afetadas:
            Vg_C = RMS[4][fault_time+fix_ms:]
            Ig_pre_C = RMS[5][130:((len(RMS[1]))-(fault_time-cycle_ms))]
            Ig_pre_avrg_C = sum(Ig_pre_C)/len(Ig_pre_C)
            Ig_C = RMS[5][fault_time+fix_ms:(len(RMS[1]) - fix_ms)]
            for i in range(len(Ig_C)):
                Ig_C[i] = Ig_C[i]-Ig_pre_avrg_C

        if fases_afetadas == 'A-T':
            for i in range(len(Ig_A)):
                m = abs(((Vg_A[i]*(Ig_A[i].conjugate())).imag)/((Z1_L*(Ig_A[i].conjugate())*(Ig_A[i]+Ig_pre_avrg_A+(K*seq_I[0][fault_time+fix_ms+i]))).imag)) #MONOFASICO
                local_falta.append(m)
        
        if fases_afetadas == 'B-T':
            for i in range(len(Ig_B)):
                m = abs(((Vg_B[i]*(Ig_B[i].conjugate())).imag)/((Z1_L*(Ig_B[i].conjugate())*(Ig_B[i]+Ig_pre_avrg_B+(K*seq_I[0][fault_time+fix_ms+i]))).imag)) #MONOFASICO
                local_falta.append(m)

        if fases_afetadas == 'C-T':
             for i in range(len(Ig_C)):
                 m = abs(((Vg_C[i]*(Ig_C[i].conjugate())).imag)/((Z1_L*(Ig_C[i].conjugate())*(Ig_C[i]+Ig_pre_avrg_C+(K*seq_I[0][fault_time+fix_ms+i]))).imag)) #MONOFASICO
                 local_falta.append(m)
        
        if fases_afetadas == 'AB' or fases_afetadas == 'AB-T' or fases_afetadas == 'ABC' or fases_afetadas == 'ABC-T':
            for i in range(len(Ig_A)):
                m = abs((((Vg_A[i]-Vg_B[i])*((Ig_A[i]-Ig_B[i]).conjugate())).imag)/((Z1_L*((Ig_A[i]+Ig_pre_avrg_A-Ig_B[i]-Ig_pre_avrg_B))*((Ig_A[i]-Ig_B[i]).conjugate())).imag)) #MONOFASICO
                local_falta.append(m)

        limiar_estabilidade = 0.1 
        fault_time_2 = fault_time + fix_ms
        final_local_falta = None
        for j in range(1, len(local_falta[fault_time_2:(len(local_falta) - cycle_ms)])):
            diferenca = abs(local_falta[fault_time_2+j] - local_falta[fault_time_2+j - 1])
            if diferenca < limiar_estabilidade:
                final_local_falta = sum(local_falta[fault_time_2+j:(len(local_falta) - cycle_ms)]) / len(local_falta[fault_time_2+j:(len(local_falta) - cycle_ms)])
                break 
        print("LOCAL DA FALTA FOI", final_local_falta)

        #print("LOCAL DA FALTA A/B/C:", m)
        if plotar_LF:
            plot_local_falta(local_falta,len(local_falta))
        
    def Tziouvaras_two_terminals(parametros, fases_afetadas, RMS, fault_time, seq_I, seq_I_2, Z_seq_r_1, Z_seq_i_1, Z_seq_r_2, Z_seq_i_2):

        m = symbols('m')
        a = []
        b = []
        c = []
        d = []
        e = []
        f = []
        Z1_L = (parametros['dadoslinha']['R1']+1j* parametros['dadoslinha']['X1'])*parametros['dadoslinha']['L']
        Z0_L = (parametros['dadoslinha']['R0']+1j* parametros['dadoslinha']['X0'])*parametros['dadoslinha']['L']
        Z2_L = Z1_L

        #METODO RUFINI
        # Z_seq_1 = (-0.096)+(1j*82.181)-(Z2_L)/2
        # Z_seq_0 = 4.032+(1j*124.611)-(Z2_L)/2

        Z_seq_1 = (-0.096)+(1j*86.481)-(Z2_L)/2
        Z_seq_0 = 4.032+(1j*127.611)-(Z2_L)/2
        #Z_seq_0 = Z_seq_1
        for j in range(len(seq_I[1])):
            a.append((seq_I[2][j]*(Z_seq_1)).real)
            b.append((seq_I[2][j]*(Z_seq_1)).imag)
            c.append((seq_I[2][j]*(Z2_L)).real)
            d.append((seq_I[2][j]*(Z2_L)).imag)
            e.append((Z_seq_0 + Z2_L).real) 
            f.append((Z_seq_0 + Z2_L).imag) 
        g = Z2_L.real
        h = Z2_L.imag

        for i in range(len(seq_I[1])):
                
            A = (abs(seq_I_2[2][i])**2) * (g**2 + h**2) - (c[i]**2 + d[i]**2)
            B = -2 *(abs(seq_I_2[2][i])**2) * (e[i]*g + f[i]*h) - 2 * (a[i]*c[i] + b[i]*d[i])
            C = (abs(seq_I_2[2][i])**2) * (e[i]**2 + f[i]**2) - (a[i]**2 + b[i]**2)
            equacao = Eq(A * m**2 + B * m + C, 0)
            solucoes = solve(equacao, m)
            modulos = [abs(solucao) for solucao in solucoes]
            
            #print(modulos)
            #print(solucoes)