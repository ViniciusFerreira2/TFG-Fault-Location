import os
import numpy as np
import matplotlib.pyplot as plt
import comtrade
from datetime import datetime
from processamentodados import process_signal, calculate_vrms, detectar_tipo_falta

def plotar_fasores(modulo, angulo, tempo_selecionado):
    """
    Plota um gráfico de fasores para os módulos e ângulos fornecidos em um tempo específico.
    """
    plt.figure(figsize=(8, 8))
    ax = plt.subplot(111, polar=True)

    for i in range(len(modulo)):
        if len(modulo[i]) > tempo_selecionado and len(angulo[i]) > tempo_selecionado:
            # Garantindo que os valores sejam reais
            magnitude = np.abs(modulo[i][tempo_selecionado])  # magnitude já é real
            theta = np.deg2rad(np.real(angulo[i][tempo_selecionado]))  # usa apenas a parte real do ângulo
            
            ax.arrow(theta, 0, 0, magnitude, 
                     head_width=0.05, head_length=0.1, fc='b', ec='b')
        else:
            print(f"Erro: Tempo selecionado {tempo_selecionado} está fora do intervalo para a série {i}")

    ax.set_ylim(0, max([max(np.abs(m)) for m in modulo]))
    plt.title("Gráfico de Fasores")
    plt.show()

def plotar_sinais(parametros, canais, label_tipo_falta, label_porcentagem_falta):
    vrms_final = []
    global modulo, angulo  # Variáveis para armazenar os valores de fasores

    arquivo1 = parametros['arquivo1']
    arquivo2 = parametros['arquivo2']
    freq_amostragem = float(parametros['freq_amostragem'])  # Frequência de amostragem fornecida na interface
    freq_corte_max = float(parametros['freq_corte_max'])    # Frequência de corte máxima fornecida na interface
    colunas_selecionadas = parametros['colunas']

    rec = comtrade.Comtrade()
    try:
        rec.load(arquivo1, arquivo2)
    except TypeError:
        print("Caminhos não foram definidos corretamente.")
        return
    signals = np.array(rec.analog)
    original_rate = 1000000  # Frequência original
    target_rate = freq_amostragem  # Frequência de amostragem alvo
    cutoff_freq = freq_corte_max  # Frequência de corte

    # Cria a pasta com data e hora no mesmo diretório onde o código está sendo executado
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")
    current_directory = os.path.dirname(os.path.abspath(__file__))  # Diretório atual do código
    pasta_nome = os.path.join(current_directory, f"Graficos_{timestamp}")
    os.makedirs(pasta_nome, exist_ok=True)

    modulo = []
    angulo = []

    for coluna in colunas_selecionadas:
        sinal = signals[coluna]
        channel_index = coluna
        scale_factor = rec.cfg.analog_channels[channel_index].a
        offset = rec.cfg.analog_channels[channel_index].b  

        # Processamento do sinal
        vrms_values, time_new, mod, ang = calculate_vrms(sinal, target_rate)
        
        if len(mod) == 0 or len(ang) == 0:
            print(f"Erro: Módulo ou ângulo vazios para a coluna {coluna}")
            continue  # Pula esta coluna se os valores estiverem vazios
        

        print(ang[16667])
        print(mod[16667])
        modulo.append(mod)
        angulo.append(ang)
        
        vrms_final.append(vrms_values)
        sinal_processado = process_signal(sinal, original_rate, target_rate, cutoff_freq)
        num_samples_processado = len(sinal_processado)
        time_processado = np.linspace(0, num_samples_processado / target_rate, num_samples_processado)

        # Plotagem e salvamento do gráfico de Vrms
        plt.figure(figsize=(12, 6))
        plt.plot(time_new, mod, label=f'Vrms - Canal {canais[coluna]}')
    
        plt.title('Vrms dos Sinais')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Vrms')
        plt.legend()
        vrms_fig_name = f"Vrms_dos_Sinais_{timestamp}.png"
        vrms_fig_path = os.path.join(pasta_nome, vrms_fig_name)
        plt.savefig(vrms_fig_path)
        plt.show()  # Exibir o gráfico
        plt.close()

        # Plotagem e salvamento do gráfico do Sinal Original e Processado
        plt.figure(figsize=(12, 6))
        plt.plot(np.arange(len(sinal)) / original_rate, sinal, label=f'Original - Canal {canais[coluna]}', alpha=0.5)
        plt.plot(time_processado, sinal_processado, label=f'Processado - Canal {canais[coluna]}')
        plt.title('Sinal Original e Processado')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Valor')
        plt.legend()
        sinal_fig_name = f"Sinal_Original_e_Processado_{timestamp}.png"
        sinal_fig_path = os.path.join(pasta_nome, sinal_fig_name)
        plt.savefig(sinal_fig_path)
        plt.show()  # Exibir o gráfico
        plt.close()

    # Detecção e exibição do tipo de falta e porcentagem
    #tipo_falta, porcentagem_falta = detectar_tipo_falta(vrms_final)
    #label_tipo_falta.config(text=f"Tipo de Falta Identificada: {tipo_falta}")
    #label_porcentagem_falta.config(text=f"Porcentagem da Linha de Transmissão: {porcentagem_falta:.2f}%")

    # Plotagem dos fasores no tempo zero
    if modulo and angulo:  # Verifica se há dados para plotar
        plotar_fasores(modulo, angulo, 0)

