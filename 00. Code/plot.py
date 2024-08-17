import os
import numpy as np
import matplotlib.pyplot as plt
import comtrade
from datetime import datetime
from processamentodados import process_signal, calculate_vrms, detectar_tipo_falta

def plotar_sinais(parametros, canais, label_tipo_falta, label_porcentagem_falta):
    vrms_final = []
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

    for coluna in colunas_selecionadas:
        sinal = signals[coluna]
        channel_index = coluna
        scale_factor = rec.cfg.analog_channels[channel_index].a
        offset = rec.cfg.analog_channels[channel_index].b  

        # Processamento do sinal
        vrms_values, time_new = calculate_vrms(sinal, target_rate, offset, scale_factor)
        vrms_final.append(vrms_values)
        sinal_processado = process_signal(sinal, original_rate, target_rate, cutoff_freq)
        num_samples_processado = len(sinal_processado)
        time_processado = np.linspace(0, num_samples_processado / target_rate, num_samples_processado)

        # Plotagem e salvamento do gráfico de Vrms
        plt.figure(figsize=(12, 6))
        plt.plot(time_new, vrms_values, label=f'Vrms - Canal {canais[coluna]}')
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
    tipo_falta, porcentagem_falta = detectar_tipo_falta(vrms_final)
    label_tipo_falta.config(text=f"Tipo de Falta Identificada: {tipo_falta}")
    label_porcentagem_falta.config(text=f"Porcentagem da Linha de Transmissão: {porcentagem_falta:.2f}%")
