import numpy as np
import matplotlib.pyplot as plt
import comtrade
from processamentodados import process_signal, calculate_vrms, detectar_tipo_falta

def plotar_sinais(parametros, canais, label_tipo_falta, label_porcentagem_falta):
    vrms_final = []
    arquivo1 = parametros['arquivo1']
    arquivo2 = parametros['arquivo2']
    freq_corte_min = float(parametros['freq_corte_min'])
    freq_corte_max = float(parametros['freq_corte_max'])
    colunas_selecionadas = parametros['colunas']

    rec = comtrade.Comtrade()
    try:
        rec.load(arquivo1, arquivo2)
    except TypeError:
        print("Caminhos não foram definidos corretamente.")
        return
    signals = np.array(rec.analog)
    original_rate = 1000000
    target_rate = 1920
    cutoff_freq = freq_corte_max

    plt.figure(figsize=(12, 8))

    for coluna in colunas_selecionadas:
        sinal = signals[coluna]
        channel_index = coluna
        scale_factor = rec.cfg.analog_channels[channel_index].a
        offset = rec.cfg.analog_channels[channel_index].b  

        # Processamento do sinal
        vrms_values, time_new = calculate_vrms(sinal, target_rate, offset, scale_factor)
        vrms_final.append(vrms_values)
        sinal_processado = process_signal(vrms_values, original_rate, target_rate, cutoff_freq)
        num_samples_processado = len(sinal_processado)
        time_processado = np.linspace(0, num_samples_processado / target_rate, num_samples_processado)

        # Plotagem dos sinais
        plt.subplot(2, 1, 1)
        plt.plot(time_new, vrms_values, label=f'Vrms - Canal {canais[coluna]}')
        plt.title('Vrms dos Sinais')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Vrms')
        plt.legend()

        plt.subplot(2, 1, 2)
        plt.plot(np.arange(len(sinal)) / original_rate, sinal, label=f'Original - Canal {canais[coluna]}', alpha=0.5)
        plt.plot(time_processado, sinal_processado, label=f'Processado - Canal {canais[coluna]}')
        plt.title('Sinal Original e Processado')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Valor')
        plt.legend()

    plt.tight_layout()
    plt.show()

    # Detecção e exibição do tipo de falta e porcentagem
    tipo_falta, porcentagem_falta = detectar_tipo_falta(vrms_final)
    label_tipo_falta.config(text=f"Tipo de Falta Identificada: {tipo_falta}")
    label_porcentagem_falta.config(text=f"Porcentagem da Linha de Transmissão: {porcentagem_falta:.2f}%")
