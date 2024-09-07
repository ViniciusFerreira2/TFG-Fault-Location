import os
import numpy as np
import matplotlib.pyplot as plt
import comtrade
from datetime import datetime
from processamentodados import processamento #process_signal, calculate_vrms, detectar_tipo_falta

def plotar_sinais(parametros, canais):
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
        print("001_____ARQUIVOS IMPORTADOS COM SUCESSO")
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
    print(f"002_____CRIAÇÃO DA PASTA: {pasta_nome}")

    modulo = []
    angulo = []

    for coluna in colunas_selecionadas:
        coluna_index = rec.analog_channel_ids.index(coluna)  # Obtém o índice da coluna selecionada
        sinal = signals[coluna_index]
        scale_factor = rec.cfg.analog_channels[coluna_index].a
        offset = rec.cfg.analog_channels[coluna_index].b 
        time = np.linspace(0, len(sinal) / original_rate, len(sinal))

        # FILTRO DOS DADOS
        print("003_____CRIANDO SINAL FILTRADO")
        sinal_filtrado, time_filtrado, sample_rate = processamento.filter_signal(sinal, original_rate, target_rate, cutoff_freq)
        plot_filter(sinal, canais, coluna_index, time_filtrado, sinal_filtrado, original_rate, timestamp, pasta_nome)
        # print(sample_rate)

        print("004_____Realizando o RMS")
        sinal_vrms, time_vrms = processamento.vrm_per_phase(sinal_filtrado, time_filtrado)
        plot_vrms(time_filtrado, sinal_filtrado, canais, coluna_index, time_vrms, sinal_vrms, original_rate, timestamp, pasta_nome)

        print("005_____Estimando o FASOR")
        signal_modulo, signal_ang = processamento.fasor(sinal_filtrado, time_filtrado)

        modulo.append(signal_modulo)
        angulo.append(signal_ang)

    plot_fasor(modulo, angulo, pasta_nome)


    # Plotagem dos fasores no tempo zero
    if modulo and angulo:  # Verifica se há dados para plotar
        plotar_fasores(modulo, angulo, 18)
    
    print(f"XXX_____FINALIZADO")

def plot_filter(sinal, canais, coluna_index, time_processado, sinal_processado, original_rate, timestamp, pasta_nome):

    pasta_new = os.path.join(pasta_nome, "00_SINAIS_FILTRADO")
    os.makedirs(pasta_new, exist_ok=True)
    plt.figure(figsize=(12, 6))
    plt.plot(np.arange(len(sinal)) / original_rate, sinal, label=f'Original - Canal {canais[coluna_index]}', alpha=0.5)
    plt.plot(time_processado, sinal_processado, label=f'Processado - Canal {canais[coluna_index]}')

    # Configurações do gráfico
    plt.title('Sinal Original e Processado')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Valor')
    plt.legend()

    # Salva o gráfico
    sinal_fig_name = f"Sinal_Original_e_Filtrado_{coluna_index}.png"
    sinal_fig_path = os.path.join(pasta_new, sinal_fig_name)
    plt.savefig(sinal_fig_path)

    # Exibe e fecha o gráfico
    plt.show()  # Exibir o gráfico
    plt.close()

def plot_vrms(time, sinal, canais, coluna, time_processado, sinal_processado, original_rate, timestamp, pasta_nome):
        
        # Plotagem e salvamento do gráfico do Sinal Original e Processado
        pasta_new = os.path.join(pasta_nome, f"01_SINAIS VRMS")
        os.makedirs(pasta_new, exist_ok=True)
        plt.figure(figsize=(12, 6))
        plt.plot(time, sinal, label=f'Original - Canal {canais[coluna]}', alpha=0.5)
        plt.plot(time_processado, sinal_processado, label=f'VMRS - Canal {canais[coluna]}')
        plt.title('Sinal Original e Processado')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Valor')
        plt.legend()
        sinal_fig_name = f"Sinal_VRMS_{coluna}.png"
        print(f"    {sinal_fig_name}")
        sinal_fig_path = os.path.join(pasta_new, sinal_fig_name)
        plt.savefig(sinal_fig_path)
        plt.show()  # Exibir o gráfico
        plt.close()

def plot_fasor(mod, ang, pasta_nome):
    pasta_new = os.path.join(pasta_nome, f"02_SINAL&ANGULO")
    os.makedirs(pasta_new, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(2, 1)


    for array in mod:
        ax1.plot(array)
    ax1.set_title('modulo')

    # Plotar os arrays do segundo grupo no segundo eixo
    for array in ang:
        ax2.plot(array)
    ax2.set_title('angulo')

    sinal_fig_name = f"Sinal_SINAL E ANGULO.png"
    print(f"    {sinal_fig_name}")
    # Ajustar o layout para evitar sobreposição
    plt.tight_layout()

    sinal_fig_path = os.path.join(pasta_new, sinal_fig_name)
    plt.savefig(sinal_fig_path)

    # Mostrar a figura
    plt.show()

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
            
            print(f"FASE {i}")
            print(modulo[i][tempo_selecionado])
            print(angulo[i][tempo_selecionado])
        else:
            print(f"Erro: Tempo selecionado {tempo_selecionado} está fora do intervalo para a série {i}")

    ax.set_ylim(0, max([max(np.abs(m)) for m in modulo]))
    plt.title("Gráfico de Fasores")
    plt.show()