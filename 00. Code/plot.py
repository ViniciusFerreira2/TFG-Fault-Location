import os
import numpy as np
import matplotlib.pyplot as plt
import comtrade
from datetime import datetime
from processamentodados import processamento #process_signal, calculate_vrms, detectar_tipo_falta

def plotar_sinais(parametros, canais):
    vrms_final = []
    global modulo, angulo

    # Novos parâmetros para verificar se os gráficos serão plotados
    plotar_filtro = parametros.get('plotar_filtro', False)
    plotar_rms = parametros.get('plotar_rms', False)
    plotar_fasores = parametros.get('plotar_fasores', False)

    arquivo1 = parametros['arquivo1']
    arquivo2 = parametros['arquivo2']
    freq_amostragem = float(parametros['freq_amostragem'])
    freq_corte_max = float(parametros['freq_corte_max'])
    colunas_selecionadas = parametros['colunas']

    rec = comtrade.Comtrade()
    try:
        rec.load(arquivo1, arquivo2)
        print("001_____ARQUIVOS IMPORTADOS COM SUCESSO")
    except TypeError:
        print("Caminhos não foram definidos corretamente.")
        return

    signals = np.array(rec.analog)
    original_rate = 1000000
    target_rate = freq_amostragem
    cutoff_freq = freq_corte_max

    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")
    current_directory = os.path.dirname(os.path.abspath(__file__))
    pasta_nome = os.path.join(current_directory, f"Graficos_{timestamp}")
    os.makedirs(pasta_nome, exist_ok=True)
    print(f"002_____CRIAÇÃO DA PASTA: {pasta_nome}")

    modulo = []
    angulo = []
    complexo = []

    for coluna in colunas_selecionadas:
        coluna_index = rec.analog_channel_ids.index(coluna)
        sinal = signals[coluna_index]
        scale_factor = rec.cfg.analog_channels[coluna_index].a
        offset = rec.cfg.analog_channels[coluna_index].b 
        time = np.linspace(0, len(sinal) / original_rate, len(sinal))

        sinal_filtrado, time_filtrado, sample_rate = processamento.filter_signal(sinal, original_rate, target_rate, cutoff_freq)

        # Plotar o sinal filtrado somente se o checkbox correspondente estiver marcado
        if plotar_filtro:
            plot_filter(sinal, canais, coluna_index, time_filtrado, sinal_filtrado, original_rate, timestamp, pasta_nome)

        # Plotar o RMS somente se o checkbox correspondente estiver marcado
        if plotar_rms:
            sinal_vrms, time_vrms = processamento.vrm_per_phase(sinal_filtrado, time_filtrado)
            plot_vrms(time_filtrado, sinal_filtrado, canais, coluna_index, time_vrms, sinal_vrms, original_rate, timestamp, pasta_nome)

        # Plotar os fasores somente se o checkbox correspondente estiver marcado
        if plotar_fasores:
            signal_modulo, signal_ang, signal_complexo = processamento.fasor(sinal_filtrado, time_filtrado)
            modulo.append(signal_modulo)
            angulo.append(signal_ang)
            complexo.append(signal_complexo)
            plot_XR(modulo, angulo, complexo)

    # print("ANGULO",angulo[5][225:295])
    # for j in range(6):  
    #     for i in range(960):  
    #         angulo[j][i] = ajustar_angulo(angulo[j][i])

    angulo0 = angulo[0]

    for i in range(960):
        #angulo[0][i] = angulo[0][i] - angulo0[i]  
        # Ajuste das demais fases com relação à Fase A
        angulo[1][i] = angulo[1][i] - angulo[0][i]
        angulo[2][i] = angulo[2][i] - angulo[0][i]  
        angulo[3][i] = angulo[3][i] - angulo[0][i]  
        angulo[4][i] = angulo[4][i] - angulo[0][i]  
        angulo[5][i] = angulo[5][i] - angulo[0][i]
        angulo[0][i] = angulo[0][i] - angulo[0][i]  

   # plot_fasor(modulo, [angulo[2], angulo[5]], pasta_nome)

# correção para valores que ultrapassam 180° ou -180°
    for j in range(6):  
        for i in range(960):  
            angulo[j][i] = ajustar_angulo(angulo[j][i])

    if plotar_fasores and modulo and angulo:
        plot_fasor([angulo[0], angulo[2], angulo[4]], [angulo[1], angulo[3], angulo[5]], pasta_nome) #[angulo[1], angulo[3], angulo[5]], pasta_nome)
        #plot_fasor(, angulo, pasta_nome)
        #plotar_fasores(modulo, angulo, 18)



    # Plotagem dos fasores no tempo zero
    if modulo and angulo:  # Verifica se há dados para plotar
        plotar_fasores(modulo, angulo, 18)
    
    print(f"XXX_____FINALIZADO")

def ajustar_angulo(angulo):
    limiteang = 190
    if angulo > limiteang:
        angulo = -180 + (angulo - 180)
    elif angulo < -limiteang:
        angulo = 180 + (angulo + 180)
    return angulo


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

    # Plotar os arrays do primeiro grupo (mod) com legendas
    for idx, array in enumerate(mod):
        ax1.plot(array, label=f'Módulo {idx+1}')  # Adiciona uma legenda para cada array
    ax1.set_title('Módulo')
    ax1.legend()  # Exibir a legenda no gráfico de módulo

    # Plotar os arrays do segundo grupo (ang) com legendas
    for idx, array in enumerate(ang):
        ax2.plot(array, label=f'Ângulo {idx+1}')  # Adiciona uma legenda para cada array
    ax2.set_title('Ângulo')
    ax2.legend()  # Exibir a legenda no gráfico de ângulo

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

def plot_XR(modulo, angulo, sinal_complexo):
   # Separar a parte real e imaginária do sinal complexo
    parte_real = np.real(sinal_complexo)
    parte_imaginaria = np.imag(sinal_complexo)

    # Criar o gráfico
    plt.figure(figsize=(8, 6))
    
    # Plotar pontos menores com scatter
    plt.scatter(parte_real, parte_imaginaria, s=20, c='b', label='Sinal Complexo')  # 's' é o tamanho dos pontos
    
    # Adicionar setas para indicar a ordem dos pontos
    for i in range(len(parte_real) - 1):
        plt.quiver(parte_real[i], parte_imaginaria[i], 
                   parte_real[i+1] - parte_real[i], parte_imaginaria[i+1] - parte_imaginaria[i],
                   angles='xy', scale_units='xy', scale=1, color='r', width=0.002, headwidth=4)
    
    # Títulos e labels
    plt.title('Gráfico da Parte Real vs Parte Imaginária do Sinal')
    plt.xlabel('Parte Real')
    plt.ylabel('Parte Imaginária')

    # Adicionar grid e legenda
    plt.grid(True)
    #plt.legend()

    # Mostrar o gráfico
    plt.show()