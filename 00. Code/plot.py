import os
import numpy as np
import matplotlib.pyplot as plt
import comtrade
import cmath
from datetime import datetime
import matplotlib

matplotlib.use("QtAgg")

def plot_filter(sinal, canais, coluna_index, time_processado, sinal_processado, original_rate, timestamp, pasta_nome):

    pasta_new = os.path.join(pasta_nome, "00_SINAIS_FILTRADO")
    os.makedirs(pasta_new, exist_ok=True)
    plt.figure(figsize=(12, 6))
    plt.plot(np.arange(len(sinal)) / original_rate, sinal, label=f'Original - Canal {canais[coluna_index]}', alpha=0.5)
    plt.plot(time_processado, sinal_processado, label=f'Processado - Canal {canais[coluna_index]}')
    plt.legend(loc='upper right')
    # Configurações do gráfico
    #plt.title('Sinal Original e Processado')
    plt.xlabel('Tempo (s)')
    if '-' in canais[coluna_index]:
        plt.ylabel('Corrente [A]')
    else: plt.ylabel('Tensão [V]')
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
        pasta_new = os.path.join(pasta_nome, f"01_SINAIS RMS")
        os.makedirs(pasta_new, exist_ok=True)
        plt.figure(figsize=(12, 6))
        plt.plot(time, sinal, label=f'Processado - Canal {canais[coluna]}', alpha=0.5)
        plt.plot(time_processado, sinal_processado, label=f'RMS - Canal {canais[coluna]}')
        #plt.title('Sinal Original e Processado')
        plt.xlabel('Tempo (s)')
        if '-' in canais[coluna]:
            plt.ylabel('Corrente [A]')
        else: plt.ylabel('Tensão [V]')
        plt.legend(loc='upper right')
        sinal_fig_name = f"Sinal_RMS_{coluna}.png"
        print(f"    {sinal_fig_name}")
        sinal_fig_path = os.path.join(pasta_new, sinal_fig_name)
        plt.savefig(sinal_fig_path)
        plt.show()  # Exibir o gráfico
        plt.close()

def plot_fasor(mod, ang, pasta_nome):
    # Criar a pasta de destino, se ainda não existir
    pasta_new = os.path.join(pasta_nome, f"02_SINAL&ANGULO")
    os.makedirs(pasta_new, exist_ok=True)

    # Dicionário de legendas para o grupo 'mod' e o grupo 'ang'
    legendas_mod = {
        1: "Módulo Tensão øA",
        2: "Módulo Corrente øA",
        3: "Módulo Tensão øB",
        4: "Módulo Corrente øB",
        5: "Módulo Tensão øC",
        6: "Módulo Corrente øC"
    }

    legendas_ang = {
        1: "Ângulo Tensão øA",
        2: "Ângulo Corrente øA",
        3: "Ângulo Tensão øB",
        4: "Ângulo Corrente øB",
        5: "Ângulo Tensão øC",
        6: "Ângulo Corrente øC"
    }

    # Figura 1: Gráficos dos Módulos
    fig1, (ax1, ax2) = plt.subplots(2, 1)

    for idx, array in enumerate(mod):
        label = legendas_mod.get(idx + 1, f'Módulo {idx + 1}')
        num_amostras = len(array)
        x = np.linspace(0, 0.2, num_amostras)  # Cria o eixo x para 0 a 0.2

        if (idx + 1) % 2 != 0:
            ax1.plot(x, array, label=label)
            ax1.set_xlim(0, 0.2)
        else:
            ax2.plot(x, array, label=label)
            ax2.set_xlim(0, 0.2)

    ax1.set_title('Módulo - Tensão')
    ax1.set_xlabel('Tempo [s]')     
    ax1.set_ylabel('Módulo [V]')
    ax1.legend(loc='upper right')  
    
    ax2.set_title('Módulo - Corrente')
    ax2.set_xlabel('Tempo [s]')  
    ax2.set_ylabel('Módulo [A]')    
    ax2.legend(loc='upper right')  

    # Ajustar layout da primeira figura
    fig1.tight_layout()
    # Salvar a primeira figura
    fig1.savefig(os.path.join(pasta_new, "Modulos.png"))
    plt.show()  # Mostrar a figura dos módulos

    # Figura 2: Gráficos dos Ângulos
    fig2, (ax3, ax4) = plt.subplots(2, 1)

    # Plotar os arrays do segundo grupo (ang) com legendas específicas
    for idx, array in enumerate(ang):
        label = legendas_ang.get(idx + 1, f'Ângulo {idx + 1}')
        num_amostras = len(array)
        x = np.linspace(0, 0.2, num_amostras)  # Cria o eixo x para 0.2

        if (idx + 1) % 2 != 0:
            ax3.plot(x, array, label=label)
            ax3.set_xlim(0, 0.2)
        else:
            ax4.plot(x, array, label=label)
            ax4.set_xlim(0, 0.2)

    ax3.set_title('Ângulo - Tensão')
    ax3.set_xlabel('Tempo [s]')  # Legenda para o eixo X do primeiro subplot
    ax3.set_ylabel('Ângulo [°]')  # Legenda para o eixo Y do primeiro subplot
    ax3.legend(loc='upper right')  # Exibir a legenda no primeiro subplot de ângulo
    
    ax4.set_title('Ângulo - Corrente')
    ax4.set_xlabel('Tempo [s]')  # Legenda para o eixo X do segundo subplot
    ax4.set_ylabel('Ângulo [°]')  # Legenda para o eixo Y do segundo subplot
    ax4.legend(loc='upper right')  # Exibir a legenda no segundo subplot de ângulo

    # Ajustar layout da segunda figura
    fig2.tight_layout()
    # Salvar a segunda figura
    fig2.savefig(os.path.join(pasta_new, "Angulos.png"))
    plt.show()  # Mostrar a figura dos ângulos

def plotar_polarformat(modulo, angulo, tempo_selecionado):
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

def plot_XR(complexo, parametros):

    R1, X1, R0, X0, L = parametros['dadoslinha']['R1'], parametros['dadoslinha']['X1'], parametros['dadoslinha']['R0'], parametros['dadoslinha']['X0'], parametros['dadoslinha']['L']
    R1=R1*L
    X1=X1*L
    R0=R0*L
    X0=X0*L

    limit_R0, limit_X0 = 3 * R0, 1 * X0
    limit_R1, limit_X1 = 3 * R1, 1 * X1
    
    # Plot para Z_seq[0]
    plt.figure(figsize=(10, 6))
    real_part = np.real(complexo[0])
    imag_part = np.imag(complexo[0])
    plt.plot(real_part, imag_part, linestyle='-', label='Z Seq. (0)')
    
    plt.plot([0, 4*R0], [0, 4*X0], linestyle='--', color='red', label='Impedância da Linha')
    plt.plot([R0, R0 + 0.2 * R0], [X0, X0], linestyle='-', color='green', label='Impedância Total da Linha')
    plt.plot([R0, R0 - 0.2 * R0], [X0, X0], linestyle='-', color='green')

    #plt.xlim(-limit_R0, limit_R0)
    #plt.ylim(0, limit_X0)
    
    # Configurar títulos e rótulos dos eixos
    plt.title("Z Seq. (0): X/R")
    plt.xlabel("X[Ω]")
    plt.ylabel("R[Ω]")
    plt.grid(True)
    plt.legend()
    plt.show()

    # Plot para Z_seq[1]
    plt.figure(figsize=(10, 6))
    real_part = np.real(complexo[1])
    imag_part = np.imag(complexo[1])
    plt.plot(real_part, imag_part, linestyle='-', label='Z Seq. (+)')
    
    # Linha da origem até 4x R1 e X1
    plt.plot([0, 4*R1], [0, 4*X1], linestyle='--', color='red', label='Impedância da Linha')
    plt.plot([R1, R1 + 1.5 * R1], [X1, X1], linestyle='-', color='green', label='Impedância Total da Linha')
    plt.plot([R1, R1 - 1.5 * R1], [X1, X1], linestyle='-', color='green')

    #plt.xlim(-limit_R1, limit_R1)
    #plt.ylim(0, limit_X1)
    
    # Configurar títulos e rótulos dos eixos
    plt.title("Z Seq. (+): X/R")
    plt.xlabel("R[Ω]")
    plt.ylabel("X[Ω]")
    plt.grid(True)
    plt.legend()
    plt.show()

    # Plot para Z_seq[2]
    plt.figure(figsize=(10, 6))
    real_part = np.real(complexo[2])
    imag_part = np.imag(complexo[2])
    plt.plot(real_part, imag_part, linestyle='-', label='Z Seq. (-)]')
    
    # Linha da origem até 4x R1 e X1
    plt.plot([0, R1], [0, X1], linestyle='--', color='red', label='Impedância da linha')
    
    # Configurar títulos e rótulos dos eixos
    plt.title("Z Seq. (-): X/R")
    plt.xlabel("R[Ω]")
    plt.ylabel("X[Ω]")
    plt.grid(True)
    plt.legend()
    plt.show()

def plot_Z_seq(parametros, Z_seq_mod, Z_seq_ang):
    legendas_Z = {
        0: "0",  # Altere a chave para começar de 0 para corresponder ao índice
        1: "+",
        2: "-",
    }

    linha_seq1 = (parametros['dadoslinha']['R1'] + 1j * parametros['dadoslinha']['X1']) * parametros['dadoslinha']['L']
    linha_seq0 = (parametros['dadoslinha']['R0'] + 1j * parametros['dadoslinha']['X0']) * parametros['dadoslinha']['L']

    linha_seq1_mod, linha_seq1_ang = cmath.polar(linha_seq1)
    linha_seq0_mod, linha_seq0_ang = cmath.polar(linha_seq0)

    # Criar a figura e os subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    num_amostras = len(Z_seq_mod[0])  # Presumindo que todas as sequências têm o mesmo tamanho
    tempo = np.linspace(0, 0.2, num_amostras)  # Eixo x de 0 a 0.2 com base no número de amostras

    # Plotar o módulo de Z_seq no primeiro subplot
    for i in range(3):
        ax1.plot(tempo, Z_seq_mod[i], label=f'Módulo Z seq. ({legendas_Z[i]})')  # Adiciona a legenda

    # Adicionar as linhas constantes para linha_seq0_mod e linha_seq1_mod
    ax1.axhline(y=linha_seq0_mod, color='r', linestyle='--', label=f'Z da linha - Seq. 0: {linha_seq0_mod:.2f} Ω')
    ax1.axhline(y=linha_seq1_mod, color='b', linestyle='--', label=f'Z da linha - Seq. (+): {linha_seq1_mod:.2f} Ω')

    ax1.set_title('Módulo de Z')
    ax1.set_xlabel('Tempo [s]')  # Alterado para refletir a nova escala do eixo x
    ax1.set_ylabel('Módulo [Ω]')
    ax1.legend(loc='upper right')
    ax1.grid(True)

    # Plotar o ângulo de Z_seq no segundo subplot
    for i in range(3):
        ax2.plot(tempo, np.degrees(Z_seq_ang[i]), label=f'Ângulo Z seq. [{i}] ({legendas_Z[i]})')  # Adiciona a legenda
    ax2.axhline(y=np.degrees(linha_seq0_ang), color='r', linestyle='--', label=f'Impedância da linha - Seq. (0): {np.degrees(linha_seq0_ang):.2f}°')
    ax2.axhline(y=np.degrees(linha_seq1_ang), color='b', linestyle='--', label=f'Impedância da linha - Seq. (+): {np.degrees(linha_seq1_ang):.2f}°')

    ax2.set_title('Ângulo de Z')
    ax2.set_xlabel('Tempo [s]')  # Alterado para refletir a nova escala do eixo x
    ax2.set_ylabel('Ângulo [°]')
    ax2.legend(loc='lower right')
    ax2.grid(True)

    plt.tight_layout()
    plt.show()

def plot_local_falta(m, len_Ig_A):

    t = np.linspace(0, 0.2, len_Ig_A)
    
    plt.figure(figsize=(10, 6))
    plt.plot(t, m, label="Local da Falta", color="blue")
    plt.xlabel("Tempo (s)")
    plt.ylabel("Local da falta")
    plt.title("Gráfico de m em função do tempo")
    plt.legend()
    plt.grid(True)
    plt.show()