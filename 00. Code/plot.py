import os
import numpy as np
import matplotlib.pyplot as plt
import comtrade
import cmath
from datetime import datetime


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

def plot_XR(complexo, R1, X1, R0, X0):
    """
    Função para plotar os dados dos fasores (parte real no eixo X e parte imaginária no eixo Y),
    criando gráficos separados para cada componente de Z_seq.
    """
    
    # Plot para Z_seq[0]
    plt.figure(figsize=(10, 6))
    real_part = np.real(complexo[0])
    imag_part = np.imag(complexo[0])
    plt.plot(real_part, imag_part, linestyle='-', label='Fasor Z_seq[0]')
    
    # Linha da origem até 4x R0 e X0
    plt.plot([0, R0], [0, X0], linestyle='--', color='red', label='Linha até (4*R0, 4*X0)')
    
    # Configurar títulos e rótulos dos eixos
    plt.title("Gráfico de Fasores Z_seq[0]: Parte Real (X) vs Parte Imaginária (Y)")
    plt.xlabel("Parte Real")
    plt.ylabel("Parte Imaginária")
    plt.grid(True)
    plt.legend()
    plt.show()

    # Plot para Z_seq[1]
    plt.figure(figsize=(10, 6))
    real_part = np.real(complexo[1])
    imag_part = np.imag(complexo[1])
    plt.plot(real_part, imag_part, linestyle='-', label='Fasor Z_seq[1]')
    
    # Linha da origem até 4x R1 e X1
    plt.plot([0, 4 * R1], [0, 4 * X1], linestyle='--', color='red', label='Linha até (4*R1, 4*X1)')
    
    # Configurar títulos e rótulos dos eixos
    plt.title("Gráfico de Fasores Z_seq[1]: Parte Real (X) vs Parte Imaginária (Y)")
    plt.xlabel("Parte Real")
    plt.ylabel("Parte Imaginária")
    plt.grid(True)
    plt.legend()
    plt.show()

    # Plot para Z_seq[2]
    plt.figure(figsize=(10, 6))
    real_part = np.real(complexo[2])
    imag_part = np.imag(complexo[2])
    plt.plot(real_part, imag_part, linestyle='-', label='Fasor Z_seq[2]')
    
    # Linha da origem até 4x R1 e X1
    plt.plot([0, 4 * R1], [0, 4 * X1], linestyle='--', color='red', label='Linha até (4*R1, 4*X1)')
    
    # Configurar títulos e rótulos dos eixos
    plt.title("Gráfico de Fasores Z_seq[2]: Parte Real (X) vs Parte Imaginária (Y)")
    plt.xlabel("Parte Real")
    plt.ylabel("Parte Imaginária")
    plt.grid(True)
    plt.legend()
    plt.show()

def plot_Z_seq(Z_seq_mod, Z_seq_ang):
    """
    Função para plotar os módulos e ângulos de Z_seq em subplots.
    Z_seq_mod: Lista com os módulos da sequência de impedâncias.
    Z_seq_ang: Lista com os ângulos da sequência de impedâncias.
    """
    linha_seq1= (0.0166*223.3) + 1j*(0.2624*223.3)
    linha_seq0 = (0.4422*223.3) + 1j*(223.3*1.3189)
    linha_seq1_mod, linha_seq1_ang = cmath.polar(linha_seq1)
    linha_seq0_mod, linha_seq0_ang = cmath.polar(linha_seq0)

    # Criar a figura e os subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # Plotar o módulo de Z_seq no primeiro subplot
    for i in range(3):
        ax1.plot(Z_seq_mod[i], label=f'Módulo Z_seq[{i}]')
    
    # Adicionar as linhas constantes para linha_seq0_mod e linha_seq1_mod
    ax1.axhline(y=linha_seq0_mod, color='r', linestyle='--', label=f'Linha Seq0 (Módulo): {linha_seq0_mod:.2f}')
    ax1.axhline(y=linha_seq1_mod, color='b', linestyle='--', label=f'Linha Seq1 (Módulo): {linha_seq1_mod:.2f}')

    ax1.set_title('Módulo de Z_seq')
    ax1.set_xlabel('Índice')
    ax1.set_ylabel('Módulo')
    ax1.legend()
    ax1.grid(True)

    # Plotar o ângulo de Z_seq no segundo subplot
    for i in range(3):
        ax2.plot(np.degrees(Z_seq_ang[i]), label=f'Ângulo Z_seq[{i}]')  # Convertendo radianos para graus
    ax2.axhline(y=np.degrees(linha_seq0_ang), color='r', linestyle='--', label=f'Linha Seq0 (Ângulo): {np.degrees(linha_seq0_ang):.2f}°')
    ax2.axhline(y=np.degrees(linha_seq1_ang), color='b', linestyle='--', label=f'Linha Seq1 (Ângulo): {np.degrees(linha_seq1_ang):.2f}°')
    
    ax2.set_title('Ângulo de Z_seq')
    ax2.set_xlabel('Índice')
    ax2.set_ylabel('Ângulo (graus)')
    ax2.legend()
    ax2.grid(True)

    # Ajustar layout
    plt.tight_layout()

    # Mostrar os gráficos
    plt.show()