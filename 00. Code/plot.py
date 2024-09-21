import os
import numpy as np
import matplotlib.pyplot as plt
import comtrade
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

def plot_XR(complexo):
    """
    Função para plotar os dados dos fasores (parte real no eixo X e parte imaginária no eixo Y),
    destacando o início, meio e fim das ligações.
    """
    # Configurar o tamanho da figura
    plt.figure(figsize=(10, 6))

    # Para cada conjunto de dados (coluna) no complexo
    for i in range(0, len(complexo)):
        # Extrair parte real (eixo X) e parte imaginária (eixo Y)
        real_part = np.real(complexo[i])
        imag_part = np.imag(complexo[i])


        # Plotar a linha conectando todos os pontos
        plt.plot(real_part, imag_part, linestyle='-', label=f'Fasor {i+1}')

        # # Destacar o início (círculo verde), meio (quadrado azul) e fim (triângulo vermelho)
        # plt.plot(real_part[inicio], imag_part[inicio], 'go', label=f'Início {i+1}')  # Círculo verde
        # plt.plot(real_part[meio], imag_part[meio], 'bs', label=f'Meio {i+1}')  # Quadrado azul
        # plt.plot(real_part[fim], imag_part[fim], 'r^', label=f'Fim {i+1}')  # Triângulo vermelho

    # Configurar títulos e rótulos dos eixos
    plt.title("Gráfico de Fasores: Parte Real (X) vs Parte Imaginária (Y)")
    plt.xlabel("Parte Real")
    plt.ylabel("Parte Imaginária")

    # Adicionar grid e legenda
    plt.grid(True)
    plt.legend()

    # Exibir o gráfico
    plt.show()