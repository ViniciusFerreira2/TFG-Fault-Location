from interface import JanelaSelecaoArquivos
from processamentodados import processamento
import tkinter as tk

def main():
    parametros = processamento.carregar_parametros()
    root = tk.Tk()
    app = JanelaSelecaoArquivos(root, parametros)
    root.mainloop()

if __name__ == "__main__":
    main()