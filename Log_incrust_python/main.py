import pandas as pd
import numpy as np
import pptx
import io
from pathlib import Path
from pptx.util import Pt
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
import matplotlib.pyplot as plt
from scipy.fft import fft # Using scipy's fft
from Graph_plot import graphs_plot, Data_processing, search, Pptx
from PIL import Image
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

def main(): 
    # %% Abertura do arquivo
    # The user has many commented-out file paths. We'll use the one that is active.
    # Original MATLAB: Table = xlsread("C:\Users\lucas\Downloads\Data Log29-05-2025 - 15h01min19s.xlsx", "Data Log");
    #file_path = r"/Users/lucasamorim/Library/CloudStorage/GoogleDrive-lucasamorim2004lt@gmail.com/My Drive/Escola/UFES VITÓRIA/Incrustacao_carbonatica/Incrustacao_log/EXCELS/Data Log31-03-2025 - 10h34min17s.xlsx"
    run_time = True
    sheet_name = "Data Log"
    while run_time:
        try:
            files_to_read, path_txt_file, path_txt_folder = search.search_data("Downloads")
        except TypeError:
            print(f"Nenhum arquivo com o prefixo Data_log foi encontrado na pasta Downloads")
            run_time = False
            return 0
            
        if files_to_read.shape[1] < 1:
            print(f"Todos os Logs foram lidos com sucesso!!!!!!!")
            run_time = False
            return 2
            
        for i in range(0, len(files_to_read)):
            figures = graphs_plot.graphs_figures(files_to_read[0,i],sheet_name)
            df_table = Data_processing.mean_std(files_to_read[0,i],sheet_name)
            Pptx.pptx_generation(figures, df_table, files_to_read[0,i].name,path_txt_folder)
            plt.close("all")
            search.write_on_text(path_txt_file, files_to_read[0,i].name)


        # Finally, show all plots
        #plt.show()

        print("Python script execution finished.")
    return 1

if __name__ == "__main__":
    a = main()
    if a == 1 or a == 2:
        print("Script runned succefully")
    else:
        print("Script was not succefully executed")
