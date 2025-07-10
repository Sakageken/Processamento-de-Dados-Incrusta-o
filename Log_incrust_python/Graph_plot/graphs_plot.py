import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft # Using scipy's fft
from matplotlib.ticker import StrMethodFormatter

def graphs_figures(filepath: str, sheetname: str):
    file_path = filepath
    sheet_name = sheetname

    # Use pandas to read the Excel file.
    # Assuming the Excel sheet does not have headers that pandas should use,
    # or that the data starts from the very first row.
    # If your data has headers in the first row, you might remove `header=None`.
    try:
        table_df = pd.read_excel(file_path, sheet_name=sheet_name, header=0)
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        exit()
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        exit()

    # %% Leitura das variáveis da planilha

    # Criação de referencia no tempo
    t_exp = 5  # seconds, sampling interval for each row in the log

    # timeref = 0: t_exp/60: t_exp*(length(Table(1:end-1,1))/60);
  
    num_rows = table_df.shape[0]
    end_time_minutes = t_exp * (num_rows - 1) / 60
    timeref = np.linspace(0, end_time_minutes, num_rows)
    
    

    #print(table_df.columns)
    # Vazão
    fit01 = table_df['FIT-01 - Flow R. [kg/h]'].astype(float) 
    fit02 = table_df['FIT-02 - Flow R. [kg/h]'].astype(float)
    fit03 = table_df['FIT-03 - Flow R. [kg/h]'].astype(float) # CO2 Gasoso
    fit04 = table_df['FIT-04 - Flow R. [kg/h]'].astype(float) # CO2 Liquido

    # Estimativa de consumo de CO2
    # var_tempo = (-Table(1:end-1,1)+Table(2:end,1))/3600;
    # Assuming column 1 in MATLAB is table_df.iloc[:, 0]
    time_column_seconds = table_df.iloc[:, 0]
    var_tempo_hours = (time_column_seconds.iloc[1:].values - time_column_seconds.iloc[:-1].values) / 3600.0

    fit03_med = (fit03[:-1].values + fit03[1:].values) / 2.0
    fit04_med = (fit04[:-1].values + fit04[1:].values) / 2.

    total_co2_g = np.sum(var_tempo_hours * fit03_med)
    total_co2_l = np.sum(var_tempo_hours * fit04_med)

    # Restrição trecho de testes
    pv05ret = table_df['PV-05 Ret. [%]'].astype(float)
    pv08ret = table_df['PV-08 Ret. [%]'].astype(float)

    # Pressão co2
    pt03 = table_df['PT-03 [bar]'].astype(float)
    pt04 = table_df['PT-04 [bar]'].astype(float)

    # Pressão de injeção das soluções
    pt01a = table_df['PT-01a [bar]'].astype(float)
    pt01b = table_df['PT-01b [bar]'].astype(float)
    pt02a = table_df['PT-02a [bar]'].astype(float)
    pt02b = table_df['PT-02b [bar]'].astype(float)

    # Pressão trecho de testes
    pt05 = table_df['PT-05 [bar]'].astype(float)
    pt07 = table_df['PT-07 [bar]'].astype(float)
    pt10 = table_df['PT-10 [bar]'].astype(float)
    pt08 = table_df['PT-08 [bar]'].astype(float)

    # pH
    ph_read1 = table_df['pH 1'].astype(float)
    ph_read2 = table_df['pH 2'].astype(float)

    # Temperatura
    tt01 = table_df['TT-01 [°C]'].astype(float)
    tt02 = table_df['TT-02 [°C]'].astype(float)
    tt05 = table_df['TT-05 [°C]'].astype(float)
    tt07 = table_df['TT-07 [°C]'].astype(float)
    tt08 = table_df['TT-08 [°C]'].astype(float)
    tt10 = table_df['TT-10 [°C]'].astype(float)

    # Funcionamento da bomba [%]
    p01 = table_df['P-01'].astype(float)
    p02 = table_df['P-02'].astype(float)

    # Vazao das solucoes nas linhas [kg/h]
    fit05 = table_df['FIT-05 - Flow R. [kg/h]'].astype(float) # Vazao linha velha
    fit06 = table_df['FIT-06 - Flow R. [kg/h]'].astype(float) # Vazao linha nova

    # Cálculo das médias moveis (using pandas rolling mean)
    # center=True to mimic MATLAB's default behavior for odd window sizes. min_periods=1 to include edges.
    window_size = 50
    ma_pt01a = pt01a.rolling(window=window_size, center=True, min_periods=1).mean()
    ma_pt01b = pt01b.rolling(window=window_size, center=True, min_periods=1).mean()
    ma_pt02a = pt02a.rolling(window=window_size, center=True, min_periods=1).mean()
    ma_pt02b = pt02b.rolling(window=window_size, center=True, min_periods=1).mean()
    ma_pt05 = pt05.rolling(window=window_size, center=True, min_periods=1).mean()
    ma_pt07 = pt07.rolling(window=window_size, center=True, min_periods=1).mean()
    ma_pt08 = pt08.rolling(window=window_size, center=True, min_periods=1).mean()
    ma_pt10 = pt10.rolling(window=window_size, center=True, min_periods=1).mean()
    ma_fit03 = fit03.rolling(window=window_size, center=True, min_periods=1).mean()
    ma_fit04 = fit04.rolling(window=window_size, center=True, min_periods=1).mean()

    # Cálculo da eficiencia temporal de teste
    k1 = 0
    while k1 < len(fit01) and fit01.iloc[k1] == 0:
        k1 += 1

    k2 = 0
    while k2 < len(fit02) and fit02.iloc[k2] == 0:
        k2 += 1

    k3 = 0 # MATLAB k3 counts from 1 up to length
    # (length(fit01)-k3) in MATLAB is (len(fit01)-1-k3_py) in Python
    while k3 < len(fit01) and p01.iloc[len(fit01) - 1 - k3] == 0:
        k3 += 1

    k4 = 0
    while k4 < len(fit02) and p02.iloc[len(fit02) - 1 - k4] == 0:
        k4 += 1

    tat = len(fit01) * t_exp / 60.0 # Total acquisition time
    # Effective start index for TTT based on when flows/pumps start
    # Effective end index for TTT based on when flows/pumps stop
    start_delay_points = (k1 + k2) / 2.0
    end_cutoff_points = (k3 + k4) / 2.0
    ttt = (len(fit01) - (start_delay_points + end_cutoff_points)) * t_exp / 60.0 # Total test time (excluding startup/shutdown)

    # tet = ((nnz(fit01)+nnz(fit02))/2)*t_exp/60;
    # nnz is number of non-zero elements
    tet = ((np.count_nonzero(fit01) + np.count_nonzero(fit02)) / 2.0) * t_exp / 60.0 # Effective experiment time
    ett = tet / ttt if ttt != 0 else 0 # Efficiency of test time

    # %% Plots Parte 1
    width = 1500/100
    height = 900/100
    fig1, axs1 = plt.subplots(4,1,figsize=(width,height))

    for axs in axs1.flat:
        axs.grid(True)
        axs.set_xlim(0, max(timeref)+1)
        axs.set_xlabel("Tempo [minutos]")

    ## Vazao solucoes
    axs1[0].plot(timeref,fit01, label='$NaHCO_{3}$')
    axs1[0].plot(timeref,fit02, label='$CaCl_{2}$')
    axs1[0].set_title('Vazão das soluções')
    axs1[0].set_ylabel('Vazão mássica [kg/h]')
    axs1[0].legend(loc='upper right')

    ## Vazão CO2
    axs1[1].plot(timeref, fit03, label='$CO_{2}$ Gasoso')
    axs1[1].plot(timeref, fit04, label='$CO_{2}$ Líquido')
    axs1[1].set_title('Vazão de CO2') 
    axs1[1].set_ylabel('Vazão mássica [kg/h]')
    axs1[1].legend(loc='upper right')
    axs1[1].yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

    ## Restrição das valvulas
    axs1[2].plot(timeref, pv05ret, label='PV-05')
    axs1[2].plot(timeref, pv08ret, label='PV-08')
    axs1[2].set_title('Restrição da Válvula') 
    axs1[2].set_ylabel('Restrição [%]')
    axs1[2].legend(loc='upper right')

    ## Média móvel da pressão pt05, pt07, pt08 e pt10
    axs1[3].plot(timeref, ma_pt05, label='PT-05')
    axs1[3].plot(timeref, ma_pt07, label='PT-07')
    axs1[3].plot(timeref, ma_pt08, label='PT-08')
    axs1[3].plot(timeref, ma_pt10, label='PT-10')
    axs1[3].set_title('Pressão no trecho de testes - Média Móvel') 
    axs1[3].set_ylabel('Pressão [bar]')
    axs1[3].legend(loc='upper right')

    fig1.subplots_adjust(left=0.055, bottom=0.06, right=0.983, top=0.949, hspace=0.496, wspace=0.151)
    fig1.tight_layout()

    # %% Plots Parte 2

    fig2, axs2 = plt.subplots(4,1,figsize=(width,height))

    for axs in axs2.flat:
        axs.grid(True)
        axs.set_xlim(0, max(timeref)+1)
        axs.set_xlabel('Tempo [minutos]')

    # Pressão na saída da bomba
    axs2[0].plot(timeref, pt01a, label='PT-01a')
    axs2[0].plot(timeref, pt02a, label='PT-02a')
    axs2[0].set_title('Pressão na saída das bombas')
    axs2[0].set_ylabel('Pressão [bar]')
    axs2[0].legend(loc='upper right')

    # Media movel da pressão na saída das bombas
    axs2[1].plot(timeref, ma_pt01a, label='PT-01a')
    axs2[1].plot(timeref, ma_pt02a, label='PT-02a')
    axs2[1].plot(timeref, ma_pt01b, label='PT-01b')
    axs2[1].plot(timeref, ma_pt02b, label='PT-02b')
    axs2[1].set_title('Pressão na saída das bombas - Média Móvel') 
    axs2[1].set_ylabel('Pressão [bar]')
    axs2[1].legend(loc='upper right')

    # Pressão pt05, pt07, pt08 e pt10
    axs2[2].plot(timeref, pt05, label='PT-05')
    axs2[2].plot(timeref, pt07, label='PT-07')
    axs2[2].plot(timeref, pt08, label='PT-08')
    axs2[2].plot(timeref, pt10, label='PT-10')
    axs2[2].set_title('Pressão no trecho de testes') 
    axs2[2].set_ylabel('Pressão [bar]')
    axs2[2].legend(loc='upper right')

    ## Média móvel da pressão pt05, pt07, pt08 e pt10
    axs2[3].plot(timeref, ma_fit03, label='Gas')
    axs2[3].plot(timeref, ma_fit04, label='Líquido')
    axs2[3].set_title('Vazão de $CO_{2}$ - Média Móvel') 
    axs2[3].set_ylabel('Vazão [kg/h]')
    axs2[3].legend(loc='upper right')
    axs2[3].yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

    fig2.subplots_adjust(left=0.055, bottom=0.06, right=0.983, top=0.949, hspace=0.496, wspace=0.151)
    fig2.tight_layout()


    # %% Plots Parte 3
    fig3, axs3 = plt.subplots(4,1,figsize=(width,height))

    for axs in axs3[:-1].flat:
        axs.grid(True)
        axs.set_xlim(0, max(timeref)+1)
        axs.legend(loc='upper right')
        axs.set_xlabel('Tempo [minutos]')

    # pH
    axs3[0].plot(timeref, ph_read1, label='pH 1')
    axs3[0].plot(timeref, ph_read2, label='pH 2')
    axs3[0].set_ylim(0,14)
    axs3[0].set_title('pH')
    axs3[0].set_ylabel('pH')

    # Temperatura
    axs3[1].plot(timeref, tt01, label='TT-01')
    axs3[1].plot(timeref, tt05, label='TT-05')
    axs3[1].plot(timeref, tt08, label='TT-08')
    axs3[1].set_title('Temperatura') 
    axs3[1].set_ylabel('T [°C]')

    # Pressão na injeção de CO2
    axs3[2].plot(timeref, pt03, label='PT-03 (Gas)')
    axs3[2].plot(timeref, pt04, label='PT-04 (Líquido)')
    axs3[2].set_title('Pressão na injeção de $CO_{2}$') 
    axs3[2].set_ylabel('Pressão [bar]')

    ## Média móvel da pressão pt05, pt07, pt08 e pt10
    labels_co2 = ['Gás', 'Líquido']
    largura_barra = 0.6
    plot_co2_values = [total_co2_g, total_co2_l]
    axs3[3].bar(labels_co2, plot_co2_values, width=largura_barra, label=labels_co2)
    axs3[3].set_title('(1) Consumo de $CO_{2}$ gasoso, (2) Consumo de $CO_{2}$ Líquido') 
    axs3[3].set_ylabel('Consumo [kg]')
    axs3[3].yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

    fig3.subplots_adjust(left=0.055, bottom=0.06, right=0.983, top=0.949, hspace=0.496, wspace=0.151)
    fig3.tight_layout()




    # %% Plots parte 4
    fig4, axs4 = plt.subplots(1,1,figsize=(width,height))

    axs4.plot(timeref, fit05, label='Vazão linha velha')
    axs4.plot(timeref, fit06, label='Vazão linha nova')
    axs4.set_title('Vazão das soluções na linha')
    axs4.set_ylabel('Vazão mássica [kg/h]')
    axs4.set_xlabel('Tempo [minutos]')
    axs4.grid(True)
    axs4.set_xlim(0, max(timeref)+1)
    axs4.legend(loc='upper right')
    fig4.subplots_adjust(left=0.055, bottom=0.06, right=0.983, top=0.949, hspace=0.496, wspace=0.151)

    plt.tight_layout()


    # %% Plot todos juntos (Tiled Layout)
    fig5, axs5 = plt.subplots(4, 3, figsize=(width,height), tight_layout=True) # tight_layout similar to "TileSpacing","tight"
    axs5 = axs5.ravel() # Flatten the 2D array of axes for easy iteration like nexttile
    fig5.subplots_adjust(left=0.055, bottom=0.06, right=0.983, top=0.949, hspace=0.496, wspace=0.151)

    for axs in axs5[:-1]:
        axs.grid(True)
        axs.set_xlim(0, max(timeref)+1)
        axs.legend(loc='upper right')
        axs.set_xlabel('Tempo [minutos]')

    # Plot 1: Vazão soluções
    axs5[0].plot(timeref, fit01, label='$NaHCO_{3}$')
    axs5[0].plot(timeref, fit02, label='$CaCl_{2}$')
    axs5[0].set_title('Vazão das soluções')
    axs5[0].set_ylim([0, 500])
    axs5[0].set_yticks(np.arange(0, 501, 100))
    axs5[0].set_ylabel('Vazão mássica [kg/h]')


    # Plot 2: Pressão na saída da bomba
    axs5[1].plot(timeref, pt01a, label='PT-01a')
    axs5[1].plot(timeref, pt02a, label='PT-02a')
    axs5[1].set_ylim([0, 200])
    axs5[1].set_title('Pressão na saída das bombas')
    axs5[1].set_ylabel('Pressão [bar]')


    # Plot 3: pH
    axs5[2].plot(timeref, ph_read1, label='pH 1')
    axs5[2].plot(timeref, ph_read2, label='pH 2')
    axs5[2].set_ylim([0, 14])
    axs5[2].set_yticks([0, 3.5, 7, 10.5, 14])
    axs5[2].set_title('pH')
    axs5[2].set_ylabel('pH')


    # Plot 4: Vazão CO2
    axs5[3].plot(timeref, fit03, label='$CO_{2}$ Gasoso')
    axs5[3].plot(timeref, fit04, label='$CO_{2}$ Líquido')
    axs5[3].set_title('Vazão de CO2') # Corrected
    axs5[3].set_ylim([0, 50])
    axs5[3].set_ylabel('Vazão mássica [kg/h]')
    axs5[3].yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))


    # Plot 5: Media movel da pressão na saída das bombas
    axs5[4].plot(timeref, ma_pt01a, label='PT-01a')
    axs5[4].plot(timeref, ma_pt02a, label='PT-02a')
    axs5[4].plot(timeref, ma_pt01b, label='PT-01b')
    axs5[4].plot(timeref, ma_pt02b, label='PT-02b')
    axs5[4].set_ylim([0, 200])
    axs5[4].set_title('Pressão na saída das bombas - Média Móvel')
    axs5[4].set_ylabel('Pressão [bar]')


    # Plot 6: Temperatura
    axs5[5].plot(timeref, tt01, label='TT-01')
    axs5[5].plot(timeref, tt05, label='TT-05')
    axs5[5].plot(timeref, tt08, label='TT-08')
    axs5[5].set_title('Temperatura')
    axs5[5].set_ylabel('°C')


    # Plot 7: Restrição
    axs5[6].plot(timeref, pv05ret, label='PV-05')
    axs5[6].plot(timeref, pv08ret, label='PV-08')
    axs5[6].set_title('Restrição da Válvula')
    axs5[6].set_ylabel('Restrição [%]')


    # Plot 8: Pressão pt05, pt07, pt08 e pt10
    axs5[7].plot(timeref, pt05, label='PT-05')
    axs5[7].plot(timeref, pt07, label='PT-07')
    axs5[7].plot(timeref, pt08, label='PT-08')
    axs5[7].plot(timeref, pt10, label='PT-10')
    axs5[7].set_ylim([0, 200])
    axs5[7].set_title('Pressão no trecho de testes')
    axs5[7].set_ylabel('Pressão [bar]')


    # Plot 9: Pressao CO2
    axs5[8].plot(timeref, pt03, label='PT-03 (Gas)')
    axs5[8].plot(timeref, pt04, label='PT-04 (Líquido)')
    axs5[8].set_ylim([0, 120])
    axs5[8].set_title('Pressão na injeção de $CO_{2}$')
    axs5[8].set_ylabel('Pressão [bar]')


    # Plot 10: Média móvel da pressão pt05, pt07, pt08 e pt10
    axs5[9].plot(timeref, ma_pt05, label='PT-05')
    axs5[9].plot(timeref, ma_pt07, label='PT-07')
    axs5[9].plot(timeref, ma_pt08, label='PT-08')
    axs5[9].plot(timeref, ma_pt10, label='PT-10')
    axs5[9].set_ylim([0, 200])
    axs5[9].set_title('Pressão no trecho de testes - Média Móvel')
    axs5[9].set_ylabel('Pressão [bar]')

    # Plot 11: Vazao CO2 Média Móvel
    axs5[10].plot(timeref, ma_fit03, label='Gas')
    axs5[10].plot(timeref, ma_fit04, label='Líquido')
    # axs[10].set_ylim([0, 200]); # As before, usually CO2 flow is lower
    axs5[10].set_title('Vazão de $CO_{2}$ - Média Móvel')
    axs5[10].set_ylabel('Vazão [kg/h]')
    axs5[10].yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

    # Plot 12: Consumo co2 (bar plot)
    axs5[11].bar(labels_co2, plot_co2_values)
    axs5[11].set_title('(1) Consumo $CO_{2}$ gasoso, (2) $CO_{2}$ Líquido')
    axs5[11].set_ylabel('Consumo [kg]')
    axs5[11].grid(True, axis='y')
    axs5[11].yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))


    # %% Uso da transformada de Fourier para análise do espectro frequencias/amplitudes
    var_fft_series = pt07.copy()  # Variável a ser analisada (use .copy() if modifying)
    var_fft_series = var_fft_series.fillna(0) # FFT doesn't like NaNs
    var_fft_np = var_fft_series.values

    Fs = 1/5  # Frequência de aquisição Hz (1 sample every 5 seconds = 1/5 Hz)
    T = 1 / Fs
    L = len(var_fft_np)

    if L > 0:
        # Perform FFT
        fft_p = fft(var_fft_np)

        # Compute the single-sided amplitude spectrum
        P2 = np.abs(fft_p / L)
        P1_1 = P2[:L // 2 + 1] # Take the first half (positive frequencies)
        P1_1[1:-1] = 2 * P1_1[1:-1] # Double the amplitude (except DC and Nyquist)

        # Define the frequency domain f
        f = Fs * np.arange(L // 2 + 1) / L # Equivalent to MATLAB: Fs*(0:(L/2))/L

        fig6 = plt.figure(6, figsize=(10, 6))
        # MATLAB plots from index 2 to end-1, which is f[1:-1] in Python for 0-indexed arrays
        if len(f) > 2 and len(P1_1) > 2:
             plt.plot(f[1:-1], P1_1[1:-1])
        elif len(f) > 1 and len(P1_1) > 1: # Plot if at least DC and one freq point exists past it
            plt.plot(f[1:], P1_1[1:]) # Plot all except DC
        else: # Not enough points to skip DC or Nyquist
            plt.plot(f,P1_1)

        plt.title('Single-Sided Amplitude Spectrum of $P_{7}(t)$')
        plt.xlabel('f (Hz)')
        plt.ylabel('$|Y [kg]|$') # Original Y label, though PT07 is pressure [bar]
        # plt.ylim([0, 3]) # Uncomment if needed
        plt.grid(True)
    else:
        print("FFT not performed: var_fft_np is empty.")

    figures = (fig1,fig2,fig3,fig4,fig5, fig6)
    
    return figures



if __name__ == "__main__":
    print("File should not be run as main")

    

