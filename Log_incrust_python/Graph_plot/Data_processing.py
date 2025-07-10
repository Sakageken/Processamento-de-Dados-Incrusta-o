import pandas as pd
import numpy as np
from scipy.fft import fft # Using scipy's fft

def mean_std(filepath: str, sheetname: str):
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
    # If N is the number of rows in Table, length(Table(1:end-1,1)) is N-1.
    # The endpoint is t_exp * (N-1) / 60.
    # The number of points is (N-1) + 1 = N.
    num_rows = table_df.shape[0]
    end_time_minutes = t_exp * (num_rows - 1) / 60
    timeref = np.linspace(0, end_time_minutes, num_rows)
    
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

    # %% Medias e desvios
    t_i_idx = 0 # Default start index
    t_f_idx = len(pt05) -1 # Default end index

    # Ensure pt05 and tt05 are not empty and fit01 is not empty for std calculation
    if not pt05.empty and not tt05.empty and not fit01.empty:
        # Condicional
        if (pt05 >= 60).sum() >= 1:
            print('Using the pressure as intervaling time')
            pressure_threshold = pt05.mean() - pt05.std()/5
            print(f"Pressure thresholed = {pressure_threshold}")

            # Find t_i: first index where pt05 >= threshold
            valid_indices_ti_pt = np.where(pt05.values >= pressure_threshold)[0]

            if len(valid_indices_ti_pt) > 0:
                t_i_idx = valid_indices_ti_pt[0]
            else:
                t_i_idx = 0 # Default if no point meets criteria

            # Find t_f: last index where pt05 > threshold
            valid_indices_tf_pt = np.where(pt05.values > pressure_threshold)[0] 
            if len(valid_indices_tf_pt) > 0:
                t_f_idx = valid_indices_tf_pt[-1]
            else:
                t_f_idx = len(pt05) -1 # Default if no point meets criteria

        elif (pt05 < 60).sum() >= 1: # Original logic: exp > 1 OR (exp < 1 AND mean >= 20)
            print('Using the flow rate as intervaling time')
            #Derivada dos vetores de vazao fit 01 e 02
            fit01_dt = np.gradient(fit01,timeref)
            fit02_dt = np.gradient(fit02,timeref)

            #Condicional de primeira derivada ser diferente de zero
            condition_ti_01 = (fit01_dt != 0)
            condition_ti_02 = (fit02_dt != 0)

            valid_indices_01 = np.where(condition_ti_01)[0]
            valid_indices_02 = np.where(condition_ti_02)[0]

            if len(valid_indices_01) > 0 and len(valid_indices_02) > 0:
                if valid_indices_01[0] < valid_indices_02[0]:
                    t_i_idx = valid_indices_02[0] # As per MATLAB logic
                elif valid_indices_01[0] > valid_indices_02[0]:
                    t_i_idx = valid_indices_01[0]
                else: 
                    t_i_idx = 0# Boundary check

            if valid_indices_01[-1] > valid_indices_02[-1]: 
                valid_indices_tf = valid_indices_02
            else:
                valid_indices_tf = valid_indices_01

            if len(valid_indices_tf) > 0:
                t_f_idx = valid_indices_tf[-1] # As per MATLAB logic

    else:
        print("Skipping Medias e desvios calculation due to empty critical series.")
        # Set defaults to avoid errors later if series were empty
        t_i_idx = 0
        t_f_idx = len(fit01)
    print(f"Valor inicial de temp(Ti)={t_i_idx}, valor final de temp(Tf)={t_f_idx}")

    # Ensure slice indices are valid
    if t_i_idx > t_f_idx or t_i_idx >= len(fit01) or t_f_idx < 0 :
        print(f"Warning: Invalid slice indices t_i_idx={t_i_idx}, t_f_idx={t_f_idx}. Using full range for calculations.")
        slice_obj = slice(None) # Use full range
    else:
        slice_obj = slice(t_i_idx, t_f_idx + 1) # Python slicing is exclusive at the end

    # Helper function for safe mean/std calculation on slices
    def safe_stat(series, slice_obj, func_name):
        if series.empty: return np.nan
        sliced_data = series.iloc[slice_obj]
        if sliced_data.empty: return np.nan
        if func_name == 'mean': return sliced_data.mean()
        if func_name == 'std': return sliced_data.std()
        if func_name == 'sum': return sliced_data.sum()
        return np.nan

    # solucoes
    media_sol1 = safe_stat(fit01, slice_obj, 'mean')
    media_sol2 = safe_stat(fit02, slice_obj, 'mean')
    media_sol = np.nanmean([media_sol1, media_sol2]) if not (np.isnan(media_sol1) and np.isnan(media_sol2)) else np.nan
    desv_sol1 = safe_stat(fit01, slice_obj, 'std')
    desv_sol2 = safe_stat(fit02, slice_obj, 'std')
    desv_sol = np.nanmean([desv_sol1, desv_sol2]) if not (np.isnan(desv_sol1) and np.isnan(desv_sol2)) else np.nan

    # Temperatura na linha dos cupons, utilizando o termopar da linha antiga
    media_temp = safe_stat(tt05, slice_obj, 'mean')
    desv_temp = safe_stat(tt05, slice_obj, 'std')

    # CO2
    media_co2 = safe_stat(fit04, slice_obj, 'mean') # Assuming fit04 (liquid CO2) is used for this calc
    co2_por_sol = (media_co2 / (2 * media_sol)) * 100 if media_sol != 0 and not np.isnan(media_co2) and not np.isnan(media_sol) else np.nan

    # Pressao
    media_p_ln = safe_stat(pt08, slice_obj, 'mean') # pt08 for linha nova
    media_p_lv = safe_stat(pt05, slice_obj, 'mean') # pt05 for linha velha
    desv_p_ln = safe_stat(pt08, slice_obj, 'std')
    desv_p_lv = safe_stat(pt05, slice_obj, 'std')

    # Sol em cada linha
    media_sol_ln = safe_stat(fit06, slice_obj, 'mean') # fit06 (nova)
    media_sol_lv = safe_stat(fit05, slice_obj, 'mean') # fit05 (velha)
    desv_sol_ln = safe_stat(fit06, slice_obj, 'std')
    desv_sol_lv = safe_stat(fit05, slice_obj, 'std')


    sol_lv_total_mass = safe_stat(fit05, slice_obj, 'sum') * (t_exp / 3600.0)
    sol_ln_total_mass = safe_stat(fit06, slice_obj, 'sum') * (t_exp / 3600.0)

    # razao sol linha por sol das bombas
    sol_total_bicarbonato_mass = safe_stat(fit01, slice_obj, 'sum') * (t_exp / 3600.0)
    sol_total_cloreto_mass = safe_stat(fit02, slice_obj, 'sum') * (t_exp / 3600.0)
    sum_sol_bombas = sol_total_bicarbonato_mass + sol_total_cloreto_mass
    razao = ((sol_lv_total_mass + sol_ln_total_mass) / sum_sol_bombas) * 100 if sum_sol_bombas != 0 and not np.isnan(sum_sol_bombas) else np.nan


    # Matriz com os valores de medias e desvios
    table_data_np = np.array([
        [np.nan, np.nan],  # ph - not calculated here, placeholder
        [media_temp, desv_temp],
        [media_co2, np.nan], # desv_co2 not calculated
        [media_sol, desv_sol],
        [co2_por_sol, np.nan],
        [media_p_ln, desv_p_ln],
        [media_p_lv, desv_p_lv],
        [2500, np.nan],  # Campo mag - hardcoded
        [np.nan, np.nan],  # DeltaV - not calculated
        [media_sol_ln, desv_sol_ln],
        [media_sol_lv, desv_sol_lv],
        [sol_ln_total_mass, np.nan],
        [sol_lv_total_mass, np.nan],
        [razao, np.nan]
    ])

    # Label da matriz a ser plotada
    columns_py = ['Média', 'Desvio Padrão']
    rows_py = [
        'pH 1 [-]', 'Temperatura [ºC]', 'Vazão CO2 [kg/h]', 'Vazão soluções [kg/h]', 'Vazão CO2/Sol [%]',
        'Pressão linha nova [bar]', 'Pressão linha velha [bar]', 'Campo mag [Gauss]', 'DeltaV [%]',
        'Vazao solucao linha nova [kg/h]', 'Vazao solucao linha velha [kg/h]',
        'Massa total linha nova [kg]', 'Massa total linha velha [kg]', 'Sol_saída/Sol_bombas [%]'
    ]

    # Plot da figura ajustando o tamanho da tabela com a janela
    df_table = pd.DataFrame(table_data_np, columns=columns_py, index=rows_py)
    df_table = df_table.round(2)
    df_table = df_table.fillna("-")
    df_table.index.name = 'Teste'
    for cols_to_format in columns_py:
        df_table[cols_to_format] = df_table[cols_to_format].astype(str).str.replace('.',',',regex=False)
    #print("Dataframe", df_table)

    return df_table

if __name__ == "__main__":
    print("File should not be run as main")
 