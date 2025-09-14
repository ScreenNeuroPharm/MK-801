import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import os

colore = 'deeppink'
# colore = 'dodgerblue'


def crea_e_o_apri(nome_nuova_cartella):
    actual_folder = os.getcwd()
    try:
        os.makedirs(actual_folder + "\\" + nome_nuova_cartella)
    except FileExistsError:
        pass
    except FileNotFoundError:
        pass
    os.chdir(nome_nuova_cartella)


def ll4(x, n, val_min, val_max, IC50):
    return val_min + (val_max - val_min) / (1 + np.power((x / IC50), n))


basale = 1

gruppo = "Cx"
# gruppo = "Hp"
compartments = gruppo
grandezza = "MFR"

starting_folder = os.getcwd()

n_fasi = 9

basale_c_iniziale = 0.01

xData = np.array([basale_c_iniziale, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0])

xFit = np.arange(xData[basale], max(xData), 0.001)

crea_e_o_apri("Results")
crea_e_o_apri("Curve dose risposta")
crea_e_o_apri(gruppo)
crea_e_o_apri(grandezza)

if basale == 1:
    text = "without"
elif basale == 0:
    text = "with"
else:
    text = " "

filetxt = open(grandezza + "_" + gruppo + "_" + text + "_basal.txt", "w")
IC50_txt = open(grandezza + "_" + gruppo + "_" + text + "_basal_IC50_" + compartments + ".txt", "w")
R2_txt = open(grandezza + "_" + gruppo + "_" + text + "_basal_R2_" + compartments + ".txt", "w")
n_txt = open(grandezza + "_" + gruppo + "_" + text + "_basal_n_" + compartments + ".txt", "w")

filetxt.write("Matrix\tIC50\t\tHill coeff.\tCoeff. Det.")

data_folder = ""

os.chdir(starting_folder)
lista_starting_folder = os.listdir()
for folder in lista_starting_folder:
    if "DATI_TXT" in folder:
        os.chdir(folder)
        data_folder = os.getcwd()

os.chdir(data_folder)
lista_gruppi = os.listdir()

n_matrici = 0
all_data = []
lista_matrici = []

for i_gruppo in lista_gruppi:
    if gruppo == i_gruppo:

        os.chdir(gruppo)
        lista_grandezze = os.listdir()

        for i_grandezza in lista_grandezze:
            if grandezza == i_grandezza:

                os.chdir(grandezza)
                lista_matrici = os.listdir()
                n_matrici = len(lista_matrici)
                all_data = [[] for _ in range(n_matrici)]

                for i in range(n_matrici):
                    all_data[i] = np.double(np.loadtxt(fname=lista_matrici[i], delimiter=','))

for i in range(n_matrici):

    n_fasi = 9
    # data management
    DATI = all_data[i]
    yData = DATI[basale:n_fasi]
    n_basal = DATI[0]
    yData = np.array(yData) / n_basal

    n_fasi = len(yData) + 1

    # fitting
    max_data = max(yData)
    min_data = min(yData)

    m = (max_data + min_data) / 2

    popt, pcov = curve_fit(ll4, xData[basale:n_fasi], yData, maxfev=16500000)

    lim_inf_n = 1
    lim_inf_IC50 = xData[1 + basale]
    lim_sup_IC50 = xData[len(xData) - 2]

    if popt[1] > popt[2]:
        param_bounds = ([lim_inf_n, m, 0, lim_inf_IC50], [np.inf, np.inf, m, lim_sup_IC50])
        initial_guess = [lim_inf_n, max_data, min_data, 1]
        IC50_condition = False
        EC50_condition = True
        condition = "EC50"

    else:
        param_bounds = ([lim_inf_n, 0, m, lim_inf_IC50], [np.inf, m, np.inf, lim_sup_IC50])
        initial_guess = [lim_inf_n, min_data, max_data, 1]
        IC50_condition = True
        EC50_condition = False
        condition = "IC50"

    popt, pcov = curve_fit(ll4, xData[basale:n_fasi], yData,
                           p0=initial_guess, bounds=param_bounds, maxfev=165000)

    y_fit = ll4(xData[basale:n_fasi], *popt)

    ss_res = np.sum((yData - y_fit) ** 2)
    ss_tot = np.sum((yData - np.mean(yData)) ** 2)
    r2 = 1 - (ss_res / ss_tot)

    # saving of results
    nome_matrice = lista_matrici[i].split("_")
    nome_matrice = nome_matrice[-1].split(".")
    nome_matrice = nome_matrice[0]

    filetxt.write("\n" + nome_matrice
                  + "\t\t" + str(round(popt[3], 2))
                  + "\t\t" + str(round(popt[0], 2))
                  + "\t\t" + str(r2))

    cifre_decimali = 4

    if i == n_matrici - 1:
        IC50_txt.write(str(round(popt[3], cifre_decimali)))
        R2_txt.write(str(round(r2, cifre_decimali)))
        n_txt.write(str(round(popt[0], cifre_decimali)))
    else:
        IC50_txt.write(str(round(popt[3], cifre_decimali)) + ",")
        R2_txt.write(str(round(r2, cifre_decimali)) + ",")
        n_txt.write(str(round(popt[0], cifre_decimali)) + ",")

    # plot
    fig, ax = plt.subplots()
    plt.figure(i)
    ax.plot(xData[basale:n_fasi], yData, 'bo')
    ax.semilogx(xFit, ll4(xFit, *popt), colore, label=condition + '=%5.3f' % popt[3])
    tratteggio = 1
    ax.plot([xData[basale], xData[8]], [tratteggio, tratteggio], 'c--')

    if grandezza == 'MFR':
        grandezza_label = 'normalized MFR'
    else:
        grandezza_label = grandezza

    ax.set(xlabel='BIC [µM]', ylabel=grandezza_label)
    delta = (max_data - min_data) / 2
    ax.set_ylim([min_data - delta, max_data + delta])

    ax.legend(bbox_to_anchor=(0.625, 0.85, 0.4, 0.2), loc='upper right', ncol=2, mode="expand", borderaxespad=2)

    fig.suptitle(gruppo + ': matrix # ' + nome_matrice + " - " + gruppo, size=15)

    os.chdir(starting_folder)
    crea_e_o_apri("Results")
    crea_e_o_apri("Curve dose risposta")
    crea_e_o_apri(gruppo)
    crea_e_o_apri(grandezza)

    # fig.savefig(grandezza+'_' + gruppo + '_' + text + '_basal_' + nome_matrice + '_' + compartments + '.tif', dpi=600)
    fig.savefig(grandezza + '_' + gruppo + '_' + text + '_basal_' + nome_matrice + '_' + compartments + '.png')

filetxt.close()
IC50_txt.close()
R2_txt.close()
n_txt.close()
