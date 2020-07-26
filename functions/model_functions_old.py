import numpy as np
import pandas as pd
from datetime import datetime
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def run_SEIR_ODE_model(covid_parameters, model_parameters) -> pd.DataFrame:
    """
    Runs the simulation
    
    output:
        dataframe for SINGLE RUN
        dataframe list for SENSITIVITY ANALYSIS AND CONFIDENCE INTERVAL
    """
    cp = covid_parameters
    mp = model_parameters
    # Variaveis apresentadas em base diaria
    # A grid of time points (in days)
    t = range(mp.t_max)
        
    # CONDICOES INICIAIS
    # Initial conditions vector
    SEIRHUM_0 = initial_conditions(mp)
    
    niveis_isolamento = len(mp.contact_reduction_elderly)
    
    if mp.IC_analysis == 2:
    
        ii = 1
        df = pd.DataFrame()
        
        # 1: without; 2: vertical; 3: horizontal isolation 
        for i in range(niveis_isolamento): # 2: paper
            omega_i = mp.contact_reduction_elderly[i]
            omega_j = mp.contact_reduction_young[i]
        
            # Integrate the SEIR equations over the time grid, t
            # PARAMETROS PARA CALCULAR DERIVADAS
            args = args_assignment(cp, mp, omega_i, omega_j, i, ii)
            for a in range(4):
                t = range(5)
                if a == 0:
                    ret = odeint(derivSEIRHUM, SEIRHUM_0, t, args)
                    Si, Sj, Ei, Ej, Ii, Ij, Ri, Rj, Hi, Hj, dHi, dHj, Ui, Uj, dUi, dUj, Mi, Mj, pHi, pHj, pUi, pUj, pMi, pMj = ret.T
                else:
                    casa_negativa = -1
                    SEIRHUM_0 = Si[casa_negativa], Sj[casa_negativa], Ei[casa_negativa], Ej[casa_negativa], Ii[casa_negativa], Ij[casa_negativa], Ri[casa_negativa], Rj[casa_negativa], Hi[casa_negativa], Hj[casa_negativa], dHi[casa_negativa], dHj[casa_negativa], Ui[casa_negativa], Uj[casa_negativa], dUi[casa_negativa], dUj[casa_negativa], Mi[casa_negativa], Mj[casa_negativa], pHi[casa_negativa], pHj[casa_negativa], pUi[casa_negativa], pUj[casa_negativa], pMi[casa_negativa], pMj[casa_negativa]
                    ret = odeint(derivSEIRHUM, SEIRHUM_0, t, args)
                    print(type(ret))
                    Si, Sj, Ei, Ej, Ii, Ij, Ri, Rj, Hi, Hj, dHi, dHj, Ui, Uj, dUi, dUj, Mi, Mj, pHi, pHj, pUi, pUj, pMi, pMj = ret.T
                # Update the variables
                #Si, Sj, Ei, Ej, Ii, Ij, Ri, Rj, Hi, Hj, dHi, dHj, Ui, Uj, dUi, dUj, Mi, Mj, pHi, pHj, pUi, pUj, pMi, pMj = ret.T
                #plt.plot(t,pHi+pHj,Hi+Hj)
                #plt.show()
                df = df.append(pd.DataFrame({'Si': Si, 'Sj': Sj, 'Ei': Ei, 'Ej': Ej,
                                 'Ii': Ii, 'Ij': Ij, 'Ri': Ri, 'Rj': Rj,
                                 'Hi': Hi, 'Hj': Hj, 'dHi': dHi, 'dHj': dHj, 'Ui': Ui, 'Uj': Uj,
                                 'dUi': dUi, 'dUj': dUj, 'Mi': Mi, 'Mj': Mj,
                                 'pHi': pHi, 'pHj': pHj, 'pUi': pUi, 'pUj': pUj, 'pMi': pMi, 'pMj': pMj}, index=t)
                                    .assign(omega_i = omega_i)
                                    .assign(omega_j = omega_j))

        t = range(mp.t_max)
        DF_list = df
        plt.plot(df.Ii)
        plt.show()
    
    else:
        DF_list = list() # list of data frames
        
        runs = len(cp.alpha)
        print('Rodando ' + str(runs) + ' casos')
        print('Para ' + str(mp.t_max) + ' dias')
        print('Para cada um dos ' + str(niveis_isolamento) + ' niveis de isolamento de entrada')
        print('')
    
        for ii in range(runs): # sweeps the data frames list
            df = pd.DataFrame()
        
            # 1: without; 2: vertical; 3: horizontal isolation 
            for i in range(niveis_isolamento): # 2: paper
                omega_i = mp.contact_reduction_elderly[i]
                omega_j = mp.contact_reduction_young[i]
            
                # Integrate the SEIR equations over the time grid, t
                # PARAMETROS PARA CALCULAR DERIVADAS
                args = args_assignment(cp, mp, omega_i, omega_j, i, ii)
                ret = odeint(derivSEIRHUM, SEIRHUM_0, t, args)
                # Update the variables
                Si, Sj, Ei, Ej, Ii, Ij, Ri, Rj, Hi, Hj, dHi, dHj, Ui, Uj, dUi, dUj, Mi, Mj, pHi, pHj, pUi, pUj, pMi, pMj = ret.T
            
                df = df.append(pd.DataFrame({'Si': Si, 'Sj': Sj, 'Ei': Ei, 'Ej': Ej,
                                 'Ii': Ii, 'Ij': Ij, 'Ri': Ri, 'Rj': Rj,
                                 'Hi': Hi, 'Hj': Hj, 'dHi': dHi, 'dHj': dHj, 'Ui': Ui, 'Uj': Uj,
                                 'dUi': dUi, 'dUj': dUj, 'Mi': Mi, 'Mj': Mj,
                                 'pHi': pHi, 'pHj': pHj, 'pUi': pUi, 'pUj': pUj, 'pMi': pMi, 'pMj': pMj}, index=t)
                                    .assign(omega_i = omega_i)
                                    .assign(omega_j = omega_j))
            DF_list.append(df)
        
    return DF_list

def initial_conditions(mp):
    """
    Assembly of the initial conditions
    input: model_parameters (namedtuple)
    output: vector SEIRHUM_0 with the variables:
    Si0, Sj0, Ei0, Ej0, Ii0, Ij0, Ri0, Rj0, Hi0, Hj0, Ui0, Uj0, Mi0, Mj0
    Suscetible, Exposed, Infected, Removed, Ward Bed demand, ICU bed demand, Death
    i: elderly (idoso, 60+); j: young (jovem, 0-59 years)
    """    
    
    Ei0 = mp.init_exposed_elderly             # Ee0
    Ej0 = mp.init_exposed_young               # Ey0
    Ii0 = mp.init_infected_elderly            # Ie0
    Ij0 = mp.init_infected_young              # Iy0
    Ri0 = mp.init_removed_elderly             # Re0
    Rj0 = mp.init_removed_young               # Ry0
    Hi0 = mp.init_hospitalized_ward_elderly # He0
    Hj0 = mp.init_hospitalized_ward_young   # Hy0
    dHi0 = mp.init_hospitalized_ward_elderly_excess
    dHj0 = mp.init_hospitalized_ward_young_excess
    Ui0 = mp.init_hospitalized_icu_elderly  # Ue0
    Uj0 = mp.init_hospitalized_icu_young    # Uy0
    dUi0 = mp.init_hospitalized_icu_elderly_excess
    dUj0 = mp.init_hospitalized_icu_young_excess
    Mi0 = mp.init_deceased_elderly            # Me0
    Mj0 = mp.init_deceased_young            # My0
    pHi0 = 0
    pHj0 = 0
    pUi0 = 0
    pUj0 = 0
    pMi0 = 0
    pMj0 = 0

    # Suscetiveis
    Si0 = mp.population * mp.population_rate_elderly - Ii0 - Ri0 - Ei0  # Suscetiveis idosos
    Sj0 = mp.population * (1 - mp.population_rate_elderly) - Ij0 - Rj0 - Ej0 # Suscetiveis jovens
    
    SEIRHUM_0 = Si0, Sj0, Ei0, Ej0, Ii0, Ij0, Ri0, Rj0, Hi0, Hj0, dHi0, dHj0, Ui0, Uj0, dUi0, dUj0, Mi0, Mj0, pHi0, pHj0, pUi0, pUj0, pMi0, pMj0
    return SEIRHUM_0


def args_assignment(cp, mp, omega_i, omega_j, i, ii):
    """
    Assembly of the derivative parameters
    input: covid_parameters, model_parameters
    output: vector args with the variables:

    N, alpha, beta, gamma,
    los_leito, los_uti, tax_int_i, tax_int_j, tax_uti_i, tax_uti_j,
    taxa_mortalidade_i, taxa_mortalidade_j,
    omega_i, omega_j
    
    Population, incubation_rate, contact_rate, infectiviy_rate,
    average_length_of_stay (regular and icu beds), internation rates (regular and icu beds, by age)
    i: elderly (idoso, 60+); j: young (jovem, 0-59 years)
    mortality_rate for young and elderly
    """    
    
    N = mp.population
    pI = mp.population_rate_elderly
    Normalization_constant = mp.Normalization_constant[0] #Because if the constant be scaled after changing the contact matrix again,
    #it should lose the effect of reducing infection rate
    if mp.IC_analysis == 2: # SINGLE RUN
        alpha = cp.alpha
        beta = cp.beta
        gamma = cp.gamma
        delta = cp.delta
    else: # CONFIDENCE INTERVAL OR SENSITIVITY ANALYSIS
        alpha = cp.alpha[ii]
        beta = cp.beta[ii]
        gamma = cp.gamma[ii]
        delta = cp.delta[ii]
    
    contact_matrix = mp.contact_matrix[i]
    taxa_mortalidade_i = cp.mortality_rate_elderly
    taxa_mortalidade_j = cp.mortality_rate_young
    pH = cp.pH
    pU = cp.pU
    los_leito = cp.los_ward
    los_uti = cp.los_icu

    infection_to_hospitalization = cp.infection_to_hospitalization
    infection_to_icu = cp.infection_to_icu
    
    tax_int_i = cp.internation_rate_ward_elderly
    tax_int_j = cp.internation_rate_ward_young
    
    tax_uti_i = cp.internation_rate_icu_elderly
    tax_uti_j = cp.internation_rate_icu_young
    
    capacidade_UTIs = mp.bed_icu
    capacidade_Ward = mp.bed_ward

    args = (N, alpha, beta, gamma,delta,
            los_leito, los_uti, tax_int_i, tax_int_j, tax_uti_i, tax_uti_j,
            taxa_mortalidade_i, taxa_mortalidade_j,omega_i, omega_j,contact_matrix,pI,
            infection_to_hospitalization,infection_to_icu,capacidade_UTIs,capacidade_Ward,Normalization_constant,pH,pU)
    return args



def derivSEIRHUM(SEIRHUM, t, N, alpha, beta, gamma, delta,
                los_leito, los_uti, tax_int_i, tax_int_j, tax_uti_i, tax_uti_j,
                taxa_mortalidade_i, taxa_mortalidade_j,omega_i, omega_j,contact_matrix,pI,
                infection_to_hospitalization,infection_to_icu,capacidade_UTIs,capacidade_Ward,Normalization_constant,pH,pU):
    """
    Computes the derivatives

    input: SEIRHUM variables for elderly (i) and young (j), 
    Suscetible, Exposed, Infected, Recovered, Hospitalized, ICU, Deacesed
    time, Brazillian population,
    incubation rate, contamination rate, infectivity rate,
    LOS, hospitalization rates for wards and icu beds,
    death rates
    attenuating factors

    output: vector with the derivatives
    """    
    
    # Vetor variaveis incognitas
    Si, Sj, Ei, Ej, Ii, Ij, Ri, Rj, Hi, Hj, dHi, dHj, Ui, Uj, dUi, dUj, Mi, Mj, pHi, pHj, pUi, pUj, pMi, pMj = SEIRHUM
    
    Iij = np.array([[Ij/((1-pI)*N)],[Ii/(pI*N)]])
    Sij = np.array([[Sj],[Si]])
    dSijdt = -(beta/Normalization_constant)*np.dot(contact_matrix,Iij)*Sij
    dSjdt = dSijdt[0]
    dSidt = dSijdt[1]
    dEidt = - dSidt - alpha * Ei
    dEjdt = - dSjdt - alpha * Ej
    dIidt = alpha * Ei - gamma * Ii
    dIjdt = alpha * Ej - gamma * Ij
    dRidt = gamma * Ii
    dRjdt = gamma * Ij

    dpHi = -tax_int_i*dSidt - pHi / infection_to_hospitalization
    dpHj = -tax_int_j*dSjdt - pHj / infection_to_hospitalization
    dpUi = -tax_uti_i*dSidt - pUi / infection_to_icu
    dpUj = -tax_uti_j*dSjdt - pUj / infection_to_icu
    dpMi = -taxa_mortalidade_i*dSidt - pMi*delta
    dpMj = -taxa_mortalidade_j*dSjdt - pMj*delta

    coisa = 1/500
    coisa2 = -coisa*(Hi+Hj-capacidade_Ward)
    coisa = 1/50
    coisa3 = -coisa*(Ui+Uj-capacidade_UTIs)

    # Leitos demandados
    dHidt = (pHi / infection_to_hospitalization )*(1-1/(1+np.exp(coisa2))) - Hi / los_leito 
    dHjdt = (pHj / infection_to_hospitalization )*(1-1/(1+np.exp(coisa2))) - Hj / los_leito 

    dUidt = (pUi / infection_to_icu)*(1-1/(1+np.exp(coisa3))) - Ui / los_uti
    dUjdt = (pUj / infection_to_icu)*(1-1/(1+np.exp(coisa3))) - Uj / los_uti
    
    #Leitos demandados em excesso
    ddHidt = (pHi / infection_to_hospitalization)*(1/(1+np.exp(coisa2)))
    ddHjdt = (pHj / infection_to_hospitalization)*(1/(1+np.exp(coisa2)))

    ddUidt = (pUi / infection_to_icu)*(1/(1+np.exp(coisa3)))
    ddUjdt = (pUj / infection_to_icu)*(1/(1+np.exp(coisa3)))

    # Obitos
    dMidt = pMi * delta + ddHidt*pH + ddUidt*pU
    dMjdt = pMj * delta + ddHjdt*pH + ddUjdt*pU
    
    return (dSidt, dSjdt, dEidt, dEjdt, dIidt, dIjdt, dRidt, dRjdt,
            dHidt, dHjdt, ddHidt, ddHjdt, dUidt, dUjdt, ddUidt, ddUjdt, dMidt, dMjdt,
            dpHi, dpHj, dpUi, dpUj, dpMi, dpMj)
