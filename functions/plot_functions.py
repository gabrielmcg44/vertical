# TO DO
# 1) Revisar fit
# 6) Incorporar plot para efeito de sobrecarga do sistema
# 7) Incorporar cenário alternativo com isolamento horizontal,
#     talvez travar o range de omega nos plots

import matplotlib.pyplot as plt
from .utils import *
from datetime import datetime
import numpy as np
import matplotlib as mpl

import pandas as pd
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter

def number_formatter(number, pos=None):
    """Convert a number into a human readable format."""
    magnitude = 0
    while abs(number) >= 1000:
        magnitude += 1
        number /= 1000.0
    return '%.1f%s' % (number, ['', 'K', 'M', 'B', 'T', 'Q'][magnitude])

def format_float(float_number, precision=2):
    return str(round(float_number,2)).replace(".", "")

def auxiliar_names(covid_parameters, model_parameters):
    """
    Provides filename with timestamp and IC_analysis type
    (2: Single Run, 1: Confidence Interval, 3: Sensitivity Analysis)
    as string for files
    """
    time = datetime.today()
    time = time.strftime('%Y%m%d%H%M')    

    if model_parameters.IC_analysis == 2: # SINGLE RUN

    	beta = covid_parameters.beta			# infectiviy_rate
    	gamma = covid_parameters.gamma			# contamination_rate

    	basic_reproduction_number = beta / gamma
    	r0 = basic_reproduction_number

    	filename = (time
            + '_single_run'
    		+ '_r' + ("%.1f" % r0)[0] + '_' + ("%.1f" % r0)[2]
    		+ '__g' + ("%.1f" % gamma)[0] + '_' + ("%.1f" % gamma)[2]
    		)
    elif model_parameters.IC_analysis == 1: # CONFIDENCE INTERVAL
        filename = (time + '_confidence_interval')
    elif model_parameters.IC_analysis == 3: # SENSITIVITY_ANALYSIS
    	filename = (time + '_sensitivity_analysis')
    else: #Rt analysis
    	filename = (time + '_Rt')
    return filename

def pos_format(title_fig,
			main_label_y,
			main_label_x,
			fsLabelTitle,leg_loc,fsPlotLegend):
    """
    põe labels e título, pós-formata fontes, eixo y, legenda
    """
    plt.title(title_fig, fontsize=fsLabelTitle)
    plt.legend(loc = leg_loc, fontsize=fsPlotLegend)
    plt.xlabel(main_label_x, fontsize=fsLabelTitle)
    plt.ylabel(main_label_y, fontsize=fsLabelTitle)		
    ax = plt.gca()
    ax.yaxis.set_major_formatter(FuncFormatter(number_formatter))

# TIPOS DE PLOT
# PLOT CONFIDENCE INTERVAL
def plot_ci(Y,cor,t_space):
    """
    plota intervalo dos percentis 5 e 95%
    """
    plt.fill_between(t_space,
		     np.quantile(Y, 0.05, axis = 0),
		     np.quantile(Y, 0.95, axis = 0).clip(Y[0,0]),
		     color = cor, alpha=0.2)
		
# PLOT CONFIDENCE INTERVAL
def plot_median(Y,cor,ls,line_label,t_space):
    """
    plota mediana, i.e. percentil 50%
    """
    plt.plot(t_space,
	     np.quantile(Y, 0.5, axis=0), # MEDIANA
	     ls,
	     color = cor,
	     label = line_label)

# PLOT TOTAL (i.e sem ser por idade)
def plot_total(Yi,Yj,name_variable,
			title_fig,
			fig_number,
			main_label_y,
			main_label_x,
			tamfig,fig_style,IC_analysis,t_space,
			ls,cor,isolation_name,i,omega_i,omega_j,
			fsLabelTitle,leg_loc,fsPlotLegend,
			plot_dir,filetype):
    """
    plota curvas totais (i.e sem ser por idade)
    com intervalo de confiança (5% percentil, mediana, 95% percentil)
    ou um único valor (SINGLE RUN)
    """
    plt.figure(fig_number, figsize = tamfig)
    plt.style.use(fig_style)
    if (IC_analysis == 2): # SINGLE RUN
    	  plt.plot(t_space,
		   (Yi+Yj),
		   ls[(i)%2], 		#0: dashed linestyle, 1: solid linestyle
		   color = cor[i],
		   label = 'Total' + isolation_name[i])
    else: # CONFIDENCE INTERVAL
        plot_median(Yi+Yj, cor[i], ls[(i)%2],'Total' + isolation_name[i], t_space)
        plot_ci(Yi+Yj, cor[i], t_space)
    		
    pos_format(title_fig, main_label_y,	main_label_x,
	       fsLabelTitle, leg_loc, fsPlotLegend)
    plt.savefig(os.path.join(plot_dir, name_variable + "_diff_isol" + filetype))

# PLOT POR IDADE - IDOSOS E JOVENS
def plot_byage(Yi,Yj,name_variable,
			title_fig,
			fig_number,
			main_label_y,
			main_label_x,
			tamfig,fig_style,IC_analysis,t_space,
			ls,cor,isolation_name,i,omega_i,omega_j,
			fsLabelTitle,leg_loc,fsPlotLegend,
			plot_dir,filetype):
    """
    plota curvas por idade
    com intervalo de confiança (5% percentil, 50% mediana, 95% percentil)
    ou um único valor (SINGLE RUN)
    """
    plt.figure(fig_number, figsize = tamfig)
    plt.style.use(fig_style)
    if (IC_analysis == 2): # SINGLE RUN
        plt.plot(t_space,
		 Yi,
		 ls[i%2],		#0: dashed linestyle, 1: solid linestyle
		 color = cor[2*i],
		 label = ('Elderly' + isolation_name[i]))
        plt.plot(t_space,
		 Yj,
		 ls[(i+1)%2],		#0: dashed linestyle, 1: solid linestyle
		 color = cor[1+2*i],
		 label = ('Young' + isolation_name[i]))
        complemento = isolation_name[i]
    else: # CONFIDENCE INTERVAL
        plot_median(Yi, cor[2*i], ls[i%2],'Elderly' + isolation_name[i], t_space)
        plot_ci(Yi,cor[2*i],t_space)
        plot_median(Yj, cor[1+2*i], ls[(i+1)%2],'Young' + isolation_name[i], t_space)
        plot_ci(Yj,cor[1+2*i],t_space)
        complemento = '_diff_isol'
    		
    pos_format(title_fig, main_label_y,	main_label_x,
		             fsLabelTitle,leg_loc,fsPlotLegend)
    plt.savefig(os.path.join(plot_dir, 
                         name_variable + "ij" + complemento + filetype))


def plots(results, covid_parameters, model_parameters, plot_dir):

	"""
	Makes plots:
        CONFIDENCE INTERVAL AND SINGLE RUN
	
	Figures numbers:

    1,4,7) Infected by age group for each degree of isolation
    
    2,5,8) Bed demand by age group for each degree of isolation
    
    3,6,9) Deceased by age group for each degree of isolation
    
    10) Infected for different degrees of isolation
    
    100) INFECTADOS - IDOSOS E JOVENS - DIFERENTES ISOLAMENTOS

    110) HOSPITALIZADOS TOTAL - DIFERENTES ISOLAMENTOS

    101) HOSPITALIZADOS UTI - IDOSOS E JOVENS - DIFERENTES ISOLAMENTOS
    
    23) Fit Infected
    
    24) Fit Deceased
    
        CONFIDENCE INTERVAL
    5% QUARTIL, MEDIANA, 95% QUARTIL
    
        SENSITIVITY ANALYSIS
    1) Infected people (r0)
    
    
    Degrees of isolation (i)
    no isolation, vertical, horizontal
        
    IC_Analysis
    1: Confidence Interval; 2: Single Run; 3: Sensitivity Analysis
    
    Age groups
    Elderly (60+); Young (0-59)
    
    Hospital Bed
    Ward; ICU
    
	"""
	
    
	mp = model_parameters    
	N = mp.population
    
	#capacidade_leitos = model_parameters.bed_ward
	#capacidade_UTIs = model_parameters.bed_icu
	IC_analysis = mp.IC_analysis

	t_max = mp.t_max
	t_space = np.arange(0, t_max)
	
# VARIÁVEIS PADRÃO DOS PLOTS
	fig_style = "ggplot" # "ggplot" # "classic" #
	# plot
	tamfig = (8,6)     # Figure Size
	fsLabelTitle = 11   # Font Size: Label and Title
	fsPlotLegend = 10   # Font Size: Plot and Legend

	main_label_x = 'Days'
	main_label_y = 'Infected people'

	cor = ['b','r','k','g','y','m']     # Line Colors
	ls = ['-.', '-']                    # Line Style
	leg_loc = 'upper left' # 'best' # 'upper right' # 
	filetype = '.pdf'      # '.pdf' # '.png' # '.svg' #
	
	isolation_name = [' (unmitigated)', ' (vertical isolation)', ' (horizontal isolation)']
	
	print('Plotando resultados')
    
	


	# 1: without; 2: vertical; 3: horizontal isolation 
	for i in range(len(mp.contact_reduction_elderly)): # 2: paper
		omega_i = mp.contact_reduction_elderly[i]
		omega_j = mp.contact_reduction_young[i]
        
		f_omega_i = format_float(omega_i, 1)
		f_omega_j = format_float(omega_j, 1)
       
		filename = f"we_{f_omega_i}_wy_{f_omega_j}"

		main_title = f'SEIR: $\omega_e={omega_i}$, $\omega_y={omega_j}$'
		
		
		if IC_analysis == 2: # SINGLE RUN
			Si = results.query('omega_i == @omega_i & omega_j == @omega_j')['Si']
			Sj = results.query('omega_i == @omega_i & omega_j == @omega_j')['Sj']
			Ei = results.query('omega_i == @omega_i & omega_j == @omega_j')['Ei']
			Ej = results.query('omega_i == @omega_i & omega_j == @omega_j')['Ej']
			Ii = results.query('omega_i == @omega_i & omega_j == @omega_j')['Ii']
			Ij = results.query('omega_i == @omega_i & omega_j == @omega_j')['Ij']
			Ri = results.query('omega_i == @omega_i & omega_j == @omega_j')['Ri']
			Rj = results.query('omega_i == @omega_i & omega_j == @omega_j')['Rj']
			Hi = results.query('omega_i == @omega_i & omega_j == @omega_j')['Hi']
			Hj = results.query('omega_i == @omega_i & omega_j == @omega_j')['Hj']
			Ui = results.query('omega_i == @omega_i & omega_j == @omega_j')['Ui']
			Uj = results.query('omega_i == @omega_i & omega_j == @omega_j')['Uj']
			Mi = results.query('omega_i == @omega_i & omega_j == @omega_j')['Mi']
			Mj = results.query('omega_i == @omega_i & omega_j == @omega_j')['Mj']
			
		else: # SENSITIVITY ANALYSIS OR CONFIDENCE INTERVAL
		
			Si = np.zeros((len(results),t_max))
			Sj = np.zeros((len(results),t_max))
			Ei = np.zeros((len(results),t_max))
			Ej = np.zeros((len(results),t_max))
			Ii = np.zeros((len(results),t_max))
			Ij = np.zeros((len(results),t_max))
			Ri = np.zeros((len(results),t_max))
			Rj = np.zeros((len(results),t_max))
			Hi = np.zeros((len(results),t_max))
			Hj = np.zeros((len(results),t_max))
			Ui = np.zeros((len(results),t_max))
			Uj = np.zeros((len(results),t_max))
			Mi = np.zeros((len(results),t_max))
			Mj = np.zeros((len(results),t_max))
			
			for ii in range(len(results)):
				Si[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Si']
				Sj[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Sj']
				Ei[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Ei']
				Ej[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Ej']
				Ii[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Ii']
				Ij[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Ij']
				Ri[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Ri']
				Rj[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Rj']
				Hi[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Hi']
				Hj[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Hj']
				Ui[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Ui']
				Uj[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Uj']
				Mi[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Mi']
				Mj[ii,] = results[ii].query('omega_i == @omega_i & omega_j == @omega_j')['Mj']
		
# SENSITIVITY ANALYSIS r0         
		if IC_analysis == 3:
			plt.figure(i, figsize = tamfig)
			plt.style.use(fig_style)
			
			r0 = covid_parameters.beta / covid_parameters.gamma
			r0min = r0[0]
			r0max = r0[len(results)-1]
         
			for ii in range(len(results)):
				a = (r0[ii] - r0min) / (r0max - r0min)
				plt.plot(t_space, Ii[ii, ] + Ij[ii, ],
						color = [a, 0, 1-a, 1],
                        label= f'$\omega_e={omega_i}$, $\omega_y={omega_j}$',
						linewidth = 0.5 )
		
			plt.title(main_title, fontsize=fsLabelTitle)
			plt.xlabel(main_label_x, fontsize=fsLabelTitle)
			plt.ylabel(main_label_y, fontsize=fsLabelTitle)
			mymap = mpl.colors.LinearSegmentedColormap.from_list(
				'mycolors',['blue','red'])
			sm = plt.cm.ScalarMappable(cmap=mymap, 
                                       norm=plt.Normalize(vmin=r0min, vmax=r0max))
			cbar = plt.colorbar(sm)
			cbar.set_label('Basic Reproduction Number',
				       rotation = 90, fontsize=fsLabelTitle)
			
			plt.savefig(os.path.join(plot_dir,
						 "I_" + filename + 'VariosR0' + filetype))
		
		else: # SINGLE RUN OR CONFIDENCE INTERVAL
# INFECTADOS TOTAL - DIFERENTES ISOLAMENTOS
			plot_total(Ii,Ij,'I',					#Yi,Yj,variable_name_Y
				   'Infected People by different isolation degrees',	#figure title
				   10,							#figure number
				   main_label_y,					#label_y
				   main_label_x,
				   tamfig,fig_style,IC_analysis,t_space,
				   ls,cor,isolation_name,i,omega_i,omega_j,
				   fsLabelTitle,leg_loc,fsPlotLegend,
				   plot_dir,filetype)

# INFECTADOS - IDOSOS E JOVENS - DIFERENTES ISOLAMENTOS
			plot_byage(Ii,Ij,'I',					#Yi,Yj,variable_name_Y
				   'Infected by age group for different isolation degrees',	#figure title
				   100,							#figure number
				   main_label_y,					#label_y
				   main_label_x,
				   tamfig,fig_style,IC_analysis,t_space,
				   ls,cor,isolation_name,i,omega_i,omega_j,
				   fsLabelTitle,leg_loc,fsPlotLegend,
				   plot_dir,filetype)

# HOSPITALIZADOS TOTAL - DIFERENTES ISOLAMENTOS
			plot_total(Hi,Hj,'H',					#Yi,Yj,variable_name_Y
				   'Hospitalized People by different isolation degrees',	#figure title
				   110,							#figure number
				   'Bed Demand',					#label_y
				   main_label_x,
				   tamfig,fig_style,IC_analysis,t_space,
				   ls,cor,isolation_name,i,omega_i,omega_j,
				   fsLabelTitle,leg_loc,fsPlotLegend,
				   plot_dir,filetype)
	 
# HOSPITALIZADOS UTI - IDOSOS E JOVENS - DIFERENTES ISOLAMENTOS
			plot_byage(Ui,Uj,'U',					#Yi,Yj,variable_name_Y
				   'ICU Bed Demand by age group for different isolation degrees',	#figure title
				   101,							#figure number
				   'ICU Bed Demand',					#label_y
				   main_label_x,
				   tamfig,fig_style,IC_analysis,t_space,
				   ls,cor,isolation_name,i,omega_i,omega_j,
				   fsLabelTitle,leg_loc,fsPlotLegend,
				   plot_dir,filetype)
		
# INFECTADOS - IDOSOS E JOVENS
			plot_byage(Ii,Ij,'I',					#Yi,Yj,variable_name_Y
				   'Infected by age group' + isolation_name[i],	#figure title
				   i,							#figure number
				   main_label_y,					#label_y
				   main_label_x,
				   tamfig,fig_style,IC_analysis,t_space,
				   ls,cor,isolation_name,i,omega_i,omega_j,
				   fsLabelTitle,leg_loc,fsPlotLegend,
				   plot_dir,filetype)
		
		
# LEITOS DEMANDADOS - IDOSOS E JOVENS
			plt.figure(3+i, figsize = tamfig)
			plt.style.use(fig_style)

			if (IC_analysis == 2): # Single run
				plt.plot(t_space, Hi, ls[0], color = cor[0], label = 'Ward for Elderly')
				plt.plot(t_space, Hj, ls[1], color = cor[1], label = 'Ward for Young')
				plt.plot(t_space, Ui, ls[1], color = cor[2], label = 'ICU for Elderly')
				plt.plot(t_space, Uj, ls[0], color = cor[3], label = 'ICU for Young')
				
			else: # Confidence interval
				plot_median(Hi, cor[0], ls[0], 'Ward for Elderly', t_space)
				plot_median(Hj, cor[1], ls[1], 'Ward for Young', t_space)
				plot_median(Ui, cor[2], ls[1], 'ICU for Elderly', t_space)
				plot_median(Uj, cor[3], ls[0], 'ICU for Young', t_space)

				plot_ci(Hi, cor[0], t_space)
				plot_ci(Hj, cor[1], t_space)
				plot_ci(Ui, cor[2], t_space)
				plot_ci(Uj, cor[3], t_space)
			pos_format(main_title, 'Bed Demand', main_label_x,
				   fsLabelTitle,leg_loc,fsPlotLegend)           
            
			plt.savefig(os.path.join(plot_dir, "HUey_" + filename  + filetype))
		
		
				# OBITOS - IDOSOS E JOVENS
			plt.figure(6+i, figsize = tamfig)
			plt.style.use(fig_style)
		
# OBITOS - IDOSOS E JOVENS
			plot_byage(Mi, Mj,'M',					#Yi,Yj,variable_name_Y
				   'Deacesed by age group' + isolation_name[i],	#figure title
				   6+i,						#figure number
				   'Deceased people',				#label_y
				   main_label_x,
				   tamfig,fig_style,IC_analysis,t_space,
				   ls,cor,isolation_name,i,omega_i,omega_j,
				   fsLabelTitle,leg_loc,fsPlotLegend,
				   plot_dir,filetype)         
            
            # ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
			plt.savefig(os.path.join(plot_dir, "Mey_" + filename  + filetype))
		
				

			
#		plt.figure(3+i)
#		plt.style.use('ggplot')

#		(results.query('omega_i == @omega_i & omega_j == @omega_j')
			# .div(1_000_000)
#			[['Hi', 'Hj', 'Ui', 'Uj']]
#			.plot(figsize=tamfig, fontsize=fsPlotLegend, logy=False)
#		)

#		plt.hlines(capacidade_leitos,
#					1, 
#					t_max, 
#					label=f'100% Leitos', colors='y', linestyles='dotted')

#		plt.hlines(capacidade_UTIs,
#					1, 
#					t_max, 
#					label=f'100% Leitos (UTI)', colors='y', linestyles='dashed')

#		plt.title(f'Demanda diaria de leitos: $N$={N}, ' + f'$\omega_i={omega_i}$, $\omega_j={omega_j}$', fontsize=fsLabelTitle)
#		plt.legend(['Leito normal idosos', 'Leito normal jovens', 'UTI idosos',
#				'UTI jovens', '100% Leitos', '100% UTIs'], fontsize=fsPlotLegend)
#		plt.xlabel('Dias', fontsize=fsLabelTitle)
#		plt.ylabel('Leitos', fontsize=fsLabelTitle)
#		plt.savefig(os.path.join(plot_dir, "HU_" + filename + ".png"))    







# FIT
			startdate = mp.startdate
			if ((IC_analysis == 4) and (not startdate == []) and (i == 0) ):
				
				for ifig in range(11):
				    plt.close(ifig)
                    
				dfMS = mp.dfMS
				state_name = mp.state_name
				data_sim = pd.to_datetime(t_space, unit='D',
				           origin=pd.Timestamp(startdate))
                # data_sim = pd.date_range(start=startdate, periods=t_max, freq='D')
                
				r0 = mp.r0_fit
				sub_report = mp.sub_report
				
				x_ticks = 14 # de quantos em quantos dias aparece tick no plot
                
                
				E0 = mp.init_exposed_elderly + mp.init_exposed_young
				I0 = mp.init_infected_elderly + mp.init_infected_young
				R0 = mp.init_removed_elderly + mp.init_removed_young
				S0 = N - E0 - I0 - R0
                
                
# OBITOS - IDOSOS E JOVENS
				plt.figure(23, figsize = tamfig)
				plt.style.use(fig_style)
		
				plot_median(Mi + Mj, cor[0], ls[0], 'Model', data_sim)

				plt.plot(dfMS['data'], dfMS['obitosAcumulado'],
    				     ls[1], color = cor[1], label = 'Reported Data')
					
				plot_ci(Mi+Mj, cor[0], data_sim)
				
				ax = plt.gca()
				ax.xaxis.set_major_locator(mdates.DayLocator(interval=x_ticks))
				ax.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
				#ax.set_xlim(pd.Timestamp('15/03/2020'), pd.Timestamp('15/04/2020'))
                
				plt.gcf().autofmt_xdate() # Rotation
				
				#ax.set_ylim(0, 1_000)
                
				title_fit = ('r0 = (' + ("%.1f" % r0[0])
					     + ', ' + ("%.1f" % r0[1]) + ') @' 
					     + pd.to_datetime((startdate), format='%Y-%m-%d')
					     .strftime('%d/%m/%Y')
					     + ' - subreport = ' + ("%d" % sub_report)
					     + ' - ' + state_name + '\n'
					     + 'SEIR(0) = ' 
					     + f"({S0:,.0f}; {E0:,.0f}; {I0:,.0f}; {R0:,.0f})")
         
				pos_format(title_fit, 'Deceased people', 'Date',
					   fsLabelTitle,leg_loc,fsPlotLegend)
		
				plt.savefig(os.path.join(plot_dir,
							 "Fit_" + state_name + "_M" + filetype))

# INFECTADOS - IDOSOS E JOVENS
				plt.figure(24, figsize = tamfig)
				plt.style.use(fig_style)
				
				plot_median(Ii+Ij, cor[0], ls[0], 'Model', data_sim)
				
				plt.plot(dfMS['data'], dfMS['casosAcumulado'] * sub_report,
					 ls[1], color = cor[1], label = 'Reported Data*sub_report')

				plot_ci(Ii+Ij, cor[0], data_sim)
				
				ax = plt.gca()
				ax.xaxis.set_major_locator(mdates.DayLocator(interval=x_ticks))
				ax.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
				#ax.set_xlim(pd.Timestamp('15/03/2020'), 
                #            pd.Timestamp('15/05/2020')) # 15/04/2020
                
				plt.gcf().autofmt_xdate() # Rotation
				
				#ax.set_ylim(0, 100_000) # 10_000
				pos_format(title_fit, 'Infected people', 'Date',
					   fsLabelTitle,leg_loc,fsPlotLegend)
				
				plt.savefig(os.path.join(plot_dir,
							 "Fit_" + state_name + "_I" + filetype))


	if ((IC_analysis == 4) and (not startdate == [])):
		for ifig in range(11):
			plt.close(ifig)
