import random
import pickle
import inspect
from scipy.optimize import minimize
from IPython.display import display, Math
import math
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import sinter
import matplotlib.scale as scale
import scipy.stats
import numpy as np
import warnings
import sys
import scipy.optimize as optimize


def objective_function(point, lines):
    """Objective function to minimize - sum of squared distances to lines."""
    return sum(distance_point_to_line(point, line)**2 for line in lines)

def distance_point_to_line(point, line):
    """Compute the distance from a point to a line."""
    a, b, c = line
    x, y = point
    return np.abs(a * x + b * y + c) / np.sqrt(a ** 2 + b ** 2)

def sinterplotthreshold(ax,mylist,order,rot,mem,num_rounds,find_pL_per,mind=0,maxd=100,maxp=1,minp=0.001,plot_inset_values = False):

    # sinter.plot_error_rate: https://github.com/quantumlib/Stim/wiki/Sinter-v1.12-Python-API-Reference#sinter.plot_error_rate:


    colors = ['cornflowerblue', '#ff7f0e', '#d62728','#2ca02c', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', 'steelblue','#673AB7', '#FF0000', '#FFA07A','c','m','gold','magenta']

    if mem == 'x' or mem == 'X':
        mem1 = 'x'
        mem2 = 'X'
    elif mem == 'z' or mem == 'Z':
        mem1 = 'z'
        mem2 = 'Z'

    # Num_rounds atm are d, 3d, 10d. So to convert this to a number need to first extract the coefficient:
    if num_rounds == 'd':
        rounds_coefficient = 1
    else:
        rounds_coefficient = int(num_rounds.strip('d'))

    # Number index to be used for plot colours:
    ds = np.arange(mind,maxd,2)
    number = len(ds)

    ax.grid(True, which='both',zorder =0)

    sinter.plot_error_rate(   
        ax=ax,
        stats=mylist,
        failure_units_per_shot_func=lambda stats: stats.json_metadata['d']*rounds_coefficient if find_pL_per == 'round' else rounds_coefficient if find_pL_per == 'd rounds' else 1, # If it's finding pL per shot then the failure_units_per_shot is just 1 :P 
        group_func=lambda stat: f"d={stat.json_metadata['d']}",
        x_func=lambda stat: stat.json_metadata['p'],
        
        filter_func=lambda s: 
            s.json_metadata['p'] <= maxp and
            s.json_metadata['p'] >= minp and

            

            (not (0.0056<=s.json_metadata['p']<=0.0058) if not plot_inset_values else True) and
            (not (0.0050<s.json_metadata['p']<0.0055) if not plot_inset_values else True) and
            (s.json_metadata['b']==mem1 or s.json_metadata['b']==mem2) and

            (rot in (s.json_metadata.get('ro'), s.json_metadata.get('rt'))) and
            (s.decoder == 'pymatching' or s.decoder == 'uncorrelat') and
            s.json_metadata['d']<=maxd and 
            s.json_metadata['d']>=mind and
            s.json_metadata['r']==num_rounds and
            # (s.json_metadata.get('noise', None) == 'SD' if 'noise' in s.json_metadata else True) #Can specify a noise model if need be, though I have separated my noise models to be different csv files.
            (s.json_metadata.get('idl',None)=='y' if 'idl' in s.json_metadata else True) and

            ## just excluding the points for high distance and low p which don't have enough samples:
            (s.json_metadata['p'] >= 0.0007 if (s.json_metadata['d']==15 and s.json_metadata['ro']=='unro') else s.json_metadata['p']<=maxp) and
            (s.json_metadata['p'] >= 0.001 if (s.json_metadata['d']>=16 and s.json_metadata['ro']=='unro') else s.json_metadata['p']<=maxp) and

            (s.json_metadata['p'] >= 0.001 if (s.json_metadata['d']>=19 and s.json_metadata['ro']=='ro') else s.json_metadata['p']<=maxp) and

            # CXSI noise 
            # (s.json_metadata['p'] >= 0.001 if (s.json_metadata['d']>=17 and s.json_metadata['ro']=='ro' and s.json_metadata['noise']=='CXSI') else s.json_metadata['p']<=maxp) and
            # (s.json_metadata['p'] >= 0.0007 if (s.json_metadata['d']==15 and s.json_metadata['ro']=='ro' and s.json_metadata['noise']=='CXSI') else s.json_metadata['p']<=maxp) and
            
            # (s.json_metadata['p'] >= 0.0007 if (s.json_metadata['d']==13 and s.json_metadata['ro']=='unro' and s.json_metadata['noise']=='CXSI') else s.json_metadata['p']<=maxp) and
            # (s.json_metadata['p'] >= 0.001 if (s.json_metadata['d']==15 and s.json_metadata['ro']=='unro' and s.json_metadata['noise']=='CXSI') else s.json_metadata['p']<=maxp) and
            # (s.json_metadata['p'] >= 0.002 if (s.json_metadata['d']==17 and s.json_metadata['ro']=='unro' and s.json_metadata['noise']=='CXSI') else s.json_metadata['p']<=maxp) and

            (not s.json_metadata['d'] == 19 if s.json_metadata.get('ro', s.json_metadata.get('rot', None)) == 'unro' else True) and
            s.json_metadata['o']==order,
    
        plot_args_func = lambda index, curve_id: {
        # 'color': ((1, 1-index/number, 0) if rot=='ro' else  (0,1-index/number, 1) if order == 21302130 else (0,1-index/number,0) ), 
        'color': colors[(index+2) % len(colors)] if rot == 'ro' else colors[index % len(colors)], # put in +2 to rotated because I'm plotting its distance from 8, whereas unro from 6
        'marker': 'x' if mem == 'x' else 'p'
        },
        # line_fits = ('log','log'), # need updated sinter for this
    )

    ax.loglog()
    # ax.set_title(f'$p_L$ vs $p$ ({rot}tated, mem {mem}, ${num_rounds}$ rounds, {order} order)')
    ax.set_title(f'{rot}tated, {order} order')
    # ax.set_title(f'Memory {mem}')
    ax.set_ylabel(f'Logical error rate ($p_L$) per $d$ rounds') if find_pL_per == 'd rounds' else ax.set_ylabel(f'Logical error rate ($p_L$) per {find_pL_per}')
    ax.set_xlabel(f'Physical Error Rate ($p$)')
    
    # ax.xaxis.set_major_locator(MaxNLocator(integer=True)) #forces integers on xaxis
    # ax.legend(loc='lower right')

    # ax.set_xlim(0.0008,0.01)
    # ax.set_ylim(1e-13,1e0)

    # Save to file and also open in a window.
    # fig.savefig(f'plots/thresholds/{rot},{mem},{num_rounds} rounds,o={order}.pdf')

    return




def bin_CNOT_orders(mylist):
    ro_orders = []
    unro_orders = []
    for stat in mylist:
        order = stat.json_metadata['o']
        rot = stat.json_metadata['ro']
        if rot == 'ro' and order not in ro_orders:
            ro_orders.append(order)
        if rot == 'unro' and order not in unro_orders:
            unro_orders.append(order)
    
    print("Rotated code CNOT orders:")
    for el in ro_orders:
        print(el)
    print("\nUnrotated code CNOT orders:")
    for el in unro_orders:
        print(el)

    return(ro_orders, unro_orders)



def plot_thresholds(mylist, roorder, unroorder, romind = 2, unromind = 2, minp = 0, maxp = 1,  output_dir='plots/thresholds',romaxd = 1000, unromaxd = 1000, find_pL_per = 'd rounds', ylims = [None, None], plot_inset_values = False, ignore_minds = False):

    orders = [unroorder, roorder]
    maxds = [unromaxd, romaxd]
    rotations = ['unro', 'ro']


    for rot, order, maxd in zip(rotations, orders, maxds):
        fig, ax = plt.subplots(1, 1, figsize=(5.4, 4.8))
        ax.grid(zorder = 0)
        plt.tight_layout(pad=3.0)

        for mem in 'xz':
            plt.gca().set_prop_cycle(None)  # reset the colour cycle
            sinterplotthreshold_v2(ax,mylist,order,rot,mem,'3d',find_pL_per,romind, unromind, maxd, minp = minp, maxp = maxp,plot_inset_values = plot_inset_values, ignore_minds = ignore_minds)
            ax.set_ylim(ylims[0], ylims[1])

        # Get the handles and labels
        handles, labels = ax.get_legend_handles_labels()

        # Create a dictionary to store unique labels
        unique_labels = {}
        unique_handles = []

        # Filter out duplicate labels
        for handle, label in zip(handles, labels):
            if label not in unique_labels:
                unique_labels[label] = True  # Mark the label as seen
                unique_handles.append(handle)  # Add corresponding handle

        # first_legend = ax.legend(handles[:len(handles)//2], labels[:len(handles)//2], title='Distance', markerscale=0, loc='lower right')
        # ax.add_artist(first_legend)
        
        first_legend = ax.legend(unique_handles, unique_labels.keys(), title='Distance', markerscale=0, loc='lower right')
        ax.add_artist(first_legend)

        for line in first_legend.get_lines():
            line.set_linewidth(3)

        marker_styles = {'X': 'x', 'Z': 'p'}
        marker_legend_handles = [plt.scatter([], [], marker=marker_styles[B], color='black', label=f'{B}') for B in 'XZ']

        second_legend = ax.legend(handles=marker_legend_handles, title='Memory:', loc='upper left')
        ax.add_artist(second_legend)
        # ax.grid(True, which='both', zorder = 0)

        ax.set_title(f'Rotated') if rot == 'ro' else ax.set_title(f'Unrotated')
        

        fig.savefig(f'{output_dir}/{rot}{str(order)}_threshold_plot.pdf', format='pdf')



def sinterplotthreshold_v2(ax,mylist,order,rot,mem,num_rounds,find_pL_per,romind = 2, unromind = 2, maxd = 100,maxp=1, minp=0, plot_inset_values = False,ignore_minds = False):
    # v2 calculates what the start colour for each graph should be based on romind (rotated code, minimum distance) and unromind (unrotated code minimum distance) to make the same colours for the same distances

    # sinter.plot_error_rate: https://github.com/quantumlib/Stim/wiki/Sinter-v1.12-Python-API-Reference#sinter.plot_error_rate:


    # find minimum plotted d's to get colours matching between both graphs
    ro_ds = []
    unro_ds = []
    for stat in mylist:
        d = stat.json_metadata['d']
        rotation = stat.json_metadata['ro']
        if rotation == 'ro' and d not in ro_ds:
            ro_ds.append(d)
        if rotation == 'unro' and d not in unro_ds:
            unro_ds.append(d)
        
    if romind < min(ro_ds):
        romind = min(ro_ds)
    if unromind < min(unro_ds):
        unromind = min(unro_ds)


    if ignore_minds == True:
        romind = 2
        unromind = 2
        colors = ['rebeccapurple', 'royalblue', 'slategray','mediumseagreen','cornflowerblue', '#ff7f0e', '#d62728','#2ca02c', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', 'steelblue','#673AB7', '#FF0000', '#FFA07A','c','m','gold','magenta']
    else:
        colors = ['cornflowerblue', '#ff7f0e', '#d62728','#2ca02c', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', 'steelblue','#673AB7', '#FF0000', '#FFA07A','c','m','gold','magenta']

    diff = romind - unromind ## if romind starts 2 above unro (for example) then its colours should start two above too.


    if mem == 'x' or mem == 'X':
        mem1 = 'x'
        mem2 = 'X'
    elif mem == 'z' or mem == 'Z':
        mem1 = 'z'
        mem2 = 'Z'

    # Num_rounds atm are d, 3d, 10d. So to convert this to a number need to first extract the coefficient:
    if num_rounds == 'd':
        rounds_coefficient = 1
    else:
        rounds_coefficient = int(num_rounds.strip('d'))

    # Number index to be used for plot colours:
    if rot == 'ro':
        mind = romind

    else:
        mind = unromind

    ds = np.arange(mind,maxd,2)
    number = len(ds)
    
    ax.grid(True, which='both',zorder = 0)

    sinter.plot_error_rate(   
        ax=ax,
        stats=mylist,
        failure_units_per_shot_func=lambda stats: stats.json_metadata['d']*rounds_coefficient if find_pL_per == 'round' else rounds_coefficient if find_pL_per == 'd rounds' else 1, # If it's finding pL per shot then the failure_units_per_shot is just 1 :P 
        group_func=lambda stat: f"d={stat.json_metadata['d']}",
        x_func=lambda stat: stat.json_metadata['p'],
        
        filter_func=lambda s: 
            s.json_metadata['p'] <= maxp and
            s.json_metadata['p'] >= minp and

            

            (not (0.0056<=s.json_metadata['p']<=0.0058) if not plot_inset_values else True) and
            (not (0.0050<s.json_metadata['p']<0.0055) if not plot_inset_values else True) and
            
            (s.json_metadata['b']==mem1 or s.json_metadata['b']==mem2) and

            (rot in (s.json_metadata.get('ro'), s.json_metadata.get('rt'))) and
            (s.decoder == 'pymatching' or s.decoder == 'uncorrelat') and
            s.json_metadata['d']<=maxd and 
            s.json_metadata['d']>=mind and
            s.json_metadata['r']==num_rounds and
            (s.json_metadata.get('idl',None)=='y' if 'idl' in s.json_metadata else True) and

        #     ## just excluding the points for high distance and low p which don't have enough samples:
            (s.json_metadata['p'] >= 0.0007 if (s.json_metadata['d']==15 and s.json_metadata['ro']=='unro') else s.json_metadata['p']<=maxp) and
            (s.json_metadata['p'] >= 0.001 if (s.json_metadata['d']>=16 and s.json_metadata['ro']=='unro') else s.json_metadata['p']<=maxp) and
            (s.json_metadata['p'] >= 0.001 if (s.json_metadata['d']>=19 and s.json_metadata['ro']=='ro') else s.json_metadata['p']<=maxp) and
            (not s.json_metadata['d'] == 19 if s.json_metadata.get('ro', s.json_metadata.get('rot', None)) == 'unro' else True) and
            
            # Finally the CNOT order should be the one specified:
            s.json_metadata['o'] == order,
    
        plot_args_func = lambda index, curve_id: {   # for each curve this function plots it gives it an index. The first curve will have the lowest index.

            'color': colors[(index + diff) % len(colors)] if rot == 'ro' else colors[index % len(colors)], # put in + diff to rotated because I'm usually plotting its distance from 8, whereas unro from 6

            'marker': 'x' if mem == 'x' else 'p'
        },
        
        # line_fits = ('log','log'), # need updated sinter for this
    )

    ax.loglog()
    ax.set_title(f'{rot}tated, {order} order')
    ax.set_ylabel(f'Logical error rate ($p_L$) per $d$ rounds') if find_pL_per == 'd rounds' else ax.set_ylabel(f'Logical error rate ($p_L$) per {find_pL_per}')
    ax.set_xlabel(f'Physical Error Rate ($p$)')

    return






def find_closest_point(lines):
    """Find the point closest to multiple lines."""
    # Initial guess for the point (you can modify this depending on your problem)
    initial_guess = [0.0, 0.0]

    # Optimize the objective function
    result = minimize(objective_function, initial_guess, args=(lines,))
    
    # Extract the optimized point
    closest_point = result.x
    
    return closest_point




def RMSE_BME(n, k): # using the analytical solution for RMSE(p_BME) for X~Bin(n,p) with n samples k hits
    MSE = (k + 1)*(n - k + 1)/ ((n + 3) * ( n + 2)**2)
    RMSE = np.sqrt(MSE)
    return RMSE

def RMSE_MLE(n, k): # using the analytical solution for RMSE(p_MLE) for X~Bin(n,p) with n samples k hits
    MSE = (6*k**2 - n * k * (k + 6) + n**2 * (k + 2)) / (n**2 * (n + 2) * (n + 3))
    RMSE = np.sqrt(MSE)    
    return RMSE

def linear_func(x, a, b):  
    return a * x + b

def calculate_pL(stats, find_pL_per = 'd rounds', pL_estimator = 'MLE', num_rounds = '3d'):

    d = stats.json_metadata['d']

    n = stats.shots
    k = stats.errors

    if pL_estimator == 'BME':
        pL = (k + 1) / (n + 2)  # Bayesian Mean Estimate for pL per shot
        RMSE = RMSE_BME(n, k)
    else:
        pL = k / n # Maximum likelihood estimate for pL per shot
        RMSE = RMSE_MLE(n, k)
    
    
    if find_pL_per != 'shot': 
        if num_rounds == 'd':
            rounds_multiple = 1
            number_of_rounds = d                    
        else: # i.e. if num_rounds is 3d, 10d, etc.
            rounds_multiple = int(num_rounds.strip('d')) # number of rounds taken as multiple of d
            number_of_rounds = d * rounds_multiple 
    
    
        if find_pL_per == 'd rounds':
            n = rounds_multiple
        elif find_pL_per == 'round':
            num_pieces = number_of_rounds

        # Updater pL and RMSE(pL) given we are not finding it per shot but per find_pL_per
        pL = sinter.shot_error_rate_to_piece_error_rate(pL, pieces = n)
        RMSE = (RMSE / n) * (1 - 2 * pL)**((1 / n) - 1) # propagating error (derivation in OneNote)

    return pL, RMSE


def extract4footprint(memtype,noise,mind,maxd,unroorder,roorder,mylist,num_rounds,find_pL_per):

    gradients = []; intercepts = []; pL = []; distances = []; qubits = []; volumes=[];

    for rotation in ['ro','unro']:  # ensure this matches csv file
        
        # Binding orders to rotation (unro could be either OG or 2130)
        if rotation == 'unro':
            order = unroorder
            def q(d):
                return (2*d-1)**2
        else:
            order = roorder
            def q(d):
                return 2*d**2 - 1
        
        xs = []
        ys = []
        log_ys = []
        qs = []
        
        for stats in mylist:
            p = stats.json_metadata['p']
            b = stats.json_metadata['b']
            rt = stats.json_metadata['ro']
            d = stats.json_metadata['d']
            o = stats.json_metadata['o']
            r = stats.json_metadata['r']

            if p != noise or rt != rotation or o != order or b!=memtype or r!=num_rounds:
                continue

            # if noise == 0.006 and memtype == 'Z':
            #     if d > 32 and (d-31)%4 != 0:   # use this for mem Z (Z_d_rounds.csv) when p=0.006 (just collected every 4d after d=31)
            #         continue

            if d > maxd or d < mind:
                continue
            if not stats.errors:
                # print(f"Didn't see any errors for d={d}")
                continue

            xs.append(d)
            qs.append(q(d))


            # Append pL values to ys. Will either be pL per shot, round or d rounds:

            # Calculate pL per shot:
            
            per_shot = stats.errors / stats.shots

            if find_pL_per == 'shot':
                ys.append(per_shot)
                log_ys.append(np.log(per_shot))
            
            else:
                if num_rounds == 'd':
                    rounds_multiple = 1
                    number_of_rounds = d
                else: # i.e. if num_rounds is 3d, 10d, etc.
                    rounds_multiple = int(num_rounds.strip('d')) # number of rounds taken as multiple of d
                    number_of_rounds = d*rounds_multiple 
                
                if find_pL_per == 'd rounds':
                    per_d_rounds = sinter.shot_error_rate_to_piece_error_rate(per_shot, pieces=rounds_multiple)
                    ys.append(per_d_rounds)
                    log_ys.append(np.log(per_d_rounds))
                elif find_pL_per == 'round':
                    per_round = sinter.shot_error_rate_to_piece_error_rate(per_shot, pieces=number_of_rounds)
                    ys.append(per_round)
                    log_ys.append(np.log(per_round))


        if len(xs) == 0:
            # print(f"No circuits matched p={noise},b={memtype},rt={rotation},o={order} in input csv file.")
            sys.exit(f"No circuits matched p={noise},b={memtype},rt={rotation},o={order} in input csv file.")
            break


        pL.append(ys)
        distances.append(xs)
        qubits.append(qs)

        fit = scipy.stats.linregress(xs, log_ys)
        gradients.append(fit.slope)
        intercepts.append(fit.intercept)

    rtqubits = [[] for _ in range(len(qubits))]
    volumes = [[] for _ in range(len(qubits))]
    
    for i in range(len(qubits)):
        rtqubits[i]=np.sqrt(qubits[i])
        if num_rounds == 'd':
            volumes[i] = [a*b for a,b in zip(distances[i],qubits[i])] # multiplies qubits by number of rounds (ran 'distance' rounds per distance)
        else:
            coefficient_string = num_rounds.strip('d') # If ran '3d' rounds this takes the number 3
            coefficient = int(coefficient_string)
            volumes[i] = [a*b for a,b in zip(coefficient * distances[i],qubits[i])] 


    # Put into ascending order of distance:
    for i in range(len(distances)): # rembember distances = [[ro data],[unro data]]. It is a list of 2 lists
        d_s, pL_s, q_s, rtq_s, v_s = map(list, zip(*sorted(zip(distances[i], pL[i], qubits[i],rtqubits[i],volumes[i]))))
        distances[i],pL[i],qubits[i],rtqubits[i],volumes[i] = d_s,pL_s,q_s,rtq_s,v_s

    return pL,distances,qubits,gradients,intercepts, rtqubits, volumes



def sinterplotpLvQ_overlaid(b,noise,mind,maxd,unroorder,roorder,mylist,num_rounds,find_pL_per,ax,sizeoffont,rotation,ps,colours,pltqubits = False):
    # Using: https://github.com/quantumlib/Stim/wiki/Sinter-v1.12-Python-API-Reference#sinter.plot_error_rate:

    if num_rounds == 'd':
        rounds_coefficient = 1
    else:
        rounds_coefficient = int(num_rounds.strip('d'))
    # ax.grid(zorder = 0)
    
    position = ps.index(noise)
    sinter.plot_error_rate(   
        ax=ax,
        stats=mylist,
        failure_units_per_shot_func=lambda stats: stats.json_metadata['d']*rounds_coefficient if find_pL_per == 'round' else rounds_coefficient if find_pL_per == 'd rounds' else 1, # If it's finding pL per shot then the failure_units_per_shot is just 1. If it's per round then it's coefficient*d. If it's per d rounds then it's just coefficient.
        
        # group_func=lambda stat: f"r={stat.json_metadata['ro']}, b={b}, o={stat.json_metadata['o']}",
        group_func=lambda stat: f"p={stat.json_metadata['p']}",
        # x_func=lambda stat: stat.json_metadata['d'],

        # x_func=lambda stat: np.sqrt((2*stat.json_metadata['d']**2)-1) if (stat.json_metadata['ro']=='ro' and pltqubits == True) else (2*stat.json_metadata['d']-1) if (stat.json_metadata['ro']=='unro' and pltqubits == True) else stat.json_metadata['d'],

        x_func=lambda stat: np.sqrt((2*stat.json_metadata['d']**2)-1) if (stat.json_metadata['ro']=='ro' and pltqubits == True) else (2*stat.json_metadata['d']-1) if (stat.json_metadata['ro']=='unro' and pltqubits == True) else stat.json_metadata['d'],

        
        filter_func=lambda s: 
        s.json_metadata['ro'] == rotation and
        s.json_metadata['b'] == b and
        s.json_metadata['p'] == noise and
        s.json_metadata['d']<=maxd and
        s.json_metadata['d']>=mind and
        s.json_metadata['idl'] == 'y' and
        (s.json_metadata['r'] == num_rounds or s.json_metadata['r'] == s.json_metadata['d']*rounds_coefficient) and
        ((s.json_metadata['ro']=='ro' and s.json_metadata['o']==roorder) or (s.json_metadata['ro']=='unro' and s.json_metadata['o']==unroorder)) 
        
        # # excluding some points without enough samples:
        and
        # (s.json_metadata['p'] > 0.001 if (s.json_metadata['d']>=16 and s.json_metadata['ro']=='unro') else s.json_metadata['p']>0) and
        (s.json_metadata['d'] <= 15 if (s.json_metadata['p']==0.0007 and s.json_metadata['ro']=='unro') else s.json_metadata['d']<=maxd) 
        and
        (s.json_metadata['d'] <= 18 if (s.json_metadata['p']==0.0007 and s.json_metadata['ro']=='ro') else s.json_metadata['d']<=maxd)
        

        ,

        plot_args_func = lambda index, curve_id: 
        {
            'color': colours[position % len(colours)],
            
            # (1-(np.log10(float(re.search(r'p=([0-9.]+)', curve_id).group(1)))+4)/2, 0,1) 
            
            # if 'unro' in rotation 
            
            # 'color': ('c' if rotation=='ro' else 'm'), 
            

            # 'color': ((1, 0.9-position/len(ps), 0) if rotation=='ro' else  (0,0.9-position/len(ps), 1)), 
            
            # else (1.34-(np.log10(float(re.search(r'p=([0-9.]+)', curve_id).group(1)))+4)/2,0, 0),

            
            'marker': '^' if rotation=='ro' else 'o',
            # 'marker': markers[position],
            
            },

        line_fits = ('linear','log'), # line_fits usually works with updated version of sinter but now isn't??
    )

    ax.semilogy()
    ax.tick_params(axis='both', which='major', labelsize=sizeoffont) 
    ax.set_ylabel(f'$p_L$ per {find_pL_per}' if find_pL_per != 'd rounds' else f'$p_L$ per $d$ rounds', fontsize=sizeoffont)  
    xlabel = 'Distance' if pltqubits == False else '$\sqrt{\mathrm{Total\ qubit\ count}}$'
    ax.set_xlabel(xlabel, fontsize=sizeoffont) 
    ax.legend(fontsize=sizeoffont)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True)) #forces integers on xaxis
    # xlim = 45 if pltqubits == True else 24


    # if mind != 2:

    #     maxd = 0
    #     mind = 40
    #     for stat in mylist:
    #         d = stat.json_metadata['d']
    #         mind = min(mind, d)
    #         maxd = max(maxd, d)


    if pltqubits == True:
        ax.set_xlim(np.sqrt(2 * mind ** 2 - 1), 1.3 * np.sqrt(4 * maxd ** 2 - 4 * maxd + 1))
        pass
    else:
        ax.set_xlim(0.4 * mind, 1.3 * maxd)
        pass
    ax.set_ylim(1e-12,1e0)
    return


def plot_pL_vs_qubit_count(mylist, b, roorder = 10231203, unroorder = 10231203, ps = None, romind = 2, unromind = 2):

    # Plot versus distance or qubits:
    pltqubits = True

    # # Load the noise model data

    # # with open(f'pickled_stats/{noise_model}_importedstats.pickle', 'rb') as file: # just odd distances
    # with open(f'pickled_stats/{noise_model}_combined_importedstats.pickle', 'rb') as file: # odd and even distances
    #     mylist = pickle.load(file)


    
    if ps == None:
        ps = []
        for stat in mylist:
            p = stat.json_metadata['p']
            d = stat.json_metadata['d']

            if p not in ps:
                ps.append(p)



    # Configuration for rotated plot

    num_rounds = '3d'
    find_pL_per = 'd rounds'
    fontsize = 10
    figure_size = (5.4, 4.8)

    fig_ro, ax_ro = plt.subplots(1,1,figsize=figure_size)
    ax_ro.set_axisbelow(True)
    ax_ro.grid(True, zorder = 0)
    rotation = 'ro'

    fig, ax = fig_ro, ax_ro


    # ps = [0.0005,0.007,0.001,0.0015, 0.002, 0.003, 0.004,0.005]
    colours = ['rebeccapurple', 'royalblue', 'magenta', 'slategray','mediumseagreen', 'indianred', 'gold', 'darkcyan', 'cornflowerblue'] * (len(ps) // 10 + 1) # Repeat colours if needed

    for rotation in ['ro','unro']:
        if rotation == 'ro':
            minds = [romind for _ in ps]
            # minds = [8 for _ in ps]
            maxds = [17,22] +  [22 for _ in range(len(ps)-2)]
        else:
            minds = [unromind for _ in ps]
            # minds = [6 for _ in ps]
            maxds = [13,17,17] + [17 for _ in range(len(ps)-3)]

        for noise, mind, maxd in zip(ps, minds, maxds):
            pL, _, _, _, _, _, _ = extract4footprint(b, noise, mind, maxd, unroorder, roorder, mylist, num_rounds, find_pL_per)
            log_pL = [math.log(x, 10) for x in pL[0]]  

            sinterplotpLvQ_overlaid(b, noise, mind, maxd, unroorder, roorder, mylist, num_rounds, find_pL_per, ax, fontsize, rotation,ps,colours,pltqubits)


    # Create custom legend handles for the lines
    legend_handles_lines = [plt.Line2D([0], [0], color=color, linewidth=2) for color in list(reversed(colours))[-len(ps):]]
    legend_labels_lines = [f'{p}' for p in reversed(ps)]

    # Add the legend for the lines
    line_legend = ax.legend(legend_handles_lines, legend_labels_lines, loc='upper right',title = '$p$')
    ax.add_artist(line_legend)

    # Create custom handles
    custom_handles = [
        plt.Line2D([], [], color='k', marker='o', linestyle='', label='Unrotated'),
        plt.Line2D([], [], color='k', marker='^', linestyle='', label='Rotated')
    ]

    # Create the legend
    first_legend = ax.legend(
        handles=custom_handles[::], 
        fontsize=fontsize, 
        loc='lower right',
        # bbox_to_anchor = (1,0.41),
        handlelength=1.1,  # 2 default
        handletextpad=0.8,  # Padding around text 0.8 default
        borderpad=0.4,  # Adjust the padding inside the legend border - default is 0.4
        labelspacing=0.5  # Adjust the spacing between labels
    )

    # Add the legend to the axes
    ax.add_artist(first_legend)


    ax.set_title(f'$p_L$ vs. Qubit Count' if pltqubits else f'$p_L$ vs. Distance')

    # custom_lines = [Line2D([0], [0], color='black', marker='o', linestyle='None') for i in range(len(ps))]
    # ax.legend(custom_lines, [f'${p}$' for p in ps[::-1]], loc='lower right',title = '$p$')

    fig_ro.savefig(f'plots/footprints/footprint_plot_roorder{roorder}_unroorder_{unroorder}_mem_{b}.pdf', format='pdf',transparent=True)
    plt.show()




def extract4teraquop(memtype,noise,mind,maxd,unroorder,roorder,mylist,num_rounds,find_pL_per, pL_estimator = 'MLE', plvsd = True,weighting_option = 1):


    pLs = []
    log_pLs = []
    uncertainties_log_pLs = []

    gradients = [] 
    gradient_errors = []
    intercepts = [] 
    intercept_errors = []
 
    distances = [] 
    qubits = [] 


    for rotation in ['ro','unro']:  # ensure this matches csv file
        
        # Binding orders to rotation (unro could be either OG or 2130)
        if rotation == 'unro':
            order = unroorder
            def q(d):
                return (2*d-1)**2
        else:
            order = roorder
            def q(d):
                return 2*d**2 - 1
        
        temp_ds = [] # distances for this rotation
        temp_qs = []
        
        temp_pLs = [] # pL_MLE estimates for this rotation
        temp_log_pLs = []
        temp_uncertainties_of_log_pLs = []
        
        
        

        for stats in mylist:
            p = stats.json_metadata['p']
            b = stats.json_metadata['b']
            rt = stats.json_metadata['ro']
            d = stats.json_metadata['d']
            o = stats.json_metadata['o']
            r = stats.json_metadata['r']

            if p != noise or rt != rotation or o != order or b!=memtype or r!=num_rounds:
                continue

            if d > maxd or d < mind:
                continue
            if not stats.errors:
                continue

            temp_ds.append(d)
            temp_qs.append(q(d))


            n = stats.shots
            k = stats.errors

            if pL_estimator == 'BME':
                pL = (k + 1) / (n + 2)  # Bayesian Mean Estimate for pL per shot
                RMSE = RMSE_BME(n, k)
            else:
                pL = k / n # Maximum likelihood estimate for pL per shot
                RMSE = RMSE_MLE(n, k)
            
            
            if find_pL_per != 'shot': 
                if num_rounds == 'd':
                    rounds_multiple = 1
                    number_of_rounds = d                    
                else: # i.e. if num_rounds is 3d, 10d, etc.
                    rounds_multiple = int(num_rounds.strip('d')) # number of rounds taken as multiple of d
                    number_of_rounds = d*rounds_multiple 
            
            
                if find_pL_per == 'd rounds':
                    n = rounds_multiple
                elif find_pL_per == 'round':
                    num_pieces = number_of_rounds

                # Updater pL and RMSE(pL) given we are not finding it per shot but per find_pL_per
                pL = sinter.shot_error_rate_to_piece_error_rate(pL, pieces = n)
                RMSE = (RMSE / n) * (1 - 2 * pL)**((1 / n) - 1) # propagating error (derivation in OneNote)

            temp_pLs.append(pL)
            temp_log_pLs.append(np.log10(pL))
            temp_uncertainties_of_log_pLs.append(RMSE / (np.log(10) * pL)) # for data y with error E_y, the uncertainty on some function u(y) is (du/dy)*E_y. In this case u is the log_10 function ⇒ du/dy = 1/(ln(10)) (see Uncertainty in OneNote for more)

        if len(temp_ds) == 0:
            sys.exit(f"No circuits matched p={noise},b={memtype},rt={rotation},o={order} in input csv file.")
            break

        # For loop above generated pL values vs. qubit counts for a particular rotation  
                
        # Fit a straight line through the log10(pLs) and ds using curve_fit.
        # curve_fit(function, xdata, ydata). Assumes ydata = f(xdata, *params) + eps
        # So we are fitting to log10(pL) = m*d + b (for a particular value of p)
        
        popt, pcov = optimize.curve_fit(linear_func, temp_ds, temp_log_pLs, sigma = temp_uncertainties_of_log_pLs, absolute_sigma = True) # sigma: This means the curve fit takes into account the sampling error on each point. Absolute sigma: you are providing absolute errors (in the same units as your data) to the curve fit, not relative errors (which are relative ratios). 
        
        m_fit, b_fit = popt


        # # Something broke here so I'm just commenting out the line fits seeing as I'm going to use line fits to the x-axis being qubits anyway (rather than to the x-axis being distance like it is here.)
        # m_fit, b_fit = None, None
        # pcov = np.array([[1.0, 0.5],  # Example values
        #          [0.5, 1.0]])


        # print(f"noise = {noise}, {rotation}, m = {m_fit}, b = {b_fit}")
        E_m = np.sqrt(pcov[0,0]) # error (standard deviation) in m_fit
        E_b = np.sqrt(pcov[1,1])

        # Append everyting we found to lists outside the loop:
        pLs.append(temp_pLs)
        log_pLs.append(temp_log_pLs)
        uncertainties_log_pLs.append(temp_uncertainties_of_log_pLs)
        
        distances.append(temp_ds)
        qubits.append(temp_qs)

        gradients.append(m_fit)  # this function does one physical error rate (noise) value at a time and returns it so the list gradients will just have two values: gradient for rotated at p = noise and gradient for unrotated at p = noise
        intercepts.append(b_fit)

        gradient_errors.append(E_m)
        intercept_errors.append(E_b)

    # Put into ascending order of distance for the two lists (rotated and unrotated) contained in distances, pL, log_pL, uncertainties_log_pL etc.
    for i in range(len(distances)): # rembember distances = [[ro data],[unro data]]. It is a list of 2 lists
        d_s, q_s, pL_s, log_pL_s, unc_log_pL_s = map(list, zip(*sorted(zip(distances[i], qubits[i], pLs[i], log_pLs[i], uncertainties_log_pLs[i]))))
        distances[i], qubits[i], pLs[i], log_pLs[i], uncertainties_log_pLs[i] = d_s, q_s, pL_s, log_pL_s, unc_log_pL_s

    return distances,qubits,gradients,intercepts, pLs, log_pLs, uncertainties_log_pLs, gradient_errors, intercept_errors




def sinterplotpLvD_forteraquop(b,noise,mind,maxd,unroorder,roorder,mylist,num_rounds,find_pL_per,ax,sizeoffont,noise_model,rotation,ps,colours,pltqubits = False):
    # Using: https://github.com/quantumlib/Stim/wiki/Sinter-v1.12-Python-API-Reference#sinter.plot_error_rate:

    if num_rounds == 'd':
        rounds_coefficient = 1
    else:
        rounds_coefficient = int(num_rounds.strip('d'))

    num_lines = len(ps)
    position = ps.index(noise)

    sinter.plot_error_rate(   
        ax=ax,
        stats=mylist,
        failure_units_per_shot_func=lambda stats: stats.json_metadata['d']*rounds_coefficient if find_pL_per == 'round' else rounds_coefficient if find_pL_per == 'd rounds' else 1, # If it's finding pL per shot then the failure_units_per_shot is just 1. If it's per round then it's coefficient*d. If it's per d rounds then it's just coefficient.
        
        # group_func=lambda stat: f"r={stat.json_metadata['ro']}, b={b}, o={stat.json_metadata['o']}",
        group_func=lambda stat: f"p={stat.json_metadata['p']}",
        # x_func=lambda stat: stat.json_metadata['d'],
        x_func=lambda stat: np.sqrt((2*stat.json_metadata['d']**2)-1) if (stat.json_metadata['ro']=='ro' and pltqubits == True) else (2*stat.json_metadata['d']-1) if (stat.json_metadata['ro']=='unro' and pltqubits == True) else stat.json_metadata['d'],
        
        filter_func=lambda s: 
        s.json_metadata['ro'] == rotation and
        s.json_metadata['b'] == b and
        s.json_metadata['p'] == noise and
        s.json_metadata['d']<=maxd and
        s.json_metadata['d']>=mind and
        s.json_metadata['idl'] == 'y' and
        (s.json_metadata['r'] == num_rounds or s.json_metadata['r'] == s.json_metadata['d']*rounds_coefficient) and
        ((s.json_metadata['ro']=='ro' and s.json_metadata['o']==roorder) or (s.json_metadata['ro']=='unro' and s.json_metadata['o']==unroorder)) 
        
        # # excluding some points without enough samples for now (12h38 13.5.24):
        # and
        # (s.json_metadata['p'] > 0.001 if (s.json_metadata['d']>=16 and s.json_metadata['ro']=='unro') else s.json_metadata['p']>0) and
        # (s.json_metadata['p'] >= 0.001 if (s.json_metadata['d']==15 and s.json_metadata['ro']=='unro') else s.json_metadata['p']>0)

        ,

        plot_args_func = lambda index, curve_id: 
        {
            'color': colours[position],
            # (1-(np.log10(float(re.search(r'p=([0-9.]+)', curve_id).group(1)))+4)/2, 0,1) 
            # if 'unro' in rotation 
            
            # else (1.34-(np.log10(float(re.search(r'p=([0-9.]+)', curve_id).group(1)))+4)/2,0, 0),

            
            'marker': 'o' if b=='z' else 'x',
            
            },

        line_fits = ('linear','log'), # line_fits  work with sinter 1.13.0
    )

    ax.semilogy()
    ax.tick_params(axis='both', which='major', labelsize=sizeoffont) 
    ax.set_ylabel(f'$p_L$ per {find_pL_per}' if find_pL_per != 'd rounds' else f'$p_L$ per $d$ rounds', fontsize=sizeoffont)  
    xlabel = 'Distance' if pltqubits == False else '$\sqrt{Qubit\ Count}$'
    ax.set_xlabel(xlabel, fontsize=sizeoffont) 
    ax.legend(fontsize=sizeoffont)
    ax.grid(zorder = 0)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True)) #forces integers on xaxis
    if pltqubits == True:
        ax.set_xlim(6,48)
        pass
    else:
        ax.set_xlim(4,24)
        pass
    ax.set_ylim(1e-12,1e0)
    return




def give_p_values(mylist,roorder, unroorder):
    
    pvalues = []

    for stat in mylist:
        if stat.json_metadata['ro'] == 'ro':
            if stat.json_metadata['o'] == roorder:
                p = stat.json_metadata['p']
                if p not in pvalues:
                    pvalues.append(p)
        elif stat.json_metadata['ro'] == 'unro':
            if stat.json_metadata['o'] == unroorder:
                p = stat.json_metadata['p']
                if p not in pvalues:
                    pvalues.append(p)

    return sorted(pvalues)





def extract_data_for_ratios(mylist, b, ps, roorder, unroorder, romind = 2, unromind = 2, noise_model = 'SD', find_pL_per = 'd rounds',):
    import com_funcs5 as funcs
    from importlib import reload
    reload(funcs)

    # Constants
    num_rounds = '3d'
    find_pL_per = 'd rounds'  # Choices are 'round', 'd rounds' or 'shot'
    rotations = ['unro', 'ro']
    orders = [unroorder, roorder]
    # ps = [0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0055]
    colors = ['b', 'r', 'g', 'salmon', 'chocolate', 'orange', 'purple']

    # ps = [0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0055]
    colors = ['b', 'r', 'magenta', 'lightseagreen', 'chocolate', 'orange', 'purple']

    is_projection_below = np.full(len(ps), -4)
    
    
    is_projection_below = np.array([-10, -9.5, -6, -4, -3, -2, -2, -1])


    is_projection_below = 10.0 ** is_projection_below

    # Load the stats
    # with open(f'pickled_stats/{noise_model}_combined_importedstats.pickle', 'rb') as file:
        # mylist = pickle.load(file)

    # Initialize lists to store data
    ro_gradients = []
    E_ro_gradients = []
    ro_intercepts = []
    E_ro_intercepts = []
    unro_gradients = []
    E_unro_gradients = []
    unro_intercepts = []
    E_unro_intercepts = []

    for rotation, order in zip(rotations, orders):
        num_lines = len(ps)

        if rotation == 'ro':
            roorder = order
            # minds = [8 for _ in ps]
            minds = [romind for _ in ps]
            maxds = [17, 18] + [22 for _ in range(len(ps) - 2)]
        else:
            unroorder = order
            # minds = [6 for _ in ps]
            minds = [unromind for _ in ps]
            maxds = [13, 15, 17] + [17 for _ in range(len(ps) - 3)]

        for i, (noise, mind, maxd) in enumerate(zip(ps, minds, maxds)):
            distances, qubits, gradients, intercepts, pLs, log_pLs, uncertainties_log_pLs, gradient_errors, intercept_errors = extract4teraquop(
                b, noise, mind, maxd, unroorder, roorder, mylist, num_rounds, find_pL_per, pL_estimator='MLE'
            )

            index = 0 if rotation == 'ro' else 1
            log_pL = log_pLs[index]

            # Calculate line fit using optimize.curve_fit
            params, pcov = optimize.curve_fit(linear_model, np.sqrt(qubits[index]), log_pL)
            gradient, intercept = params
            E_gradient = np.sqrt(pcov[0, 0])  # Standard error in gradient fit
            E_intercept = np.sqrt(pcov[1, 1])

            if rotation == 'ro':
                ro_gradients.append(gradient)
                ro_intercepts.append(intercept)
                E_ro_intercepts.append(E_intercept)
                E_ro_gradients.append(E_gradient)
            else:
                unro_gradients.append(gradient)
                unro_intercepts.append(intercept)
                E_unro_gradients.append(E_gradient)
                E_unro_intercepts.append(E_intercept)

    return {
        'ro_gradients': ro_gradients,
        'E_ro_gradients': E_ro_gradients,
        'ro_intercepts': ro_intercepts,
        'E_ro_intercepts': E_ro_intercepts,
        'unro_gradients': unro_gradients,
        'E_unro_gradients': E_unro_gradients,
        'unro_intercepts': unro_intercepts,
        'E_unro_intercepts': E_unro_intercepts,
        'ps': ps,
        'is_projection_below': is_projection_below,
        'colors': colors
    }

def linear_model(x, m, b):
    return m * x + b


def plot_ratio(mylist, b, roorder, unroorder, romind = 2, unromind = 2, ps = None, paper_ylims = False):

    if ps == None:
        ps = []
        for stat in mylist:
            p = stat.json_metadata['p']
            if p not in ps:
                ps.append(p)


    data = extract_data_for_ratios(mylist, b, ps, roorder, unroorder, romind, unromind)
    
    ro_gradients = data['ro_gradients']
    E_ro_gradients = data['E_ro_gradients']
    ro_intercepts = data['ro_intercepts']
    E_ro_intercepts = data['E_ro_intercepts']
    unro_gradients = data['unro_gradients']
    E_unro_gradients = data['E_unro_gradients']
    unro_intercepts = data['unro_intercepts']
    E_unro_intercepts = data['E_unro_intercepts']
    ps = data['ps']
    is_projection_below = data['is_projection_below']
    
    # Usual colours:
    some_ps = [0.0005,0.007,0.001,0.0015, 0.002, 0.003, 0.004,0.005,0.0055]
    some_colours  = ['rebeccapurple', 'royalblue', 'magenta', 'slategray','mediumseagreen', 'indianred', 'gold', 'darkcyan', 'cornflowerblue']

    colours_dict = dict(zip(some_ps, some_colours))


    # List to store the colors
    colors = []

    # List of potential random colors
    other_colours = ['mediumvioletred', 'darkcyan', 'gold', 'deeppink', 'royalblue', 
          'darkgoldenrod', 'indigo', 'darkmagenta', 'steelblue', 'saddlebrown']

    # Loop through pvalues in ps
    for p in ps:
        if p in colours_dict:
            colors.append(colours_dict[p])
        else:
            colors.append(random.choice(other_colours))  # Append a random color if pvalue not found





    def calc_uncertainty_in_ratio(pL, m_ro, b_ro, m_unro, b_unro, Em_ro, Em_unro, Eb_ro, Eb_unro):

        log_pL = np.log10(pL)
        
        a1 = (log_pL - b_ro) / m_ro
        a2 = (log_pL - b_unro) / m_unro
        
        y = (a1 / a2)**2
        
        # Partial derivatives
        dy_dm_ro = -2 * y / m_ro
        dy_dm_unro = 2 * y / m_unro

        dy_db_ro = - 2 * y / (log_pL - b_ro)
        dy_db_unro = 2 * y / (log_pL - b_unro)
        
        # Uncertainty in y
        Ey = np.sqrt((dy_dm_ro * Em_ro)**2 +
                    (dy_db_ro * Eb_ro)**2 +
                    (dy_dm_unro * Em_unro)**2 +
                    (dy_db_unro * Eb_unro)**2)
        
        return Ey


    # Divide the lines to find the ratio of qubit counts for a given pL:

    def ratio_function(x, m_ro, b_ro, m_unro, b_unro):
        return ((x-b_ro)/m_ro / ((x-b_unro)/m_unro))**2


    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    def ratio_function(pL, m_ro, b_ro, m_unro, b_unro):
        return ((np.log10(pL) - b_ro) / m_ro / ((np.log10(pL) - b_unro) / m_unro)) ** 2


    plt.figure(figsize=(5.4,4.8))

    minpL = 10**(-20)
    maxpL = 10**(-4)
    x = np.logspace(-20, -4, num = 400, base = 10) # these are pL values.

    limits = []
    Elimits = []

    # Plotting
    for i in reversed(range(len(ps))):  # can do len(ps) -1 here to not plot p = 0.0055 
        m_ro = ro_gradients[i]
        b_ro = ro_intercepts[i]

        m_unro = unro_gradients[i]
        b_unro = unro_intercepts[i]

        Em_ro = E_ro_gradients[i]
        Em_unro = E_unro_gradients[i]
        Eb_ro = E_ro_intercepts[i]
        Eb_unro = E_unro_intercepts[i]

        limit = (m_unro / m_ro)**2
        Emunrooverro = (m_unro / m_ro) * np.sqrt((Em_unro / m_unro)**2 + (Em_ro / m_ro)**2)
        Elimit = 2*(m_unro / m_ro) * Emunrooverro
        limits.append(limit)
        Elimits.append(Elimit)

        y = ratio_function(x, m_ro, b_ro, m_unro, b_unro)

        Ey = calc_uncertainty_in_ratio(x, m_ro, b_ro, m_unro, b_unro, Em_ro, Em_unro, Eb_ro, Eb_unro)

        # print(f"p = {ps[i]}, y = {y}, Ey = {Ey}")
        
        # for i in range(len(y)):
            # print(f"{y[i]} ± {Ey[i]}") # Ok .... those uncertainties are way too large lol.

        project_below = is_projection_below[i]
        x_dash = x[x < project_below]
        y_dash = y[x < project_below]
        x_solid = x[x >= project_below]
        y_solid = y[x >= project_below]

        # solid lines for ratios calculated from data:
        plt.plot(x_solid, y_solid, color=colors[i % len(colors)], linestyle='-', label=f'p = {ps[i]}', linewidth = 2)
        
        # dashed lines for ratios calculated from line-fit projections:
        plt.plot(x_dash, y_dash, color=colors[i % len(colors)], linestyle='--', linewidth = 2)

        # X marks the limit as pL -> 0 
        plt.scatter(minpL, limit, color=colors[i % len(colors)], marker='x', s = 100, zorder = 3)

        plt.fill_between(x, y + Ey, y - Ey , alpha = 0.3, color = colors[i % len(colors)]) # uncertainty

    x_ticks = np.arange(maxpL, minpL - 1, -4)
    plt.xticks(x_ticks)

    # def log_tick_formatter(val, pos):
        # return f'$10^{{{int(val)}}}$'

    # formatter = FuncFormatter(log_tick_formatter)
    # plt.gca().xaxis.set_major_formatter(formatter)

    legend = plt.legend(loc='upper right', fontsize=10)#, bbox_to_anchor=(0.03,0.99))
    for line in legend.get_lines():
        line.set_linewidth(2)  # Adjust the linewidth here

    plt.ylabel('$n_r\ /\ n_u$')
    plt.xlabel(f'$p_L$')
    plt.title(f'Qubit Count Ratio')
    # plt.grid()
    
    if paper_ylims == True:
        plt.ylim(0.68,0.93) # just making it the same as the paper for these p values
    
    plt.xlim(minpL, maxpL)
    plt.xscale('log')
    plt.grid()

    plt.savefig('plots/ratio_plot/ratio_plot.pdf', format='pdf')
    plt.show()





















def plot_memory_times(mylist, b, roorder, unroorder, ps = None, romind = 2, unromind = 2, plotagainst = 'rtqubits', plot_choice = 'Time'):

    # plot_choice = 'Time'  # 'Rounds' or 'Time'
    # plotagainst = 'rtqubits'  # 'rtqubits', 'qubits' or 'distance'


    if ps == None:
        ps = []
        for stat in mylist:
            p = stat.json_metadata['p']
            
            if p not in ps:
                ps.append(p)


    # Define colors and markers

    # Usual colours:
    some_ps = [0.0005,0.007,0.001,0.0015, 0.002, 0.003, 0.004,0.005,0.0055]
    some_colours  = ['rebeccapurple', 'royalblue', 'magenta', 'slategray','mediumseagreen', 'indianred', 'gold', 'darkcyan', 'cornflowerblue']
    colours_dict = dict(zip(some_ps, some_colours))
    # List to store the colors
    colours = []
    # List of potential random colors
    other_colours = ['mediumvioletred', 'darkcyan', 'gold', 'deeppink', 'royalblue', 
          'darkgoldenrod', 'indigo', 'darkmagenta', 'steelblue', 'saddlebrown']
    # Loop through pvalues in ps
    for p in ps:
        if p in colours_dict:
            colours.append(colours_dict[p])
        else:
            colours.append(random.choice(other_colours))  # Append a random color if pvalue not found


    noise_model = 'SD'

    time_per_round = 10**(-6)  # Characteristic times

    # with open(f'pickled_stats/{noise_model}_combined_importedstats.pickle', 'rb') as file:
    #     mylist = pickle.load(file)

    fig, ax = plt.subplots(figsize=(5.4, 4.8))

    # Define and max distances

    # ps = [0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0055]

    unromaxds = [17, 17, 17, 17, 17, 17, 17, 17]
    romaxds = [17, 22, 22, 22, 22, 22, 22, 22]

    for p, unromaxd, romaxd, color in zip(ps, unromaxds, romaxds, colours):
        for rot, order in zip(['ro', 'unro'], [roorder, unroorder]):
            if rot == 'ro':
                mind = romind # 8
                maxd = romaxd
            else:
                mind = unromind # 6
                maxd = unromaxd

            nrs = []
            qubitss = []
            memorytimes = []
            memorytimeuncertainties = []

            for stat in mylist:
                if stat.json_metadata['p'] != p: continue
                if stat.json_metadata['o'] != order: continue
                if stat.json_metadata['ro'] != rot: continue
                if stat.json_metadata['b'] != b: continue
                if stat.json_metadata['d'] < mind: continue
                if stat.json_metadata['d'] > maxd: continue
                if not stat.errors or not stat.shots: continue

                num_rounds = stat.json_metadata['r']
                num_rounds_as_number = int(num_rounds.strip('d')) * stat.json_metadata['d'] if num_rounds != 'd' else stat.json_metadata['d']
                
                n = stat.shots
                k = stat.errors
                per_shot = k / n
                RMSE = RMSE_MLE(n, k)

                n = num_rounds_as_number
                pL_per_round = sinter.shot_error_rate_to_piece_error_rate(per_shot, pieces = n)
                RMSE = (RMSE / n) * (1 - 2 * pL_per_round)**((1 / n) - 1) # propagating error (derivation in OneNote)

                p = stat.json_metadata['p']
                
                nr = p / pL_per_round
                nr_uncertainty = nr * (RMSE / pL_per_round)

                
                memorytime = time_per_round * nr
                memorytime_uncertainty = memorytime * (nr_uncertainty / nr)



                d = stat.json_metadata['d']
                if plotagainst == 'rtqubits':
                    qubits = np.sqrt(2 * d**2 - 1) if rot == 'ro' else (2 * d - 1)
                elif plotagainst == 'qubits':
                    qubits = 2 * d**2 - 1 if rot == 'ro' else (2 * d - 1)**2
                elif plotagainst == 'distance':
                    qubits = d

                memorytimes.append(memorytime)
                memorytimeuncertainties.append(memorytime_uncertainty)
                nrs.append(nr)
                qubitss.append(qubits)

            if nrs and qubitss and memorytimes and memorytimeuncertainties:
                nrs, qubitss, memorytimes, memorytimeuncertainties = zip(*sorted(zip(nrs, qubitss, memorytimes, memorytimeuncertainties)))
                yvalues = nrs if plot_choice == 'Rounds' else memorytimes if plot_choice == 'Time' else None

                lower_bounds = np.array(memorytimes) - np.array(memorytimeuncertainties)
                upper_bounds = np.array(memorytimes) + np.array(memorytimeuncertainties)

                marker = 'o' if rot == 'unro' else '^'
                ax.scatter(qubitss, yvalues, color=color, marker=marker, linestyle='-',zorder = 3)
                ax.fill_between(qubitss, lower_bounds, upper_bounds, color = color, alpha = 0.3, zorder = 2)


    ax.set_yscale('log')
    ax.set_axisbelow(True)
    ax.grid(zorder = 0)
    ax.set_title(f'Memory Time vs. Qubit Count' if plotagainst == 'rtqubits' else f'Memory Time vs. Distance' if plotagainst == 'distance' else f'')

    ax.set_xlabel('$\sqrt{\mathrm{Total\ qubit\ count}}$' if plotagainst == 'rtqubits' else 'Qubit Count' if plotagainst == 'qubits' else 'Distance' if plotagainst == 'distance' else '')
    ax.set_ylabel(f"{plot_choice}{' (s)' if plot_choice == 'Time' else ''} before $p_L = p$")

    # Create custom legend handles for the lines
    # legend_handles_lines = [plt.Line2D([0], [0], color=color, linewidth=2) for color in list(reversed(colours))[-len(ps):]]
    # legend_labels_lines = [f'p = {p}' for p in reversed(ps)]
    legend_handles_lines = [plt.Line2D([0], [0], color=color, linewidth=2) for color in colours[:len(ps)]]
    legend_labels_lines = [f'p = {p}' for p in ps]


    # Add the legend for the lines
    line_legend = ax.legend(legend_handles_lines, legend_labels_lines, loc='upper left', bbox_to_anchor = (0,0.85))
    ax.add_artist(line_legend)

    # Create custom handles for rotated and unrotated
    custom_handles = [
        plt.Line2D([], [], color='k', marker='^', linestyle='', label='Rotated'),
        plt.Line2D([], [], color='k', marker='o', linestyle='', label='Unrotated')
    ]

    # Create the legend
    first_legend = ax.legend(
        handles=custom_handles[::], 
        fontsize=10, 
        loc='upper left',
        handlelength=1.1,  # 2 default
        handletextpad=0.8,  # Padding around text 0.8 default
        borderpad=0.4,  # Adjust the padding inside the legend border - default is 0.4
        labelspacing=0.5  # Adjust the spacing between labels
    )

    # Add the legend to the axes
    ax.add_artist(first_legend)

    ax.yaxis.set_minor_locator(plt.NullLocator())

    plt.savefig(f'plots/memory_times/Memory_Times_{"X" if b == "x" else "Z"}_ro{roorder}_unro{unroorder}.pdf', format='pdf')
    plt.show()








def plot_teraquop(mylist, b, roorder, unroorder, ps = None,  noise_model = 'SD', romind = 2, unromind = 2, optional_plot = False, teraquop_inset = False):


# # Load the stats:
# if noise_model == 'SD':
#     ps = [0.0005,0.0007,0.001,0.0015, 0.002, 0.003,0.004,0.005,0.0055]
#     with open(f'pickled_stats/SD_combined_importedstats.pickle', 'rb') as file:
#         mylist = pickle.load(file)
# elif noise_model == 'CXSI':
#     ps = [0.0005,0.001, 0.002, 0.003,0.004,0.005]
#     with open(f'pickled_stats/CXSI_importedstats.pickle', 'rb') as file:
#         mylist = pickle.load(file)

        
    if ps == None:
        ps = []
        for stat in mylist:
            p = stat.json_metadata['p']
            
            if p not in ps:
                ps.append(p)


    

    def linear_model(x, m, b):
        return m * x + b



    # Specify which stats to plot:
    num_rounds = '3d'
    find_pL_per = 'd rounds'  # Choices are 'round', 'd rounds' or 'shot'
    rotations = ['unro_1203','ro']

    pL_estimator = 'MLE'



    # for optional plot of pL vs. d for each p (i.e. 'footprint' plots)
    if optional_plot:
        fig, axs = plt.subplots(1, len(rotations), figsize=(18, 5))


    if pL_estimator == 'MLE':
        estimator_display = r"\hat{p}_{\text{MLE}}"
    elif pL_estimator == 'BME':
        estimator_display = r"\hat{p}_{\text{BME}}"

    displayb = 'Z' if b == 'z' else 'X'

    # Display the LaTeX string
    display(Math(r"\text{Memory  }  " + f"{displayb} : {estimator_display}"))
    ## Plot the footprints (pL versus d) with line fits to find and save the teraquop distances:


    ## Set up arrays to store values:
    d_teras = [[] for _ in rotations] 
    d_teras_uncertainties = [[] for _ in rotations] 
    ceil_d_teras = [[] for _ in rotations] 
    q_teras = [[] for _ in rotations]
    Eq_teras = [[] for _ in rotations]
    q_teras_upper_bounds =  [[] for _ in rotations]
    q_teras_lower_bounds =  [[] for _ in rotations]
    ceil_d_to_q_teras = [[] for _ in rotations]
    ceil_d_to_q_teras_upper_bounds =  [[] for _ in rotations]
    ceil_d_to_q_teras_lower_bounds =  [[] for _ in rotations]

    fontsize = 10

    i = 0 
    for j in range(len(rotations)):
        
        if optional_plot:
            ax = axs[j]
        
        rotation = rotations[j]
        

        # index = 1
        # unroorder = 10231203 if noise_model == 'SD' else 12031203
        # if rotation == 'unro_1203':
        #     rotationtitle = 'Unrot. 1203'
        #     r = 'unro'
        #     unroorder = 10231203 if noise_model == 'SD' else 12031203
        # elif rotation == 'unro_2130':
        #     rotationtitle = 'Unrot. 2130'
        #     r = 'unro'
        #     unroorder = 23102130
        # else: 
        #     r = 'ro'
        #     rotationtitle = 'Rotated'
        #     index = 0

        if 'unro' in rotation:
            rotationtitle = 'Unrotated'
            r = 'unro'
        else:
            r = 'ro'
            rotationtitle = 'Rotated'


        num_lines = len(ps)

        if rotation == 'ro':
            # minds = [8 for _ in ps]
            minds = [romind for _ in ps]
            maxds = [17,18] + [22 for _ in range(len(ps)-2)]
            maxds = maxds[:len(ps)]
            color_start = np.array([255,165,0]) / 255 
            color_end = np.array([255,0,0]) / 255

        else:
            minds = [unromind for _ in ps]
            # minds = [6 for _ in ps]
            maxds = [13,15] + [17 for _ in range(len(ps)-2)]
            maxds = maxds[:len(ps)]
            color_end = np.array([112,8,112]) / 255      # purple
            color_start = np.array([254,1,154]) / 255    # neon pink

        colours = [color_start + (color_end - color_start) * (i / (num_lines - 1)) for i in range(num_lines)]
        
        # print(f"\n{rotation}")

        for noise, mind, maxd in zip(ps, minds, maxds):

            # GRADIENTS AND INTERCEPTS HERE ARE FITS FOR log10(p_L) = m * d + b 
            distances, qubits, gradients, intercepts, pL, log_pL, uncertainties_log_pL, gradient_errors, intercept_errors     =     extract4teraquop(b, noise, mind, maxd, unroorder, roorder, mylist, num_rounds, find_pL_per, pL_estimator, plvsd = True)

            # LET'S INSTEAD GET GRADIENTS AND INTERCEPTS FOR log10(p_L) = m * rt(q) + b WHERE Q IS QUBIT COUNT
            
            index = 0 if rotation == 'ro' else 1 # pL is a list of list. The first list (at index 0) is for rotated, the second is for unrotated:


            params, pcov = optimize.curve_fit(linear_model, np.sqrt(qubits[index]), log_pL[index])
            gradient, intercept = params
            E_gradient = np.sqrt(pcov[0,0]) # standard error in m_fit
            E_intercept = np.sqrt(pcov[1,1])

            ## Optional plot 
            if optional_plot:
                sinterplotpLvD_forteraquop(b, noise, mind, maxd, unroorder, roorder, mylist, num_rounds, find_pL_per, ax, fontsize, noise_model, r, ps, colours, pltqubits = True)


            # Fits based on qubit counts:
            b_fit = intercept
            E_b = E_intercept
            m_fit = gradient
            E_m = E_gradient

            yval = -12 # on a log10 plot this corresponds to a pL of 10^-12

            # d = (yval - b_fit) / m_fit # the teraquop distance
            rtq = (yval - b_fit) / m_fit # the teraquop qubit count

            q = rtq ** 2

            Ertq = rtq * np.sqrt(( E_b / (yval - b_fit) )**2 + (E_m / m_fit)**2)
            Eq = 2 * rtq * Ertq


            if noise == 0.001:
                print(rotation)
                print(q)
                print(Eq)



            if rotation == 'ro':
                d = np.sqrt((1 + q) / 2)
                Ed = 0.5 * Eq * 1/(np.sqrt(1 + q))
            else:
                d = (1 + rtq) / 2
                Ed = 0.5 * Ertq

            ceil_d = math.ceil(d) # Ceiling of the distance for teraquop: (rounding up to the next achievable distance)

            # convert ceil d to qubit count:
            if 'unro' in rotation:
                ceil_d_to_q = round((2 * ceil_d - 1) ** 2)
            else:
                ceil_d_to_q = round(2 * ceil_d ** 2 - 1)

            # Save what you found:
            d_teras[i].append(d)
            d_teras_uncertainties[i].append(Ed)
            
            ceil_d_teras[i].append(ceil_d)

            q_teras[i].append(q)
            Eq_teras[i].append(Eq)

            q_teras_upper_bounds[i].append(q + Eq)
            q_teras_lower_bounds[i].append(q - Eq)

            # if noise > 0: # 0.001:
                # print(f"p = {noise}")
                # print(f"d = {round(d,1)} ± {round(Ed,1)}")
                # print(f"q = {round(q,2)} ± {round(Eq,2)}")
                # print("ceil to next d:")
                # print(f"d = {math.ceil(d)} ± {round(Ed)}")
                # print(f"q = {ceil_d_to_q} ± {round(Eq)}")

            ceil_d_to_q_teras[i].append(ceil_d_to_q)
            ceil_d_to_q_teras_upper_bounds[i].append(math.ceil(ceil_d_to_q + Eq))
            ceil_d_to_q_teras_lower_bounds[i].append(math.ceil(ceil_d_to_q - Eq))

        if optional_plot:
            ax.set_title(rotationtitle) # Set title for each subplot
            handles, labels = ax.get_legend_handles_labels() # Get the legend handles and labels from the first subplot
            ax.legend(handles[::-1], labels[::-1], fontsize=fontsize) # Reverse the order of the handles and labels and set the legend for each subplot
            ax.grid(zorder=0) # Set the grid zorder to ensure it is beneath other plot elements
        
        i += 1



    # Plot required distance or qubit count to achieve teraquop regime (pL = 10^-12):
    plt.figure(figsize=(5.4,4.8))
    plt.grid(zorder=0) # Ensure the grid is beneath all other plot elements

    # loop over the rotations:
    for rotation,j in zip(rotations,range(len(rotations))):

        thelabel = 'Unrotated' if 'unro' in rotation else 'Rotated'
        colour = 'c' if rotation == 'ro' else 'm'

        plt.scatter(ps , list((q_teras[j])), marker='o' if 'unro' in rotation else '^', label = f'{thelabel}', color = colour, zorder=2)

        plt.fill_between(ps, q_teras_upper_bounds[j], q_teras_lower_bounds[j], alpha = 0.3, color = colour, zorder=1) # uncertainty

    B = 'Z' if b == 'z' else 'X'


    plt.title(f'Teraquop Qubit Count vs. $p$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Physical Error Rate ($p$)')
    plt.legend(loc = 'lower right')
    plt.ylabel('Teraquop qubit count')
    




    if teraquop_inset == True:

        ## NOW DO INSET GRAPH:

        ax_inset = plt.axes([0.22, 0.38, 0.45, 0.4]) # [left, bottom, width, height]
        figure_size = (5.4, 4.8)
        ax_inset.grid(zorder=0) # Ensure the grid is beneath all other plot elements




        for rotation,j in zip(rotations,range(0,len(rotations))):

            colour = 'c' if rotation == 'ro' else 'm'

            ax_inset.scatter(ps , list((ceil_d_to_q_teras[j])), linewidth=1, marker='o' if 'unro' in rotation else '^', s=10, color=colour, zorder=2)
            ax_inset.fill_between(ps, ceil_d_to_q_teras_lower_bounds[j], ceil_d_to_q_teras_upper_bounds[j], alpha=0.3, color=colour, zorder=1) # uncertainty

        ax_inset.set_xscale('log')
        ax_inset.set_yscale('log')
        ax_inset.tick_params(axis='both', which='major', labelsize=fontsize)
        ax_inset.tick_params(axis='both', which='minor', labelsize=fontsize)
        ax_inset.set_title(f'Rounded up to next $d$', fontsize=fontsize)
        plt.savefig(f'plots/teraquops/teraquop_plot_mem_{B}_ro{roorder}_unro{unroorder}.pdf', format='pdf')
        plt.show()

        # if optional_plot:
            # display(Math(r"\text{(Sinter footprint plots are always }\hat{p}_{\text{MLE}} \text{)}"))

        # display(Math(r"\text{Teraquop plot is memory  }  " + f"{displayb} : {estimator_display}"))











def fit_scaling_and_plot(combined_list, distances, basis, roorder, unroorder, minp = 0, maxp = 1, romind = 2, romaxd = 1000, unromind = 2, unromaxd = 1000, weighting_option = 2, optional_plot = True, plot_fits = True, ylims = [None, None], ignore_minds = False):

    # In paper used romind = 8, unromind = 6, (distances below this didn't fit scaling)
    # romaxd = 22, unromaxd = 17  (didn't have enough data for distances above these)

    # equation = r"\text{Fitting to:}\ \ \ p_L = \alpha \left( p/\beta \right)^{\gamma d\ -\ \delta}"
    # display(Math(equation))

    # weighting_option (for uncertainties)
    # option 1: weight each point based on its uncertainty in subsequent calculations

    # Other option, considering that using option 1 means your function fits will be dragged to fit low distances and low error rates because these are the points you have the most data for with the lowest uncertainty. This biases the fit to not be as relevant to higher distances and lower p values:
    # 2: fit lines as if all points have equal weight but still propagate the uncertainty through into the fit parameters, meaning they have larger uncertainty. This fits better to the points at lower p and with high d, points which we care about more, but the sacrifice is that it increases the uncertainty.



    # Split into even and odd distances to compare:
    my_even_list = [] # will contain whole element
    my_odd_list = []
    
    unro_combined_ds = []
    ro_combined_ds = []
    
    unro_odd_ds = []
    ro_odd_ds = []
    
    unro_even_ds = []
    ro_even_ds = []

    for el in combined_list:
        d = el.json_metadata['d']
        p = el.json_metadata['p']
        rot = el.json_metadata['ro']

        if rot == 'ro':
            if (d%2) == 0:
                my_even_list.append(el)
                if d not in ro_even_ds:
                    ro_even_ds.append(d)
            else:
                my_odd_list.append(el)
                if d not in ro_odd_ds:
                    ro_odd_ds.append(d)
            
            if d not in ro_combined_ds:
                ro_combined_ds.append(d)


        if rot == 'unro':
            if d != 19:  # don't have enough data for this distance
                
                if (d%2) == 0:
                    my_even_list.append(el)
                    if d not in unro_even_ds:
                        unro_even_ds.append(d)
                else:
                    my_odd_list.append(el)
                    if d not in unro_odd_ds:
                        unro_odd_ds.append(d)
                
                if d not in unro_combined_ds:
                    unro_combined_ds.append(d)



    print(f"\ndistances = {distances}")

    
    
    num_rounds = '3d'; 
    find_pL_per = 'd rounds'

    # colors = colors # different colour for every distance
    if ignore_minds == True:
        romind = 2
        unromind = 2
        colors = ['rebeccapurple', 'royalblue', 'slategray','mediumseagreen','cornflowerblue', '#ff7f0e', '#d62728','#2ca02c', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', 'steelblue','#673AB7', '#FF0000', '#FFA07A','c','m','gold','magenta']
    else:
        colors = ['cornflowerblue', '#ff7f0e', '#d62728','#2ca02c', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', 'steelblue','#673AB7', '#FF0000', '#FFA07A','c','m','gold','magenta']



    for rot, order, mem in zip(['ro','unro'],[roorder,unroorder],[basis,basis]):

        if distances == 'even':
            mylist = my_even_list

            if rot == 'ro':
                ds = ro_even_ds
            elif rot == 'unro':
                ds = unro_even_ds
        
        
        elif distances == 'odd':
            mylist = my_odd_list
            
            if rot == 'ro':
                ds = ro_odd_ds
            elif rot == 'unro':
                ds = unro_odd_ds
            
        elif distances == 'combined':
            mylist = combined_list
            
            if rot == 'ro':
                ds = ro_combined_ds
            elif rot == 'unro':
                ds = unro_combined_ds
        else:
            print("choose distances = 'odd', 'even' or 'combined'")


        # Order the ds (only makes a difference to plot legend)
        ds = np.array(ds)
        ds = np.sort(ds)

        # print(f"ds = {ds}")
        
        if rot == 'ro':

            mind = max(romind, min(ds))
            maxd = min(romaxd, max(ds))


        if rot == 'unro':
            mind = max(unromind, min(ds))
            maxd = min(unromaxd, max(ds))
        


        print(f"\n{rot} {mem} {order}\n")
        # print(f"mind = {mind}, maxd = {maxd}")

        # For a given d, I need to run over all the stats to fit a curve for log10(pL) = m * log10(p) + b

        lines = []
        Ems = []
        Ebs = []
        these_ds = []
        these_ds_ps = []

        # Below is based on extract4teraquop, but that fit to pL vs. d. Here we fit to pL vs. p
        index = 0 



        for d in ds: # find m and b for each distance's log10(pL) = m * log10(p) + b 

            if d < mind or d > maxd:
                continue
            
            these_ds.append(d)

            this_ds_pLs = []
            this_ds_log10_pLs = []

            this_ds_ps = []
            this_ds_log10_ps = []

            this_ds_pLs_RMSE = []
            
            this_ds_uncertainties_in_log10_pLs = []

            for stat in mylist:
                p = stat.json_metadata['p']
                b = stat.json_metadata['b']
                rt = stat.json_metadata['ro']
                stat_d = stat.json_metadata['d']
                o = stat.json_metadata['o']
                r = stat.json_metadata['r']

                if stat_d != d:
                    continue
                if not stat.errors:
                    continue
                if rt != rot or o != order or b!=mem or r!=num_rounds or p > maxp:
                    continue
                # excluding a couple points which don't have enough samples:
                if d == 19 and rt == 'ro' and p < 0.001: 
                    continue
                if d == 16 and rt == 'unro' and p < 0.001:
                    continue

                this_ds_ps.append(p)
                this_ds_log10_ps.append(np.log10(p))


                pL, RMSE = calculate_pL(stat)

                this_ds_pLs.append(pL)
                this_ds_pLs_RMSE.append(RMSE)

                this_ds_log10_pLs.append(np.log10(pL))
                this_ds_uncertainties_in_log10_pLs.append(RMSE / (np.log(10) * pL)) 

            # print(this_ds_ps)

            if len(this_ds_ps) == 0:
                line_number = inspect.currentframe().f_lineno
                # sys.exit(f"Line {line_number + 1}: no circuits matched d = {d}, b = {mem}, rt = {rot} , o = {order} in input csv file.")
                print(f"no circuits matched d = {d}, b = {mem}, rt = {rot} , o = {order} in input csv file.")
                continue
                # break

            # Now, for a given d, we have log10(pL) and log10(p). Let's fit to them

            xvalues = np.array(this_ds_log10_ps)
            yvalues = np.array(this_ds_log10_pLs)
            yvalues_uncertainties = np.array(this_ds_uncertainties_in_log10_pLs)
            
            if weighting_option == 1:
                popt, pcov = optimize.curve_fit(linear_func, xvalues, yvalues, sigma = this_ds_uncertainties_in_log10_pLs, absolute_sigma = True) 

                m, b = popt

                Em = np.sqrt(pcov[0,0]) # error (standard deviation) in m_fit
                Eb = np.sqrt(pcov[1,1])
            
            else:

                N = len(xvalues)

                # Step 1: Fit the line (Ordinary Least Squares)
                X = np.vstack((xvalues, np.ones(N))).T

                # Perform the least-squares fit to get m and b
                # (X^T * X)^(-1) * X^T * y
                coeffs, residuals, _, _ = np.linalg.lstsq(X, yvalues, rcond=None)
                m, b = coeffs

                # Step 2: Calculate the residuals
                residuals = yvalues - (m * xvalues + b)

                # Step 3: Estimate the variance of the residuals
                residual_variance = np.sum(residuals**2) / (N - 2)

                # Step 4: Include the known uncertainties in y-values
                effective_variance = residual_variance + np.mean(yvalues_uncertainties**2)

                # Step 5: Calculate the covariance matrix

                cov_matrix = effective_variance * np.linalg.inv(X.T @ X)

                # Step 6: Extract uncertainties in m and b
                Em = np.sqrt(cov_matrix[0, 0])
                Eb = np.sqrt(cov_matrix[1, 1])


            line = [] # will be of the form ax + by + c = 0. i.e. mx - y + b = 0

            line.append(m) 
            line.append(-1)
            line.append(b)

            lines.append(line)

            Ems.append(Em)
            Ebs.append(Eb)

            these_ds_ps.append(this_ds_ps)

            if optional_plot:

                diff = romind - unromind ## if romind starts 2 above unro (for example) then its colours should start two above too.
                color = colors[(index + diff) % len(colors)] if rot == 'ro' else colors[index % len(colors)], # put in + diff to rotated because I'm usually plotting its distance from 8, whereas unro from 6
                
                scatterplot = plt.scatter(this_ds_ps, this_ds_pLs, marker = '.', label = f"d = {d}", color = color,zorder = 3)
                color = scatterplot.get_facecolor()[0]

                # Plot the uncertainty (RMSE(p_MLE)):
                upper_log_pL = np.array(this_ds_pLs) + np.array(this_ds_pLs_RMSE)
                lower_log_pL = np.array(this_ds_pLs) - np.array(this_ds_pLs_RMSE)
                plt.fill_between(this_ds_ps, upper_log_pL, lower_log_pL, color = color, alpha = 0.3, zorder=2)
                index += 1


        # For each line (distance), we now have the m and b paramaters for the curve log10(pL) = m * log10(p) + b

        # Threshold is where these curves intersect (or, at least, the region around the closest point to all of them as we know they don't perfectly intersect)

        # for a given perp. dist., what is its error? 


        closest_point = find_closest_point(lines)
        x, y = closest_point # NOTICE THIS IS NOT LOG

        pth = 10 ** x
        pLth = 10 ** y


        # Calculate uncertainty:
        use_max_perp_dist = True

        if use_max_perp_dist == False:
            Ex = find_error_in_point(lines, Ems, Ebs, closest_point)
            Ey = Ex

            # print(f"Point and error before converting: \n {x} ± {Ex}\n")

            Epth = np.log(10) * (10 ** x) * Ex
            EpLth = np.log(10) * (10 ** y) * Ey

            # print(f"Using uncertainty in position of closest point to all lines:\n {pth} ± {Epth}\n") # this is uncertainty in the position of the closest point to all the lines. The uncertainty in this is very small due to small sampling errors. So if you want to report uncertainty on the closest point to all the lines you can.
        
        else:
            ## Alternatively the maximum perpindicular distance as uncertainty would encapsulate all the lines, which I would say is a better indication of where the threshold would be.
            largest_perp_dist = 0
            for line in lines:
                temp_dist = distance_point_to_line(closest_point, line) 
                
                if temp_dist > largest_perp_dist:
                    largest_perp_dist = temp_dist

            largest_perp_dist_on_log_graph = largest_perp_dist * (10 ** x) # when converting distances on a logarithmic scale you multiply by the scaling at the point you are at
            Epth = largest_perp_dist_on_log_graph


        ## NEXT STEP: find alpha.
        # We have each line in the form log10(pL) = m * log10(p) + b
        # This rearranges to pL = 10^b * p^m.
        # To make equivalent to pL = α (p/p_th)^m ⇒  α = 10^b * (p_th)^m

        # (If you just tried to make equivalent to pL = α p^m ⇒  α = 10^b, this would result in a different alpha for every distance.)

        # I'm just going to find each curve's α and then get a weighted mean, propagating the uncertainty through.

        alphas = []
        Ealphas = []

        for j, line in enumerate(lines):
            m = line[0]
            b = line[2]

            
            Em = Ems[j]
            Eb = Ebs[j]

            alpha = 10 ** b * pth ** m

            Etenb = 10 ** b * np.log(10) * Eb
            u = pth ** m
            dudp = m * (pth ** (m - 1))
            dudm = (pth ** m) * np.log(pth)
            Eu = np.sqrt((dudp * Epth) ** 2 + (dudm * Em) ** 2)
            
            Ealpha = alpha * np.sqrt((Etenb / (10 ** b)) ** 2 + (Eu / u) ** 2)
            alphas.append(alpha)
            Ealphas.append(Ealpha)

            # print(f"α = {alpha} ± {Ealpha}")

        alphas = np.array(alphas)
        Ealphas = np.array(Ealphas)

        weights = 1 / np.square(Ealphas)
        weighted_mean = np.sum(weights * alphas) / np.sum(weights)
        weighted_mean_error =  np.sqrt(1 / np.sum(weights))

        alpha = weighted_mean

        print(f"    α = {weighted_mean:.3f} ± {weighted_mean_error:.3f}" )
        print(f"    β = {pth:.5f} ± {Epth:.5f}")  # printing pth after α simply because it appears after in eqn.

        # NEXT STEP: find polynomial in d which describes m. m(d) = β * d + γ 

        ms = []
        for line in lines:
            ms.append(line[0])

        Ems = np.array(Ems)

        # # Using np.polyfit:
        # coefficients, coeff_errors = weighted_polyfit_with_errors(these_ds, ms, Ems, degree = 1)
        # gamma = coefficients[0]
        # delta = coefficients[1]
        # Egamma = coeff_errors[0]
        # Edelta = coeff_errors[1]
        
        if weighting_option == 1:
            # Using scipy optimize:
            popt, pcov = optimize.curve_fit(linear_func, these_ds, ms, sigma = Ems, absolute_sigma = True) 
            gamma, delta = popt
            Egamma = np.sqrt(pcov[0,0]) # error (standard deviation) in m_fit
            Edelta = np.sqrt(pcov[1,1])
        
        else:

            xvalues = np.array(these_ds)
            yvalues = np.array(ms)
            yvalues_uncertainties = np.array(Ems)
            
            N = len(xvalues)

            # Step 1: Fit the line (Ordinary Least Squares)
            X = np.vstack((xvalues, np.ones(N))).T

            # Perform the least-squares fit to get m and b
            # (X^T * X)^(-1) * X^T * y
            coeffs, residuals, _, _ = np.linalg.lstsq(X, yvalues, rcond=None)
            gamma, delta = coeffs

            # Step 2: Calculate the residuals
            residuals = yvalues - (gamma * xvalues + delta)

            # Step 3: Estimate the variance of the residuals
            residual_variance = np.sum(residuals**2) / (N - 2)

            # Step 4: Include the known uncertainties in y-values
            effective_variance = residual_variance + np.mean(yvalues_uncertainties**2)

            # Step 5: Calculate the covariance matrix
            cov_matrix = effective_variance * np.linalg.inv(X.T @ X)


            # Step 6: Extract uncertainties in gamma and delta
            Egamma = np.sqrt(cov_matrix[0, 0])
            Edelta = np.sqrt(cov_matrix[1, 1])

        
        print(f"    γ = {gamma:.3f} ± {Egamma:.3f}")
        print(f"    δ = {-delta:.2f} ± {Edelta:.2f}")

        # print(f"\n⇒ pL = {alpha:.3f}(p/{pth:.4f})^({gamma:.2f}d - {-delta:.2f})")


        if optional_plot: # plots the line fits as well


            plt.grid(zorder = 0)

            
            def pL(p,d):
                return alpha*(p/pth)**(gamma*d + delta)

            index = 0
            for j in range(len(these_ds)): # for d in list(range(6,30)): # if want to extrapolate
                
                # color = colors[(index+2) % len(colors)] if rot == 'ro' else colors[index % len(colors)]
                color = colors[index % len(colors)]
                
                d = these_ds[j]
                # ps = these_ds_ps[j]
                # ps = np.sort(ps)
                
                ps = np.linspace(minp, maxp, len(these_ds))
                
                if plot_fits == True:
                    linestyle = '--'
                    plt.plot(ps, pL(ps, d), linestyle = linestyle, linewidth = 1, color = 'grey', label = 'Curve fit' if j == 0 else None, zorder = 4)
            
                index += 1

            plt.xscale('log')
            plt.yscale('log')
            ax = plt.gca()

            ax.set_ylim(ylims[0],ylims[1])

            # ax.set_ylim(1e-12,2e-1)

            plt.legend(loc = 'lower right')
            title = "Rotated" if rot == 'ro' else "Unrotated"
            plt.title(title)
            plt.ylabel('$p_L$ per $d$ rounds')
            plt.xlabel('$p$')
            plt.savefig(f'plots/thresholds/{rot}{str(order)}_threshold_plot_WITH_FIT.pdf', format='pdf')
            plt.show()

        # # Plot fits to each individual line:
        # # plotlines(lines, maxd, mem) 
















