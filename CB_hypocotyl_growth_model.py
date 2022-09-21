#!/usr/bin/env python
# coding: utf-8

# # <center> <font color='lightseagreen'> Circadian Clock Linked Hypocotyl Elongation Model </font>

# <center> <b> Chloe Bantock</b>  <br>
# Master in Bioinformatics and Computational Biology, UAM.<br>  
# Course 2021-2022 </center> 
# 

# ***

# 
# **Modelling the effects of circadian clock regulation of ELF3 and COP1 expression on hypocotyl elongation in *Arabidopsis thaliana***

# Recent advances in computational and experimental studies of the circadian clock in  Arabidopsis thaliana have led to new models of the clock gene network that include a three-component repressilator circuit (Pokhilko et al., 2012). This new circadian clock model of Pokhilko et al. (2012) suggests that the plant clock functions as an integrated, multi-feedback system that maintains robust oscillations and entrainment under multiple perturbations. <br>
# 
# The main objective of this master’s thesis will be to integrate the model developed by Pokhilko et al. (2013) with one developed by Cristina Nieto et al. (2021) for modeling light and temperature cues in hypocotyl elongation. This is possible as hypocotyl elongation is subject to photoperiodic regulation by the circadian clock through PIF4 and PIF5 (phytochrome interacting factors) transcription (Seaton et al., 2015). The integration of these two models would allow us to take a holistic approach to understanding the circadian clock through its interaction in output pathways. Further simulation and analysis of this coupled model would highlight the underlying molecular mechanisms that coordinate plant growth across changing conditions and may provide insights into the overall contribution of the clock to plant fitness. Future advancements in this area may in turn enable one to engineer aspects of plant physiology for improved growth in both existing and new environments (Seaton et al., 2015). 

# <font color="coral"><h3> **Code**</font>
# 
# The following sections include: 
#   
# 1.   [Replicating results from hypocotyl elongation model](#replicating)
#       -   Mutation Simulations
#       -   Protein Dynamics in Short Day 
#       -   Protein Dynamics in Long Day
# 2.   [Calling Matlab Code from Python](#calling_matlab) 
# 3.   [Clock Linked Growth Model (MATLAB)](#joint_matlab_model)
#       -   Comparison to original growth model 
#       -   Normalizing values of ELF3, COP1, and the EC 
# 4. [Original Clock Model](#original_clock_model)
# 5. [Final Linked Matlab Model](#final_linked_matlab_model) 
#     - Short Day Protein Dynamics
#     - Long Day Protein Dynamics
#     - Growth Curves 
#     - Mutation studies <br>
# 

# In[1]:


#imports

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as mticker
prop = fm.FontProperties(fname='/usr/share/fonts/truetype/Mark Simonson - Proxima Nova Alt Regular-webfont.ttf')
import pylab
import numpy as np
#!pip install seaborn
import seaborn as sns
from scipy.integrate import odeint
import itertools
import subprocess
import pandas as pd
import requests
import io
import openpyxl

import matlab.engine

from platform import python_version
print(python_version())
import sys
sys.executable


# <a id='replicating'></a>
# ## <font color='lightseagreen'> REPLICATING RESULTS FROM HYPOCOTYL ELONGATION MODEL </font>
# 
# Code from [Pablo Catalan](https://github.com/PabloCatalan/hypocotyl)

# In[2]:


# retrieve experimental data
day_length_data_url = 'https://raw.githubusercontent.com/PabloCatalan/hypocotyl/main/data/daylength_def.csv'
day_length_data = requests.get(day_length_data_url).content
DF = pd.read_csv(io.StringIO(day_length_data.decode('utf-8')))


# In[3]:


def get_key(df, key):
    '''
    given a dataframe and key (ie column name), returns a list containing 
    all values of the specified key. For example, calling get_key(DF, 'Genotype') would 
    return a list containing all Genotypes found in the experimental data
    
    :param a: dataframe 
    :param b: column name/key 
    
    :return: sorted list of values 
    '''
    kL=df[key].tolist()
    kM={}
    for k in kL:
        kM[k]=1
    kL=[]
    for k in kM:
        kL.append(k)
    return sorted(kL)

def read_data():
    '''
    returns the average growth and standard deviation in growth for each growth condition in 
    each mutant type in the experimental data (where the growth condition refers to the
    temperature and light conditions)
    '''
    #DF=pd.read_csv('data/daylength_def.csv')
    mut=get_key(DF,'Genotype')
    Length=get_key(DF,'Daylength')
    Temp=get_key(DF,'Temperature')    
    avgdata={}
    stddata={}
    for m in mut:
        DFm=DF[DF.Genotype==m]
        m=m.replace(' ','_')
        avgdata[m]={}
        stddata[m]={}
        # what does itertools.prodcut do???
        for d,t in itertools.product(Length,Temp):
            key=str(t)+'_'+str(d)
            DFnew=DFm[(DFm.Daylength==d) & (DFm.Temperature==t)]
            avgdata[m][key]=DFnew.Growth.mean()
            stddata[m][key]=DFnew.Growth.std()
    return avgdata, stddata

def is_day(t,Daylength):
    '''
    the light parameter can take one of two values – 0 at night, or 1 during daytime.
    Thus, the function returns 1 - soln to the heavyside function, which is given by: 
            t1 < 0,  = 0 (daytime)
            t1 == 0, = 0 (daytime)
            t1 > 0,  = 1 (nightime)
    '''
    t1=t%24-Daylength
    return 1-np.heaviside(t1,1)

# elf3 protein dynamics
def elf3p(t,Daylength,pE1,pE2):  
    k0=5
    t2=t-Daylength
    t3=t-24.0
    
    if Daylength==0:
        return pE1+pE2
    elif Daylength==24:
        return pE1-pE2    
    else:
        SigT=2.0/(1.0+np.exp(-k0*t))
        SigT2=2.0/(1.0+np.exp(-k0*t2))
        SigT3=2.0/(1.0+np.exp(-k0*t3))
        return pE1-pE2*(-1.0+SigT-SigT2+SigT3)

    
# where t a timepoint in parameter 'time' as defined when calling the odeint solver (120 hours, 500 timepoints)    
    
def growth(y, t, Temp, Day, mut, pB28, kr22, kr28, pE122, pE128,
           pE222, pE228, dE, pPE22, pPE28, dP, kPC, dPB, pCL28,
           pCD, dC, pG, kG, pGP, pGE, pGB, pGH, pHC, mutBox,
           mutEox, mutPox, mutPko1, mutPko2, mutCox, mutCko1,
           mutCko2):
    '''
    y is an array containing [PHYB], [ELF3], [PIF], [COP1], and growth (m) - protein concentrations are non-dimensional 
    parameters that are a function of T will have different values at 22C and 28C 
    For a given species K, dk = decay rate 
    mutk is a multiplier that alters molecule K's production to accomodate knout and over-expressor lines
    '''
        
    #Variables
    B=y[0]#PHYB
    E=y[1]#ELF3
    P=y[2]#PIF
    C=y[3]#COP1
    G=y[4]#Hypocotyl
    
    #Parameters
    # t1 is the hour in the day 
    t1=t%24
    L=is_day(t1,Day)
    pB=10.0
    kr=kr22#0.232 datos de Casal
    pE1=pE122#adimensional
    pE2=pE222
    pP=1.0#adimensional
    pPE=pPE22
    pCL=1.0#adimensional
    mB=1.0#maximum PHYB value
    if Temp==28:
        pB=pB28
        kr=kr28#0.411 datos Casal
        pE1=pE128
        pE2=pE228
        pPE=pPE28
        pCL=pCL28
    if 'PHYBox' in mut:
        mB*=mutBox
    if 'ELF3ox' in mut:
        pE1*=mutEox
    if 'PIF4ox' in mut:
        pP*=mutPox
    if 'pif4' in mut:
        pP*=mutPko1
    if 'pifq' in mut:
        pP*=mutPko2
    if 'COP1' in mut:
        pCL*=mutCox
        pCD*=mutCox
    if 'cop1-4' in mut:
        pCL*=mutCko1
        pCD*=mutCko1
    if 'cop1-6' in mut:
        pCL*=mutCko2
        pCD*=mutCko2
    if 'hy5' in mut:
        pGH=0
        
    #Equations (all are nuclear proteins)
    
    # active form of phyB, Pfr
    # where pB is rate of activation and translocation to nuclues 
    # kr is the rate of dark reversion (assume happens during day and night)
    dBdt=pB*L*(mB-B)-kr*B
    
    # ELF3p concentration 
    # where t1 is time (120h in 500 timesteps) % 24 --> specifies hour in day ( to determine if L vs D)

    dEdt=elf3p(t1,Day,pE1,pE2)-dE*E

    
    # PIF1, PIF3, PIF4, and PIF5 concentration 
    # pP = PIF production rate 
    # pPE = intensity of ELF3 inhibition of PIF expression 
    # kPC intensity of COP1s inhibition of PIF degradation 
    # dPB = intensity of phyBs promotion of PIF degradation and inativation 
    dPdt=pP/(1+pPE*E)-dP*P/(1+kPC*C)-dPB*P*B
    
    # COP1 concentration 
    # pCL, pCD = COP1 produciton rates during day and night 
    dCdt=pCL*L+pCD*(1-L)-dC*C
    
    # hypocotyl growth in mm 
    # pG = basal rate of hypocotyl growth 
    # kG = conversion between PIFs targets (PIL1, XTH7, ATHB2) and growth 
    # pGK is molecule K's intensity of its effect on growth, where pGH is related to HY5
    # pHC is the intensity of COP1s inhibition on HY5 
    dGdt=pG+kG*pGP*P/(1+pGP*P+pGE*E+pGB*B+pGH/(1+pHC*C))
    
    if 'elf3-8' in mut:
        dEdt=0
    if 'phyB' in mut:
        dBdt=0
    dydt=[dBdt, dEdt, dPdt, dCdt, dGdt]
    return dydt  
   
def model_results(Temp, Daylength, params, mutants):
    hypo={}
    tot={}
    for mut in mutants:
        hypo[mut]={}
        tot[mut]={}
        for T,D in itertools.product(Temp,Daylength):
                key=str(T)+'_'+str(D)
                # start, stop, number 
                time=np.linspace(0,120,500)
                y0=[0,0,0,0,0]
                # *params unpacks the list of imported parameters as the growth function 
                # has y, time, T, D, mut, and all params as arguments 
                sol=odeint(growth, y0, time, args=(T, D, mut, *params))
                # R is returning the growth at the last time step – last row and lost col of soln matrix
                R=sol[-1,4]
                hypo[mut][key]=R
                tot[mut][key]=sol
    return hypo, tot

## still unclear what arrows and decorrelate functions do exactly 

#magnitude of effect of expression - can delete
def arrows(Temp, Daylength, params, mutants):
    S={}
    mut='Col'
    for T,D in itertools.product(Temp,Daylength):
        key=str(T)+'_'+str(D)
        time=np.linspace(0,120,500)
        y0=[0,0,0,0,0]
        sol=odeint(growth, y0, time, args=(T, D, mut, *params))
        S[key]=pd.DataFrame(sol, columns=['phyB', 'ELF3', 'PIF4', 'COP1', 'Growth'])
    return S

#unfinished
def decorrelate(Temp, Daylength, params, mutants):
    hypo={}
    tot={}
    for mut in mutants:
        hypo[mut]={}
        tot[mut]={}
        for T,D in itertools.product(Temp,Daylength):
                key=str(T)+'_'+str(D)
                time=np.linspace(0,120,500)
                y0=[0,0,0,0,0]
                sol=odeint(growth, y0, time, args=(T, D, mut, *params))
                R=sol[-1,4]
                hypo[mut][key]=R
                tot[mut][key]=sol
    return hypo, tot


# In[4]:


#READ DATA
avgdata, stddata=read_data()

#SUFFIX FOR SIMULATIONS
suffix='01_paper'

#PARAMETERS
Temp=[22,28]
# experimental data daylength 
Daylength=[0,8,12,16,24]
# model simulation daylengths, shouldnt this be range (6, 25, 2)
Daylength2=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4]+list(range(4,25,2))
params=[]
with open('./parameters.txt', 'r') as f:
    for line in f:
        w=line.split()
        params.append(float(w[0]))


# In[5]:


Daylength2


# ### FIGURE 1 – MUTATION SIMULATIONS
# 
# #### SIMULATIONS

# In[6]:


mutants=['Col', 'cop1-4', 'COP1-OE', 'pif4', 'elf3-8', 'elf3-8_cop1-4',
         'ELF3ox', 'ELF3ox_cop1-4', 'hy5', 'phyB-9', 'PHYBox', 
         'pif4', 'PIF4ox', 'pifq']
# where hypopython is a dictionary of dictionaries containing the predicted growth 
# of each mutant for each temp/day length (keys = mutant, temp_daylength), values = predicted growth 
hypo_python, tot_python=model_results(Temp,Daylength2,params,mutants)


# #### WRITE DATA FROM SIMULATIONS

# In[ ]:





# In[7]:


with pd.ExcelWriter('./hypo_results_fig1_data.xlsx') as writer:
    #writer = './hypo_results/fig1_data.xlsx'
    for mut in mutants:
        W=[]
        for (D,T) in itertools.product(Daylength2,Temp):
            key=str(T)+'_'+str(D)
            try:
                W.append([D,T,avgdata[mut][key],stddata[mut][key],hypo_python[mut][key]])
            except:
                W.append([D,T,'NaN','NaN',hypo_python[mut][key]])
        DFdata=pd.DataFrame.from_records(W, columns=['Daylength', 'Temperature', 'Average Growth (experimental)', 'Standard Deviation Growth (experimental)', 'Model prediction'])
        DFdata.to_excel(writer, sheet_name=mut, float_format='%.3f',index=False)


# In[8]:


DFdata


# In[9]:


hypo_python


# #### PLOT SIMULATIONS

# In[10]:


for i1,mut in enumerate(mutants):
    print(i1, mut)


# In[11]:


get_ipython().run_line_magic('matplotlib', 'inline')
fig=plt.figure(figsize=(7,6))
ncols=4
nrows=4
for i1,mut in enumerate(mutants):
    ax=fig.add_subplot(nrows,ncols,i1+1)
    #PLOT SIMULATIONS PYTHON
    hp22=[]
    hp28=[]
    for D in Daylength2:
        key22='22_'+str(D)
        key28='28_'+str(D)
        hp22.append(hypo_python[mut][key22])
        hp28.append(hypo_python[mut][key28])  
    ax.plot(Daylength2, hp22, 'k', label='22 ºC')
    ax.plot(Daylength2, hp28, 'r', label='28 ºC')
    #DATA
    d22=[]
    d28=[]
    s22=[]
    s28=[]
    for D in Daylength:
        key22='22_'+str(D)
        key28='28_'+str(D)
        if mut in avgdata:
            d22.append(avgdata[mut][key22])
            d28.append(avgdata[mut][key28])
            s22.append(stddata[mut][key22])
            s28.append(stddata[mut][key28])
    if mut in avgdata:
        ax.errorbar(Daylength, d22, yerr=s22, fmt='o', color='k')
        ax.errorbar(Daylength, d28, yerr=s28, fmt='o', color='r')
    if mut=='Col':
        ax.set_title(mut, size=10)
    else:
        ax.set_title(mut, style='italic', size=10)
    ax.set_ylim([0,20])
    if i1==0:
        ax.legend(loc='upper right', frameon=False)
    if i1>7:
        ax.set_xlabel('Daylength (hours)', size=10)
    if i1%4==0:
        ax.set_ylabel('growth (mm)', size=10)
    ax.set_xticks([0,4,8,12,16,20,24])
    ax.set_xticklabels([0,4,8,12,16,20,24], size=5)
    # setting label format to integer with 0 folating decimals
    label_format = '{:,.0f}'
    ticks_loc = ax.get_yticks().tolist()
    ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    ax.set_yticklabels([label_format.format(x) for x in ticks_loc], size=5)
    
fig.tight_layout()


# ### FIGURE 2 – PROTEIN DYNAMICS IN SHORT DAY

# In[12]:


time=np.linspace(0,120,500)
mut='Col'
B22=tot_python[mut]['22_8'][:,0]
B28=tot_python[mut]['28_8'][:,0]
E22=tot_python[mut]['22_8'][:,1]
E28=tot_python[mut]['28_8'][:,1]
P22=tot_python[mut]['22_8'][:,2]
P28=tot_python[mut]['28_8'][:,2]
C22=tot_python[mut]['22_8'][:,3]
C28=tot_python[mut]['28_8'][:,3]
G22=tot_python[mut]['22_8'][:,4]
G28=tot_python[mut]['28_8'][:,4]


# #### WRITE RESULTS

# In[13]:


DF=pd.DataFrame({'Time':time, 'phyb22':B22, 'phyb28':B28,
                 'ELF322':E22, 'ELF328':E28, 'PIF422':P22,
                 'PIF428':P28, 'COP122':C22, 'COP128':C28})
DF.to_csv('./hypo_results_short_day_proteins.csv', float_format='%.3f', index=False)

with pd.ExcelWriter('./hypo_results_short_day_proteins.xlsx') as writer:
    DF.to_excel(writer, float_format='%.3f', index=False)


# #### PLOT RESULTS

# In[14]:


get_ipython().run_line_magic('matplotlib', 'inline')
def modsavefigs(time, b22, b28, Daylength, Title):
    xmin=24
    xmax=120
    mask=(time>=xmin) & (time<=xmax)
    Time=time[mask]
    B22=b22[mask]
    B28=b28[mask]
    fig=plt.figure()
    ax=plt.gca()
    ax.plot(Time,B22,'k', label='22 ºC')
    ax.plot(Time,B28,'r', label='28 ºC')
    #
    for days in range(0,5):
        time_day=np.linspace(days*24+Daylength,(days*24)+24, 100)
        ax.fill_between(time_day, 0, 10000, facecolor='grey', alpha=0.5)
    ax.set_ylim([0.9*min(min(B22),min(B28)), 1.10*max(max(B22),max(B28))])
    ax.set_xlabel('ZT (h)', size=15)
    ax.set_ylabel(Title, size=15)
    
    label_format_x = '{:,.0f}'
    label_format_y = '{:,.1f}'
    ticks_loc_x = ax.get_xticks().tolist()
    ticks_loc_y = ax.get_yticks().tolist()
    ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc_x))
    ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc_y))
    ax.set_xticklabels([label_format_x.format(x) for x in ticks_loc_x], size=10)
    ax.set_yticklabels([label_format_y.format(x) for x in ticks_loc_y], size=10)
    ax.legend(loc='upper right',frameon=False)
    ax.set_xlim([xmin, xmax])
    if Title=='COP1': 
        ax.set_yscale('log')
    #fig.savefig('figures/fig2_'+Title+'_'+suffix+'.pdf', bbox_inches='tight')
modsavefigs(time, B22, B28, 8, 'PHYB')
modsavefigs(time, E22, E28, 8, 'ELF3')
modsavefigs(time, P22, P28, 8, 'PIF4')
modsavefigs(time, C22, C28, 8, 'COP1')
modsavefigs(time, G22, G28, 8, 'growth (mm)')


# ### SUPPFIG 2 - PROTEIN DYNAMICS IN LONG DAY
# 
# For WT plants

# In[15]:


time=np.linspace(0,120,500)
B22=tot_python['Col']['22_16'][:,0]
B28=tot_python['Col']['28_16'][:,0]
E22=tot_python['Col']['22_16'][:,1]
E28=tot_python['Col']['28_16'][:,1]
P22=tot_python['Col']['22_16'][:,2]
P28=tot_python['Col']['28_16'][:,2]
C22=tot_python['Col']['22_16'][:,3]
C28=tot_python['Col']['28_16'][:,3]
G22=tot_python['Col']['22_16'][:,4]
G28=tot_python['Col']['28_16'][:,4]


# #### WRITE RESULTS

# In[16]:


DF=pd.DataFrame({'Time':time, 'phyb22':B22, 'phyb28':B28,
                 'ELF322':E22, 'ELF328':E28, 'PIF422':P22,
                 'PIF428':P28, 'COP122':C22, 'COP128':C28})
DF.to_csv('./hypo_results_long_day_proteins.csv', float_format='%.3f', index=False)

with pd.ExcelWriter('./hypo_results_long_day_proteins.xlsx') as writer:
    DF.to_excel(writer, float_format='%.3f', index=False)


# In[17]:


DF


# #### PLOT

# In[18]:


modsavefigs(time, B22, B28, 16, 'PHYB')
modsavefigs(time, E22, E28, 16, 'ELF3')
modsavefigs(time, P22, P28, 16, 'PIF4')
modsavefigs(time, C22, C28, 16, 'COP1')
modsavefigs(time, G22, G28, 16, 'growth (mm)')


# <a id='calling_matlab' ></a>
# ## <font color='lightseagreen'>CALLING MATLAB CODE FROM PYTHON</font>

# Initially, I had attempted to import the outputs from the clock model and feed them directly into the hypocotyl elongation model. This method proved to be unsuccessful as the the ode23 solver from matlab uses flexible timesteps in order to solve the differential equations, while the solver used in the hyptoctyl elongation model uses fixed timesteps. Thus could be worked around as definitive time steps can be defined for the ode23 solver, however, this method is not the most robust. The results of this attempt are shown below. 
# <br> 
# 
# Instead, the hypocotyl elongation model as defined above was combined with the clock model (clock_linked_growth_model.m and all_params.m files as found on github) and the joint system was solved using the ode23 solver. For more details and results, continue to [Joint Matlab Model](#joint_matlab_model) 

# ### Flexible timesteps
# https://uk.mathworks.com/matlabcentral/answers/92961-how-do-i-use-a-fixed-step-size-with-ode23-and-ode45-in-matlab

# In[19]:


eng = matlab.engine.start_matlab()
path = '/Users/chloe/Desktop/thesis'
eng.cd(path, nargout=0)
clock_model_test = eng.clock_model_test(nargout=0)
time_steps = pd.Series(eng.workspace['T'])
num_days = eng.workspace['t'][0][1]
time_steps = time_steps.apply(lambda x: x[0])
clock_protein_levels = eng.workspace['Y']
df_cpl = pd.DataFrame.from_records(clock_protein_levels)

clock_outputs = 'LHY mRNA,P,GI-ZTL,GI-ELF3 cytoplasm,LHY prot,TOC1 mRNA,PRR9 prot,PRR5 (NI) mRNA,PRR5 (NI) prot,GI prot cytoplasm,TOC1 prot,ZTL,EC,GI mRNA,PRR9 mRNA,PRR7 mRNA,PRR7 prot,ELF4 mRNA,ELF4 prot,LHY prot modif,ABAR mRNA,COP1 cytoplasm,ELF3 mRNA,ELF3 cytoplasm,ELF3 nuclear,COP1 nuclear night,COP1 nuclear day,LUX mRNA,LUX prot,PP2C prot,SnRK2 prot,stomata'

def add_clock_header(df, clock_outputs): 
    clock_outputs = clock_outputs.split(',')
    unknown_outputs = ['?', '?', '?']
    clock_outputs = clock_outputs + unknown_outputs
    header = []
    for item in clock_outputs: 
        item = item.replace(" ", '_')
        header.append(item)
    df.columns = header

add_clock_header(df_cpl, clock_outputs)
    
df_cpl['time_steps'] = time_steps
df_cpl = df_cpl[['time_steps']+[x for x in df_cpl.columns if x != 'time_steps']]

ELF3_nuc = df_cpl[['time_steps', 'ELF3_nuclear']]
# plotting against index 
ELF3_nuc.ELF3_nuclear.plot()
plt.title('Flexible Time Steps')


# In[20]:


df_cpl


# ### Re-ran with fixed time intervals

# In[21]:


eng = matlab.engine.start_matlab()
path = '/Users/chloe/Desktop/thesis'
eng.cd(path, nargout=0)
clock_model_test = eng.clock_model_test_fixed_t(nargout=0)
time_steps = pd.Series(eng.workspace['T'])
num_days = eng.workspace['t'][0][1]
time_steps = time_steps.apply(lambda x: x[0])
clock_protein_levels = eng.workspace['Y']
df_cpl = pd.DataFrame.from_records(clock_protein_levels)

clock_outputs = 'LHY mRNA,P,GI-ZTL,GI-ELF3 cytoplasm,LHY prot,TOC1 mRNA,PRR9 prot,PRR5 (NI) mRNA,PRR5 (NI) prot,GI prot cytoplasm,TOC1 prot,ZTL,EC,GI mRNA,PRR9 mRNA,PRR7 mRNA,PRR7 prot,ELF4 mRNA,ELF4 prot,LHY prot modif,ABAR mRNA,COP1 cytoplasm,ELF3 mRNA,ELF3 cytoplasm,ELF3 nuclear,COP1 nuclear night,COP1 nuclear day,LUX mRNA,LUX prot,PP2C prot,SnRK2 prot,stomata'

def add_clock_header(df, clock_outputs): 
    clock_outputs = clock_outputs.split(',')
    unknown_outputs = ['?', '?', '?']
    clock_outputs = clock_outputs + unknown_outputs
    header = []
    for item in clock_outputs: 
        item = item.replace(" ", '_')
        header.append(item)
    df.columns = header

add_clock_header(df_cpl, clock_outputs)
    
df_cpl['time_steps'] = time_steps
df_cpl = df_cpl[['time_steps']+[x for x in df_cpl.columns if x != 'time_steps']]

ELF3_nuc = df_cpl[['time_steps', 'ELF3_nuclear']]
# plotting against index 
ELF3_nuc.ELF3_nuclear.plot()
plt.title('Fixed Time Steps')


# In[22]:


ELF3_nuc


# We need to retrieve the expression values of ELF3 at 500 different timepoints that are equally spaced over a 120 hour period. The only problem with this, when trying to integrate the models, is that ode23 uses a flexible solver that does not use a fixed timestep – 

# <a id='joint_matlab_model'></a>
# ## <font color='lightseagreen'>JOINT MATLAB MODEL</font>
# 
# The above hypocotyl elongation model was linked to the plant circadian clock model developed by Pokhilko et al., 2013 and solved using ode23 solver in MATLAB. The ad hoc expression profiles of ELF3 and COP1 were replaced by those generated via the clock model. Specifically, ELF3 dynamics were set to equal the ELF3 nuclear concentrations as defined by y(25) and COP1 was set to be the sum of the light and dark nuclear proteins given by y(26) and y(27). The starting concentrations of these, while 0 in the hypocotyl elongation model, were set to equal the corresponding initial values as defined in Pokhilko et al – (0.2234 and 1.2513 respectively). 
# 
# This model can be found in the **clock_linked_growth_model.m** and **all_params.m** files on my [github](https://github.com/cbantock/arabidopsis_clock_growth_model)
# 
# 
# The Matlab engine for python is used in order to call and run the model in matlab and then subsequently import the outputs of the model to a jupyter notebook for further analysis and comparison to the unlinked hypocotyl elongation model. 
# 

# In[23]:


#starting engine
eng = matlab.engine.start_matlab()
path = '/Users/chloe/Desktop/thesis'
eng.cd(path, nargout=0)

#calling the linked growth model 
linked_growth_model = eng.all_params(nargout=0)


# ### **Extracting MATLAB outputs** 
# 
# Extracting the timesteps and protein levels for each simulation run. Simulations were run for the following different daylengths: 
# 
# * 0L/24D
# * 4L/20D
# * 8L/16D
# * 12L/12L
# * 16L/8D
# * 20L/4D
# * 24L
# 

# In[24]:


# consistent for all runs 
num_days = eng.workspace['t'][0][1]
daylength = np.arange(0,28,4).tolist()
sim_time_points = []
daylength_sims = []
for i in daylength: 
    daylength_sims.append('Y'+str(i))
    sim_time_points.append('T'+str(i))
sim_time_points


# In[25]:


def get_outputs(model, sim_names, sim_time_points): 
    '''
    Returns protein dynamics and timesteps used to solve ODE system of
    model simulations. 
    
    :param a: matlab model run through matlab.eng
    :param b: names of simulations (ie Y0, Y4 etc)
    :type b: list 
    :param c: names of T time steps used to solve ODE system from sim (T0, T4 etc)
    :type c: list 
    
    
    :rtype: dict of dictionaries 
    :return: dictionary of simulations run and their corresponding outputs and time steps. 
    '''
    sim_outputs = {}
    for (name, time) in zip(sim_names, sim_time_points): 
        sim_outputs[name] = {}
        sim_outputs[name]['time_steps'] = pd.Series(eng.workspace[time]).apply(lambda x: x[0])
        sim_outputs[name]['prot_dyn'] = eng.workspace[name]
    return sim_outputs


# In[26]:


sim_outputs = get_outputs(linked_growth_model, daylength_sims, sim_time_points)


# In[27]:


clock_outputs = 'LHY mRNA,P,GI-ZTL,GI-ELF3 cytoplasm,LHY prot,TOC1 mRNA,PRR9 prot,PRR5 (NI) mRNA,PRR5 (NI) prot,GI prot cytoplasm,TOC1 prot,ZTL,EC,GI mRNA,PRR9 mRNA,PRR7 mRNA,PRR7 prot,ELF4 mRNA,ELF4 prot,LHY prot modif,ABAR mRNA,COP1 cytoplasm,ELF3 mRNA,ELF3 cytoplasm,ELF3 nuclear,COP1 nuclear night,COP1 nuclear day,LUX mRNA,LUX prot,PP2C prot,SnRK2 prot,stomata,?,?,?,PHYB,ELF3,PIF,COP1,growth'

clock_outputs = clock_outputs.split(',')
clock_outputs
header = ['time_steps']
for item in clock_outputs: 
    item = item.replace(" ", '_')
    header.append(item)


# In[28]:


# def create_df(output_dict, sim): 
#     '''
#     Returns a dataframe of protein dynamics over time of the system being simulated.
    
#     :param a: dictionary of simulations containing outputs Y, T
#     :param b: name of simulation
#     :type b: str
    
#     rtype: dataframe
#     '''
#     sim_df = pd.concat([pd.DataFrame(output_dict[sim]['time_steps']), pd.DataFrame(output_dict[sim]['prot_dyn'])],axis=1)
#     sim_df.columns=(header)
#     return sim_df

def create_df(output_dict, sim): 
    '''
    Returns a dataframe of protein dynamics over time of the system being simulated.
    
    :param a: dictionary of simulations containing outputs Y, T
    :param b: name of simulation
    :type b: str
    
    rtype: dataframe
    '''
    if pd.DataFrame(output_dict[sim]).shape[1] <=40: 
        sim_df = pd.concat([pd.DataFrame(output_dict[sim]['time_steps']), pd.DataFrame(output_dict[sim]['prot_dyn'])],axis=1)
        sim_df.columns=(header)
    else: 
        sim_df = pd.DataFrame(output_dict[sim])
        sim_df.columns=(header)
    return sim_df
    


# In[29]:


Y0 = create_df(sim_outputs,'Y0')
Y4 = create_df(sim_outputs,'Y4')
Y8 = create_df(sim_outputs,'Y8')
Y12 = create_df(sim_outputs,'Y12')
Y16 = create_df(sim_outputs,'Y16')
Y20 = create_df(sim_outputs,'Y20')
Y24 = create_df(sim_outputs,'Y24')


# In[30]:


Y0


# In[31]:


Y12_COP1 = pd.concat([Y12['time_steps'], Y12['COP1_nuclear_night'], Y12['COP1_nuclear_day'], Y12['COP1']], axis=1)
Y12_COP1


# In[32]:


fig, ax = plt.subplots(1, sharex=True, figsize=(8,5))


# COP1
ax.plot(Y12_COP1['time_steps'], Y12_COP1['COP1_nuclear_night'], label = 'night', color='blue' )
ax.plot(Y12_COP1['time_steps'], Y12_COP1['COP1_nuclear_day'], label = 'day', color='green')
ax.plot(Y12_COP1['time_steps'], Y12_COP1['COP1'], label = 'combined', color='orange', linestyle='--')
#ax[1].legend(loc='best')
ax.set_ylabel('Expression')
ax.set_xlabel('Hours')
ax.set_title('COP1 Dynamics')
for days in range(0,5):
    time_day=np.linspace(days*24+12,(days*24)+24, 100)
    ax.fill_between(time_day, 0, 1.5, facecolor='grey', alpha=0.5)
ax.legend()


# In[34]:


linked_growth_12 = Y12[['time_steps', 'PHYB', 'ELF3', 'PIF', 'COP1', 'growth']]
linked_growth_12


# In[35]:


# results for 12L/12D
fig, ax = plt.subplots(figsize=(8,5))

ax.plot(linked_growth_12['time_steps'], linked_growth_12['ELF3'], label='ELF3')
ax.plot(linked_growth_12['time_steps'], linked_growth_12['COP1'], label='COP1')
ax.plot(linked_growth_12['time_steps'], linked_growth_12['growth'], label='growth')
plt.title('12L Sim dynamics')
plt.xlabel('Hours')
plt.ylabel('Expression levels')

for days in range(0,5):
    time_day=np.linspace(days*24+12,(days*24)+24, 100)
    ax.fill_between(time_day, 0, 2, facecolor='grey', alpha=0.5)
plt.legend()


# In[36]:


#plt.plot(linked_growth_12['time_steps'], linked_growth_12['ELF3'])


# ### **Comparison to original growth model** 

# #### Import data for WT simultions 
# (mut = 'Col')

# In[37]:


# retrieve experimental data
hypo_results_long_day_url = 'https://raw.githubusercontent.com/cbantock/arabidopsis_clock_growth_model/main/hypo_results_long_day_proteins.csv'
hypo_results_long_day_data = requests.get(hypo_results_long_day_url).content
hypo_results_long_day_df = pd.read_csv(io.StringIO(hypo_results_long_day_data.decode('utf-8')))


# In[38]:


hypo_results_long_day_df


# In[39]:


hypo_results_short_day_url='https://raw.githubusercontent.com/cbantock/arabidopsis_clock_growth_model/main/hypo_results_short_day_proteins.csv'
hypo_results_short_day_data = requests.get(hypo_results_short_day_url).content
hypo_results_short_day_df = pd.read_csv(io.StringIO(hypo_results_short_day_data.decode('utf-8')))


# In[40]:


hypo_results_short_day_df


# In[41]:


DFdata_22 = DFdata[DFdata['Temperature']==22]
DFdata_22                   


# #### **Plots for short day** 
# 
# dusk = 8 
# <br>
# 8L, 16D conditions 

# In[42]:


linked_growth_8 = Y8[['time_steps', 'PHYB', 'ELF3', 'PIF', 'COP1', 'growth']]


# In[43]:


fig, ax = plt.subplots(4, sharex=True, figsize=(10,15))

# ELF3
ax[0].plot(linked_growth_8['time_steps'], linked_growth_8['ELF3'], label = 'joint')
ax[0].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['ELF322'], label = 'original')
#ax[0].legend(loc='best')
ax[0].set_title('ELF3')

# COP1
ax[1].plot(linked_growth_8['time_steps'], linked_growth_8['COP1'], label = 'joint' )
ax[1].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['COP122'], label = 'original')
#ax[1].legend(loc='best')
ax[1].set_title('COP1')

# PHYB
ax[2].plot(linked_growth_8['time_steps'], linked_growth_8['PHYB'], label = 'joint' )
ax[2].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['phyb22'], label = 'original')
#ax[2].legend(loc='best')
ax[2].set_title('PHYB')

# PIF
ax[3].plot(linked_growth_8['time_steps'], linked_growth_8['PIF'], label = 'joint' )
ax[3].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['PIF422'], label = 'original')
#ax[3].legend(loc='best')
ax[3].set_title('PIF')

fig.legend(labels = ['linked', 'original'])


# ##### **LUX, ELF4, ELF3 and the EC**

# It seems weird that ELF3 peaks at the end of the night, it should be the other way around (see Fig S1b and c in our paper, attached)
# 
# % y(13) EC <br>
# % y(18) ELF4 mRNA<br>
# % y(19) ELF4 prot<br>
# % y(23) ELF3 mRNA<br>
# % y(24) ELF3 cytoplasm<br>
# % y(25) ELF3 nuclear<br>
# % y(28) LUX mRNA<br>
# % y(29) LUX prot<br>

# In[44]:


fig, ax = plt.subplots(4, sharex=True, figsize=(10,15))

# ELF3
ax[0].plot(Y8['time_steps'], Y8['ELF3_mRNA'], label = 'ELF3_mRNA', color='k', linestyle='--')
ax[0].plot(Y8['time_steps'], Y8['ELF3_nuclear'], label = 'ELF3_nuclear', color='red')
ax[0].plot(Y8['time_steps'], Y8['ELF3_cytoplasm'], label = 'ELF3_cytoplasm', color='k',linestyle='dotted' )
ax[0].legend(loc='best')
ax[0].set_ylabel('Expression')
ax[0].set_xlabel('Hours')
ax[0].set_title('ELF3 Dynamics')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[0].fill_between(time_day, 0, 0.6, facecolor='grey', alpha=0.5)
    
# ELF4
ax[1].plot(Y8['time_steps'], Y8['ELF4_mRNA'], label = 'ELF4_mRNA', color='k', linestyle='--')
ax[1].plot(Y8['time_steps'], Y8['ELF4_prot'], label = 'ELF4_prot', color='red')
ax[1].legend(loc='best')
ax[1].set_ylabel('Expression')
ax[1].set_xlabel('Hours')
ax[1].set_title('ELF4 Dynamics')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[1].fill_between(time_day, 0, 1.6, facecolor='grey', alpha=0.5)
    
# LUX
ax[2].plot(Y8['time_steps'], Y8['LUX_mRNA'], label = 'LUX_mRNA', color='k', linestyle='--')
ax[2].plot(Y8['time_steps'], Y8['LUX_prot'], label = 'LUX_prot', color='red')
ax[2].legend(loc='best')
ax[2].set_ylabel('Expression')
ax[2].set_xlabel('Hours')
ax[2].set_title('LUX Dynamics')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[2].fill_between(time_day, 0, 2.6, facecolor='grey', alpha=0.5)
    
# EC
ax[3].plot(Y8['time_steps'], Y8['EC'], label = 'EC', color='red')
ax[3].legend(loc='best')
ax[3].set_ylabel('Expression')
ax[3].set_xlabel('Hours')
ax[3].set_title('EC Dynamics')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[3].fill_between(time_day, 0, 0.2, facecolor='grey', alpha=0.5)


# ##### **Rescaled**

# In[45]:


fig, ax = plt.subplots(4, sharex=True, figsize=(10,15))

# ELF3
ax[0].plot(linked_growth_8['time_steps'], linked_growth_8['ELF3'], label = 'joint', color='red')
ax[0].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['ELF322'], label = 'original', color='k')
#ax[0].legend(loc='best')
ax[0].set_yscale('log')
ax[0].set_ylabel('Expression, logscale')
ax[0].set_xlabel('Hours')
ax[0].set_title('ELF3')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[0].fill_between(time_day, 0, 6, facecolor='grey', alpha=0.5)

# COP1
ax[1].plot(linked_growth_8['time_steps'], linked_growth_8['COP1'], label = 'joint', color='red' )
ax[1].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['COP122'], label = 'original', color='k')
#ax[1].legend(loc='best')
ax[1].set_ylabel('Expression, logscale')
ax[1].set_xlabel('Hours')
ax[1].set_yscale('log')
ax[1].set_title('COP1')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[1].fill_between(time_day, 0, 60, facecolor='grey', alpha=0.5)

# PHYB
ax[2].plot(linked_growth_8['time_steps'], linked_growth_8['PHYB'], label = 'joint', color='red')
ax[2].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['phyb22'], label = 'original', color='k')
#ax[2].legend(loc='best')
ax[2].set_ylabel('Expression')
ax[2].set_xlabel('Hours')
ax[2].set_title('PHYB')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[2].fill_between(time_day, 0, 1, facecolor='grey', alpha=0.5)

# PIF
ax[3].plot(linked_growth_8['time_steps'], linked_growth_8['PIF'], label = 'joint', color='red' )
ax[3].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['PIF422'], label = 'original', color='k')
#ax[3].legend(loc='best')
ax[3].set_ylabel('Expression')
ax[3].set_xlabel('Hours')
ax[3].set_title('PIF')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[3].fill_between(time_day, 0, 4, facecolor='grey', alpha=0.5)

fig.legend(labels = ['linked', 'original'])

plt.suptitle('Protein dynamics in short day – 8L/16D', size=20)


# #### **Plots for long day** 

# In[46]:


linked_growth_16 = Y16[['time_steps', 'PHYB', 'ELF3', 'PIF', 'COP1', 'growth']]


# In[47]:


fig, ax = plt.subplots(4, sharex=True, figsize=(10,15))

# ELF3
ax[0].plot(linked_growth_16['time_steps'], linked_growth_16['ELF3'], label = 'joint')
ax[0].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['ELF322'], label = 'original')
#ax[0].legend(loc='best')
ax[0].set_title('ELF3')

# COP1
ax[1].plot(linked_growth_16['time_steps'], linked_growth_16['COP1'], label = 'joint' )
ax[1].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['COP122'], label = 'original')
#ax[1].legend(loc='best')
ax[1].set_title('COP1')

# PHYB
ax[2].plot(linked_growth_16['time_steps'], linked_growth_16['PHYB'], label = 'joint' )
ax[2].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['phyb22'], label = 'original')
#ax[2].legend(loc='best')
ax[2].set_title('PHYB')

# PIF
ax[3].plot(linked_growth_16['time_steps'], linked_growth_16['PIF'], label = 'joint' )
ax[3].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['PIF422'], label = 'original')
#ax[3].legend(loc='best')
ax[3].set_title('PIF')

fig.legend(labels = ['linked', 'original'])


# ##### **Rescaled** 

# In[48]:


fig, ax = plt.subplots(4, sharex=True, figsize=(10,15))

# ELF3
ax[0].plot(linked_growth_16['time_steps'], linked_growth_16['ELF3'], label = 'joint', color='red')
ax[0].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['ELF322'], label = 'original', color='k')
#ax[0].legend(loc='best')
ax[0].set_yscale('log')
ax[0].set_ylabel('Expression, logscale')
ax[0].set_xlabel('Hours')
ax[0].set_title('ELF3')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[0].fill_between(time_day, 0, 6, facecolor='grey', alpha=0.5)

# COP1
ax[1].plot(linked_growth_16['time_steps'], linked_growth_16['COP1'], label = 'joint', color='red' )
ax[1].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['COP122'], label = 'original', color='k')
#ax[1].legend(loc='best')
ax[1].set_ylabel('Expression, logscale')
ax[1].set_xlabel('Hours')
ax[1].set_yscale('log')
ax[1].set_title('COP1')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[1].fill_between(time_day, 0, 60, facecolor='grey', alpha=0.5)

# PHYB
ax[2].plot(linked_growth_16['time_steps'], linked_growth_16['PHYB'], label = 'joint', color='red')
ax[2].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['phyb22'], label = 'original', color='k')
#ax[2].legend(loc='best')
ax[2].set_ylabel('Expression')
ax[2].set_xlabel('Hours')
ax[2].set_title('PHYB')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[2].fill_between(time_day, 0, 1, facecolor='grey', alpha=0.5)

# PIF
ax[3].plot(linked_growth_16['time_steps'], linked_growth_16['PIF'], label = 'joint', color='red' )
ax[3].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['PIF422'], label = 'original', color='k')
#ax[3].legend(loc='best')
ax[3].set_ylabel('Expression')
ax[3].set_xlabel('Hours')
ax[3].set_title('PIF')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[3].fill_between(time_day, 0, 4, facecolor='grey', alpha=0.5)

fig.legend(labels = ['linked', 'original'])

plt.suptitle('Protein dynamics in long day – 16L/8D', size=20)


# #### Growth vs Daylength

# In[49]:


sim_growth = []
sim_outputs = [Y0, Y4, Y8, Y12, Y16, Y20, Y24]
for sim in sim_outputs:
    sim_growth.append(sim['growth'].iloc[-1])
linked_growth_df = pd.DataFrame(sim_growth, daylength, columns=['growth'])



# In[50]:


linked_growth_df


# In[51]:


fig, ax = plt.subplots()

ax.plot(daylength,linked_growth_df['growth'], color='red', label = 'linked', marker='.')
hp22=[]
for D in Daylength2:
    key22='22_'+str(D)
    hp22.append(hypo_python['Col'][key22])  
ax.plot(Daylength2, hp22, 'k', label='original', marker='.')
ax.legend(loc='best')
plt.title("Growth vs Daylength")
plt.xlabel('Daylength')
plt.ylabel('Growth (mm)')


# #### Comments 
# 
# – large discrepancies in expression values for ELF3 and COP1
# 
# – ELF3 exhibits softer changes in expression during transition from day to night/night to day. In the hypocotyl model, its expression is modeled as a square wave that peaks at the onset of dusk. In the linked model, ELF3 peaks at the onset of dawn and then steadily declines.
# 
# – the COP1 expression profile differs in that expression initially decreases at the onset of dusk, and slowly increases until midday before it begins to slowly decrease until dusk, where it rapidly decreases  whereas in the hypocotyl elongation model, expression peaks during the night and then decreases with the onset of dawn 
# 
# – PHYB expression is the same as unlinked model for short day sim as well as the long day simulation 
# 
# – the PIF expression profile is consistent in both long and short day simulations, but the absolute value is slightly higher in the linked model
# 
# – growth response as a function of daylength not as strong as seen in the original model 

# ### **Normalizing values of ELF3 and COP1** 

# Here we aim to set the average levels of ELF3 and COP1 to be the same as in the hypocotyl model. The hypocotyl elongation parameters had previously been found fitting the model to experimental data, thus, using different absolute levels of ELF3 and COP1 could effect the model outputs. We can do this by multiplying instances of ELF3 and COP1 in the growth model by constants:
# 
# **cE = (average level E in hypocotyl model) / (average level of E in linked model)**
# 
# **cC = (average level C in hypocotyl model) / (average level of C in linked model)**
# 
# **cEC=(average ELF3 levels in hypocotyl model) / (average EC levels in linked model)**
# 

# First lets find the average level of ELF3 and COP1 in each model: 

# In[52]:


daylength_sims


# In[53]:


Y0['ELF3_nuclear']


# In[54]:


sims = [Y0, Y4, Y8, Y12, Y16, Y20, Y24]


# In[55]:


avg_levels_dict = {}
avg_levels_dict['linked_model'] = {}
avg_levels_dict['growth_model'] = {}
avg_levels_dict['linked_model']['ELF3_linked'] = []
avg_levels_dict['linked_model']['COP1_linked'] = []
avg_levels_dict['linked_model']['EC_linked'] = []
avg_levels_dict['growth_model']['ELF3_growth'] = []
avg_levels_dict['growth_model']['COP1_growth'] = []

for sim in sims: 
    avg_levels_dict['linked_model']['ELF3_linked'].append(np.mean(sim['ELF3_nuclear']))
    avg_levels_dict['linked_model']['COP1_linked'].append(np.mean(sim['COP1']))
    avg_levels_dict['linked_model']['EC_linked'].append(np.mean(sim['EC']))


# In[56]:


daylengths_growth_model = ['22_0', '22_4', '22_8', '22_12', '22_16', '22_20', '22_24']

for key in daylengths_growth_model: 
    temp = pd.DataFrame(tot_python['Col'][key])
    avg_levels_dict['growth_model']['ELF3_growth'].append(np.mean(temp[1]))
    avg_levels_dict['growth_model']['COP1_growth'].append(np.mean(temp[3]))


# In[57]:


avg_levels_dict


# In[58]:


linked_model_avgs = pd.DataFrame(avg_levels_dict['linked_model'], index=daylength_sims)
growth_model_avgs = pd.DataFrame(avg_levels_dict['growth_model'], index=daylength_sims)


# In[59]:


model_avgs = pd.concat([linked_model_avgs, growth_model_avgs], axis=1)
model_avgs['cE'] = model_avgs['ELF3_growth']/model_avgs['ELF3_linked']
model_avgs['cC'] = model_avgs['COP1_growth']/model_avgs['COP1_linked']
model_avgs['cEC'] = model_avgs['ELF3_growth']/model_avgs['EC_linked']


# In[60]:


model_avgs


# #### **Running Normalized Linked Model**

# In[61]:


#starting engine
eng = matlab.engine.start_matlab()
path = '/Users/chloe/Desktop/thesis'
eng.cd(path, nargout=0)

#calling the linked growth model 
linked_growth_model_normalized = eng.all_params_normalized(nargout=0)


# In[62]:


sim_outputs_normalized = get_outputs(linked_growth_model_normalized, daylength_sims, sim_time_points)


# In[63]:


Y0_normalized = create_df(sim_outputs_normalized,'Y0')
Y4_normalized = create_df(sim_outputs_normalized,'Y4')
Y8_normalized = create_df(sim_outputs_normalized,'Y8')
Y12_normalized = create_df(sim_outputs_normalized,'Y12')
Y16_normalized = create_df(sim_outputs_normalized,'Y16')
Y20_normalized = create_df(sim_outputs_normalized,'Y20')
Y24_normalized = create_df(sim_outputs_normalized,'Y24')


# In[64]:


model_avgs['cE']['Y16']


# In[65]:


model_avgs['cC']['Y16']


# In[66]:


linked_growth_16_normalized = Y16_normalized[['time_steps', 'PHYB', 'ELF3', 'PIF', 'COP1', 'growth']]

pd.options.mode.chained_assignment = None  # default='warn'
linked_growth_16_normalized['ELF3_normalized'] = linked_growth_16_normalized['ELF3']*(35.36030532800484)
linked_growth_16_normalized['COP1_normalized'] = linked_growth_16_normalized['COP1']*model_avgs['cC']['Y16']


# In[67]:


linked_growth_16_normalized


# In[68]:


fig, ax = plt.subplots(4, sharex=True, figsize=(10,15))

# ELF3
ax[0].plot(linked_growth_16_normalized['time_steps'], linked_growth_16_normalized['ELF3'], label = 'joint', color='red')
ax[0].plot(linked_growth_16_normalized['time_steps'], linked_growth_16_normalized['ELF3_normalized'], label = 'joint', color='red', linestyle='--')
ax[0].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['ELF322'], label = 'original', color='k')
#ax[0].legend(loc='best')
#ax[0].set_yscale('log')
ax[0].set_ylabel('Expression')
ax[0].set_xlabel('Hours')
ax[0].set_title('ELF3')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[0].fill_between(time_day, 0, 8, facecolor='grey', alpha=0.5)
ax[0].legend(labels = ['linked', 'linked_normalized', 'original'], loc='best')

# COP1
ax[1].plot(linked_growth_16_normalized['time_steps'], linked_growth_16_normalized['COP1'], label = 'joint', color='red' )
ax[1].plot(linked_growth_16_normalized['time_steps'], linked_growth_16_normalized['COP1_normalized'], label = 'joint', color='red', linestyle ='--' )
ax[1].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['COP122'], label = 'original', color='k')
#ax[1].legend(loc='best')
ax[1].set_ylabel('Expression')
ax[1].set_xlabel('Hours')
#ax[1].set_yscale('log')
ax[1].set_title('COP1')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[1].fill_between(time_day, 0, 60, facecolor='grey', alpha=0.5)
ax[1].legend(labels = ['linked', 'linked_normalized', 'original'], loc='best')

# PHYB
ax[2].plot(linked_growth_16['time_steps'], linked_growth_16['PHYB'], label = 'joint', color='red' )
ax[2].plot(linked_growth_16_normalized['time_steps'], linked_growth_16_normalized['PHYB'], label = 'joint_normalized', color='red', linestyle = '--')
ax[2].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['phyb22'], label = 'original', color='k')
#ax[2].legend(loc='best')
ax[2].set_ylabel('Expression')
ax[2].set_xlabel('Hours')
ax[2].set_title('PHYB')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[2].fill_between(time_day, 0, 1, facecolor='grey', alpha=0.5)
ax[2].legend(labels = ['linked', 'linked_normalized', 'original'], loc='best')


# PIF
ax[3].plot(linked_growth_16['time_steps'], linked_growth_16['PIF'], label = 'joint', color='red' )
ax[3].plot(linked_growth_16_normalized['time_steps'], linked_growth_16_normalized['PIF'], label = 'linked_normalized', color='red', linestyle = '--' )
ax[3].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['PIF422'], label = 'original', color='k')
#ax[3].legend(loc='best')
ax[3].set_ylabel('Expression')
ax[3].set_xlabel('Hours')
ax[3].set_title('PIF')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[3].fill_between(time_day, 0, 5, facecolor='grey', alpha=0.5)
ax[3].legend(labels = ['linked', 'linked_normalized', 'original'], loc='best')


#fig.legend(labels = ['linked', 'original'])

plt.suptitle('Protein dynamics in long day – 16L/8D', size=20)


# In[69]:


sim_growth_normalized = []
sim_outputs_normalized_list = [Y0_normalized, Y4_normalized, Y8_normalized, Y12_normalized, Y16_normalized, Y20_normalized, Y24_normalized]
for sim in range(len(sim_outputs_normalized_list)):
    sim_growth_normalized.append(sim_outputs_normalized_list[sim]['growth'].iloc[-1])
linked_growth_normalized_df = pd.DataFrame(sim_growth_normalized, daylength, columns=['growth'])


fig, ax = plt.subplots()

ax.plot(daylength,linked_growth_normalized_df['growth'], color='red', label = 'linked', marker='.')
hp22=[]
for D in Daylength2:
    key22='22_'+str(D)
    hp22.append(hypo_python['Col'][key22])  
ax.plot(Daylength2, hp22, 'k', label='original', marker='.')
ax.legend(loc='best')
plt.title("Growth vs Daylength")
plt.xlabel('Daylength')
plt.ylabel('Growth (mm)')


# <a id='original_clock_model'></a>
# ## <font color='lightseagreen'>Original Clock Model – P2013 </font>

# ### As Published by Pokhilko et al. (2013)

# In[70]:


eng = matlab.engine.start_matlab()
path = '/Users/chloe/Desktop/thesis'
eng.cd(path, nargout=0)

#calling original model 
clock_model = eng.ABA_art_final(nargout=0)


# In[71]:


clock_model_sim = get_outputs(clock_model, ['Y'], ['T'])

clock_outputs_og = 'LHY mRNA,P,GI-ZTL,GI-ELF3 cytoplasm,LHY prot,TOC1 mRNA,PRR9 prot,PRR5 (NI) mRNA,PRR5 (NI) prot,GI prot cytoplasm,TOC1 prot,ZTL,EC,GI mRNA,PRR9 mRNA,PRR7 mRNA,PRR7 prot,ELF4 mRNA,ELF4 prot,LHY prot modif,ABAR mRNA,COP1 cytoplasm,ELF3 mRNA,ELF3 cytoplasm,ELF3 nuclear,COP1 nuclear night,COP1 nuclear day,LUX mRNA,LUX prot,PP2C prot,SnRK2 prot,stomata,?,?,?'

clock_outputs_og = clock_outputs_og.split(',')
header = ['time_steps']
for item in clock_outputs_og: 
    item = item.replace(" ", '_')
    header.append(item)
clock_outputs_df = create_df(clock_model_sim, 'Y')


# In[72]:


clock_outputs_df


# In[73]:


fig, ax = plt.subplots(4, sharex=True, figsize=(10,15))

# ELF3 joint sim
ax[0].plot(Y12['time_steps'], Y12['ELF3_mRNA'], label = 'ELF3_mRNA linked', color='red', linestyle='--')
ax[0].plot(Y12['time_steps'], Y12['ELF3_nuclear'], label = 'ELF3_nuclear linked', color='red')
ax[0].plot(Y12['time_steps'], Y12['ELF3_cytoplasm'], label = 'ELF3_cytoplasm linked', color='red',linestyle='dotted' )
ax[0].legend(loc='best')
ax[0].set_ylabel('Expression')
ax[0].set_xlabel('Hours')
ax[0].set_title('ELF3 Dynamics')
for days in range(0,5):
    time_day=np.linspace(days*24+12,(days*24)+24, 100)
    ax[0].fill_between(time_day, 0, 0.6, facecolor='grey', alpha=0.5)

# ELF3 original clock model
ax[0].plot(clock_outputs_df['time_steps'], clock_outputs_df['ELF3_mRNA'], label = 'ELF3_mRNA', color='k', linestyle='--')
ax[0].plot(clock_outputs_df['time_steps'], clock_outputs_df['ELF3_nuclear'], label = 'ELF3_nuclear', color='k')
ax[0].plot(clock_outputs_df['time_steps'], clock_outputs_df['ELF3_cytoplasm'], label = 'ELF3_cytoplasm', color='k',linestyle='dotted' )
ax[0].legend(loc='best')
ax[0].set_ylabel('Expression')
ax[0].set_xlabel('Hours')
ax[0].set_title('ELF3 Dynamics')

    
# ELF4 linked model 
ax[1].plot(Y12['time_steps'], Y12['ELF4_mRNA'], label = 'ELF4_mRNA linked model', color='red', linestyle='--')
ax[1].plot(Y12['time_steps'], Y12['ELF4_prot'], label = 'ELF4_prot linked model', color='red')
ax[1].legend(loc='best')
ax[1].set_ylabel('Expression')
ax[1].set_xlabel('Hours')
ax[1].set_title('ELF4 Dynamics')
for days in range(0,5):
    time_day=np.linspace(days*24+12,(days*24)+24, 100)
    ax[1].fill_between(time_day, 0, 1.6, facecolor='grey', alpha=0.5)
    
# ELF4 original clock model
ax[1].plot(clock_outputs_df['time_steps'], clock_outputs_df['ELF4_mRNA'], label = 'ELF4_mRNA', color='k', linestyle='--')
ax[1].plot(clock_outputs_df['time_steps'], clock_outputs_df['ELF4_prot'], label = 'ELF4_prot', color='k')
ax[1].legend(loc='best')
ax[1].set_ylabel('Expression')
ax[1].set_xlabel('Hours')
ax[1].set_title('ELF4 Dynamics')

    
# LUX linked model
ax[2].plot(Y12['time_steps'], Y12['LUX_mRNA'], label = 'LUX_mRNA linked model', color='red', linestyle='--')
ax[2].plot(Y12['time_steps'], Y12['LUX_prot'], label = 'LUX_prot linked model', color='red')
ax[2].legend(loc='best')
ax[2].set_ylabel('Expression')
ax[2].set_xlabel('Hours')
ax[2].set_title('LUX Dynamics')
for days in range(0,5):
    time_day=np.linspace(days*24+12,(days*24)+24, 100)
    ax[2].fill_between(time_day, 0, 2.6, facecolor='grey', alpha=0.5)
    
# LUX original clock model
ax[2].plot(clock_outputs_df['time_steps'], clock_outputs_df['LUX_mRNA'], label = 'LUX_mRNA', color='k', linestyle='--')
ax[2].plot(clock_outputs_df['time_steps'], clock_outputs_df['LUX_prot'], label = 'LUX_prot', color='k')
ax[2].legend(loc='best')
ax[2].set_ylabel('Expression')
ax[2].set_xlabel('Hours')
ax[2].set_title('LUX Dynamics')

    
# EC linked model
ax[3].plot(Y12['time_steps'], Y12['EC'], label = 'EC linked model', color='red')
ax[3].legend(loc='best')
ax[3].set_ylabel('Expression')
ax[3].set_xlabel('Hours')
ax[3].set_title('EC Dynamics')
for days in range(0,5):
    time_day=np.linspace(days*24+12,(days*24)+24, 100)
    ax[3].fill_between(time_day, 0, 0.2, facecolor='grey', alpha=0.5)
    
# EC original clock model
ax[3].plot(clock_outputs_df['time_steps'], clock_outputs_df['EC'], label = 'EC', color='k')
ax[3].legend(loc='best')
ax[3].set_ylabel('Expression')
ax[3].set_xlabel('Hours')
ax[3].set_title('EC Dynamics')


# ### Original clock model with 2012 parameter values 

# 
# The parameter values for the 2012 and 2013 Pokhilko models differ and thus affect some of the
# expression dynamics of some clock genes that may affect the growth model, namely ELF3 and
# the EC. For a full comparison of the parameter values, refer to [this google sheet](https://docs.google.com/spreadsheets/d/1kb4nD8FqJTmONcxI_lefmHqHUeXDIWFlTbQCZqWdLeA/edit?usp=sharing). Pokhilko et al. (2013) discusses changes to the parameter space as compared to the P2012 model – in short, parameter values were fitted to multiple time-series datasets or constrained based on available experimental data, as in the 2012 model. Most parameters have the same value, however, the differences are attributed to changes in the structure of the model
# itself. It is highlighted that the new model run under these optimal parameter values retains most
# of the properties observed in the 2012 model (showing good fit to the data) while also explaining
# new experimental data that was not described in the 2012 model - ie, inhibition of clock genes
# by TOC1, changing of the clock period by ABA, as well as the dynamics of stomata. They do note, however, that while the new parameter values are robust to perturbations, an in depth study of this expanded parameter space falls outside the scope of the paper due to the increased complexity of the model. Thus, they cannot conclude that other parameter sets may equally describe the data.
# <br>
# 
# Looking at the parameters from the 2012 paper, 28 parameter values were measured directly or
# estimated from the experimental data and therefore had higher confidence than the rest of the
# parameters. Of these, m25, m28, m38, and g5 were changed, where:
# <br>
# 
#     – m25 and m28 are the degradation rates of HY5,COP1n and HFR21,COP1n proteins respectively
#     – m38 is the degradation rate of HY5,COP1d protein
#     – the equations for HY5/HFR21 were only used for the optimization of COP1 parameters
#     – g5 is the michaelis menten constant for TOC1 transcription
# <br>
# The 2013 parameter values were used for all simulations presented in the Final Linked Matlab Model.  Refer to the plot below to see results of the linked model run with 2012 parameter values. The protein dynamics are from 12L/12D simulations.
# 

# In[74]:


eng = matlab.engine.start_matlab()
path = '/Users/chloe/Desktop/thesis'
eng.cd(path, nargout=0)

#calling original model 
clock_model_old_params = eng.ABA_old_param_values(nargout=0)


# In[75]:


clock_model_old_params_sim = get_outputs(clock_model_old_params, ['Y_old'], ['T_old'])

clock_outputs_og = 'LHY mRNA,P,GI-ZTL,GI-ELF3 cytoplasm,LHY prot,TOC1 mRNA,PRR9 prot,PRR5 (NI) mRNA,PRR5 (NI) prot,GI prot cytoplasm,TOC1 prot,ZTL,EC,GI mRNA,PRR9 mRNA,PRR7 mRNA,PRR7 prot,ELF4 mRNA,ELF4 prot,LHY prot modif,ABAR mRNA,COP1 cytoplasm,ELF3 mRNA,ELF3 cytoplasm,ELF3 nuclear,COP1 nuclear night,COP1 nuclear day,LUX mRNA,LUX prot,PP2C prot,SnRK2 prot,stomata,?,?,?'

clock_outputs_og = clock_outputs_og.split(',')
header = ['time_steps']
for item in clock_outputs_og: 
    item = item.replace(" ", '_')
    header.append(item)
clock_model_old_params_df = create_df(clock_model_old_params_sim, 'Y_old')


# In[76]:


clock_model_old_params_df


# In[77]:


fig, ax = plt.subplots(4, sharex=True, figsize=(10,15))

# ELF3 joint sim
ax[0].plot(clock_model_old_params_df['time_steps'], clock_model_old_params_df['ELF3_mRNA'], label = 'mRNA, 2012', color='red', linestyle='--')
ax[0].plot(clock_model_old_params_df['time_steps'], clock_model_old_params_df['ELF3_nuclear'], label = 'nuclear prot, 2012', color='red')
ax[0].plot(clock_model_old_params_df['time_steps'], clock_model_old_params_df['ELF3_cytoplasm'], label = 'cyto_prot, 2012', color='red',linestyle='dotted' )
ax[0].legend(loc='best')
ax[0].set_ylabel('Expression')
ax[0].set_xlabel('Hours')
ax[0].set_title('ELF3 Dynamics')
for days in range(0,3):
    time_day=np.linspace(days*24+12,(days*24)+24, 100)
    ax[0].fill_between(time_day, 0, 0.6, facecolor='grey', alpha=0.5)

# ELF3 original clock model
ax[0].plot(clock_outputs_df['time_steps'], clock_outputs_df['ELF3_mRNA'], label = 'mRNA', color='k', linestyle='--')
ax[0].plot(clock_outputs_df['time_steps'], clock_outputs_df['ELF3_nuclear'], label = 'nuclear prot', color='k')
ax[0].plot(clock_outputs_df['time_steps'], clock_outputs_df['ELF3_cytoplasm'], label = 'cyto_prot', color='k',linestyle='dotted' )
ax[0].legend(loc='best')
ax[0].set_ylabel('Expression')
ax[0].set_xlabel('Hours')
ax[0].set_title('ELF3 Dynamics')

    
# ELF4 linked model 
ax[1].plot(clock_model_old_params_df['time_steps'], clock_model_old_params_df['ELF4_mRNA'], label = 'mRNA, 2012', color='red', linestyle='--')
ax[1].plot(clock_model_old_params_df['time_steps'], clock_model_old_params_df['ELF4_prot'], label = 'nuclear prot, 2012', color='red')
ax[1].legend(loc='best')
ax[1].set_ylabel('Expression')
ax[1].set_xlabel('Hours')
ax[1].set_title('ELF4 Dynamics')
for days in range(0,3):
    time_day=np.linspace(days*24+12,(days*24)+24, 100)
    ax[1].fill_between(time_day, 0, 1.7, facecolor='grey', alpha=0.5)
    
# ELF4 original clock model
ax[1].plot(clock_outputs_df['time_steps'], clock_outputs_df['ELF4_mRNA'], label = 'mRNA', color='k', linestyle='--')
ax[1].plot(clock_outputs_df['time_steps'], clock_outputs_df['ELF4_prot'], label = 'nuclear prot', color='k')
ax[1].legend(loc='best')
ax[1].set_ylabel('Expression')
ax[1].set_xlabel('Hours')
ax[1].set_title('ELF4 Dynamics')

    
# LUX linked model
ax[2].plot(clock_model_old_params_df['time_steps'], clock_model_old_params_df['LUX_mRNA'], label = 'mRNA, 2012', color='red', linestyle='--')
ax[2].plot(clock_model_old_params_df['time_steps'], clock_model_old_params_df['LUX_prot'], label = 'nuclear prot, 2012', color='red')
ax[2].legend(loc='best')
ax[2].set_ylabel('Expression')
ax[2].set_xlabel('Hours')
ax[2].set_title('LUX Dynamics')
for days in range(0,3):
    time_day=np.linspace(days*24+12,(days*24)+24, 100)
    ax[2].fill_between(time_day, 0, 3.65, facecolor='grey', alpha=0.5)
    
# LUX original clock model
ax[2].plot(clock_outputs_df['time_steps'], clock_outputs_df['LUX_mRNA'], label = 'mRNA', color='k', linestyle='--')
ax[2].plot(clock_outputs_df['time_steps'], clock_outputs_df['LUX_prot'], label = 'nuclear prot', color='k')
ax[2].legend(loc='best')
ax[2].set_ylabel('Expression')
ax[2].set_xlabel('Hours')
ax[2].set_title('LUX Dynamics')

    
# EC linked model
ax[3].plot(clock_model_old_params_df['time_steps'], clock_model_old_params_df['EC'], label = 'EC, 2012', color='red')
ax[3].legend(loc='best')
ax[3].set_ylabel('Expression')
ax[3].set_xlabel('Hours')
ax[3].set_title('EC Dynamics')
for days in range(0,3):
    time_day=np.linspace(days*24+12,(days*24)+24, 100)
    ax[3].fill_between(time_day, 0, 0.16, facecolor='grey', alpha=0.5)
    
# EC original clock model
ax[3].plot(clock_outputs_df['time_steps'], clock_outputs_df['EC'], label = 'EC', color='k')
ax[3].legend(loc='best')
ax[3].set_ylabel('Expression')
ax[3].set_xlabel('Hours')
ax[3].set_title('EC Dynamics')

for item in [int(x) for x in (np.linspace(0,3,4))]:
    ax[item].set_xticks([0, 24, 48, 72])
    ax[item].set_xticklabels([0, 24, 48, 72], size=10)


# **It seems as though the EC dynamics revert to what we expected – EC expression peaks at midnight when the model is run with the 2012 param values.**
# 
# 

# <a id='final_linked_matlab_model'></a>
# ## <font color='lightseagreen'>Final Linked Matlab Model</font>

# Now that we have finalized the model, for ease of running different simulations for different mutant types and daylight hours, I have created the function `run_model` in matlab that runs the `clock linked growth model` taking *mut, daylight_hours, cE, cC, cEC* as paramters. Thus, instead of running a script and exporting the matlab workspace variables to python, we can now call the function `run_model` from python which returns the timesteps (T) and outputs (Y) from the ode23 solver for each simulation. 
# 
# In order to be able to feed integers to the matlab model from python, one first needs to convert all integer and float types into recognizable matlab types. The constants used in the model are summarized in the table `model_avgs`, and the daylight hours in the list *daylight_hours*. Finally, we can call `run_model` and save the output to a pandas dataframe for each simulation. 
# 

# ### Mutation Studies 

# In[78]:


model_avgs


# In[79]:


daylight_hours = np.linspace(0, 24, 7)
daylight_hours


# In[80]:


eng = matlab.engine.start_matlab()
path = '/Users/chloe/Desktop/thesis'
eng.cd(path, nargout=0)


# In[81]:


final_sim_outputs = {}
mutants = ['Col', 'COP1-OE', 'cop1-4', 'ELF3ox', 'elf3-8', 'phyB', 'lhy', 'toc1', 'gi']
daylengths = ['0L', '4L', '8L', '12L', '16L', '20L', '24L']


for mut in mutants:
    count = 0
    final_sim_outputs[mut] = {}
    for sim,hour in  zip(daylengths, daylight_hours): 
        constants = [float(x) for x in list(model_avgs.iloc[count, 5:])]
        #print(constants)
        final_sim_outputs[mut][sim] = eng.run_model(mut, float(hour), constants[0], constants[1], constants[2])
        count +=1


# In[82]:


clock_outputs = 'LHY mRNA,P,GI-ZTL,GI-ELF3 cytoplasm,LHY prot,TOC1 mRNA,PRR9 prot,PRR5 (NI) mRNA,PRR5 (NI) prot,GI prot cytoplasm,TOC1 prot,ZTL,EC,GI mRNA,PRR9 mRNA,PRR7 mRNA,PRR7 prot,ELF4 mRNA,ELF4 prot,LHY prot modif,ABAR mRNA,COP1 cytoplasm,ELF3 mRNA,ELF3 cytoplasm,ELF3 nuclear,COP1 nuclear night,COP1 nuclear day,LUX mRNA,LUX prot,PP2C prot,SnRK2 prot,stomata,?,?,?,PHYB,ELF3,PIF,COP1,growth'

clock_outputs = clock_outputs.split(',')
clock_outputs
header = ['time_steps']
for item in clock_outputs: 
    item = item.replace(" ", '_')
    header.append(item)


WT_0_df = create_df(final_sim_outputs['Col'], '0L')
WT_4_df = create_df(final_sim_outputs['Col'], '4L')
WT_8_df = create_df(final_sim_outputs['Col'], '8L')
WT_12_df = create_df(final_sim_outputs['Col'], '12L')
WT_16_df = create_df(final_sim_outputs['Col'], '16L')
WT_20_df = create_df(final_sim_outputs['Col'], '20L')
WT_24_df = create_df(final_sim_outputs['Col'], '24L')

COP1ox_0_df = create_df(final_sim_outputs['COP1-OE'], '0L')
COP1ox_4_df = create_df(final_sim_outputs['COP1-OE'], '4L')
COP1ox_8_df = create_df(final_sim_outputs['COP1-OE'], '8L')
COP1ox_12_df = create_df(final_sim_outputs['COP1-OE'], '12L')
COP1ox_16_df = create_df(final_sim_outputs['COP1-OE'], '16L')
COP1ox_20_df = create_df(final_sim_outputs['COP1-OE'], '20L')
COP1ox_24_df = create_df(final_sim_outputs['COP1-OE'], '24L')

ELF3ox_0_df = create_df(final_sim_outputs['ELF3ox'], '0L')
ELF3ox_4_df = create_df(final_sim_outputs['ELF3ox'], '4L')
ELF3ox_8_df = create_df(final_sim_outputs['ELF3ox'], '8L')
ELF3ox_12_df = create_df(final_sim_outputs['ELF3ox'], '12L')
ELF3ox_16_df = create_df(final_sim_outputs['ELF3ox'], '16L')
ELF3ox_20_df = create_df(final_sim_outputs['ELF3ox'], '20L')
ELF3ox_24_df = create_df(final_sim_outputs['ELF3ox'], '24L')

cop1_0_df = create_df(final_sim_outputs['cop1-4'], '0L')
cop1_4_df = create_df(final_sim_outputs['cop1-4'], '4L')
cop1_8_df = create_df(final_sim_outputs['cop1-4'], '8L')
cop1_12_df = create_df(final_sim_outputs['cop1-4'], '12L')
cop1_16_df = create_df(final_sim_outputs['cop1-4'], '16L')
cop1_20_df = create_df(final_sim_outputs['cop1-4'], '20L')
cop1_24_df = create_df(final_sim_outputs['cop1-4'], '24L')

elf3_0_df = create_df(final_sim_outputs['elf3-8'], '0L')
elf3_4_df = create_df(final_sim_outputs['elf3-8'], '4L')
elf3_8_df = create_df(final_sim_outputs['elf3-8'], '8L')
elf3_12_df = create_df(final_sim_outputs['elf3-8'], '12L')
elf3_16_df = create_df(final_sim_outputs['elf3-8'], '16L')
elf3_20_df = create_df(final_sim_outputs['elf3-8'], '20L')
elf3_24_df = create_df(final_sim_outputs['elf3-8'], '24L')

lhy_0_df = create_df(final_sim_outputs['lhy'], '0L')
lhy_4_df = create_df(final_sim_outputs['lhy'], '4L')
lhy_8_df = create_df(final_sim_outputs['lhy'], '8L')
lhy_12_df = create_df(final_sim_outputs['lhy'], '12L')
lhy_16_df = create_df(final_sim_outputs['lhy'], '16L')
lhy_20_df = create_df(final_sim_outputs['lhy'], '20L')
lhy_24_df = create_df(final_sim_outputs['lhy'], '24L')

toc1_0_df = create_df(final_sim_outputs['toc1'], '0L')
toc1_4_df = create_df(final_sim_outputs['toc1'], '4L')
toc1_8_df = create_df(final_sim_outputs['toc1'], '8L')
toc1_12_df = create_df(final_sim_outputs['toc1'], '12L')
toc1_16_df = create_df(final_sim_outputs['toc1'], '16L')
toc1_20_df = create_df(final_sim_outputs['toc1'], '20L')
toc1_24_df = create_df(final_sim_outputs['toc1'], '24L')

phyb_0_df = create_df(final_sim_outputs['phyB'], '0L')
phyb_4_df = create_df(final_sim_outputs['phyB'], '4L')
phyb_8_df = create_df(final_sim_outputs['phyB'], '8L')
phyb_12_df = create_df(final_sim_outputs['phyB'], '12L')
phyb_16_df = create_df(final_sim_outputs['phyB'], '16L')
phyb_20_df = create_df(final_sim_outputs['phyB'], '20L')
phyb_24_df = create_df(final_sim_outputs['phyB'], '24L')

gi_0_df = create_df(final_sim_outputs['gi'], '0L')
gi_4_df = create_df(final_sim_outputs['gi'], '4L')
gi_8_df = create_df(final_sim_outputs['gi'], '8L')
gi_12_df = create_df(final_sim_outputs['gi'], '12L')
gi_16_df = create_df(final_sim_outputs['gi'], '16L')
gi_20_df = create_df(final_sim_outputs['gi'], '20L')
gi_24_df = create_df(final_sim_outputs['gi'], '24L')


# In[83]:


elf3_8_df


# #### Short Day Protein Dynamics 

# In[84]:


fig, ax = plt.subplots(4, sharex=True, figsize=(10,15))
shortday_dfs = [WT_8_df, COP1ox_8_df, cop1_8_df, ELF3ox_8_df, elf3_8_df, phyb_8_df, lhy_8_df, toc1_8_df, gi_8_df]


max_value_ELF3 = 0 
max_value_COP1 = 0 
max_value_PHYB = 0 
max_value_PIF = 0 

for df, mut in zip(shortday_dfs, mutants): 
    df['ELF3_normalized'] = df['ELF3']*model_avgs['cE']['Y8']
    df['COP1_normalized'] = df['COP1']*model_avgs['cC']['Y8']


    # ELF3
    if df['ELF3_normalized'].max() > max_value_ELF3: 
        max_value_ELF3 = df['ELF3_normalized'].max()
    ax[0].plot(df['time_steps'], df['ELF3_normalized'])



    # COP1
    #if mut =='COP1-OE': 
        #pass
    #else: 
    if df['COP1_normalized'].max() > max_value_COP1: 
        max_value_COP1 = df['COP1_normalized'].max()
    ax[1].plot(df['time_steps'], df['COP1_normalized'])


    # PHYB
    if df['PHYB'].max() > max_value_PHYB: 
        max_value_PHYB = df['PHYB'].max()
    ax[2].plot(df['time_steps'], df['PHYB'])


    # PIF
    if df['PIF'].max() > max_value_PIF: 
        max_value_PIF = df['PIF'].max()
    ax[3].plot(df['time_steps'], df['PIF'])


    # ELF3
#ax[0].plot(WT_8_df['time_steps'], WT_8_df['ELF3_normalized'], label='Col', color='k',  linestyle=':')
ax[0].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['ELF322'], label = 'original model', color='k', linestyle=':')
ax[0].set_ylabel('Expression')
ax[0].set_xlabel('Hours')
ax[0].set_title('ELF3')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[0].fill_between(time_day, 0, max_value_ELF3, facecolor='grey', alpha=0.5)
#ax[0].legend(labels = mutants + ['original_model'], loc='best')

# COP1
if hypo_results_short_day_df['COP122'].max() > max_value_COP1: 
        max_value_COP1 = hypo_results_short_day_df['COP122'].max()
ax[1].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['COP122'], label = 'original model', color='k', linestyle=':')
#ax[1].legend(loc='best')
ax[1].set_ylabel('Expression')
ax[1].set_xlabel('Hours')
#ax[1].set_yscale('log')
ax[1].set_title('COP1')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[1].fill_between(time_day, 0, max_value_COP1, facecolor='grey', alpha=0.5)
#ax[1].legend(labels = mutants + ['original_model'], loc='best')

# PHYB
ax[2].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['phyb22'], label = 'original model', color='k', linestyle='-.')
#ax[2].legend(loc='best')
ax[2].set_ylabel('Expression')
ax[2].set_xlabel('Hours')
ax[2].set_title('PHYB')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[2].fill_between(time_day, 0, max_value_PHYB, facecolor='grey', alpha=0.5)
#ax[2].legend(labels = mutants + ['original_model'], loc='best')


# PIF
ax[3].plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['PIF422'], label = 'original model', color='k', linestyle=':')
#ax[3].legend(loc='best')
ax[3].set_ylabel('Expression')
ax[3].set_xlabel('Hours')
ax[3].set_title('PIF')
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    ax[3].fill_between(time_day, 0, max_value_PIF, facecolor='grey', alpha=0.5)
#ax[3].legend(labels = mutants + ['original_model'], loc='best')

for item in [int(x) for x in (np.linspace(0,3,4))]:
    ax[item].set_xticks([0, 24, 48, 72, 96, 120])
    ax[item].set_xticklabels([0, 24, 48, 72, 96, 120], size=10)

ax[0].legend(labels = mutants + ['original_model'], bbox_to_anchor = (0.3, 1.07), ncol=3)
    #fig.legend(labels = ['linked', 'original'])

plt.suptitle('Protein dynamics in short day – 8L/16D', size=20, y=1.0)
plt.tight_layout()


# In[85]:


shortday_dfs_no_COP1OE = [WT_8_df, cop1_8_df, ELF3ox_8_df,  elf3_8_df, phyb_8_df, lhy_8_df, toc1_8_df, gi_8_df]
mutants_no_COP1OE = ['Col', 'cop1-4', 'ELF3ox', 'elf3-8', 'phyB', 'lhy', 'toc1', 'gi']
max_value_COP1 = 0 
for df in shortday_dfs_no_COP1OE: 
    # COP1
    if df['COP1_normalized'].max() > max_value_COP1: 
        max_value_COP1 = df['COP1_normalized'].max()
    plt.plot(df['time_steps'], df['COP1_normalized'])
    
if hypo_results_short_day_df['COP122'].max() > max_value_COP1: 
        max_value_COP1 = hypo_results_short_day_df['COP122'].max()
plt.plot(hypo_results_short_day_df['Time'], hypo_results_short_day_df['COP122'], label = 'original model', color='k', linestyle=':')
#ax[1].legend(loc='best')
plt.ylabel('Expression')
plt.xlabel('Hours')
#ax[1].set_yscale('log')
plt.title('COP1 Dynamics', y=1.3)
for days in range(0,5):
    time_day=np.linspace(days*24+8,(days*24)+24, 100)
    plt.fill_between(time_day, 0, max_value_COP1, facecolor='grey', alpha=0.5)
plt.legend(labels = mutants_no_COP1OE + ['original_model'], bbox_to_anchor = (0.15, 1), ncol = 3)


# In[86]:


mutants


# In[87]:


mutants
mutants_no_COP1OE = ['Col', 'cop1-4', 'ELF3ox', 'elf3-8', 'phyB', 'lhy', 'toc1', 'gi']


# In[89]:


# SHORT DAY DYNAMICS 
fig=plt.figure(figsize=(18, 12))
WT_8_df['ELF3_normalized'] = WT_8_df['ELF3']*model_avgs['cE']['Y8']
WT_8_df['COP1_normalized'] = WT_8_df['COP1']*model_avgs['cC']['Y8']

proteins = ['ELF3_normalized', 'COP1_normalized', 'PHYB', 'PIF']
proteins_orig = ['ELF322', 'COP122', 'phyb22', 'PIF422']
titles=['ELF3', 'COP1', 'PHYB', 'PIF']

for i1, prot in enumerate(proteins):
    ax=fig.add_subplot(2,2,i1+1)
    max_value = 0 
    if WT_8_df[prot].max() > max_value:
        max_value = WT_8_df[prot].max()
        if hypo_results_short_day_df[proteins_orig[i1]].max() > max_value:
            max_value = hypo_results_short_day_df[proteins_orig[i1]].max()
    for days in range(0,3):
        time_day=np.linspace(days*24+8,(days*24)+24, 100)
        ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
    ax.plot(WT_8_df['time_steps'][WT_8_df['time_steps'] <= 72], WT_8_df[prot][WT_8_df['time_steps'] <=72], label = 'clock linked model', color='red')
    ax.plot(hypo_results_short_day_df['Time'][hypo_results_short_day_df['Time'] <=72], hypo_results_short_day_df[proteins_orig[i1]][hypo_results_short_day_df['Time']<=72], label = 'original model', color='k')
    ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
    ax.set_xticks([0, 24, 48, 72])
    ax.set_xticklabels([0, 24, 48, 72], size=10)
    if i1>1: 
        ax.set_xlabel('Hours', fontsize=11)
    if i1%2 == 0: 
        ax.set_ylabel('Concentration', fontsize=11)
plt.legend(loc='upper center', bbox_to_anchor=(-0.1, 1.15), ncol=4, fancybox=True, shadow=True, fontsize=12)
#plt.tight_layout() 


# In[90]:


# SHORT DAY DYNAMICS FOR ALL MUTANTS 

shortday_dfs = [WT_8_df, COP1ox_8_df, cop1_8_df, ELF3ox_8_df, elf3_8_df, phyb_8_df, lhy_8_df, toc1_8_df, gi_8_df]
proteins = ['ELF3_normalized', 'COP1_normalized', 'PHYB', 'PIF']
proteins_orig = ['ELF322', 'COP122', 'phyb22', 'PIF422']
titles=['ELF3', 'COP1', 'PHYB', 'PIF']

    
for df, mut in zip(shortday_dfs, mutants): 
    df['ELF3_normalized'] = df['ELF3']*model_avgs['cE']['Y8']
    df['COP1_normalized'] = df['COP1']*model_avgs['cC']['Y8']

    fig=plt.figure(figsize=(18, 12))
    if mut == 'Col': 
        plt.suptitle(f'Short Day Protein Dynamics in {mut} lines ', size=18, y=0.93)
    else: 
        plt.suptitle(f'Short Day Protein Dynamics in  $\it{mut}$ mutant lines ', size=18, y=0.93)

    for i1, prot in enumerate(proteins):
        ax=fig.add_subplot(2,2,i1+1)
        max_value = np.max([WT_8_df[prot].max(),df[prot].max(), hypo_results_short_day_df[proteins_orig[i1]].max()])
        for days in range(0,3):
            time_day=np.linspace(days*24+8,(days*24)+24, 100)
            ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
        ax.plot(WT_8_df['time_steps'][WT_8_df['time_steps'] <= 72], WT_8_df[prot][WT_8_df['time_steps'] <=72], label = 'WT', color='k')
        ax.plot(df['time_steps'][df['time_steps'] <= 72], df[prot][df['time_steps'] <=72], label = mut, color='red')
        ax.plot(hypo_results_short_day_df['Time'][hypo_results_short_day_df['Time'] <=72], hypo_results_short_day_df[proteins_orig[i1]][hypo_results_short_day_df['Time']<=72], label = 'original model', color='k', linestyle=':')
        ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
        ax.set_xticks([0, 24, 48, 72])
        ax.set_xticklabels([0, 24, 48, 72], size=10)
        if i1>1: 
            ax.set_xlabel('Hours', fontsize=11)
        if i1%2 == 0: 
            ax.set_ylabel('Concentration', fontsize=11)
    plt.legend(loc='upper center', bbox_to_anchor=(-0.1, 1.15), ncol=4, fancybox=True, shadow=True, fontsize=12)
    #plt.tight_layout() 


# #### 12L/12D Protein Dynamics 

# In[91]:


WT_8_df.columns


# In[92]:


# 12L/12D DYNAMICS FOR ALL MUTANTS 

dfs = [WT_12_df, COP1ox_12_df, cop1_12_df, ELF3ox_12_df, elf3_12_df, phyb_12_df, lhy_12_df, toc1_12_df, gi_12_df]
proteins = ['ELF3_normalized', 'COP1_normalized', 'PHYB', 'PIF', 'LUX_prot', 'ELF4_prot', 'GI_prot_cytoplasm', 'TOC1_prot']
proteins_orig = ['ELF322', 'COP122', 'phyb22', 'PIF422']
titles=['ELF3', 'COP1', 'PHYB', 'PIF', 'LUX', 'ELF4', 'GI', 'TOC1']

    
for df, mut in zip(dfs, mutants): 
    df['ELF3_normalized'] = df['ELF3']*model_avgs['cE']['Y12']
    df['COP1_normalized'] = df['COP1']*model_avgs['cC']['Y12']

    fig=plt.figure(figsize=(18, 12))
    if mut == 'Col': 
        plt.suptitle(f'12L/12D Protein Dyanmics in {mut} lines ', size=18, y=0.93)
    else: 
        plt.suptitle(f'12L/12D Protein Dyanmics in $\it{mut}$ mutant lines ', size=18, y=0.93)
    for i1, prot in enumerate(proteins):
        ax=fig.add_subplot(4,2,i1+1)
        max_value = np.max([WT_12_df[prot].max(),df[prot].max()])
        for days in range(0,3):
            time_day=np.linspace(days*24+12,(days*24)+24, 100)
            ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
        ax.plot(WT_12_df['time_steps'][WT_12_df['time_steps'] <= 72], WT_12_df[prot][WT_12_df['time_steps'] <=72], label = 'WT', color='k')
        ax.plot(df['time_steps'][df['time_steps'] <= 72], df[prot][df['time_steps'] <=72], label = mut, color='red')
        ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
        ax.set_xticks([0, 24, 48, 72])
        ax.set_xticklabels([0, 24, 48, 72], size=10)
        if i1>5: 
            ax.set_xlabel('Hours', fontsize=11)
        if i1%2 == 0: 
            ax.set_ylabel('Concentration', fontsize=11)
    plt.legend(loc='upper center', ncol=4, fancybox=True, shadow=True, fontsize=12)
    #plt.tight_layout() 


# #### Long Day Protein Dynamics  

# In[94]:


# LONG DAY PROTEIN DYNAMICS, WT
fig=plt.figure(figsize=(18, 12))
WT_16_df['ELF3_normalized'] = WT_16_df['ELF3']*model_avgs['cE']['Y16']
WT_16_df['COP1_normalized'] = WT_16_df['COP1']*model_avgs['cC']['Y16']

proteins = ['ELF3_normalized', 'COP1_normalized', 'PHYB', 'PIF']
proteins_orig = ['ELF322', 'COP122', 'phyb22', 'PIF422']
titles=['ELF3', 'COP1', 'PHYB', 'PIF']

for i1, prot in enumerate(proteins):
    ax=fig.add_subplot(2,2,i1+1)
    max_value = 0 
    if WT_16_df[prot].max() > max_value:
        max_value = WT_16_df[prot].max()
        if hypo_results_long_day_df[proteins_orig[i1]].max() > max_value:
            max_value = hypo_results_long_day_df[proteins_orig[i1]].max()
    for days in range(0,3):
        time_day=np.linspace(days*24+16,(days*24)+24, 100)
        ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
    ax.plot(WT_16_df['time_steps'][WT_16_df['time_steps'] <= 72], WT_16_df[prot][WT_16_df['time_steps'] <=72], label = 'clock linked model', color='red')
    ax.plot(hypo_results_long_day_df['Time'][hypo_results_long_day_df['Time'] <=72], hypo_results_long_day_df[proteins_orig[i1]][hypo_results_long_day_df['Time']<=72], label = 'original model', color='k')
    ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
    ax.set_xticks([0, 24, 48, 72])
    ax.set_xticklabels([0, 24, 48, 72], size=10)
    if i1>1: 
        ax.set_xlabel('Hours', fontsize=11)
    if i1%2 == 0: 
        ax.set_ylabel('Concentration', fontsize=11)
plt.legend(loc='upper center', bbox_to_anchor=(0.0, 2.4), ncol=4, fancybox=True, shadow=True, fontsize=12)
#plt.tight_layout() 


# In[95]:


# LONG DAY DYNAMICS FOR ALL MUTANTS, plotted together
fig, ax = plt.subplots(4, sharex=True, figsize=(10,15))
longday_dfs = [WT_16_df, COP1ox_16_df, cop1_16_df, ELF3ox_16_df,  elf3_16_df, phyb_16_df, lhy_16_df, toc1_16_df, gi_16_df]

max_value_ELF3 = 0 
max_value_COP1 = 0 
max_value_PHYB = 0 
max_value_PIF = 0 

for df in longday_dfs: 
    df['ELF3_normalized'] = df['ELF3']*model_avgs['cE']['Y16']
    df['COP1_normalized'] = df['COP1']*model_avgs['cC']['Y16']

    # ELF3
    if df['ELF3_normalized'].max() > max_value_ELF3: 
        max_value_ELF3 = df['ELF3_normalized'].max()
    ax[0].plot(df['time_steps'], df['ELF3_normalized'])


    # COP1
    if df['COP1_normalized'].max() > max_value_COP1: 
        max_value_COP1 = df['COP1_normalized'].max()
    ax[1].plot(df['time_steps'], df['COP1_normalized'])


    # PHYB
    if df['PHYB'].max() > max_value_PHYB: 
        max_value_PHYB = df['PHYB'].max()
    ax[2].plot(df['time_steps'], df['PHYB'])


    # PIF
    if df['PIF'].max() > max_value_PIF: 
        max_value_PIF = df['PIF'].max()
    ax[3].plot(df['time_steps'], df['PIF'])


    # ELF3
ax[0].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['ELF322'], label = 'original model', color='k')
ax[0].set_ylabel('Expression')
ax[0].set_xlabel('Hours')
ax[0].set_title('ELF3')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[0].fill_between(time_day, 0, max_value_ELF3, facecolor='grey', alpha=0.5)
ax[0].legend(labels = mutants + ['original_model'], loc='best')

# COP1
ax[1].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['COP122'], label = 'original model', color='k')
#ax[1].legend(loc='best')
ax[1].set_ylabel('Expression')
ax[1].set_xlabel('Hours')
#ax[1].set_yscale('log')
ax[1].set_title('COP1')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[1].fill_between(time_day, 0, max_value_COP1, facecolor='grey', alpha=0.5)
ax[1].legend(labels = mutants + ['original_model'], loc='best')

# PHYB
ax[2].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['phyb22'], label = 'original model', color='k')
#ax[2].legend(loc='best')
ax[2].set_ylabel('Expression')
ax[2].set_xlabel('Hours')
ax[2].set_title('PHYB')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[2].fill_between(time_day, 0, max_value_PHYB, facecolor='grey', alpha=0.5)
ax[2].legend(labels = mutants + ['original_model'], loc='best')


# PIF
ax[3].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['PIF422'], label = 'original model', color='k')
#ax[3].legend(loc='best')
ax[3].set_ylabel('Expression')
ax[3].set_xlabel('Hours')
ax[3].set_title('PIF')
for days in range(0,5):
    time_day=np.linspace(days*24+16,(days*24)+24, 100)
    ax[3].fill_between(time_day, 0, max_value_PIF, facecolor='grey', alpha=0.5)
ax[3].legend(labels = mutants + ['original_model'], loc='best')

for item in [int(x) for x in (np.linspace(0,3,4))]:
    ax[item].set_xticks([0, 24, 48, 72, 96, 120])
    ax[item].set_xticklabels([0, 24, 48, 72, 96, 120], size=10)


    #fig.legend(labels = ['linked', 'original'])

plt.suptitle('Protein dynamics in long day – 16L/8D', size=20)
plt.tight_layout()


# In[96]:


# LONG DAY DYNAMICS FOR ALL MUTANTS, individual graphs

longday_dfs = [WT_16_df, COP1ox_16_df, cop1_16_df, ELF3ox_16_df, elf3_16_df, phyb_16_df, lhy_16_df, toc1_16_df, gi_16_df]
proteins = ['ELF3_normalized', 'COP1_normalized', 'PHYB', 'PIF']
proteins_orig = ['ELF322', 'COP122', 'phyb22', 'PIF422']
titles=['ELF3', 'COP1', 'PHYB', 'PIF']

    
for df, mut in zip(longday_dfs, mutants): 
    df['ELF3_normalized'] = df['ELF3']*model_avgs['cE']['Y16']
    df['COP1_normalized'] = df['COP1']*model_avgs['cC']['Y16']

    fig=plt.figure(figsize=(18, 12))
    if mut == 'Col': 
        plt.suptitle(f'Long Day Protein Dynamics in {mut} lines ', size=18, y=0.93)
    else: 
        plt.suptitle(f'Long Day Protein Dynamics in  $\it{mut}$ mutant lines ', size=18, y=0.93)

    for i1, prot in enumerate(proteins):
        ax=fig.add_subplot(2,2,i1+1)
        max_value = np.max([WT_16_df[prot].max(),df[prot].max(), hypo_results_long_day_df[proteins_orig[i1]].max()])
        for days in range(0,3):
            time_day=np.linspace(days*24+16,(days*24)+24, 100)
            ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
        ax.plot(WT_16_df['time_steps'][WT_16_df['time_steps'] <= 72], WT_16_df[prot][WT_16_df['time_steps'] <=72], label = 'WT', color='k')
        ax.plot(df['time_steps'][df['time_steps'] <= 72], df[prot][df['time_steps'] <=72], label = mut, color='red')
        ax.plot(hypo_results_long_day_df['Time'][hypo_results_long_day_df['Time'] <=72], hypo_results_long_day_df[proteins_orig[i1]][hypo_results_short_day_df['Time']<=72], label = 'original model', color='k', linestyle=':')
        ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
        ax.set_xticks([0, 24, 48, 72])
        ax.set_xticklabels([0, 24, 48, 72], size=10)
        if i1>1: 
            ax.set_xlabel('Hours', fontsize=11)
        if i1%2 == 0: 
            ax.set_ylabel('Concentration', fontsize=11)
    plt.legend(loc='upper center', bbox_to_anchor=(-0.1, 1.15), ncol=4, fancybox=True, shadow=True, fontsize=12)
    #plt.tight_layout() 


# In[97]:


# Long day Protein DYNAMICS FOR ALL MUTANTS 

longday_dfs = [WT_16_df, COP1ox_16_df, cop1_16_df, ELF3ox_16_df, elf3_16_df, phyb_16_df, lhy_16_df, toc1_16_df, gi_16_df]
proteins = ['ELF3_normalized', 'COP1_normalized', 'PHYB', 'PIF', 'LUX_prot', 'LHY_prot', 'GI_prot_cytoplasm', 'TOC1_prot']
proteins_orig = ['ELF322', 'COP122', 'phyb22', 'PIF422']
titles=['ELF3', 'COP1', 'PHYB', 'PIF', 'LUX', 'LHY/CCA1', 'GI', 'TOC1']

    
for df, mut in zip(longday_dfs, mutants): 
    df['ELF3_normalized'] = df['ELF3']*model_avgs['cE']['Y16']
    df['COP1_normalized'] = df['COP1']*model_avgs['cC']['Y16']

    fig=plt.figure(figsize=(18, 12))
    if mut == 'Col': 
        plt.suptitle(f'Long Day Protein Dynamics in {mut} lines ', size=18, y=0.95)
    else: 
        plt.suptitle(f'Long Day Protein Dynamics in $\it{mut}$ mutant lines ', size=18, y=0.95)

    for i1, prot in enumerate(proteins):
        ax=fig.add_subplot(4,2,i1+1)
        if i1 <= 3: 
            max_value = np.max([WT_16_df[prot].max(),df[prot].max(), hypo_results_long_day_df[proteins_orig[i1]].max()])
            ax.plot(hypo_results_long_day_df['Time'][hypo_results_long_day_df['Time'] <=72], hypo_results_long_day_df[proteins_orig[i1]][hypo_results_short_day_df['Time']<=72], label = 'original model', color='k', linestyle=':')
            if i1 == 0: 
                ax.legend(bbox_to_anchor=(1.24, 0))
        else: 
            max_value = np.max([WT_16_df[prot].max(),df[prot].max()])
        for days in range(0,3):
            time_day=np.linspace(days*24+16,(days*24)+24, 100)
            ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
        ax.plot(WT_16_df['time_steps'][WT_16_df['time_steps'] <= 72], WT_16_df[prot][WT_16_df['time_steps'] <=72], label = 'WT', color='k')
        ax.plot(df['time_steps'][df['time_steps'] <= 72], df[prot][df['time_steps'] <=72], label = mut, color='red')
        ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
        ax.set_xticks([0, 24, 48, 72])
        ax.set_xticklabels([0, 24, 48, 72], size=10)
        if i1>5: 
            ax.set_xlabel('Hours', fontsize=11)
        if i1%2 == 0: 
            ax.set_ylabel('Concentration', fontsize=11)
    plt.legend(loc='upper center', bbox_to_anchor=(-0.1, 4.9), ncol=4, fancybox=True, shadow=True, fontsize=12)
    #plt.tight_layout() 


# ##### COP1 Mutant Dynamics in Long Day 

# In[98]:


# Long day Protein DYNAMICS FOR ALL MUTANTS 

longday_dfs = [WT_16_df, COP1ox_16_df, cop1_16_df]
proteins = ['ELF3_normalized', 'COP1_normalized', 'PHYB', 'PIF', 'LUX_prot', 'LHY_prot', 'GI_prot_cytoplasm', 'TOC1_prot']
proteins_orig = ['ELF322', 'COP122', 'phyb22', 'PIF422']
titles=['ELF3', 'COP1', 'PHYB', 'PIF', 'LUX', 'LHY/CCA1', 'GI', 'TOC1']
fig=plt.figure(figsize=(18, 12))

    
        

for i1, prot in enumerate(proteins):
    ax=fig.add_subplot(4,2, i1+1)
    if i1 <= 3: 
        max_value = np.max([WT_16_df[prot].max(),COP1ox_16_df[prot].max(),cop1_16_df[prot].max(), hypo_results_long_day_df[proteins_orig[i1]].max()])
        ax.plot(hypo_results_long_day_df['Time'][hypo_results_long_day_df['Time'] <=72], hypo_results_long_day_df[proteins_orig[i1]][hypo_results_short_day_df['Time']<=72], label = 'original model', color='k', linestyle=':')
        if i1 == 0: 
            ax.legend(bbox_to_anchor=(1.24, 0))
    else: 
        max_value = np.max([WT_16_df[prot].max(),COP1ox_16_df[prot].max(),cop1_16_df[prot].max()])
    for days in range(0,3):
        time_day=np.linspace(days*24+16,(days*24)+24, 100)
        ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
    ax.plot(WT_16_df['time_steps'][WT_16_df['time_steps'] <= 72], WT_16_df[prot][WT_16_df['time_steps'] <=72], label = 'WT', color='k')
    ax.plot(COP1ox_16_df['time_steps'][COP1ox_16_df['time_steps'] <= 72], COP1ox_16_df[prot][COP1ox_16_df['time_steps'] <=72], label = 'COP1-OE', color='red')
    ax.plot(cop1_16_df['time_steps'][cop1_16_df['time_steps'] <= 72], cop1_16_df[prot][cop1_16_df['time_steps'] <=72], label = 'cop1-4', color='red', linestyle='--')
    ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
    ax.set_xticks([0, 24, 48, 72])
    ax.set_xticklabels([0, 24, 48, 72], size=10)
    if i1>5: 
        ax.set_xlabel('Hours', fontsize=11)
    if i1%2 == 0: 
        ax.set_ylabel('Concentration', fontsize=11)
plt.legend(loc='upper center', bbox_to_anchor=(-0.1, 4.9), ncol=4, fancybox=True, shadow=True, fontsize=12)
plt.suptitle(f'Long Day Protein Dynamics in $\itcop1$ mutant lines ', size=18, y=0.95)
#plt.tight_layout() 


# ##### gi Mutant Dynamics in Long Day 

# In[157]:


# Long day Protein DYNAMICS FOR gi mutant

longday_dfs = [WT_16_df, gi_16_df]
proteins = ['ELF3_normalized', 'COP1_normalized', 'PIF', 'LHY_prot' ,'LUX_prot', 'EC', 'GI_prot_cytoplasm', 'TOC1_prot']
proteins_orig = ['ELF322', 'COP122', 'PIF422']
titles=['ELF3', 'COP1','PIF', 'LHY/CCA1', 'LUX', 'EC', 'GI', 'TOC1']
fig=plt.figure(figsize=(18, 12))

    
        

for i1, prot in enumerate(proteins):
    ax=fig.add_subplot(5,2, i1+1)
    if i1 <= 2: 
        max_value = np.max([WT_16_df[prot].max(),gi_16_df[prot].max(),gi_16_df[prot].max(), hypo_results_long_day_df[proteins_orig[i1]].max()])
        ax.plot(hypo_results_long_day_df['Time'][hypo_results_long_day_df['Time'] <=72], hypo_results_long_day_df[proteins_orig[i1]][hypo_results_short_day_df['Time']<=72], label = 'original model', color='k', linestyle=':')
        if i1 == 0: 
            ax.legend(bbox_to_anchor=(1.24, 0))
    else: 
        max_value = np.max([WT_16_df[prot].max(),gi_16_df[prot].max(),gi_16_df[prot].max()])
    for days in range(0,3):
        time_day=np.linspace(days*24+16,(days*24)+24, 100)
        ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
    ax.plot(WT_16_df['time_steps'][WT_16_df['time_steps'] <= 72], WT_16_df[prot][WT_16_df['time_steps'] <=72], label = 'WT', color='k')
    ax.plot(gi_16_df['time_steps'][gi_16_df['time_steps'] <= 72], gi_16_df[prot][gi_16_df['time_steps'] <=72], label = '$\itgi$', color='red')
    ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
    ax.set_xticks([0, 24, 48, 72])
    ax.set_xticklabels([0, 24, 48, 72], size=10)
    if i1>5: 
        ax.set_xlabel('Hours', fontsize=11)
    if i1%2 == 0: 
        ax.set_ylabel('Concentration', fontsize=11)
plt.legend(loc='upper center', bbox_to_anchor=(-0.1, 5), ncol=4, fancybox=True, shadow=True, fontsize=12)
plt.suptitle(f'Long Day Protein Dynamics in gi mutant lines ', size=18, y=0.95)
#plt.tight_layout() 


# ##### Long day Protein dynamics for lhy mutant

# In[100]:


# Long day Protein DYNAMICS FOR lhy mutant

longday_dfs = [WT_16_df, lhy_16_df]
proteins = ['ELF3_normalized', 'PIF', 'LHY_prot', 'ELF4_prot', 'LUX_prot', 'EC', 'GI_prot_cytoplasm', 'TOC1_prot']
proteins_orig = ['ELF322', 'PIF422']
titles=['ELF3','PIF', 'LHY/CCA1', 'ELF4', 'LUX', 'EC', 'GI', 'TOC1']
fig=plt.figure(figsize=(18, 12))

    
        

for i1, prot in enumerate(proteins):
    ax=fig.add_subplot(4,2, i1+1)
    if i1 <= 1: 
        max_value = np.max([WT_16_df[prot].max(),lhy_16_df[prot].max(),lhy_16_df[prot].max(), hypo_results_long_day_df[proteins_orig[i1]].max()])
        ax.plot(hypo_results_long_day_df['Time'][hypo_results_long_day_df['Time'] <=72], hypo_results_long_day_df[proteins_orig[i1]][hypo_results_short_day_df['Time']<=72], label = 'original model', color='k', linestyle=':')
        if i1 == 0: 
            ax.legend(bbox_to_anchor=(1.24, 0))
    else: 
        max_value = np.max([WT_16_df[prot].max(),lhy_16_df[prot].max(),lhy_16_df[prot].max()])
    for days in range(0,3):
        time_day=np.linspace(days*24+16,(days*24)+24, 100)
        ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
    ax.plot(WT_16_df['time_steps'][WT_16_df['time_steps'] <= 72], WT_16_df[prot][WT_16_df['time_steps'] <=72], label = 'WT', color='k')
    ax.plot(lhy_16_df['time_steps'][lhy_16_df['time_steps'] <= 72], lhy_16_df[prot][lhy_16_df['time_steps'] <=72], label = 'lhy', color='red')
    ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
    ax.set_xticks([0, 24, 48, 72])
    ax.set_xticklabels([0, 24, 48, 72], size=10)
    if i1>5: 
        ax.set_xlabel('Hours', fontsize=11)
    if i1%2 == 0: 
        ax.set_ylabel('Concentration', fontsize=11)
plt.legend(loc='upper center', bbox_to_anchor=(-0.1, 4.9), ncol=4, fancybox=True, shadow=True, fontsize=12)
plt.suptitle(f'Long Day Protein Dynamics in $\itlhy$ mutant lines ', size=18, y=0.95)
#plt.tight_layout() 


# #### **PIF4 Expression in long and short days** 

# In[101]:


# PIF4 dynamics in short and long days 
fig=plt.figure(figsize=(18, 12))

pif_dfs = [WT_8_df, WT_16_df]
proteins = ['PIF']
proteins_orig = ['PIF422']
titles=['PIF']

for i1, df in enumerate(pif_dfs):
    ax=fig.add_subplot(2,2,i1+1)
    max_value = np.max([WT_8_df['PIF'].max(), WT_16_df['PIF'].max(), hypo_results_short_day_df['PIF422'].max(), hypo_results_long_day_df['PIF422'].max()])
    if i1 ==0: 
        max_value = np.max([WT_8_df['PIF'].max(), WT_16_df['PIF'].max(), hypo_results_short_day_df['PIF422'].max()])
        for days in range(0,3):
            time_day=np.linspace(days*24+8,(days*24)+24, 100)
            ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
        ax.plot(hypo_results_short_day_df['Time'][hypo_results_short_day_df['Time'] <=72], hypo_results_short_day_df[proteins_orig[0]][hypo_results_short_day_df['Time']<=72], label = 'original model', color='k')
        ax.plot(df['time_steps'][df['time_steps'] <= 72], df['PIF'][df['time_steps'] <=72], label = 'clock linked model', color='red')
        ax.set_title(f'{titles[0]}, 8L/16D', fontsize=14, fontweight="bold")
        ax.text(-7, 7, 'A', fontweight="bold", fontsize=18)
    if i1==1: 
        max_value = np.max([WT_16_df['PIF'].max(), hypo_results_long_day_df['PIF422'].max()])
        for days in range(0,3):
            time_day=np.linspace(days*24+16,(days*24)+24, 100)
            ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
        ax.plot(hypo_results_long_day_df['Time'][hypo_results_long_day_df['Time'] <=72], hypo_results_long_day_df[proteins_orig[0]][hypo_results_long_day_df['Time']<=72], label = 'original model', color='k')
        ax.plot(df['time_steps'][df['time_steps'] <= 72], df['PIF'][df['time_steps'] <=72], label = 'clock linked model', color='red')
        ax.set_title(f'{titles[0]}, 16L/8D', fontsize=12, fontweight="bold")
        ax.text(-7, 3.2, 'B', fontweight="bold", fontsize=18)
    ax.set_xticks([0, 24, 48, 72])
    ax.set_xticklabels([0, 24, 48, 72], size=14)
    ax.set_xlabel('Hours', fontsize=11)
    
    if i1%2 == 0: 
        ax.set_ylabel('Concentration', fontsize=11)
plt.legend(loc='upper center', bbox_to_anchor=(-0.1, 1.2), ncol=4, fancybox=True, shadow=True, fontsize=12)
#plt.tight_layout() 


# In[144]:


# PIF4 dynamics in long days, elf3-8 mutant 
fig=plt.figure(figsize=(18, 12))

pif_dfs = [elf3_16_df]
proteins = ['PIF']
proteins_orig = ['PIF422']
titles=['PIF']

for i1, df in enumerate(pif_dfs):
    ax=fig.add_subplot(2,2,i1+1)
    max_value = np.max([elf3_16_df['PIF'].max(), WT_16_df['PIF'].max(), hypo_results_long_day_df['PIF422'].max()])
    for days in range(0,3):
        time_day=np.linspace(days*24+16,(days*24)+24, 100)
        ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
    ax.plot(df['time_steps'][df['time_steps'] <= 72], df['PIF'][df['time_steps'] <=72], label = '$\itelf3-8$', color='red')
    ax.plot(WT_16_df['time_steps'][WT_16_df['time_steps'] <= 72], WT_16_df['PIF'][WT_16_df['time_steps'] <=72], label = 'WT', color='black')
    ax.set_title(f'{titles[0]}, 16L/8D', fontsize=12, fontweight="bold")
    #ax.text(-7, 3.2, 'B', fontweight="bold", fontsize=18)
    ax.set_xticks([0, 24, 48, 72])
    ax.set_xticklabels([0, 24, 48, 72], size=14)
    ax.set_xlabel('Hours', fontsize=11)
    ax.set_ylabel('Concentration', fontsize=11)
plt.legend()    
#plt.legend(loc='upper center', bbox_to_anchor=(-0.1, 1.2), ncol=4, fancybox=True, shadow=True, fontsize=12)
#plt.tight_layout() 


# #### **EC levels in long and short days** 

# In[102]:


# PIF4 dynamics in short and long days 
fig=plt.figure(figsize=(18, 12))

pif_dfs = [WT_8_df, WT_16_df]
proteins = ['EC']
titles = ['EC']

for i1, df in enumerate(pif_dfs):
    ax=fig.add_subplot(2,2,i1+1)
    if i1 ==0: 
        max_value = np.max([WT_8_df['EC'].max()])
        for days in range(0,3):
            time_day=np.linspace(days*24+8,(days*24)+24, 100)
            ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
        ax.plot(df['time_steps'][df['time_steps'] <= 72], df['EC'][df['time_steps'] <=72], label = 'clock linked model', color='red')
        ax.set_title(f'{titles[0]}, 8L/16D', fontsize=14, fontweight="bold")
        ax.text(-7, 0.155, 'C', fontweight="bold", fontsize=18)
    if i1==1: 
        max_value = np.max([WT_16_df['EC'].max()])
        for days in range(0,3):
            time_day=np.linspace(days*24+16,(days*24)+24, 100)
            ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.5)
        ax.plot(df['time_steps'][df['time_steps'] <= 72], df['EC'][df['time_steps'] <=72], label = 'clock linked model', color='red')
        ax.set_title(f'{titles[0]}, 16L/8D', fontsize=12, fontweight="bold")
        ax.text(-7, 0.19, 'D', fontweight="bold", fontsize=18)
    ax.set_xticks([0, 24, 48, 72])
    ax.set_xticklabels([0, 24, 48, 72], size=14)
    ax.set_xlabel('Hours', fontsize=11)
    
    if i1%2 == 0: 
        ax.set_ylabel('Concentration', fontsize=11)
#plt.legend(loc='upper center', bbox_to_anchor=(0, 1.25), ncol=4, fancybox=True, shadow=True, fontsize=12)
#plt.tight_layout() 


# #### Constant light Mutant Protein Dynamics 

# In[103]:


# CONSTANT LIGHT AND DARK DYNAMICS FOR ALL MUTANTS 

dfs = [WT_24_df, COP1ox_24_df, cop1_24_df, ELF3ox_24_df, elf3_24_df, phyb_24_df, lhy_24_df, toc1_24_df, gi_24_df]
dark_dfs = [WT_0_df, COP1ox_0_df, cop1_0_df, ELF3ox_0_df, elf3_0_df, phyb_0_df, lhy_0_df, toc1_0_df, gi_0_df]
proteins = ['ELF3_normalized', 'COP1_normalized', 'PHYB', 'PIF', 'LUX_prot', 'LHY_prot', 'GI_prot_cytoplasm', 'TOC1_prot']
proteins_orig = ['ELF322', 'COP122', 'phyb22', 'PIF422']
titles=['ELF3', 'COP1', 'PHYB', 'PIF', 'LUX', 'LHY/CCA1', 'GI', 'TOC1']
count = 0
    
for df, mut in zip(dfs, mutants): 
    df['ELF3_normalized'] = df['ELF3']*model_avgs['cE']['Y24']
    df['COP1_normalized'] = df['COP1']*model_avgs['cC']['Y24']
    dark_dfs[count]['ELF3_normalized'] = dark_dfs[count]['ELF3']*model_avgs['cE']['Y0']
    dark_dfs[count]['COP1_normalized'] = dark_dfs[count]['COP1']*model_avgs['cC']['Y0']

    fig=plt.figure(figsize=(18, 12))
    if mut == 'Col': 
        plt.suptitle(f'Constant Light Protein Dynamics in {mut} lines ', size=18, y=0.95)
    else: 
        plt.suptitle(f'Constant Light Protein Dynamics in $\it{mut}$ mutant lines ', size=18, y=0.95)

    for i1, prot in enumerate(proteins):
        ax=fig.add_subplot(4,2,i1+1)
        if i1 <= 3: 
            max_value = np.max([WT_24_df[prot].max(),df[prot].max()])
        else: 
            max_value = np.max([WT_24_df[prot].max(),df[prot].max()])
        ax.plot(WT_24_df['time_steps'][WT_24_df['time_steps'] <= 72], WT_24_df[prot][WT_24_df['time_steps'] <=72], label = 'WT_LL', color='k')
        ax.plot(WT_0_df['time_steps'][WT_0_df['time_steps'] <= 72], WT_0_df[prot][WT_0_df['time_steps'] <=72], label = 'WT_DD', color='k', linestyle = '--')
        ax.plot(df['time_steps'][df['time_steps'] <= 72], df[prot][df['time_steps'] <=72], label = mut+', LL', color='red')
        ax.plot(dark_dfs[count]['time_steps'][dark_dfs[count]['time_steps'] <= 72], dark_dfs[count][prot][dark_dfs[count]['time_steps'] <=72], label = mut+', DD', color='red', linestyle = '--')
        ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
        ax.set_xticks([0, 24, 48, 72])
        ax.set_xticklabels([0, 24, 48, 72], size=10)
        if i1>5: 
            ax.set_xlabel('Hours', fontsize=11)
        if i1%2 == 0: 
            ax.set_ylabel('Concentration', fontsize=11)
    count +=1
    plt.legend(loc='upper center', bbox_to_anchor=(-0.1, 4.9), ncol=4, fancybox=True, shadow=True, fontsize=12)
    #plt.tight_layout() 


# In[ ]:





# #### COP1 PROTEIN DYNAMICS 

# In[104]:


WT_12_df.columns


# In[105]:


# 12L/12D dynamics 
n = 3
colors = plt.cm.inferno(np.linspace(0,0.8,n))
fig=plt.figure(figsize=(10, 6))
ax=fig.add_subplot()

proteins = ['COP1_cytoplasm', 'COP1_nuclear_night', 'COP1_nuclear_day' ]
titles=['cystosolic COP1', 'COP1n', 'COP1d']

for i1, prot in enumerate(proteins):
    max_value = np.max([WT_12_df['COP1_cytoplasm'].max(), WT_12_df['COP1_nuclear_night'].max(), WT_12_df['COP1_nuclear_day'].max()])
    for days in range(0,3):
        time_day=np.linspace(days*24+12,(days*24)+24, 100)
        ax.fill_between(time_day, 0, max_value, facecolor='grey', alpha=0.3)
    ax.plot(WT_12_df['time_steps'][WT_12_df['time_steps'] <= 72], WT_12_df[prot][WT_12_df['time_steps'] <=72], label = titles[i1], color=colors[i1], linewidth=2)
    ax.set_xticks([0, 24, 48, 72])
    ax.set_xticklabels([0, 24, 48, 72], size=10)
    if i1>1: 
        ax.set_xlabel('Hours', fontsize=11)
    if i1%2 == 0: 
        ax.set_ylabel('Concentration', fontsize=11)
plt.legend(loc='upper center', bbox_to_anchor = (0.5, 1.1), ncol=4, fancybox=True, shadow=True, fontsize=12)
#plt.tight_layout() 


# **Total Abundance of ELF3 protein is differently regulated by COP1 activity in constant light and dark conditions**

# In[106]:


WT_0_df.columns


# In[107]:


n = 2
colors = plt.cm.inferno(np.linspace(0,0.8,n))
styles = '--', '-'


light_cond = [WT_0_df, WT_24_df]


fig=plt.figure(figsize=(10, 8))
ax=fig.add_subplot()

proteins = ['COP1_nuclear_night', 'ELF3_nuclear']
labels = ['COP1n, DD', 'COP1n, LL', 'ELF3, LL', 'ELF3, DD']
plot_labels = ['A', 'B', 'C', 'D']
count = 0 
for i1, prot in enumerate(proteins):
    for i2, cond in enumerate(light_cond): 
        ax.plot(cond['time_steps'][cond['time_steps'] <=72], cond[prot][cond['time_steps'] <=72], color=colors[i1], label=labels[count], linestyle = styles[i2])
        count +=1
    ax.set_xticks([0, 24, 48, 72])
    ax.set_xticklabels([0, 24, 48, 72], size=10)

    ax.set_xlabel('Hours', fontsize=11)
    ax.set_ylabel('Expression', fontsize=11)

    
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=4, fancybox=True, shadow=True, fontsize=12)
#plt.tight_layout()
        


# #### WT Protein Dynamics in different Light Conditions

# In[108]:


fig, ax = plt.subplots(4, sharex=True, figsize=(10,15))
WT_0_df['ELF3_normalized'] = WT_0_df['ELF3']*model_avgs['cE']['Y0']
WT_0_df['COP1_normalized'] = WT_0_df['COP1']*model_avgs['cC']['Y0']
WT_24_df['ELF3_normalized'] = WT_24_df['ELF3']*model_avgs['cE']['Y24']
WT_24_df['COP1_normalized'] = WT_24_df['COP1']*model_avgs['cC']['Y24']

# ELF3
ax[0].plot(WT_0_df['time_steps'], WT_0_df['ELF3_normalized'], label = 'constant darkness', color='purple')
ax[0].plot(WT_24_df['time_steps'], WT_24_df['ELF3_normalized'], label = 'constant light', color='red')
ax[0].plot(WT_8_df['time_steps'], WT_8_df['ELF3_normalized'], label = 'short day ', color='k', linestyle=':')
ax[0].plot(WT_16_df['time_steps'], WT_16_df['ELF3_normalized'], label = 'long day ', color='blue', linestyle=':')
#ax[0].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['ELF322'], label = 'original model', color='k')
ax[0].set_ylabel('Expression')
ax[0].set_xlabel('Hours')
ax[0].set_title('ELF3')

#ax[0].legend(labels = ['clock linked', 'original'], loc='best')

# COP1
ax[1].plot(WT_0_df['time_steps'], WT_0_df['COP1_normalized'], label = 'constant darkness', color='purple')
ax[1].plot(WT_24_df['time_steps'], WT_24_df['COP1_normalized'], label = 'constant light', color='red')
ax[1].plot(WT_8_df['time_steps'], WT_8_df['COP1_normalized'], label = 'short day ', color='k', linestyle=':')
ax[1].plot(WT_16_df['time_steps'], WT_16_df['COP1_normalized'], label = 'long day ', color='blue', linestyle=':')
#ax[1].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['COP122'], label = 'original model', color='k')
#ax[1].legend(loc='best')
ax[1].set_ylabel('Expression')
ax[1].set_xlabel('Hours')
#ax[1].set_yscale('log')
ax[1].set_title('COP1')

#ax[1].legend(labels = ['linked', 'linked_normalized', 'original'], loc='best')

# PHYB
ax[2].plot(WT_0_df['time_steps'], WT_0_df['PHYB'], label = 'constant darkness', color='purple')
ax[2].plot(WT_24_df['time_steps'], WT_24_df['PHYB'], label = 'constant light', color='red')
ax[2].plot(WT_8_df['time_steps'], WT_8_df['PHYB'], label = 'short day ', color='k', linestyle=':')
ax[2].plot(WT_16_df['time_steps'], WT_16_df['PHYB'], label = 'long day ', color='blue', linestyle=':')
#ax[2].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['phyb22'], label = 'original model', color='k')
#ax[2].legend(loc='best')
ax[2].set_ylabel('Expression')
ax[2].set_xlabel('Hours')
ax[2].set_title('PHYB')

#ax[2].legend(labels = ['linked', 'linked_normalized', 'original'], loc='best')


# PIF
ax[3].plot(WT_0_df['time_steps'], WT_0_df['PIF'], label = 'constant darkness', color='purple')
ax[3].plot(WT_24_df['time_steps'], WT_24_df['PIF'], label = 'constant light', color='red')
ax[3].plot(WT_8_df['time_steps'], WT_8_df['PIF'], label = 'short day', color='k', linestyle=':')
ax[3].plot(WT_16_df['time_steps'], WT_16_df['PIF'], label = 'long day', color='blue', linestyle=':')
#ax[3].plot(hypo_results_long_day_df['Time'], hypo_results_long_day_df['PIF422'], label = 'original model', color='k')
#ax[3].legend(loc='best')
ax[3].set_ylabel('Expression')
ax[3].set_xlabel('Hours')
ax[3].set_title('PIF')

#ax[3].legend(labels = ['linked', 'linked_normalized', 'original'], loc='best')
for item in [int(x) for x in (np.linspace(0,3,4))]:
    ax[item].set_xticks([0, 24, 48, 72, 96, 120])
    ax[item].set_xticklabels([0, 24, 48, 72, 96, 120], size=10)

#fig.legend(labels = ['linked', 'original'])
ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.0),
          ncol=4, fancybox=True, shadow=True)

plt.suptitle('WT Protein Dynamics in Different Light Conditions', size=20, y=0.99)
plt.tight_layout()


# In[109]:


n = 4
colors = plt.cm.inferno(np.linspace(0,0.8,n))

WT_0_df['ELF3_normalized'] = WT_0_df['ELF3']*model_avgs['cE']['Y0']
WT_0_df['COP1_normalized'] = WT_0_df['COP1']*model_avgs['cC']['Y0']
WT_24_df['ELF3_normalized'] = WT_24_df['ELF3']*model_avgs['cE']['Y24']
WT_24_df['COP1_normalized'] = WT_24_df['COP1']*model_avgs['cC']['Y24']

light_cond = [WT_0_df, WT_8_df, WT_16_df, WT_24_df]
lines8 = [8, 32, 56, 64]
lines16 = [16, 40, 64]

fig=plt.figure(figsize=(18, 12))

proteins = ['ELF3_normalized', 'COP1_normalized', 'PHYB', 'PIF']
titles=['ELF3', 'COP1', 'PHYB', 'PIF']
labels = ['constant darkness', 'short day', 'long day', 'constant light']
plot_labels = ['A', 'B', 'C', 'D']
for i1, prot in enumerate(proteins):
    ax=fig.add_subplot(2,2,i1+1)
    for i2, cond in enumerate(light_cond): 
        if i2 == 3: 
            weight = 0.5
        elif i2 == 2: 
            weight = 0.8
        else: 
            weight = 1
        ax.plot(cond['time_steps'][cond['time_steps'] <=72], cond[prot][cond['time_steps'] <=72], color=colors[i2], label=labels[i2], linewidth=3, alpha=weight)
        ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
        ax.set_xticks([0, 24, 48, 72])
        ax.set_xticklabels([0, 24, 48, 72], size=10)
        for line8, line16 in zip(lines8, lines16):
            ax.axvline(x=line8, color = colors[1], linestyle=':', alpha=0.2)
            ax.axvline(x=line16, color = colors[2], linestyle=':', alpha=0.2)
            
    if i1 == 1: 
        ax.axhline(y=0, color = 'grey', linestyle='--', alpha=0.5)
    if i1>1: 
        ax.set_xlabel('Hours', fontsize=11)
    if i1%2 == 0: 
        ax.set_ylabel('Concentration', fontsize=11)
    if i1==0: 
        ax.text(-10.3,43, plot_labels[i1], fontweight="bold", size='14')
    if i1==1: 
        ax.text(-10.3,70, plot_labels[i1], fontweight="bold", size='14')
    if i1==2: 
        ax.text(-10.3,1.05, plot_labels[i1], fontweight="bold", size='14')
    if i1==3: 
        ax.text(-10.3,24.0, plot_labels[i1], fontweight="bold", size='14')
    
plt.legend(loc='upper center', bbox_to_anchor=(0.0, 2.4), ncol=4, fancybox=True, shadow=True, fontsize=12)
#plt.tight_layout()
        


# In[110]:


n = 4
colors = plt.cm.inferno(np.linspace(0,0.8,n))

WT_0_df['ELF3_normalized'] = WT_0_df['ELF3']*model_avgs['cE']['Y0']
WT_0_df['COP1_normalized'] = WT_0_df['COP1']*model_avgs['cC']['Y0']
WT_24_df['ELF3_normalized'] = WT_24_df['ELF3']*model_avgs['cE']['Y24']
WT_24_df['COP1_normalized'] = WT_24_df['COP1']*model_avgs['cC']['Y24']

light_cond = [WT_0_df, WT_8_df, WT_16_df, WT_24_df]
lines8 = [8, 32, 56, 64]
lines16 = [16, 40, 64]

fig=plt.figure(figsize=(18, 12))

proteins = ['ELF3_normalized', 'COP1_normalized', 'PHYB', 'PIF']
titles=['ELF3', 'COP1', 'PHYB', 'PIF']
labels = ['constant darkness', 'short day', 'long day', 'constant light']
plot_labels = ['A', 'B', 'C', 'D']
for i1, prot in enumerate(proteins):
    ax=fig.add_subplot(2,2,i1+1)
    for i2, cond in enumerate(light_cond): 
        if i2 == 3: 
            weight = 0.5
        elif i2 == 2: 
            weight = 0.8
        else: 
            weight = 1
        ax.plot(cond['time_steps'], cond[prot], color=colors[i2], label=labels[i2], linewidth=3, alpha=weight)
        ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
        ax.set_xticks([0, 24, 48, 72])
        ax.set_xticklabels([0, 24, 48, 72], size=10)
        for line8, line16 in zip(lines8, lines16):
            ax.axvline(x=line8, color = colors[1], linestyle=':', alpha=0.2)
            ax.axvline(x=line16, color = colors[2], linestyle=':', alpha=0.2)
    if i1>1: 
        ax.set_xlabel('Hours', fontsize=11)
    if i1%2 == 0: 
        ax.set_ylabel('Expression', fontsize=11)
    if i1==0: 
        ax.text(-10.3,43, plot_labels[i1], fontweight="bold", size='14')
    if i1==1: 
        ax.text(-10.3,70, plot_labels[i1], fontweight="bold", size='14')
    if i1==2: 
        ax.text(-10.3,1.05, plot_labels[i1], fontweight="bold", size='14')
    if i1==3: 
        ax.text(-10.3,24.0, plot_labels[i1], fontweight="bold", size='14')
    
plt.legend(loc='upper center', bbox_to_anchor=(0.0, 2.4), ncol=4, fancybox=True, shadow=True, fontsize=12)
#plt.tight_layout()
        


# In[112]:


n = 4
colors = plt.cm.inferno(np.linspace(0,0.8,n))

WT_0_df['ELF3_normalized'] = WT_0_df['ELF3']*model_avgs['cE']['Y0']
WT_0_df['COP1_normalized'] = WT_0_df['COP1']*model_avgs['cC']['Y0']
WT_24_df['ELF3_normalized'] = WT_24_df['ELF3']*model_avgs['cE']['Y24']
WT_24_df['COP1_normalized'] = WT_24_df['COP1']*model_avgs['cC']['Y24']

light_cond = [WT_0_df, WT_8_df, WT_16_df, WT_24_df]
lines8 = [8, 32, 56, 64]
lines16 = [16, 40, 64]

fig=plt.figure(figsize=(18, 12))

proteins = ['ELF3_normalized', 'COP1_normalized', 'PHYB', 'PIF', 'LHY_prot', 'GI_prot_cytoplasm']
titles=['ELF3', 'COP1', 'PHYB', 'PIF', 'LHY', 'GI']
labels = ['constant darkness', 'short day', 'long day', 'constant light']
plot_labels = ['A', 'B', 'C', 'D']
for i1, prot in enumerate(proteins):
    ax=fig.add_subplot(2,3,i1+1)
    for i2, cond in enumerate(light_cond): 
        if i2 == 3: 
            weight = 0.5
        elif i2 == 2: 
            weight = 0.8
        else: 
            weight = 1
        ax.plot(cond['time_steps'][cond['time_steps'] <=72], cond[prot][cond['time_steps'] <=72], color=colors[i2], label=labels[i2], linewidth=3, alpha=weight)
        ax.set_title(f'{titles[i1]}', fontsize=12, fontweight="bold")
        ax.set_xticks([0, 24, 48, 72])
        ax.set_xticklabels([0, 24, 48, 72], size=10)
        for line8, line16 in zip(lines8, lines16):
            ax.axvline(x=line8, color = colors[1], linestyle=':', alpha=0.2)
            ax.axvline(x=line16, color = colors[2], linestyle=':', alpha=0.2)
    if i1>1: 
        ax.set_xlabel('Hours', fontsize=11)
    if i1>= 0: 
        ax.set_ylabel('Expression', fontsize=11)
    if i1==0: 
        ax.text(-10.3,43, plot_labels[i1], fontweight="bold", size='14')
    if i1==1: 
        ax.text(-10.3,70, plot_labels[i1], fontweight="bold", size='14')
    if i1==2: 
        ax.text(-10.3,1.05, plot_labels[i1], fontweight="bold", size='14')
    if i1==3: 
        ax.text(-10.3,24.0, plot_labels[i1], fontweight="bold", size='14')
    
plt.legend(loc='upper center', bbox_to_anchor=(0.0, 2.4), ncol=4, fancybox=True, shadow=True, fontsize=12)
#plt.tight_layout()
        


# #### **Growth Curves**

# In[113]:


daylengths


# In[114]:


growth_dict = {}
for mut in mutants:
    growth_dict[mut] = []
    for hour in daylengths: 
        growth_dict[mut].append(final_sim_outputs[mut][hour][-1][-1])
    


# In[115]:


mutants


# In[116]:


growth_dict


# In[117]:


get_ipython().run_line_magic('matplotlib', 'inline')
fig=plt.figure(figsize=(10,13))
fig.subplots_adjust(hspace=0.5)
ncols=3
nrows=3
for i1,mut in enumerate(mutants):
    ax=fig.add_subplot(nrows,ncols,i1+1)
    #PLOT SIMULATIONS PYTHON
    hp22=[]
    ax.plot(daylength, growth_dict[mut], 'red', label='clock_linked model', marker='.')
    #if mut != 'Col': 
    ax.plot(daylength, growth_dict['Col'], 'red', label='WT clock_linked model', linestyle='--')
    for D in Daylength2:
        if mut == 'phyB': 
            mut = mut+'-9'
        if mut not in hypo_python.keys(): 
            continue
        key22='22_'+str(D)
        hp22.append(hypo_python[mut][key22])
    if mut in hypo_python.keys(): 
        ax.plot(Daylength2, hp22, 'k', label='original model', marker='.')
    #DATA
    d22=[]
    s22=[]
    for D in Daylength:
        key22='22_'+str(D)
        if mut in avgdata:
            d22.append(avgdata[mut][key22])
            s22.append(stddata[mut][key22])
    if mut in avgdata:
        ax.errorbar(Daylength, d22, yerr=s22, fmt='o', color='grey', label='observed')
    if mut=='Col':
        ax.set_title(mut+'/WT', size=12)
    else:
        ax.set_title(mut, style='italic', size=12)
    ax.set_ylim([0,16])
    #if i1==0:
        #ax.legend(loc='upper right', frameon=False)
    if i1 == 0: 
        #ax.legend(loc='best', frameon=False)
        fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.02),
          ncol=4, fancybox=True, shadow=True)
    if i1>5:
        ax.set_xlabel('Daylength (hours)', size=12)
    if i1%1==0:
        ax.set_ylabel('growth (mm)', size=12)
    ax.set_xticks([0,4,8,12,16,20,24])
    ax.set_xticklabels([0,4,8,12,16,20,24], size=10)
    # setting label format to integer with 0 floating decimals
    label_format = '{:,.0f}'
    ticks_loc = ax.get_yticks().tolist()
    ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    ax.set_yticklabels([label_format.format(x) for x in ticks_loc], size=8)

#plt.suptitle('Growth as a Function of Daylength', size=20, y=1)
    
fig.tight_layout()


# In[118]:


fig, ax = plt.subplots(1, sharex=True, figsize=(5,5))
for mut in ['cop1-4']:
    hp22=[]
    #PLOT SIMULATIONS PYTHON
    for D in Daylength2:
        key22='22_'+str(D)
        hp22.append(hypo_python[mut][key22])
    ax.plot(Daylength2, hp22, 'k', label='original model', marker='.')
    ax.plot(daylength, growth_dict[mut], 'red', label='clock_linked model', marker='.')
    d22 = []
    s22=[]
    for D in Daylength:
        key22='22_'+str(D)
        if mut in avgdata:
            d22.append(avgdata[mut][key22])
            s22.append(stddata[mut][key22])
    if mut in avgdata:
        ax.errorbar(Daylength, d22, yerr=s22, fmt='o', color='grey', label='observed')
    ax.legend()
    ax.set_xlabel('Daylength (hours)', size=10)
    ax.set_ylabel('Growth (mm)', size=10)
    plt.title('cop1-4 dynamics', style='italic') 


# ##### **Changing mutC values for COP1OE and cop1-4 lines.** 

# In[119]:


eng = matlab.engine.start_matlab()
path = '/Users/chloe/Desktop/thesis'
eng.cd(path, nargout=0)


# In[120]:


cop1_mut_sims = {}
mutants_cop1 = ['COP1-OE', 'cop1-4']
daylengths = ['0L', '4L', '8L', '12L', '16L', '20L', '24L']


for mut in mutants_cop1:
    count = 0
    cop1_mut_sims[mut] = {}
    for sim,hour in  zip(daylengths, daylight_hours): 
        constants = [float(x) for x in list(model_avgs.iloc[count, 5:])]
        #print(constants)
        cop1_mut_sims[mut][sim] = eng.run_model_cop1_mutc(mut, float(hour), constants[0], constants[1], constants[2])
        count +=1


# In[121]:


clock_outputs = 'LHY mRNA,P,GI-ZTL,GI-ELF3 cytoplasm,LHY prot,TOC1 mRNA,PRR9 prot,PRR5 (NI) mRNA,PRR5 (NI) prot,GI prot cytoplasm,TOC1 prot,ZTL,EC,GI mRNA,PRR9 mRNA,PRR7 mRNA,PRR7 prot,ELF4 mRNA,ELF4 prot,LHY prot modif,ABAR mRNA,COP1 cytoplasm,ELF3 mRNA,ELF3 cytoplasm,ELF3 nuclear,COP1 nuclear night,COP1 nuclear day,LUX mRNA,LUX prot,PP2C prot,SnRK2 prot,stomata,?,?,?,PHYB,ELF3,PIF,COP1,growth'

clock_outputs = clock_outputs.split(',')
clock_outputs
header = ['time_steps']
for item in clock_outputs: 
    item = item.replace(" ", '_')
    header.append(item)


COP1ox_0_df_100 = create_df(final_sim_outputs['COP1-OE'], '0L')
COP1ox_4_df_100 = create_df(final_sim_outputs['COP1-OE'], '4L')
COP1ox_8_df_100 = create_df(final_sim_outputs['COP1-OE'], '8L')
COP1ox_12_df_100 = create_df(final_sim_outputs['COP1-OE'], '12L')
COP1ox_16_df_100 = create_df(final_sim_outputs['COP1-OE'], '16L')
COP1ox_20_df_100 = create_df(final_sim_outputs['COP1-OE'], '20L')
COP1ox_24_df_100 = create_df(final_sim_outputs['COP1-OE'], '24L')


cop1_0_df_01 = create_df(final_sim_outputs['cop1-4'], '0L')
cop1_4_df_01 = create_df(final_sim_outputs['cop1-4'], '4L')
cop1_8_df_01 = create_df(final_sim_outputs['cop1-4'], '8L')
cop1_12_df_01 = create_df(final_sim_outputs['cop1-4'], '12L')
cop1_16_df_01 = create_df(final_sim_outputs['cop1-4'], '16L')
cop1_20_df_01 = create_df(final_sim_outputs['cop1-4'], '20L')
cop1_24_df_01 = create_df(final_sim_outputs['cop1-4'], '24L')


# In[122]:


growth_dict_cop1_mut = {}
for mut in mutants_cop1:
    growth_dict_cop1_mut[mut] = []
    for hour in daylengths: 
        growth_dict_cop1_mut[mut].append(cop1_mut_sims[mut][hour][-1][-1])


# In[123]:


growth_dict_cop1_mut


# **COP1 OE lines** 

# In[124]:


fig, ax = plt.subplots(1, sharex=True, figsize=(5,5))
for mut in ['COP1-OE']:
    hp22=[]
    #PLOT SIMULATIONS PYTHON
    for D in Daylength2:
        key22='22_'+str(D)
        hp22.append(hypo_python[mut][key22])
    ax.plot(Daylength2, hp22, 'k', label='original model', marker='.')
    ax.plot(daylength, growth_dict[mut], 'red', label='clock_linked model, mutC=498', marker='.')
    ax.plot(daylength, growth_dict_cop1_mut[mut], 'red', label='clock_linked model, mutC=100', marker='.', linestyle='--')
    d22 = []
    s22=[]
    for D in Daylength:
        key22='22_'+str(D)
        if mut in avgdata:
            d22.append(avgdata[mut][key22])
            s22.append(stddata[mut][key22])
    if mut in avgdata:
        ax.errorbar(Daylength, d22, yerr=s22, fmt='o', color='grey', label='observed')
    ax.legend()
    ax.set_xlabel('Daylength (hours)', size=10)
    ax.set_ylabel('Growth (mm)', size=10)
    plt.title('COP-OE', style='italic') 


# cop1-4 mutant lines 

# In[125]:


fig, ax = plt.subplots(1, sharex=True, figsize=(5,5))
for mut in ['cop1-4']:
    hp22=[]
    #PLOT SIMULATIONS PYTHON
    for D in Daylength2:
        key22='22_'+str(D)
        hp22.append(hypo_python[mut][key22])
    ax.plot(Daylength2, hp22, 'k', label='original model', marker='.')
    ax.plot(daylength, growth_dict[mut], 'red', label='clock_linked model, mutC=0.032', marker='.')
    ax.plot(daylength, growth_dict_cop1_mut[mut], 'red', label='clock_linked model, mutC=0.1', marker='.', linestyle='--')
    d22 = []
    s22=[]
    for D in Daylength:
        key22='22_'+str(D)
        if mut in avgdata:
            d22.append(avgdata[mut][key22])
            s22.append(stddata[mut][key22])
    if mut in avgdata:
        ax.errorbar(Daylength, d22, yerr=s22, fmt='o', color='grey', label='observed')
    ax.legend()
    ax.set_xlabel('Daylength (hours)', size=10)
    ax.set_ylabel('Growth (mm)', size=10)
    ax.text(-5, 6.5, 'A', fontweight="bold", size='14')
    plt.title('cop1-4', style='italic') 


# In[126]:


#!pip install scikit-learn
#from sklearn.metrics import mean_squared_error
true_values = d22
pred_values_01 = growth_dict_cop1_mut['cop1-4']
pred_values = growth_dict['cop1-4']
pred_values.pop(1)
pred_values.pop(5)
pred_values_01.pop(1)
pred_values_01.pop(5)
print(pred_values)
print(true_values)
print(np.square(np.subtract(true_values, pred_values_01)).mean())
print(np.square(np.subtract(true_values, pred_values)).mean())


# In[127]:


fig, ax = plt.subplots(1, sharex=True, figsize=(5,5))
d22_norm = []
for item in d22: 
    d22_norm.append(item/d22[0])

cop1_norm = []   
for item in growth_dict_cop1_mut['cop1-4']: 
    cop1_norm.append(item/growth_dict_cop1_mut['cop1-4'][0])
    
cop1_norm_03 = []   
for item in growth_dict['cop1-4']: 
    cop1_norm_03.append(item/growth_dict['cop1-4'][0])
    
ax.plot(d22_norm, label = 'observed', color = 'gray')
ax.plot(cop1_norm, label = 'mutC = 0.1, sim', color = 'red', linestyle = '--')
ax.plot(cop1_norm_03, label = 'mutC = 0.032, sim', color = 'red')
plt.legend()
ax.set_xlabel('Daylength (hours)', size=10)
ax.set_ylabel('Relative Hypocotyl Length', size=10)
ax.set_xticks([0, 1, 2, 3, 4])
ax.set_xticklabels([0, 8, 12, 16, 24])

ax.text(-1, 1.1, 'B', fontweight="bold", size='14')
plt.title(r'Growth Response in $cop1-4$ mutants')


# In[128]:


cop1_growth_data = pd.DataFrame()
cop1_growth_data['observed'] = d22
cop1_growth_data['mutc = 0.1'] = pred_values_01
cop1_growth_data['mutc = 0.032'] = pred_values 


# In[129]:


cop1_growth_data


# #### saving data 

# In[ ]:


WT_data = [WT_0_df, WT_4_df, WT_8_df , WT_12_df, WT_16_df, WT_20_df, WT_24_df] 
WT_data_names = ['WT_0_df', 'WT_4_df', 'WT_8_df' ,'WT_12_df', 'WT_16_df', 'WT_20_df', 'WT_24_df'] 


# In[ ]:


for name,data in zip(WT_data_names,WT_data): 
    data.to_csv(name, encoding = "utf-8", index=False)


# In[ ]:


model_avgs.to_csv('model_avgs', encoding = "utf-8", index=False)


# In[ ]:


hypo_results_long_day_df.to_csv('hypo_results_long_day_df', encoding = "utf-8", index=False)


# In[135]:


dict_url = 'https://raw.githubusercontent.com/cbantock/arabidopsis_clock_growth_model/main/WT_0_df'
pd.read_csv(dict_url, engine='python')

