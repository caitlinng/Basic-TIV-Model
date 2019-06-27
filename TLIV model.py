#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

print(1)
# Parameters
pV = 210; # Rate of viral production from infected cell
beta = 5e-7; # Rate of free viral particle consumption by target cell
betap = 3e-8; # Rate of target cell infection by free viral particle (less efficient that beta)
gamma = 5e-7 # Rate of latent cells becoming productively infectious
V0 = 1e+4; # Initial number of free viral particles
I0 = 0 # Initial number of infected cells
L0 = 0 # Initial number of latent cells
T0 = 7e+7; # Initial number of target clles
gT = 0.8; # Target cell growth rate
deltaV = 5; # Rate of viral decay
deltaI = 2; # Rate of infected cell death

time = np.linspace(0,22,100)
y_init = [4e+8, 0, 9.3e-2]

# TIV differential equations
def TLIV(t, y): 
    T,L,I,V = y
    return [gT * T * (1 - (T+1)/T0) - (betap * V * T), 
            (betap * T * V) - (gamma * L),
            (gamma * L) - (deltaI * I)
            (pV * I) - (deltaV * V) - (beta * V * T)]

# Solve TIV
sol = solve_ivp(TIV, [time[0], time[-1]], y_init, method = 'RK45', t_eval = time)

# Plot
fig1, ax1 = plt.subplots()
fig1.suptitle("TILV Model")

ax1.plot(sol.t,sol.y[0],sol.t,sol.y[1])
ax1.set_title("target and infected cells")
ax1.set_xlabel("Days")
ax1.set_ylabel("Cells")
ax1.legend(('T','I'))

fig1, ax2 = plt.subplots()
ax2.plot(sol.t,sol.y[2])
ax2.set_title('Viral Load')
ax2.set_xlabel("Days")
ax2.set_ylabel("Viral load (TCID_{50})")
ax2.set_yscale('log')


plt.show()


# In[ ]:


''' 
Varying beta values (rate of viral consumption by binding to target cells) 
e.g. If viral entry enhanced
'''

# Varying beta values
beta_values = np.linspace(1e-7, 5e-7, 4)
time = np.linspace(0,8,100)

# beta_sltns = [beta, x, y]
beta_sltns = []

for i in beta_values:
    beta = i
    # Solve TIV
    sol = solve_ivp(TIV, [time[0], time[-1]], y_init, method = 'RK45', t_eval = time)
    
    # Storing values
    beta_sltn = (beta, sol.t, sol.y[2])   
    beta_sltns.append(beta_sltn)

# Plotting viral load for dif beta values
# Create plot
fig3, ax3 = plt.subplots()
legend = []

for i in beta_sltns:
    print (i[0])
    legend.append(i[0])
    x = i[1]
    y = i[2]
    ax3.plot(x, y)
    
ax3.set_title('Viral Load')
ax3.set_xlabel("Days")
ax3.set_ylabel("Viral load (TCID_{50})")
ax3.set_yscale('log')
ax3.legend((legend), title = '\u03B2 values')

plt.show()

'''
Curve shifts right with increased target cell consumption of virus. 
Time to reach Vmax prolonged -- Because initialy less free viral particles? 
But that would be represented by initial sharp decrease
Prolonged duration of infection -- As overall increase in virus production?
'''


# In[343]:


'''Extract Vmax'''

# Create storage for Vmax values
Vmax_beta = []

# Extract Vmax and store in Vmax_vals as (beta, Vmax)
for i in beta_sltns:
    beta = i[0]
    Vmax = np.amax(i[2])
    Vmax_beta.append((beta, Vmax))
    
# Plot Vmax against beta

fig3, ax4 = plt.subplots()
x = []
y = []
for i in Vmax_beta:   
    x.append(i[0])
    y.append(i[1])

ax4.plot(x, y)

ax4.set_title('Vmax over \u03B2')
ax4.set_xlabel('\u03B2')
ax4.set_ylabel('Vmax')
plt.show()

'''
Decreasing Vmax as increased target cell consumption of virus. 
Less free viral particles since enhanced entry (though possibly more virus -- can check infected cell population.
'''


# In[500]:


''' 
Varying pV values (viral production rate)
e.g. Increased viral translation/replication/budding
'''

pV_values = np.linspace(100, 200, 5)
time = np.linspace(0,28,100)

# pV_sltns is list in which each item = [pV parameter, x values (time), y values (viral load)]
pV_sltns = []

for i in pV_values:
    pV = i
    # Solve TIV
    sol = solve_ivp(TIV, [time[0], time[-1]], y_init, method = 'RK45', t_eval = time)
    
    # Storing values in pV_sltn
    pV_sltn = (pV, sol.t, sol.y[2])   
    pV_sltns.append(pV_sltn)  

    
# Plotting viral load for dif beta values
# Create plot
fig5, ax5 = plt.subplots()
legend = []

for i in pV_sltns:
    print (i[0])
    legend.append(i[0])
    x = i[1]
    y = i[2]
    ax5.plot(x, y)
    
ax5.set_title('Viral Load')
ax5.set_xlabel("Days")
ax5.set_ylabel("Viral load (TCID_{50})")
ax5.set_yscale('log')
ax5.legend((legend), title = 'pV values')

plt.show()


# In[501]:


'''Extract Vmax coordinates and plot time'''

# pV_sltns is list in which each item = [pV parameter, x values (time), y values (viral load)]
pV_sltns = []

for i in pV_values:
    pV = i
    # Solve TIV
    sol = solve_ivp(TIV, [time[0], time[-1]], y_init, method = 'RK45', t_eval = time)
    
    # Storing values in pV_sltn
    pV_sltn = (pV, sol.t, sol.y[2])   
    pV_sltns.append(pV_sltn)  
    
# Create storage for Vmax values
Vmax_pV = []

# Extract Vmax and store in Vmax_vals as (beta, Vmax)
for i in pV_sltns:
    pV = i[0]
    Vmax_x = np.amax(i[2])
    index = np.argmax(i[2]) # of Vmax_x
    Vmax_y = i[1][index]
    Vmax_pV.append((pV, Vmax_x, Vmax_y))

# Plot time to reach Vmax (Vmax_y) against pV

fig6, ax6 = plt.subplots()
x = [] # pV
y = [] # time to reach Vmax (Vmax_y)
for i in Vmax_pV:   
    x.append(i[0])
    y.append(i[2])

ax6.plot(x, y)

ax6.set_title('Time to reach Vmax')
ax6.set_xlabel('pV')
ax6.set_ylabel('Days')
plt.show()

'''
Vmax reached exponentially faster as viral production rate increases
'''


# In[503]:


# Duration of infection
min_thresh = 1e-02 # Infection over when viral load (y) reaches 1e-5 TCID50

# Where y = 1e-5, retrieve x (time) and keep corresponding y value
val = []

for i in pV_sltns:
    pV = i[0]
    print (pV)
    x_time = i[1]
    y_viral = i[2]
    
    index = np.where(y_viral <= min_thresh)
    index = index[0]
    try:
        x_time = x_time[index[0]]
        y_viral = y_viral[index[0]]
    
    except IndexError:
        x_time = 0
        y_viral = 0
        print (("pV = " + str(pV) + " is out of bounds"))

    val.append((pV, x_time, y_viral))

# Plot time to reach Vmax (Vmax_y) against pV
fig7, ax7 = plt.subplots()

x = [] # pV
y = [] # time to reach Vmax (Vmax_y)

for i in val:   
    x.append(i[0])
    y.append(i[1])

ax7.plot(x, y)

ax7.set_title('Duration of infection')
ax7.set_xlabel('pV')
ax7.set_ylabel('Days')
plt.show()

