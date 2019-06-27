#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters
pV = 210; # Rate of viral production from infected cell
beta = 5e-7; # Rate of free viral particle consumption by target cell
betap = 3e-8; # Rate of target cell infection by free viral particle (less efficient that beta)
gamma = 2 # Rate of latent cells becoming productively infectious (eclipse phase = 5-12 hours = ~5/24 to 1/2 days)
V0 = 1e+4; # Initial number of free viral particles
I0 = 0 # Initial number of infected cells
L0 = 0 # Initial number of latent cells
T0 = 7e+7; # Initial number of target clles
gT = 0.8; # Target cell growth rate
deltaV = 5; # Rate of viral decay
deltaI = 2; # Rate of infected cell death

time = np.linspace(0,2,120)
#time = np.linspace(0,22,100)
y_init = [T0, L0, I0, V0]

# TIV differential equations
def TLIV(t, y): 
    T,L,I,V = y
    return [gT * T * (1 - (T+1)/T0) - (betap * V * T), 
            (betap * T * V) - (gamma * L),
            (gamma * L) - (deltaI * I),
            (pV * I) - (deltaV * V) - (beta * V * T)]

# Solve TIV
sol = solve_ivp(TLIV, [time[0], time[-1]], y_init, method = 'RK45', t_eval = time)

# Plot
fig1, ax1 = plt.subplots()
fig1.suptitle("TLIV Model")

ax1.plot(sol.t,sol.y[0],sol.t,sol.y[1],sol.t,sol.y[2])
ax1.set_title("target and infected cells")
ax1.set_xlabel("Days")
ax1.set_ylabel("Cells")
ax1.legend(('T', 'L', 'I'))

fig1, ax2 = plt.subplots()
ax2.plot(sol.t,sol.y[3])
ax2.set_title('Viral Load')
ax2.set_xlabel("Days")
ax2.set_ylabel("Viral load (TCID_{50})")
ax2.set_yscale('log')


plt.show()


# In[ ]:


''' 
Varying gamma values
'''

# Varying gamma values
gamma_values = np.linspace(0.01, 5, 10)
time = np.linspace(0,2,100)

# gamma_sltns = [gamma, x, y]
gamma_sltns = []

for i in gamma_values:
    gamma = i
    # Solve TLIV
    sol = solve_ivp(TLIV, [time[0], time[-1]], y_init, method = 'RK45', t_eval = time)
    
    # Storing values
    gamma_sltn = (gamma, sol.t, sol.y[3])   
    gamma_sltns.append(gamma_sltn)

# Plotting viral load for dif beta values
# Create plot
fig3, ax3 = plt.subplots()
legend = []

for i in gamma_sltns:
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


# In[ ]:


'''Find turning point (start of detectable viral production)'''
print (1)
# Create storage for Vmax values
vals = []

# Extract Vmax and store in Vmax_vals as (beta, Vmax)
for i in gamma_sltns:
    print(1)
    gamma = i[0]
    val = np.amin(i[2])
    vals.append((gamma, vals))
    
# Plot Vmax against beta

fig3, ax4 = plt.subplots()
x = []
y = []
for i in vals:   
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


# In[ ]:




