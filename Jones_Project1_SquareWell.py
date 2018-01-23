
# coding: utf-8

# In[16]:

import numpy as np

#declare constants:
em = 9.109*10**-31 #mass of electron in kilogram
eq = 1.602*10**-19 #charge of electron
l = 1*10**-10 #width of well in meter
hbar = 1.055*10**-34 #Reduced Planck constant
vx = 100*eq #potential inside well.
n = 1000 
h = l/n


#I have transformed the second order differential form of Schroedinger
#into first order differential and store in arrays.
def f(r,x,E):
    psi = r[0]
    phi = r[1]
    dpsi = phi
    dphi = (2*em/hbar**2)*(vx-E)*psi
    return np.array([dpsi,dphi],float)


#Solve wave equation for a specific energy using RK4.
def solve(E):
    psi = 0.0
    phi = 1.0
    r = np.array([psi,phi],float)
    for x in np.arange(0,l,h):
        k1 = h*f(r,x,E)
        k2 = h*f(r+0.5*k1,x+0.5*h,E)
        k3 = h*f(r+0.5*k2,x+0.5*h,E)
        k4 = h*f(r+k3,x+h,E)
        r += (k1+2*k2+2*k3+k4)/6
        
    return r[0]


#Execute Secant method to find the energy / varying energy until psi=0 at L
E1=0.0
E2 = eq
psi2 = solve(E1)

target = eq/1000
while abs(E1-E2)>target:
    psi1,psi2 = psi2,solve(E2)
    E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)
    

print ("The ground state energy for this particle in this potential energy well is equal to", E2,"Joules")



# In[ ]:



