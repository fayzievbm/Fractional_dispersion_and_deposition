# -*- coding: utf-8 -*-
"""
Created on Thu May 23 18:26:55 2019

@author: Komiljon
"""
import numpy as np
import modul
from scipy.special import gamma
#boshlang'ich ma'lumotlar
m0=0.3
v=1e-4
c00=0.05

ra0=0.1
ra1=0.005
betar=32
betaa=42
betad=0.002

alpha = 2
D=1e-5

nn=41
tm=4801 

tmax=range(1,tm) #vaqt
tau=1#0.5
h=0.01
#Procedure Nach_Dan;{Boshlang'ich va chegaraviy shartlar}
roa=np.zeros((nn,tm), dtype=float)
c=np.zeros((nn,tm), dtype=float)
c[0,:]=c00

#  /// bettaa ni fisoblash funksiyasi
def bettaa(ra,cij):
   if  ra<=ra1: 
     return betar*v*cij
   elif ra<ra0:
     return betaa*v*cij-betad*ra
   else: 
     return 0.0

#Boshlandi
for j in range(tm-1):
     for i in range(nn):
         roa[i,j+1]=roa[i,j]+tau*bettaa(roa[i,j],c[i,j])
     for i in range(1,nn-1):
         if alpha == 2:
             cc = c[i-1,j]-2*c[i,j]+c[i,j+1]
         else:
             cc = 0
             for k in range(i):
                 cc+=(c[i-(k+1),j]-2*c[i-k,j]+c[i-(k-1),j+1])*((k+1)**(2-alpha)-k**(2-alpha))
         c[i,j+1]=(m0*h*c[i,j]+v*tau*c[i-1,j+1]-h*(roa[i,j+1]-roa[i,j]))/(m0*h+v*tau)+(tau*D*h*cc)/(gamma(3-alpha)*h**alpha)
     c[nn-1,j+1]=c[nn-2,j+1]    
#grafiklar
c=c/c00
x = np.linspace(0, h*(nn-1), num=nn, endpoint=True)
modul.grafik(x,c,roa,tm)

