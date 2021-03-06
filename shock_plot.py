import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def analytic_solution(rL,rR,pL,pR,Ms,aL,y,x00,tt,xx):
    u1=(2/(y+1))*(Ms-(1/Ms))
    p1=pR*((2*y/(y+1))*Ms*Ms-((y-1)/(y+1)))
    a2=aL*((p1/pL)**((y-1)/(2*y)))
    rLDE=(2/(Ms*Ms*(y+1)))+((y-1)/(y+1))
    r1=rR/rLDE

    x1=np.around(x00-aL*tt,decimals=3)
    x2=np.around(x00+(u1-a2)*tt,decimals=3)
    x3=np.around(x00+u1*tt,decimals=3)
    x4=np.around(x00+Ms*tt,decimals=3)

    rho_T=np.zeros(len(xx))
    p_T=np.zeros(len(xx))
    v_T=np.zeros(len(xx))

    for i in range(len(xx)):
        if xx[i]<=x1:
            rho_T[i]=rL
            p_T[i]=pL
            v_T[i]=0
        elif(xx[i]>x1 and xx[i]<=x2):
            v_T[i]=(2/(y+1))*(aL+(1/tt)*(xx[i]-x00))
            a_Ti=aL-(y-1)*0.5*v_T[i]
            p_T[i]=pL*(a_Ti/aL)**(2*y/(y-1))
            rho_T[i]=y*p_T[i]/(a_Ti**2)
        elif(xx[i]>x2 and xx[i]<=x3):
            rho_T[i]=rL*(p1/pL)**(1/y)
            p_T[i]=p1
            v_T[i]=(2/(y-1))*(aL-a2)
        elif(xx[i]>x3 and xx[i]<=x4):
            rho_T[i]=r1
            p_T[i]=p1
            v_T[i]=(2/(y-1))*(aL-a2)
        elif xx[i]>x4:
            rho_T[i]=rR
            p_T[i]=pR
            v_T[i]=0

    return rho_T,p_T,v_T

rho=np.genfromtxt("density.txt")
p=np.genfromtxt("pressure.txt")
v=np.genfromtxt("velocity.txt")
t=np.genfromtxt("time_shock.txt")

yy=1.4
rho_L=rho[0]
rho_R=rho[-1]
p_L=p[0]
p_R=p[-1] 
M_s=1.6 #Calculated and rounded from the Compatibility Equation
a_L=(yy*p_L/rho_L)**0.5

x=np.linspace(0,1,1001)
x0=0.5

rho_Teo,p_Teo,v_Teo=analytic_solution(rho_L,rho_R,p_L,p_R,M_s,a_L,yy,x0,t,x)


fig1=plt.figure()
plt.plot(x,rho,label=r'Lax-Wendroff')
plt.plot(x,rho_Teo,label=r'Exact Solution')
plt.xlabel(r'$x/L$')
plt.ylabel(r'$\rho/\rho_R$')
plt.ylim(0,11)
plt.legend()
plt.savefig("density_shock.pdf")
plt.close()

fig2=plt.figure()
plt.plot(x,p,label=r'Lax-Wendroff')
plt.plot(x,p_Teo,label=r'Exact Solution')
plt.xlabel(r"$x/L$")
plt.ylabel(r"$P/(\gamma P_R)$")
plt.ylim(0,8)
plt.legend()
plt.savefig("pressure_shock.pdf")
plt.close()

fig3=plt.figure()
plt.plot(x,v,label='Lax-Wendroff')
plt.plot(x,v_Teo,label='Exact Solution')
plt.xlabel(r"$x/L$")
plt.ylabel(r"$U/a_R$")
plt.ylim(0,1)
plt.legend(loc=2)
plt.savefig("velocity_shock.pdf")
plt.close()
