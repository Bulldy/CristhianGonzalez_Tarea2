import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rho10=np.genfromtxt("r10.txt")
rho60=np.genfromtxt("r60.txt")
rho120=np.genfromtxt("r120.txt")

t10=np.genfromtxt("t10.txt")
t60=np.genfromtxt("t60.txt")
t120=np.genfromtxt("t120.txt")

print(len(rho10))

r=np.linspace(0,128,65)

fig1=plt.figure()
plt.plot(r,rho10)
plt.title(r'$t= %f $ ' % t10)
plt.xlabel(r'$r (m)$')
plt.ylabel(r'$\rho_{10} / \rho_{atm}$')
plt.savefig("sedov10.pdf")
plt.close()

fig1=plt.figure()
plt.plot(r,rho60)
plt.title(r'$t= %f $' % t60)
plt.xlabel(r'$r (m)$')
plt.ylabel(r'$\rho_{60} / \rho_{atm}$')
plt.savefig("sedov60.pdf")
plt.close()

fig1=plt.figure()
plt.plot(r,rho120)
plt.title(r'$t= %f $' % t120)
plt.xlabel(r'$r (m)$')
plt.ylabel(r'$\rho_{120} / \rho_{atm}$')
plt.savefig("sedov120.pdf")
plt.close()

