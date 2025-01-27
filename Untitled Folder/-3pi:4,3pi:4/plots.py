import numpy as np
import matplotlib.pyplot as plt

t=np.load("t_60 _0.515.npy")
Ep5=np.load("Eparti_5 _0.515.npy")
Ep10=np.load("Eparti_10 _0.515.npy")
Ep14=np.load("Eparti_14 _0.515.npy")
Ep18=np.load("Eparti_18 _0.515.npy")
Ep20=np.load("Eparti_20 _0.515.npy")
Ep26=np.load("Eparti_26 _0.515.npy")
Ep30=np.load("Eparti_30 _0.515.npy")
Ep60=np.load("Eparti_60 _0.515.npy")

n=np.load("numeros0.515.npy")
c=np.load("counts0.515.npy")

plt.plot(t,Ep5+1,label=r"$N=5$",linewidth=0.5)
plt.plot(t,Ep10+1,label=r"$N=10$",linewidth=0.5)
plt.plot(t,Ep14+1,label=r"$N=14$",linewidth=0.5)
plt.plot(t,Ep18+1,label=r"$N=18$",linewidth=0.5)
plt.plot(t,Ep20+1,label=r"$N=20$",linewidth=0.5)
plt.plot(t,Ep26+1,label=r"$N=26$",linewidth=0.5)
plt.plot(t,Ep30+1,label=r"$N=30$",linewidth=0.5)
plt.plot(t,Ep60+1,label=r"$N=60$",linewidth=0.5)
plt.xlabel(r"$t$")
plt.ylabel(r"$<E>$")
"""
plt.xscale("log")
plt.yscale("log")
"""
plt.tight_layout()
plt.title(r"$(-\pi/2, \pi)$")
leg = plt.legend()
leg_lines = leg.get_lines()
plt.setp(leg_lines, linewidth=2)



plt.savefig('myfile.png', bbox_inches = "tight")
plt.show()
plt.plot(n,c)
plt.xlabel("Number of oscilators")
plt.ylabel("Average number of jumps")
plt.show()
print(c)
