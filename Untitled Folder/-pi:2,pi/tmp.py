import numpy as np
import matplotlib.pyplot as plt

def wn(n):
    return (n/N)**(1./(1.+s))*np.exp((n/N)/2.)

N=1000000
s=0.6
w=[]
for i in range(100000,2500000):
    i=i+1
    w.append(wn(i))
#print(w)
weights = np.ones_like(w)/float(len(w))

a,b,c=plt.hist(w,weights=weights,bins=1000)
#plt.show()
np.save("values.npy",b[1:])
np.save("probs.npy",a)
#plt.plot(a,b[1:])
#plt.show()
