import numpy as np
import matplotlib.pyplot as plt

def wn(n):
    return (n/N)**(1./(1.+s))*np.exp((n/N)/2.)

N=1000000
lista=[0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8]
s=0.8
w=[]
for j in lista:
    s=j
    w=[]
    for i in range(1000,3000000):
        i=i+1
        w.append(wn(i))
    #print(w)
    w=np.asarray(w)
    w=w[w<8]
    w=w[w>0.05]
    weights = np.ones_like(w)/float(len(w))

    a,b,c=plt.hist(w,weights=weights,bins=1000)
    #plt.show()
    np.save("values{}.npy".format(s),b[1:])
    np.save("probs{}.npy".format(s),a)
    #plt.plot(b[1:],a)
    #plt.show()
