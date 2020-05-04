#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:07:25 2020

@author: santiagopena
"""


import numpy as np
import matplotlib.pyplot as plt

n=np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
cs03=np.load("counts0.3.npy")
cs04=np.load("counts0.4.npy")
cs05=np.load("counts0.5.npy")

cs08=np.load("counts0.8.npy")
cs09=np.load("counts0.9.npy")
cs10=np.load("counts1.0.npy")


s=3.
a=0.6
plt.plot(n,cs03,marker="o",label=r"$s=0.3$",markersize=s,alpha=a)
plt.plot(n,cs04,marker="o",label=r"$s=0.4$",markersize=s,alpha=a)
plt.plot(n,cs05,marker="o",label=r"$s=0.5$",markersize=s,alpha=a)
plt.plot(n,cs08,marker="o",label=r"$s=0.8$",markersize=s,alpha=a)
plt.plot(n,cs09,marker="o",label=r"$s=0.9$",markersize=s,alpha=a)
plt.plot(n,cs10,marker="o",label=r"$s=1.0$",markersize=s,alpha=a)
plt.xlabel("Number of Oscillators")
plt.ylabel("Jumps between wells")
plt.title("Jumps between one well to the other in function \n of the oscillator number for t_{max}=10000 for 100 average cases \n for E_{bath}=1.")
plt.legend()
plt.show()

