#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:07:25 2020

@author: santiagopena
"""


import numpy as np
import matplotlib.pyplot as plt

n=np.load("ntotales.npy")
cs03=np.load("cs03.npy")
cs04=np.load("cs04.npy")
cs05=np.load("cs05.npy")
cs06=np.load("cs06.npy")
cs07=np.load("cs07.npy")
cs08=np.load("cs08.npy")
cs09=np.load("cs09.npy")
cs10=np.load("cs10.npy")
cs15=np.load("cs15.npy")

s=3.
a=0.6
plt.plot(n,cs03,marker="o",label=r"$s=0.3$",markersize=s,alpha=a)
plt.plot(n,cs04,marker="o",label=r"$s=0.4$",markersize=s,alpha=a)
plt.plot(n,cs05,marker="o",label=r"$s=0.5$",markersize=s,alpha=a)
plt.plot(n,cs06,marker="o",label=r"$s=0.6$",markersize=s,alpha=a)
plt.plot(n,cs07,marker="o",label=r"$s=0.7$",markersize=s,alpha=a)
plt.plot(n,cs08,marker="o",label=r"$s=0.8$",markersize=s,alpha=a)
plt.plot(n,cs09,marker="o",label=r"$s=0.9$",markersize=s,alpha=a)
plt.plot(n,cs10,marker="o",label=r"$s=1.0$",markersize=s,alpha=a)
plt.plot(n,cs15,marker="o",label=r"$s=1.5$",markersize=s,alpha=a)
plt.xlim([0,16])
plt.xlabel("Number of Oscillators")
plt.ylabel("Jumps between wells")
plt.title("Jumps between one well to the other in function \n of the oscillator number for t_{max}=10000 for 100 average cases \n for E_{bath}=0.1")
plt.legend()
plt.show()

