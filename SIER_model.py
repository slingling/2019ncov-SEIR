# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

class SIER:
    """
    'eons' (number of time points to model, default 1000)
    'Susceptible' (number of susceptible individuals at time 0, default 950)
    'Exposed' (number of individuals during incubation period)
    'Infected' (number of infected individuals at time 0, default 50)
    'Resistant' (number of resistant individuals at time 0, default 0)
    'rateSI' (base rate 'beta' from S to E, default 0.05)
    'rateIR' (base rate 'gamma' from I to R, default 0.01)
    'rateAl' (base rate of isolation 'altha', from E to I, default 0.1)
    """
    def __init__(self, eons=1000, Susceptible=950, Exposed = 100, Infected=50,Sum = 0, Resistant=0, rateSI=0.05,rateSI2=0.05, rateIR=0.01, rateAl = 0.1):
        self.eons = eons
        self.I0 = Infected
        self.Susceptible = Susceptible
        self.Exposed = Exposed
        self.Infected = Infected
        self.Sum = Infected
        self.Resistant = Resistant
        self.rateSI = rateSI
        self.rateSI2 = rateSI2
        self.rateIR = rateIR
        self.rateAl = rateAl
        self.numIndividuals = Susceptible + Infected + Resistant + Exposed # total population
        self.results = None
        self.modelRun = False

    def run(self, death_rate):
        Susceptible = [self.Susceptible]
        Exposed = [self.Exposed]
        Infected = [self.Infected]
        Resistant = [self.Resistant]
        Sum = [self.Sum]

        for step in range(1, self.eons):
            #有修改 考虑了潜伏期也存在传染性
            S_to_E = (Susceptible[-1] *(self.rateSI * Infected[-1] + Exposed[-1]*self.rateSI2)) / self.numIndividuals 
            E_to_I = (self.rateAl * Exposed[-1])
            I_to_R = (Infected[-1] * self.rateIR)
            
            Susceptible.append(Susceptible[-1] - S_to_E)
            Exposed.append(Exposed[-1] + S_to_E - E_to_I)
            Infected.append(Infected[-1] + E_to_I - I_to_R)
            Sum.append(Sum[-1]+E_to_I)
            Resistant.append(Resistant[-1] + I_to_R)
        
        # Death is death_rate* recovery group
        Death = list(map(lambda x: (x * death_rate), Resistant))
        # Heal is removed - Death
        Heal = list(map(lambda x: (x * (1-death_rate)), Resistant))
        self.results = pd.DataFrame.from_dict({'Time':list(range(len(Susceptible))),
            'Susceptible':Susceptible, 'Exposed': Exposed, 'Infected':Infected, 'Sum':Sum,'Resistant':Resistant,
                                               'Death':Death, 'Heal': Heal},
            orient='index').transpose()
        self.modelRun = True
        return self.results

        
    def plot_show(self, title, ylabel, xlabel):
        if self.modelRun == False:
            print('Error: Model has not run. Please call SIR.run()')
            return
        print("预计最大感染人数 ",format(int(max(self.results['Infected']))))
        print("预计累计感染人数 ",format(int(max(self.results['Sum']))+self.I0))
        print("预计死亡人数 ",format(int(max(self.results['Death']))))
        fig, ax = plt.subplots(figsize=(10,6))
        plt.plot(self.results['Time'], self.results['Infected'], color='red')
        plt.plot(self.results['Time'], self.results['Resistant'], color='palegreen')
        plt.plot(self.results['Time'], self.results['Exposed'], color='orange')
        plt.plot(self.results['Time'], self.results['Heal'], color='green')
        plt.plot(self.results['Time'], self.results['Death'], color='grey')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend(['Infected','Removed','Exposed','Heal','Death'], prop={'size': 12}, bbox_to_anchor=(0.5, 1.02), ncol=5, fancybox=True, shadow=True)
        plt.title(title, fontsize = 20)
        plt.show()