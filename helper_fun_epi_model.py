# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 19:44:27 2020

@author: 15810
"""

import numpy as np
import scipy.optimize as optimization
import pandas as pd
import pandas
from SIER_model import SIER
import re

class Estimate_parameter:
    """
    Estimate other parameters: 
            beta in SIR model and R0(basic reproduction number)
    """
    def __init__(self, nu: float, k: int, t: np.ndarray, I: np.ndarray):
        self.nu = nu
        self.k = k # the number of people a confirmed case contacts daily
        self.t = t # time step
        self.I = I  # observation
        # Estimated beta
        self.beta = self._estimate_beta()
        # self.R0
        self.R0 = self._estimate_R0()
        
    def func(self, t: np.ndarray, b):
        """
        K is the mean number of people a confirmed case contacts daily
        """
        return np.exp((self.k*b-self.nu)*t)

    def _estimate_transmission_probablity(self):
        """
        Estimate the transmission probablity by non-linear OLS
        """
        return optimization.curve_fit(self.func, self.t, 
                                      self.I, maxfev=20000)[0][0]
    def _estimate_beta(self):
        """
        Estimate beta of SIR model
        """
        return self.k * self._estimate_transmission_probablity()
    
    def _estimate_R0(self):
        """
        Estimate R0(basic reproduction number)
        """
        return (self.beta)/self.nu
    
    def __str__(self): 
        """
        Representation
        """
        return f"""估计的传染率: {round(self._estimate_transmission_probablity(), 2)} 
估计的 R0(基本再生数): {round(self.R0,1)}
"""



class Estimate_Wuhan_Outbreak(Estimate_parameter):
    
    def __init__(self, Est: Estimate_parameter, k: int, ke: int, N: int,
                E0: int, I0: int, R0: int, T: int, econ: int):
        self.Est = Est 
        #print(self.Est) # print R0 
        self.k = k # the number of people one case contacts daily, which influnenced by government force
        self.ke = ke #the number of incubation people one case contacts daily
        self.N = N
        self.E0 = E0 # initial number of Enfective cases
        self.I0 = I0 # initial number of recovered cases
        self.R0 = R0
        self.S0 = N - E0 - I0 - R0# initial number of suspective cases
        self.alpha = 1/T # T is the mean of incubation period
        self.econ = econ
        self.model = None
        
        
    def _run_SIER(self,title:str, ylabel:str, xlabel:str, death_rate: float, show = True) -> pandas.core.frame.DataFrame:
        """
        Run SIER model
        """
        Est_beta = self.Est._estimate_transmission_probablity()*self.k
        Est_beta2 = self.Est._estimate_transmission_probablity()*self.ke
        sier = SIER(eons=self.econ, Susceptible=self.N-self.E0-self.I0-self.R0, Exposed = self.E0, 
                    Infected=self.I0, Resistant=self.R0, rateSI=Est_beta,rateSI2=Est_beta2, rateIR=self.Est.nu, 
                    rateAl = self.alpha)
        result = sier.run(death_rate)
        # Draw plot
        if show:
            sier.plot_show(title, ylabel, xlabel)
        
        self.model = sier
        return result


class tra_data:
    def __init__(self, in_path, out_path, population_path):
        self.in_path = in_path #各城市进入武汉
        self.out_path = out_path # 武汉离开到各城市
        self.population_path = population_path #常驻人口数据
        #出行强度,非具体数值 而是占比
        self.in_level = [2.85,3.09,4.22,4.45,5.08,4.31,4.25,4.47,4.81,4.60,4.64,4.37,4.83,4.08,4.06,4.00,4.40,4.23,4.15,4.18,4.24,2.90,1.75,0.88,0.63]
        #出行强度
        self.out_level = [3.46,3.52,5.52,6.10,5.32,5.60,6.41,7.34,8.14,6.62,7.56,6.22,5.76,5.46,5.91,6.00,6.44,7.71,7.41,8.31,10.74,11.84,11.14,3.89,1.30]
        self.Wuhan_in = {}
        self.Wuhan_out = {}
        
    def wuhan_in(self):
    #各城市进入武汉的人员数据
        with open(self.in_path,'rb') as f:  # 打开文件
            data = f.readlines()  # 读取文件
            for line in data[1:]:
                self.Wuhan_in[line.decode('UTF-8').strip('\r\n').split('\t')[0]]=[float(i) for i in line.decode('UTF-8').strip('\r\n').split('\t')[1:]]   
        return self.Wuhan_in
    
    #武汉离开到各城市的数据
    def wuhan_out(self):
        with open(self.out_path,'rb') as f:  # 打开文件
            data = f.readlines()  # 读取文件
            for line in data:
                self.Wuhan_out[line.decode('UTF-8').strip('\r\n').split('\t')[0]]=[float(i) for i in line.decode('UTF-8').strip('\r\n').split('\t')[1:]]
        return self.Wuhan_out
    
    #获取各个城市的常驻人口信息,单独考虑武汉市
    def population(self):
        lis = []
        with open(self.population_path,'rb')as f:
            for i in range(31):
                line = f.readline()
                temp = re.findall(r"\d+\.?\d*",line.decode('UTF-8'))[0]
                lis.append(float(temp)*10000)
        self.wuhan_in()
        key = [i for i in self.Wuhan_in.keys()]
        Province_n = dict(zip(key,lis))
        Province_n['武汉市'] = float(1100.0)*10000
        
        
        return Province_n
    if __name__ == '__main__':
        tra.wuhan_in()
        tra.wuhan_out()
        tra.population()
    