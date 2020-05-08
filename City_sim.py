from helper_fun_epi_model import *
from matplotlib import pyplot as plt
class city(tra_data):
    
    def __init__(self, data: tra_data, I0, city_num, city='武汉市'):
        self.data = data
        self.I0 = I0
        self.city = city
        key = [i for i in self.data.population().keys()]
        ini = [pd.DataFrame.from_dict({'Time':[0],
            'Susceptible':[data.population()[ci]], 'Exposed': [0], 'Infected':[0], 'Sum':[0],'Resistant':[0],
                                               'Death':[0], 'Heal': [0]},orient='index').transpose() for ci in self.data.population()]
        self.all_result = dict(zip(key,ini))
        self.all_result['武汉市']['Infected'] = I0
        self.all_result['武汉市']['Sum'] = I0
        #print(self.all_result)
    
    def city_sim(self,time,Est):
        for timestep in range(time):
            for key1 in self.all_result:
                city_data = Estimate_Wuhan_Outbreak(Est, k = 5,ke = 10, N=self.data.population()[key1],
                        E0 = float(self.all_result[key1]['Exposed'][-1:]), 
                        I0 = float(self.all_result[key1]['Infected'][-1:]), 
                        R0 = float(self.all_result[key1]['Resistant'][-1:]), T=10, econ = 2)
            
                subsub_result = city_data._run_SIER('Estimated 2019-nCoV Outbreak in Hubei',
                    'Hubei Population', 'Days from 2020.1.1 to 2020.1.25', death_rate = 0.04,show = False)
                #self.all_result[key] = subsub_result
            
                #人口转移
                day_infected_in = 0
                day_infected_out = 0
                day_exposed_in = 0
                day_exposed_out = 0
                if key1 =='武汉市':
                    day_in = self.data.in_level[timestep]*3*10000
                    day_out = self.data.out_level[timestep]*3*10000
                    for key2 in self.data.wuhan_in():
                        day_infected_in += float(self.data.wuhan_in()[key2][timestep]*3*self.data.in_level[timestep]*100/self.data.population()[key2]*float(self.all_result[key2]['Infected'][-1:]))
                        day_exposed_in += float(self.data.wuhan_in()[key2][timestep]*3*self.data.in_level[timestep]*100/self.data.population()[key2]*self.all_result[key2]['Exposed'][-1:])
                    day_infected_out = float(self.data.out_level[timestep]*10000*3/self.data.population()['武汉市'] *self.all_result['武汉市']['Infected'][-1:])
                    day_exposed_out = float(self.data.out_level[timestep]*10000/self.data.population()['武汉市'] *self.all_result['武汉市']['Exposed'][-1:])
                    #print(float(self.all_result['武汉市']['Infected'][-1:]))

                elif key1 in self.data.wuhan_in():
                    day_in = float(self.data.out_level[timestep]*3*100*self.data.wuhan_out()[key1][timestep])
                    day_out = float(self.data.in_level[timestep]*3*100*self.data.wuhan_in()[key1][timestep])
                    day_infected_in = float(self.all_result['武汉市']['Infected'][-1:]*self.data.wuhan_out()[key1][timestep]/100)
                    day_infected_out = float(self.all_result[key1]['Infected'][-1:]*self.data.wuhan_in()[key1][timestep]/100)
                    day_exposed_out = float(self.all_result[key1]['Exposed'][-1:]*self.data.wuhan_in()[key1][timestep]/100)
                    day_exposed_in = float(self.all_result['武汉市']['Exposed'][-1:]*self.data.wuhan_out()[key1][timestep]/100)
                    #print(day_exposed_in)
                else:
                    print('城市名称输入错误')
                    return

                subsub_result[-1:]['Exposed'] = subsub_result[-1:]['Exposed']+day_exposed_in-day_exposed_out
                subsub_result[-1:]['Infected'] = subsub_result[-1:]['Infected']+day_infected_in-day_infected_out
                subsub_result[-1:]['Susceptible'] = subsub_result[-1:]['Susceptible']+ day_in - day_out - day_infected_in - day_exposed_out -day_exposed_in - day_infected_out
                #subsub_result[-1:]['Infected'] = subsub_result[-1:]['Infected']+1

                self.all_result[key1] = pd.concat([self.all_result[key1],subsub_result[-1:]],axis = 0)
                #print(self.sub_result)
                self.all_result[key1]['Time'] = [i for i in range(len(self.all_result[key1]))]
                #self.all_result[self.city] = self.sub_result
        return self.all_result
    
    def plot_all(self,result):
        plt.rcParams['font.sans-serif']=['SimHei']
        plt.rcParams['axes.unicode_minus']=False
        plt.figure(figsize=(30,60))
        plt.subplot(11,1,1)
        #print("预计最大感染人数 ",format(int(max(result['武汉市']['Infected']))))
        #print("预计累计感染人数 ",format(int(max(self.results['Sum']))+self.I0))
        #print("预计死亡人数 ",format(int(max(result['武汉市']['Death']))))
        #fig, ax = plt.subplots(figsize=(10,6))
        plt.plot(result['武汉市']['Time'], result['武汉市']['Infected'], color='red')
        plt.plot(result['武汉市']['Time'], result['武汉市']['Resistant'], color='palegreen')
        plt.plot(result['武汉市']['Time'], result['武汉市']['Exposed'], color='orange')
        plt.plot(result['武汉市']['Time'], result['武汉市']['Heal'], color='green')
        plt.plot(result['武汉市']['Time'], result['武汉市']['Death'], color='grey')
        plt.xlabel('Days from 2020.1.1 to 2020.1.25')
        plt.ylabel('Population')
        plt.legend(['Infected','Removed','Exposed','Heal','Death'], prop={'size': 12}, bbox_to_anchor=(0.5, 1.02), ncol=5, fancybox=True, shadow=True)
        plt.title('2019-nCoV Outbreak in Wuhan', fontsize = 30)
        #plt.plot([0,1],[0,1])
        plt.tight_layout()

        for t in range(4,34):
            plt.subplot(11,3,t)
            plt.plot(result[list(result.keys())[t-4]]['Time'], result[list(result.keys())[t-4]]['Infected'], color='red')
            plt.plot(result[list(result.keys())[t-4]]['Time'], result[list(result.keys())[t-4]]['Resistant'], color='palegreen')
            plt.plot(result[list(result.keys())[t-4]]['Time'], result[list(result.keys())[t-4]]['Exposed'], color='orange')
            plt.plot(result[list(result.keys())[t-4]]['Time'], result[list(result.keys())[t-4]]['Heal'], color='green')
            plt.plot(result[list(result.keys())[t-4]]['Time'], result[list(result.keys())[t-4]]['Death'], color='grey')
            plt.xlabel('Days from 2020.1.1 to 2020.1.25')
            plt.ylabel('Population')
            #plt.legend(['Infected','Removed','Exposed','Heal','Death'], prop={'size': 12}, bbox_to_anchor=(0.5, 1.02), ncol=5, fancybox=True, shadow=True)
            plt.title('2019-nCoV Outbreak in '+ list(result.keys())[t-4] , fontsize = 30)
            plt.tight_layout()
            
            
    def plot(self,result,city):
        plt.rcParams['font.sans-serif']=['SimHei']
        plt.rcParams['axes.unicode_minus']=False
        print("预计最大感染人数 ",format(int(max(result[city]['Infected']))))
        #print("预计累计感染人数 ",format(int(max(results['Sum']))+self.I0))
        print("预计死亡人数 ",format(int(max(result[city]['Death']))))
        fig, ax = plt.subplots(figsize=(10,6))
        plt.plot(result[city]['Time'], result[city]['Infected'], color='red')
        plt.plot(result[city]['Time'], result[city]['Resistant'], color='palegreen')
        plt.plot(result[city]['Time'], result[city]['Exposed'], color='orange')
        plt.plot(result[city]['Time'], result[city]['Heal'], color='green')
        plt.plot(result[city]['Time'], result[city]['Death'], color='grey')
        plt.xlabel('Days from 2020.1.1 to 2020.1.25')
        plt.ylabel('Population')
        plt.legend(['Infected','Removed','Exposed','Heal','Death'], prop={'size': 12}, bbox_to_anchor=(0.5, 1.02), ncol=5, fancybox=True, shadow=True)
        plt.title('2019-nCoV Outbreak in '+ city , fontsize = 30)
        plt.show()