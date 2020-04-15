# 2019ncov-SEIR
 武汉肺炎SEIR模型模拟
武汉新型冠状病毒建模&预测（SEIR）

数据来源：卫健委

附：新型肺炎疫情跟踪石墨文档 https://shimo.im/sheets/tyWrrrqppYVwQtCW/gVSL1

SEIR模型:

S（Susceptible）：易感者
I（Infectious）：已被感染者
R（Removed = recovered + dead）：因治愈或死亡而不再参与感染和被感染过程的人
E（Exposed）：处于潜伏期的感染者，在潜伏期后以一定概率被确诊
N（Population）：常驻人口数量
β：疾病的有效接触率
γ：平均恢复率：I -> R 的转换率，1/γ 表示被感染者从确诊被感染到康复或者死亡的平均日期
βI：S -> I 的转换率
Baseline:

参数估计与假设：

γ = 1/14，新型肺炎从确诊到康复或死亡的日期约为14天
S = N = 11000000，武汉常驻人口大约11000000，初始时全部人口统一归为易感人群
潜伏期 T = 10，潜伏期7天+确诊所需时间3天
感染者每天密切接触的人数 K = 5
模型修正：

考虑潜伏期也具有传染性 传染率与确诊患者相同
考虑潜伏期由于未发病，患者每日接触人员数目为确诊患者的两倍
import numpy as np
import pandas as pd
import pandas
from SIER_model import SIER
from matplotlib import pyplot as plt
from helper_fun_epi_model import *
import scipy.optimize as optimization
import warnings
warnings.filterwarnings('ignore')
# 2019-12-08 为首例确诊病例出现日期. 
t = np.asarray([0]+[x for x in range(34,69)])  # 时间，官方数据从1.11日开始发布，因此T2=34
'''I = np.asarray([1,38,34,33,33,33,27,28,41,94,169,227,326,
                380,441,502,533,593,1460,1726,2050,2377,2884,3714,4653,5768,7616,9272,10606,12360
])# 官方发布的确诊人数,不含死亡和治愈'''
I = np.asarray([1,41,41,41,41,41,41,45,62,121,198,258,363,
                425,495,572,618,698,1590,1905,2261,2639,3215,
                4109,5142,6384,8351,10117,11618,13603,14982,
                16902,18454,19558,32994,35991])
# 官方发布的全部确诊人数,截至2.16日
Est = Estimate_parameter(nu = 1/14, k = 5, t = t, I = I)
print(Est)
估计的传染率: 0.04 
估计的 R0(基本再生数): 3.1

没有政府干预的情况

γ = 1/14，新型肺炎从确诊到康复或死亡的日期约为14天
S = N = 11000000，武汉常驻人口大约11000000，初始时全部人口统一归为易感人群
潜伏期 T = 10，潜伏期7天+确诊所需时间3天
感染者每天接触的人数 K = 5
潜伏者每天接触的人数 = 2*K = 10
武汉当前死亡率0.04
baseline = Estimate_Wuhan_Outbreak(Est, k = 5,ke = 10, N=11000000,
                E0 = 0, I0 = 1, R0 = 0, T=10, econ = 200)
result = baseline._run_SIER('Estimated 2019-nCoV Outbreak in Wuhan Since 2019-12-08 without Government Control',
        'Wuhan Population', 'Days after 2019-12-08', death_rate = 0.04,show_Sus = False)
估计的传染率: 0.04 
估计的 R0(基本再生数): 3.1

预计最大感染人数  4532735
预计累计感染人数  10997637
预计死亡人数  439892

CASE1: 2020.1.23日之后武汉封城

γ = 1/14，新型肺炎从确诊到康复或死亡的日期约为14天
N = 9000000，武汉封城之后留驻居民约为9000000
潜伏期 T = 8，潜伏期7天+确诊所需时间1天(确诊时间加快)
感染者每天接触的人数 K = 1
潜伏者每天接触的人数 = 2*K = 2
武汉当前死亡率0.04，治愈率0.05
1月23日已有累计确诊380例：I0 = 380，正在接受观察的密切接触者1693人，E0 = 1693，治愈28例，死亡17例：R0 = 28+17 = 45
case1 = Estimate_Wuhan_Outbreak(Est, k =1,ke = 2, N=9000000,
                E0 = 1693, I0 = 380, R0 = 45, T=8, econ = 450)
result1 = case1._run_SIER('Forecat 2019-nCoV Outbreak in Wuhan after 2020-01-23 using offical number',
        'Wuhan Population','Days after 2020-01-23', death_rate = 0.04, show_Sus = False)
估计的传染率: 0.04 
估计的 R0(基本再生数): 3.1

预计最大感染人数  311558
预计累计感染人数  4040249
预计死亡人数  159644

CASE2: 2020.2.1日起（1.23 T天后）政府实施严格控制

感染者确诊之后将被严格监控隔离，每天接触的人数 K = 0
潜伏者每天接触的人数仍为2
武汉当前死亡率0.04，治愈率0.05
2月1日已有累计确诊3215例：I0 = 3215，正在接受观察的密切接触者16557人，E0 = 16557，治愈139例，死亡192例：R0 = 331
case2 = Estimate_Wuhan_Outbreak(Est, k =0, ke = 2,N=9000000,
                E0 = 16557, I0 = 3215, R0 = 331, T=8, econ = 150)
result2 = case2._run_SIER('Forecat 2019-nCoV Outbreak in Wuhan after 2020-02-01 using offical number',
        'Wuhan Population','Days after 2020-02-01', death_rate = 0.04, show_Sus = False)
估计的传染率: 0.04 
估计的 R0(基本再生数): 3.1

预计最大感染人数  15693
预计累计感染人数  57421
预计死亡人数  2428

结论： 在政府严格管控下，武汉疫情有望在2.18日左右达到峰值

CASE3:根据交通数据模拟全国范围内的传播走向

迁移数据：武汉迁出，迁入强度，武汉迁出到各个省份的占比，各个省份人员迁入武汉的占比。
数据示例（武汉人员迁出占比）：{'安徽省': [1.94, 2.58, 2.58, 2.14, 2.3, 2.28, 2.33, 2.42, 2.44, 2.52, 2.26, 2.12, 2.32, 2.41, 2.53, 2.43, 2.52, 2.41, 2.27, 2.27, 2.27, 2.1, 1.92, 1.97, 1.81], '北京市': [1.28, 1.9, 1.9, 1.41, 2.21, 2.11, 1.94, 1.83, 1.49, 1.63, 1.16, 1.21, 1.44, 1.42, 1.3, 1.28, 1.07, 0.73, 0.68, 0.56, 0.47, 0.39, 0.36, 0.39, 0.42], '福建省': [0.61, 1.03, 1.03, 0.84, 0.84, 0.91, 0.9, 1.01, 1.02, 0.98, 0.87, 0.93, 1.02, 1.01, 1.0, 1.07, 1.04, 0.92, 1.01, 0.9, 0.84, 0.74, 0.72, 0.87, 0.68], '甘肃省': [0.21, 0.46, 0.46, 0.39, 0.44, 0.53, 0.55, 0.57, 0.58, 0.52, 0.4, 0.49, 0.52, 0.48, 0.41, 0.35, 0.33, 0.28, 0.31, 0.29, 0.28, 0.27, 0.22, 0.33, 0.29], '广东省': [2.61, 3.76, 3.76, 2.76, 3.48, 3.11, 2.95, 2.92, 2.77, 2.83, 2.22, 2.27, 2.45, 2.18, 1.93, 1.95, 1.85, 1.62, 1.77, 1.66, 1.69, 1.56, 1.55, 2.9, 3.12], '广西壮族自治区': [0.52, 0.93, 0.93, 0.84, 0.89, 0.88, 1.0, 1.17, 1.26, 1.17, 0.98, 1.2, 1.2, 1.01, 0.86, 0.81, 0.78, 0.6, 0.76, 0.72, 0.62, 0.55, 0.52, 0.66, 0.71], '贵州省': [0.64, 0.99, 0.99, 0.94, 1.02, 1.17, 1.16, 1.3, 1.24, 1.02, 0.76, 0.87, 0.98, 0.8, 0.59, 0.59, 0.51, 0.45, 0.47, 0.45, 0.38, 0.3, 0.29, 0.37, 0.4], '海南省': [0.45, 0.67, 0.67, 0.47, 0.52, 0.61, 0.66, 0.54, 0.53, 0.46, 0.36, 0.46, 0.4, 0.39, 0.38, 0.42, 0.38, 0.32, 0.39, 0.37, 0.36, 0.33, 0.32, 0.5, 0.68], '河北省': [0.68, 0.98, 0.98, 0.84, 1.01, 1.04, 1.12, 1.21, 1.27, 1.25, 0.99, 1.04, 1.13, 1.01, 1.16, 1.07, 1.12, 0.97, 1.0, 0.94, 0.82, 0.71, 0.56, 0.59, 0.53], '河南省': [4.13, 4.93, 4.93, 4.4, 4.82, 4.85, 5.0, 5.03, 5.31, 5.23, 4.92, 4.93, 5.34, 5.77, 6.03, 6.46, 6.37, 5.6, 5.86, 6.22, 6.18, 5.67, 5.32, 4.98, 4.33], '黑龙江省': [0.22, 0.32, 0.32, 0.26, 0.28, 0.3, 0.34, 0.34, 0.35, 0.31, 0.29, 0.36, 0.34, 0.33, 0.33, 0.34, 0.33, 0.28, 0.32, 0.27, 0.23, 0.2, 0.18, 0.24, 0.26], '湖北省': [72.23, 61.33, 61.33, 67.74, 62.99, 62.33, 62.13, 61.38, 61.52, 61.81, 67.43, 66.38, 63.9, 64.8, 66.2, 65.75, 66.96, 69.96, 69.5, 70.67, 71.84, 74.77, 76.57, 72.46, 74.56], '湖南省': [3.09, 3.74, 3.74, 3.25, 3.45, 3.64, 3.76, 3.77, 3.61, 3.9, 3.59, 3.45, 3.77, 3.82, 3.57, 3.63, 3.59, 3.75, 3.43, 3.36, 3.4, 3.24, 3.07, 3.12, 2.36], '吉林省': [0.13, 0.18, 0.18, 0.18, 0.19, 0.2, 0.2, 0.25, 0.23, 0.19, 0.16, 0.2, 0.23, 0.22, 0.19, 0.18, 0.22, 0.17, 0.19, 0.16, 0.15, 0.12, 0.11, 0.18, 0.16], '江苏省': [1.58, 2.28, 2.28, 1.88, 2.24, 2.19, 2.17, 2.14, 2.1, 2.21, 1.78, 1.68, 1.78, 1.75, 1.67, 1.65, 1.7, 1.55, 1.39, 1.26, 1.16, 1.03, 0.95, 1.26, 1.31], '江西省': [1.88, 2.38, 2.38, 2.0, 2.14, 2.2, 2.15, 2.13, 2.28, 2.22, 2.15, 1.98, 2.25, 2.24, 2.27, 2.4, 2.4, 2.24, 2.23, 2.09, 2.04, 1.95, 1.84, 1.73, 1.36], '辽宁省': [0.33, 0.43, 0.43, 0.38, 0.43, 0.48, 0.5, 0.43, 0.43, 0.43, 0.36, 0.37, 0.42, 0.41, 0.36, 0.42, 0.39, 0.34, 0.38, 0.34, 0.29, 0.21, 0.17, 0.19, 0.23], '内蒙古自治区': [0.12, 0.21, 0.21, 0.22, 0.2, 0.26, 0.28, 0.3, 0.33, 0.31, 0.23, 0.3, 0.28, 0.22, 0.18, 0.21, 0.18, 0.12, 0.14, 0.14, 0.14, 0.12, 0.09, 0.15, 0.17], '宁夏回族自治区': [0.06, 0.14, 0.14, 0.13, 0.14, 0.13, 0.13, 0.17, 0.22, 0.17, 0.14, 0.16, 0.15, 0.11, 0.08, 0.09, 0.05, 0.07, 0.08, 0.06, 0.05, 0.04, 0.02, 0.05, 0.06], '青海省': [0.05, 0.11, 0.11, 0.1, 0.1, 0.14, 0.14, 0.16, 0.11, 0.13, 0.07, 0.13, 0.1, 0.1, 0.07, 0.06, 0.04, 0.04, 0.03, 0.03, 0.02, 0.03, 0.02, 0.04, 0.05], '山东省': [0.84, 1.35, 1.35, 1.15, 1.27, 1.41, 1.51, 1.45, 1.58, 1.56, 1.12, 1.32, 1.3, 1.26, 1.2, 1.24, 1.25, 1.13, 1.13, 1.03, 1.0, 0.85, 0.69, 0.92, 0.85], '山西省': [0.41, 0.66, 0.66, 0.59, 0.51, 0.64, 0.72, 0.83, 0.9, 0.83, 0.67, 0.77, 0.8, 0.69, 0.65, 0.51, 0.57, 0.53, 0.57, 0.55, 0.52, 0.48, 0.47, 0.53, 0.39], '陕西省': [0.67, 0.99, 0.99, 0.74, 0.91, 0.95, 0.93, 0.91, 0.87, 0.84, 0.77, 0.79, 0.84, 0.82, 0.75, 0.8, 0.72, 0.7, 0.83, 0.73, 0.71, 0.57, 0.54, 0.73, 0.65], '上海市': [1.04, 1.35, 1.35, 1.14, 1.59, 1.33, 1.3, 1.26, 1.11, 1.18, 0.97, 0.85, 1.12, 1.04, 0.91, 0.92, 0.8, 0.6, 0.56, 0.46, 0.33, 0.28, 0.31, 0.53, 0.49], '四川省': [0.83, 1.36, 1.36, 1.14, 1.34, 1.36, 1.47, 1.51, 1.65, 1.48, 1.28, 1.54, 1.51, 1.52, 1.41, 1.44, 1.28, 1.18, 1.34, 1.21, 1.13, 0.97, 0.83, 1.2, 1.17], '天津市': [0.2, 0.25, 0.25, 0.27, 0.32, 0.31, 0.27, 0.3, 0.27, 0.29, 0.22, 0.22, 0.22, 0.17, 0.22, 0.2, 0.15, 0.14, 0.14, 0.11, 0.09, 0.07, 0.06, 0.09, 0.08], '西藏自治区': [0.01, 0.03, 0.03, 0.03, 0.03, 0.05, 0.04, 0.04, 0.05, 0.04, 0.02, 0.05, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.0, 0.01, 0.0, 0.01, 0.0], '新疆维吾尔自治区': [0.18, 0.34, 0.34, 0.29, 0.32, 0.43, 0.41, 0.46, 0.4, 0.36, 0.32, 0.37, 0.38, 0.38, 0.25, 0.24, 0.15, 0.19, 0.12, 0.13, 0.11, 0.09, 0.09, 0.13, 0.13], '云南省': [0.54, 0.8, 0.8, 0.64, 0.72, 0.85, 0.79, 0.87, 0.75, 0.74, 0.65, 0.69, 0.8, 0.65, 0.55, 0.58, 0.54, 0.45, 0.49, 0.49, 0.46, 0.39, 0.33, 0.43, 0.54], '浙江省': [1.39, 2.0, 2.0, 1.47, 1.84, 1.74, 1.62, 1.71, 1.57, 1.62, 1.38, 1.28, 1.37, 1.25, 1.14, 1.19, 1.16, 1.14, 1.04, 0.99, 0.89, 0.71, 0.66, 0.9, 0.94], '重庆市': [0.85, 1.15, 1.15, 1.04, 1.14, 1.23, 1.24, 1.27, 1.46, 1.48, 1.2, 1.34, 1.33, 1.43, 1.54, 1.45, 1.32, 1.27, 1.29, 1.27, 1.25, 1.04, 1.0, 1.29, 0.92]}

数据来自：百度慧眼

时间跨度：2020.01.01-2020.01.25
考虑的等级为省级，每一天人员动态移动，以武汉作为传播起点，假设当前武汉有感染人数20人
具体步骤

每一步，所有人员以一定概率进行选择全局出行（省间出行，只考虑与武汉的交互），局部出行（省内出行）
计算人员的转移情况（由于计算所有人员的活动较为复杂，并且已知只有感染者和潜伏者对实验结果密切相关，故只按照迁移比例计算这两种人的转移）
计算接触感染情况
每一步转移和感染之后统计每个省的人员感染情况
1.25日之后默认为城市封闭，均为省内传播

from City_sim import *
#获取出行数据
tra = tra_data("WUHAN_INdata.txt","Wuhan_out.txt","Population.txt")
#tra.wuhan_in()
#tra.wuhan_out()
#tra.population()
#武汉初始感染人是20，共有31个省级城市
city1 = city(tra,20,31)
#时间为25天
result = city1.city_sim(25,Est)
city1.plot_all(result)

#单独查看某城市
city1.plot(result,'四川省')
预计最大感染人数  1598
预计死亡人数  10

#传播过程
result['四川省']
Time	Susceptible	Exposed	Infected	Sum	Resistant	Death	Heal
0	0	8.341000e+07	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000
1	1	8.341017e+07	0.000000	0.166000	0.000000	0.000000	0.000000	0.000000
1	2	8.341061e+07	0.098339	0.402873	0.166000	0.011857	0.000474	0.011383
1	3	8.341116e+07	0.360718	0.616725	0.412707	0.040634	0.001625	0.039008
1	4	8.341086e+07	0.822433	0.795844	0.652797	0.084686	0.003387	0.081298
1	5	8.341075e+07	1.648959	1.044576	0.878088	0.141532	0.005661	0.135870
1	6	8.341097e+07	2.998263	1.376032	1.209472	0.216144	0.008646	0.207498
1	7	8.341135e+07	5.189040	1.870771	1.675858	0.314432	0.012577	0.301855
1	8	8.341188e+07	8.637866	2.615133	2.389675	0.448059	0.017922	0.430136
1	9	8.341224e+07	14.097079	3.776226	3.478919	0.634854	0.025394	0.609460
1	10	8.341149e+07	22.190769	5.474118	5.185933	0.904584	0.036183	0.868401
1	11	8.341150e+07	33.938073	7.933299	7.693195	1.295593	0.051824	1.243769
1	12	8.341150e+07	52.247607	11.790798	11.327106	1.862257	0.074490	1.787767
1	13	8.341120e+07	79.531302	17.553979	17.015559	2.704457	0.108178	2.596279
1	14	8.341114e+07	120.098754	26.142675	25.507110	3.958312	0.158332	3.799980
1	15	8.341096e+07	179.064308	38.648826	38.152551	5.825646	0.233026	5.592621
1	16	8.341083e+07	265.860617	57.079934	56.555257	8.586277	0.343451	8.242826
1	17	8.341046e+07	390.390616	83.553938	83.665995	12.663415	0.506537	12.156878
1	18	8.341051e+07	569.017272	121.556178	122.592999	18.631553	0.745262	17.886291
1	19	8.341041e+07	832.547812	177.610098	178.457905	27.314138	1.092566	26.221572
1	20	8.340998e+07	1209.021275	257.734671	260.864879	40.000573	1.600023	38.400550
1	21	8.340990e+07	1747.080094	372.278258	378.636799	58.410192	2.336408	56.073785
1	22	8.340940e+07	2509.612218	534.636802	546.986267	85.001497	3.400060	81.601437
1	23	8.340781e+07	3587.230615	764.137682	785.598024	123.189840	4.927594	118.262246
1	24	8.340467e+07	5190.822926	1105.224959	1122.860744	177.771103	7.110844	170.660258
1	25	8.340095e+07	7506.487240	1598.476891	1624.307252	256.715743	10.268630	246.447113