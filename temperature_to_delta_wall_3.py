# -*- coding: utf-8 -*-
"""
Created on 20200523
    温度->厚度
@author: weililai@foxmail.com

"""

from scipy.optimize import root,fsolve
import numpy as np
import math
from matplotlib import pyplot as plt

#设定（除温度单位为摄氏度外，其他单位都按SI国际单位制）
H = 1 #平壁面的高度（垂直壁面）或宽度（水平壁面）
T_0 = 500 #介质温度
T_a = 20 #环境温度
epsilon = 0.25 #外护层的黑度
v = 1 #风速
material_1 = 8 #内层材料种类
material_2 = 9 #中层材料种类
material_3 = 1 #中层材料种类
Q_max = 236 #允许热损
T_1 = 408.8 #界面1温度
T_2 = 257.5 #界面2温度
T_s = 57.7 #表面温度


#通过指定材料确定该材料的导热系数，参照标准《GB50264-2013》附录A 常用绝热材料性能 ,下面第8、9种为CAS产品
def lamb(x,T1,T2) :
    if x == 1: #170kg/m3 硅酸钙,  推荐使用温度 <= 550℃
        return 0.0479+0.00010185*(T1/2+T2/2)+9.65015*10**-11*(T1/2+T2/2)**3
    if x == 2: #60~80kg/m3 复合硅酸盐,  <= 450℃
        return 0.043+0.00015*((T1/2+T2/2)-70)
    if x == 3: #60~100kg/m3 岩棉,  <= 400℃
        if -20 <= (T1/2+T2/2) <=100:
            return 0.0337+0.000151*(T1/2+T2/2)
        else:
            return 0.0395+4.71*10**-5*(T1/2+T2/2)+5.03*10**-7*(T1/2+T2/2)**2
    if x == 4: #80~100kg/m3 矿渣棉,  <= 300℃
        if -20 <= (T1/2+T2/2) <=100:
            return 0.0337+0.000151*(T1/2+T2/2)
        else:
            return 0.0395+4.71*10**-5*(T1/2+T2/2)+5.03*10**-7*(T1/2+T2/2)**2
    if x == 5: #48kg/m3 玻璃棉,  <= 350℃
        return 0.041 +0.00017*((T1/2+T2/2)-70)
    if x == 6: #各密度的硅酸铝毡， <= 800℃或<= 1000℃，树脂结合毡<= 350℃
        if (T1/2+T2/2) <=400:
            return 0.044+0.0002*((T1/2+T2/2)-70)
        else:
            return 0.044+0.0002*(400-70)+0.00036*((T1/2+T2/2)-400)
    if x == 7: #硅酸镁纤维毯, <= 700℃
        return 0.0397-2.741*10**-6*(T1/2+T2/2)+4.526*10**-7*(T1/2+T2/2)**2
    if x == 8: #CAS-A0,<= 650℃
        return 0.0275 + (T1/2+T2/2)*1.010*10**-4 -(T1/2+T2/2)**2*1.729*10**-8 + (T1/2+T2/2)**3*2.506*10**-10
    if x == 9: #CAS-YN450,<= 450℃
        return 0.0275 + (T1/2+T2/2)*9.766*10**-5 + (T1/2+T2/2)**2*1.247*10**-7
    if x == 10: #气凝胶毡，<= 650℃，标准：《GB/T 34336-2017 纳米孔气凝胶复合绝热制品》
        return 2.9665*10**-7*(T1/2+T2/2)**2 -0.00002732*(T1/2+T2/2)+0.0234976076555024


#表面换热系数选取，公式取自标准《ISO 12241 2008》4.1.3节
def alpha_s(theta_se,theta_a,v,H):
    sigma = 5.67*10**-8 
    C_r = epsilon*sigma
    alpha_r = ((theta_se+273.15)**4 - (theta_a+273.15)**4)/(theta_se - theta_a)
    h_r = alpha_r*C_r
    if v == 0 :
        if H**3*(theta_se - theta_a) <= 10 :
            h_cv = 1.32*((theta_se - theta_a)/H)**0.25
        else :
            h_cv = 1.74*(theta_se - theta_a)**(1/3)
    elif v*H <= 8:
        h_cv = 3.96*(v/H)**0.5
    else:
        h_cv = 5.76*(v**4/H)**0.2
    h_se = h_r + h_cv
    alpha_s = h_se
    return alpha_s  
    
#传热方程组设立
def GB50264_delta_to_T_s(x):
    Q,delta_1,delta_2,delta_3 = x[0],x[1],x[2],x[3]
    return np.array([( T_0 -T_1 )/(delta_1/lamb(material_1,T_0,T_1))-Q,
                     ( T_1 -T_2 )/(delta_2/lamb(material_2,T_1,T_2))-Q,
                     ( T_2 -T_s )/(delta_3/lamb(material_3,T_2,T_s))-Q,
                     ( T_s -T_a )/(1/alpha_s(T_s,T_a,v,H))-Q
                     ])

#传热方程组求解
sol_root = root(GB50264_delta_to_T_s,[Q_max,0.1,0.1,0.1]) 
sol_fsolve = fsolve(GB50264_delta_to_T_s,[Q_max,0.1,0.1,0.1])
Q = sol_fsolve[0]
delta_1 = sol_fsolve[1]
delta_2 = sol_fsolve[2]
delta_3 = sol_fsolve[3]

#打印结果
print(">>>>>>>>>>> When the T_1 T_2 T_s is ",T_1,T_2,T_s," :")
print("      T_1 = ",T_1)
print("      T_2 = ",T_2)
print("      T_s = ",T_s)
print("     lamb_1 = ",lamb(material_1,T_0,T_1))
print("     lamb_2 = ",lamb(material_2,T_1,T_2))
print("     lamb_3 = ",lamb(material_3,T_2,T_s))
print("         Q  = ",Q)
print("    delta_1 = ",delta_1)
print("    delta_2 = ",delta_2)
print("    delta_3 = ",delta_3)
print("    alpha_s = ",alpha_s(T_s,T_a,v,H))
