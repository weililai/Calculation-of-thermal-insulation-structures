# -*- coding: utf-8 -*-
"""
Created on 20200523

@author: weililai@foxmail.com

#一、外表面换热系数
alpha_r = 5.669*epsilon/(T_s-T_a)*(((273+T_s)/100)**4-((273+T_a)/100)**4) #辐射换热系数
if W == 0:
    alpha_c = 26.4/(297+0.5*(T_s+T_a))**0.5*((T_s-T_a)/D_2)**0.25 #无风时的对流换热系数,20200408版中有*号没表示出来，20200513版更正全部相关错误！
elif W*D_2 > 0.8: 20200405版中错写为W > 0.8，20200408版更正全部相关错误！
    alpha_c = 4.53*W**0.805/D_2**0.195 
else:
    alpha_c = 0.008/D_2+4.2*W**0.618/D_2**0.382 
alpha_s = alpha_r + alpha_c #外表面换热系数应为辐射换热系数与对流换热系数之和
#二、一堆传热方程:
#1.保温层外表面与环境 的传热
q = (T_s-T_a)*math.pi*D_2*alpha_s
#2.内外保温层界面处-保温层外表面 的传热
q = (T_1-T_s)/np.log(D_2/D_1)*2*math.pi*lamb_2
#3.管道壁面-内外保温层界面处 的传热
q = (T_0-T_1)/np.log(D_1/D_0)*2*math.pi*lamb_1
#三、两种保温材料的导热系数方程
 ...
#四、线热损与面热损的转换
q=math.pi*D_2*Q

"""

from scipy.optimize import root,fsolve
import numpy as np
import math
from matplotlib import pyplot as plt

#设定（除温度单位为摄氏度外，其他单位都按SI国际单位制）
D_0 = 0.410 #管径
T_0 = 540 #介质温度
T_a = 25 #环境温度
T_1_max = 300 #界面温度
Q_max = 236 #允许热损
T_s_max = 50 #最高表面温度
epsilon = 0.75 #外护层的黑度
W = 0 #风速
material_1 = 8 #指定内层材料
material_2 = 5 #指定外层材料

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
        return 0.0397-2.741*10**-6(T1/2+T2/2)+4.526*10**-7*(T1/2+T2/2)**2
    if x == 8: #CAS-A0,<= 650℃
        return 0.0275 + (T1/2+T2/2)*1.010*10**-4 -(T1/2+T2/2)**2*1.729*10**-8 + (T1/2+T2/2)**3*2.506*10**-10
    if x == 9: #CAS-YN450,<= 450℃
        return 0.0275 + (T1/2+T2/2)*9.766*10**-5 + (T1/2+T2/2)**2*1.247*10**-7
    if x == 10: #气凝胶毡，<= 650℃，标准：《GB/T 34336-2017 纳米孔气凝胶复合绝热制品》
        return 2.9665*10**-7*(T1/2+T2/2)**2 -0.00002732*(T1/2+T2/2)+0.0234976076555024

		
#传热方程组设立，由热损Q求解 T_s、D_1、D_2、q，外表面换热系数的选取方法参照标准《GB50264-2013》,与以线热损q为要求的程序不同，增加了q=pi*D_2*Q这样一个等式
def GB50264_q_to_T_s(x):
    T_s,D_1,D_2,q = x[0],x[1],x[2],x[3] #方程中哪些是未知数未知数
    alpha_r = 5.669*epsilon/(T_s-T_a)*(((273+T_s)/100)**4-((273+T_a)/100)**4) #辐射换热系数
    if W == 0:
        alpha_c = 26.4/(297+0.5*(T_s+T_a))**0.5*((T_s-T_a)/D_2)**0.25 #无风时的对流换热系数
    elif W*D_2 > 0.8:
        alpha_c = 4.53*W**0.805/D_2**0.195 
    else:
        alpha_c = 0.008/D_2+4.2*W**0.618/D_2**0.382 
    alpha_s = alpha_r + alpha_c #外表面换热系数应为辐射换热系数与对流换热系数之和
    return np.array([(T_s-T_a)*math.pi*D_2*alpha_s-q, #表面换热方程
                     (T_1-T_s)/np.log(D_2/D_1)*2*math.pi*lamb(material_2,T_1,T_s)-q, #外层材料传热方程
                     (T_0-T_1)/np.log(D_1/D_0)*2*math.pi*lamb(material_1,T_0,T_1)-q, #内层材料传热方程
                     math.pi*D_2*Q-q #线热损与面热损换算
                     ])

#表面换热系数，公式取自标准《GB50264-2013》公式5.8.4
def alpha_s(T_s,T_a,W,D_1):
    alpha_r = 5.669*epsilon/(T_s-T_a)*(((273+T_s)/100)**4-((273+T_a)/100)**4) #辐射换热系数
    if W == 0:
        alpha_c = 26.4/(297+0.5*(T_s+T_a))**0.5*((T_s-T_a)/D_1)**0.25 #无风时的对流换热系数
    elif W*D_1 > 0.8:
        alpha_c = 4.53*W**0.805/D_1**0.195 
    else:
        alpha_c = 0.008/D_1+4.2*W**0.618/D_1**0.382 
    alpha_s = alpha_r + alpha_c #外表面换热系数应为辐射换热系数与对流换热系数之和
    return alpha_s

#传热方程组设立，由表面温度T_s求解 q、D_1、D_2，外表面换热系数的选取方法参照标准《GB50264-2013》
def GB50264_T_s_to_q(x):
    q,D_1,D_2 = x[0],x[1],x[2]
    return np.array([(T_s-T_a)*math.pi*D_2*alpha_s(T_s,T_a,W,D_2)-q,
                     (T_1-T_s)/np.log(D_2/D_1)*2*math.pi*lamb(material_2,T_1,T_s)-q, #外层材料传热方程
                     (T_0-T_1)/np.log(D_1/D_0)*2*math.pi*lamb(material_1,T_0,T_1)-q, #内层材料传热方程
                     ])

#传热方程组设立，由D_1、D_2求解 q、T_1、T_s，外表面换热系数的选取方法参照标准《GB50264-2013》
def GB50264_D_to_T_s(x):
    q,T_1,T_s = x[0],x[1],x[2]
    return np.array([(T_s-T_a)*math.pi*D_2*alpha_s(T_s,T_a,W,D_1)-q,
                     (T_1-T_s)/np.log(D_2/D_1)*2*math.pi*lamb(material_2,T_1,T_s)-q, #外层材料传热方程
                     (T_0-T_1)/np.log(D_1/D_0)*2*math.pi*lamb(material_1,T_0,T_1)-q, #内层材料传热方程
                     ])



#传热方程组求解，由表面温度T_s_max求解 q、D_1、D_2
Q = Q_max
T_1=T_1_max
T_s=T_s_max
sol_root2 = root(GB50264_T_s_to_q,[math.pi*(D_0+0.20)*Q,D_0+0.04,D_0+0.20]) 
sol_fsolve2 = fsolve(GB50264_T_s_to_q,[math.pi*(D_0+0.20)*Q,D_0+0.04,D_0+0.20])
q = sol_fsolve2[0]
D_1 = sol_fsolve2[1]
D_2 = sol_fsolve2[2]
Q=q/math.pi/D_2
print(" ")
print("solution :",sol_root2)
print(">>>>>>>>>>> When the T_1 T_s is ",T_1,T_s," :")
print("  q = ",q)
print("  Q = ",Q)
print("  D_1 = ",D_1)
print("  D_2 = ",D_2)
print("  alpha_s = ",alpha_s(T_s,T_a,W,D_2))

#毡材厚度为1cm时，可行的实际保温层厚度
delta = round(D_2/2-D_0/2+0.005,2)
delta_1 = round(D_1/2-D_0/2+0.005,2)
delta_2 = round(delta - delta_1,2)
print(" ")
print(">>>>>>>>>>> When the blanket is 1 cm thick :")
print("  delta = ",delta)
print("  delta_1 = ",delta_1)
print("  delta_2 = ",delta_2)

#传热方程组求解，由D_1、D_2求解 q、T_1、T_s
D_1 = D_0 + 2*delta_1
D_2 = D_0 + 2*delta
sol_root3 = root(GB50264_D_to_T_s,[q,T_1,T_s]) 
sol_fsolve3 = fsolve(GB50264_D_to_T_s,[q,T_1,T_s])
q = sol_fsolve3[0]
T_1 = sol_fsolve3[1]
T_s = sol_fsolve3[2]
print(" ")
print("solution :",sol_root3)
print(">>>>>>>>>>> When the delta_1 and delta_2 are ",delta_1,delta_2," :")
print("  delta_1 = ",delta_1)
print("  delta_2 = ",delta_2)
print("       q  = ",q)
print("      T_1 = ",T_1)
print("      T_s = ",T_s)
print("  alpha_s = ",alpha_s(T_s,T_a,W,D_2))

