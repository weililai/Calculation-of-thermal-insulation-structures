# -*- coding: utf-8 -*-
"""
Created on 20200814
    从各保温层厚度计算表面温度和界面温度
@author: weililai@foxmail.com

#一、外表面换热系数
alpha_r = 5.669*epsilon/(T_s-T_a)*(((273+T_s)/100)**4-((273+T_a)/100)**4) #辐射换热系数
if W == 0:
    alpha_c = 26.4/(297+0.5*(T_s+T_a))**0.5*((T_s-T_a)/D_3)**0.25 #无风时的对流换热系数,20200408版中有*号没表示出来，20200513版更正全部相关错误！
elif W*D_3 > 0.8: 20200405版中错写为W > 0.8，20200408版更正全部相关错误！
    alpha_c = 4.53*W**0.805/D_3**0.195 
else:
    alpha_c = 0.008/D_3+4.2*W**0.618/D_3**0.382 
alpha_s = alpha_r + alpha_c #外表面换热系数应为辐射换热系数与对流换热系数之和
#二、一堆传热方程:
#1.保温层外表面与环境 的传热
q = (T_s-T_a)*math.pi*D_3*alpha_s
#2.外层材料的传热
q = (T_2-T_s)/np.log(D_3/D_2)*2*math.pi*lamb_3
#3.中层材料的传热
q = (T_1-T_2)/np.log(D_2/D_1)*2*math.pi*lamb_2
#4.内层材料的传热
q = (T_0-T_1)/np.log(D_1/D_0)*2*math.pi*lamb_1
#三、线热损与面热损的转换
q=math.pi*D_3*Q

"""

from scipy.optimize import root,fsolve
import numpy as np
import math
from matplotlib import pyplot as plt
from iapws import IAPWS97

#设定（除温度单位为摄氏度外，其他单位都按SI国际单位制）
D_0 = 0.377 #管径
l   = 1500 #管长
T_0 = 530 #介质温度（近似为管壁温度）
P_0 = 4  #入口压力（MPa）
P_f = 3.8 #出口压力（MPa）
G = 50 #质量流量（T/h）
K = 1.2 #支架修正系数
T_a = 25 #环境温度
T_f_assump = 480 #预计出口温度
T_a = 20 #环境温度
T_1_assump = 400 #第一个界面温度猜测值
T_2_assump = 300 #第二个界面温度猜测值
T_s_assump = 50 #表面温度猜测值
Q_assump = 236 #GB50264-2013附录B要求500摄氏度管道的允许热损（W/m2）
epsilon = 0.25 #镀锌钢板的黑度，其他外护层材料参考GB50264-2013的5.8.9
W = 1 #风速
K = 1.2 #支架修正系数
material_1 = 8 #指定内层材料
material_2 = 5 #指定中层材料
material_3 = 1 #指定外层材料
delta_1 = 0.04 #内层厚度
delta_2 = 0.06 #中层厚度
delta_3 = 0.06 #外层厚度

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
        return 0.03 + (T1/2+T2/2)*5*10**-5 + (T1/2+T2/2)**2*2*10**-7
    if x == 9: #CAS-YN450,<= 400℃
        return 0.035 + (T1/2+T2/2)*5*10**-5 + (T1/2+T2/2)**2*3*10**-7
    
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

#传热方程组设立，由D_1、D_2求解 q、T_1、T_s，外表面换热系数的选取方法参照标准《GB50264-2013》
def GB50264_D_to_T_s(x):
    q,T_1,T_2,T_s = x[0],x[1],x[2],x[3]
    return np.array([(T_s-T_a)*math.pi*D_3*alpha_s(T_s,T_a,W,D_3)-q,
                     (T_2-T_s)/np.log(D_3/D_2)*2*math.pi*lamb(material_3,T_2,T_s)-q,
                     (T_1-T_2)/np.log(D_2/D_1)*2*math.pi*lamb(material_2,T_1,T_2)-q,
                     (T_0-T_1)/np.log(D_1/D_0)*2*math.pi*lamb(material_1,T_0,T_1)-q
                     ])



#传热方程组求解，由D_1、D_2求解 q、T_1、T_s
delta = delta_1 + delta_2 + delta_2
D_1 = D_0 + 2*delta_1
D_2 = D_0 + 2*(delta_1 + delta_2)
D_3 = D_0 + 2*delta
sol_root3 = root(GB50264_D_to_T_s,[math.pi*D_3*Q_assump,T_1_assump,T_2_assump,T_s_assump]) 
sol_fsolve3 = fsolve(GB50264_D_to_T_s,[math.pi*D_3*Q_assump,T_1_assump,T_2_assump,T_s_assump])
q = sol_fsolve3[0]
Q = q/(math.pi*D_3)
T_1 = sol_fsolve3[1]
T_2 = sol_fsolve3[2]
T_s = sol_fsolve3[3]
print(" ")
print("solution :",sol_root3)
print(">>>>>>>>>>> When the delta_1 delta_2 delta_3 are ",delta_1,delta_2,delta_3," :")
print("  delta_1 = ",delta_1)
print("  delta_2 = ",delta_2)
print("  delta_3 = ",delta_3)
print("   lamb_1 = ",lamb(material_1,T_0,T_1))
print("   lamb_2 = ",lamb(material_2,T_1,T_2))
print("   lamb_3 = ",lamb(material_3,T_2,T_s))
print("       q  = ",q)
print("       Q  = ",Q)
print("      T_1 = ",T_1)
print("      T_2 = ",T_2)
print("      T_s = ",T_s)
print("  alpha_s = ",alpha_s(T_s,T_a,W,D_3))

#末端温度，公式来自 《长输蒸汽管道的温降和压降的计算方法研究-薛永明》、《蒸汽输热管道的温降设计-赵光显》
steam_0 = IAPWS97(P=P_0, T=T_0+273.15)  #入口蒸汽物性
H_0 = steam_0.h #入口蒸汽比焓
deltaH = q*l*K/((1000/3600)*G)/1000 #出入口蒸汽焓变,国际单位制公式：deltaH=q*l*K*G
H_f = H_0 - deltaH #出口蒸汽比焓
steam_f = IAPWS97(P=P_f, h=H_f) #出口蒸汽物性
T_f = steam_f.T-273.15 #出口蒸汽温度
print("      H_0 = ",H_0)
print("   deltaH = ",deltaH)
print("      H_f = ",H_f)
print("      T_f = ",T_f)
