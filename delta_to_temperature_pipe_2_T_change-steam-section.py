# -*- coding: utf-8 -*-
"""
Created on 20200919
    从各保温层厚度计算表面温度和界面温度，且分段计算蒸汽管道内的温降
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
#三、线热损与面热损的转换
q=math.pi*D_2*Q

"""

from scipy.optimize import root,fsolve
import numpy as np
import math
from matplotlib import pyplot as plt
from iapws import IAPWS97

#设定（无标明时，除温度单位为摄氏度外，其他单位都按SI国际单位制）
D_0 = 0.325 #管外径
D_i = 0.300 #管内径
l   = 5100 #管长
l_z = 100  #局部阻力当量长度
K_p = 0.0002 #管道内壁绝对粗糙度（m）
pipe_number = 50 #计算温降时，管道分段数量(多分段使计算更准确)
T_0 = 285.4 #入口介质温度（近似为管壁温度）
P_0 = 1.13  #入口压力（MPa）
G = 20 #质量流量（T/h）
K = 1.2 #支架修正系数
T_a = 25 #环境温度
T_f_assump = 480 #预计出口温度
T_1_assump = 200 #界面温度猜测值
Q_assump = 236 #GB50264-2013附录B要求500摄氏度管道的允许热损（W/m2）
T_s_assump = 45 #表面温度猜测值
epsilon = 0.25 #镀锌钢板的黑度，其他外护层材料参考GB50264-2013的5.8.9
W = 2 #风速
material_1 = 2 #指定内层材料，各数值代表的材料类型参看下面的 def lamb(x,T1,T2) 一段
material_2 = 5 #指定外层材料
delta_1 = 0.12 #内保温层厚度
delta_2 = 0.0001 #外保温层厚度

delta = delta_1 + delta_2
D_1 = D_0 + 2*delta_1
D_2 = D_0 + 2*delta

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
    q,T_1,T_s = x[0],x[1],x[2]
    return np.array([(T_s-T_a)*math.pi*D_2*alpha_s(T_s,T_a,W,D_2)-q,
                     (T_1-T_s)/np.log(D_2/D_1)*2*math.pi*lamb(material_2,T_1,T_s)-q,
                     (T_0-T_1)/np.log(D_1/D_0)*2*math.pi*lamb(material_1,T_0,T_1)-q
                     ])


#传热方程组求解，由T_0求解 q 。 （不同管段内介质温度有差异，所以热损也有差异）
def D_to_q(T):
    global T_0 #声明为全局变量，因为本函数需要引用“GB50264_D_to_T_s”
    T_0 = T
    sol_root3 = root(GB50264_D_to_T_s,[math.pi*D_2*Q_assump,T_1_assump,T_s_assump]) 
    sol_fsolve3 = fsolve(GB50264_D_to_T_s,[math.pi*D_2*Q_assump,T_1_assump,T_s_assump])
    q = sol_fsolve3[0]
    Q = q/(math.pi*D_2)
    T_1 = sol_fsolve3[1]
    T_s = sol_fsolve3[2]
    return q,T_s

#压降计算函数,单位为MPa，公式参考：《工业长距离输送蒸汽管道温降压降计算_张伟》
def deltaP(L,D,K_p,rho,v): #压降是管段长度(管段长度+局部阻力当量长度)、管内径、管道内壁绝对粗糙度、介质密度、介质流速的函数
    lambR = 0.11*(K_p/D)**0.25 #摩擦阻力系数的希弗林公式
    return 1.15*lambR*L/D*(rho*v**2)/2*10**-6

#“入口附近管道保温结构”的传热方程组求解，由D_1、D_2求解 q、T_1、T_s
sol_root3 = root(GB50264_D_to_T_s,[math.pi*D_2*Q_assump,T_1_assump,T_s_assump]) 
sol_fsolve3 = fsolve(GB50264_D_to_T_s,[math.pi*D_2*Q_assump,T_1_assump,T_s_assump])
q = sol_fsolve3[0]
Q = q/(math.pi*D_2)
T_1 = sol_fsolve3[1]
T_s = sol_fsolve3[2]
print(" ")
print("solution :",sol_root3)
print(">>>>>>>>>>> When the delta_1 and delta_2 are ",delta_1,delta_2," :")
print("  delta_1 = ",delta_1)
print("  delta_2 = ",delta_2)
print("   lamb_1 = ",lamb(material_1,T_0,T_1))
print("   lamb_2 = ",lamb(material_2,T_1,T_s))
print("       q  = ",q)
print("       Q  = ",Q)
print("      T_1 = ",T_1)
print("      T_s = ",T_s)
print("  alpha_s = ",alpha_s(T_s,T_a,W,D_2))

#末端温度，公式来自 《长输蒸汽管道的温降和压降的计算方法研究-薛永明》、《蒸汽输热管道的温降设计-赵光显》、《工业长距离输送蒸汽管道温降压降计算_张伟》
i = 0 #用于分段计算压降、温降的循环变量，可以当作管段标号（标号从0开始）
H_sec      = np.zeros(pipe_number+1) #定义变量来储存各管段入口蒸汽比焓
deltaH_sec = np.zeros(pipe_number+1) #定义变量来储存各管段蒸汽焓变
P_sec      = np.zeros(pipe_number+1) #定义变量来储存各管段入口蒸汽压力
T_sec      = np.zeros(pipe_number+1) #定义变量来储存各管段入口蒸汽温度
rho_sec    = np.zeros(pipe_number+1) #定义变量来储存各管段入口蒸汽密度
q_sec      = np.zeros(pipe_number+1) #定义变量来储存各管段热损（W/m）
T_s_sec      = np.zeros(pipe_number+1) #定义变量来储存各管段表面温度
steam    = IAPWS97(P=P_0, T=T_0+273.15)  #入口蒸汽物性，由入口蒸汽压力及温度来确定
while i < pipe_number :
    H_sec[i]       = steam.h #第i段入口蒸汽比焓(kJ/kg)
    T_sec[i]       = steam.T-273.15 #第i段入口蒸汽温度
    P_sec[i]       = steam.P #第i段入口蒸汽压力
    rho_sec[i]     = steam.rho #第i段入口蒸汽密度
    q_sec[i]       = D_to_q(T_sec[i])[0] #第i段热损（W/m）
    T_s_sec[i]     = D_to_q(T_sec[i])[1] #第i段表面温度
    deltaH_sec[i]  = q_sec[i]*(l/pipe_number)*K/((1000/3600)*G)/1000 #第i段出入口蒸汽焓变,国际单位制公式：deltaH=q*l*K/G
    H_sec[i+1]     = H_sec[i] - deltaH_sec[i] #第i+1段入口蒸汽比焓(kJ/kg)
    P_sec[i+1]     = P_sec[i] - deltaP((l+l_z)/pipe_number,D_i,K_p,rho_sec[i],G*1000/3600/rho_sec[i]/(math.pi*D_i**2/4)) #第i+1段入口蒸汽压力(MPa)
    steam          = IAPWS97(P=P_sec[i+1], h=H_sec[i+1]) #第i+1段入口蒸汽物性
    T_sec[i+1]     = steam.T-273.15 #第i+1段入口蒸汽温度
    H_f            = H_sec[i+1] #出口蒸汽比焓
    T_f            = T_sec[i+1] #出口蒸汽温度
    P_f            = P_sec[i+1] #出口蒸汽压力
    i = i+1
print("      H_0 = ",H_sec[0])
print("      H_f = ",H_f)
print("      T_f = ",T_f)
print("      P_f = ",P_f)
