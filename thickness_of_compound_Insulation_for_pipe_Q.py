# -*- coding: utf-8 -*-
"""
Created on 20200523

@author: weililai

#一、外表面换热系数
W = 0.8
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
lamb_1 = (   k1_0 + (T_0/2+T_1/2)*k1_1 + (T_0/2+T_1/2)**2*k1_2      )
lamb_2 = (   k2_0 + (T_1/2+T_s/2)*k2_1 + (T_1/2+T_s/2)**2*k2_2      )
#四、线热损与面热损的转换
q=math.pi*D_2*Q

"""

from scipy.optimize import root,fsolve
import numpy as np
import math
from matplotlib import pyplot as plt

#设定（除温度单位为摄氏度外，其他单位都按SI国际单位制）
D_0 = 0.3 #管径
T_0 = 400 #介质温度（近似为管壁温度）
T_a = 20 #环境温度
T_1 = 350 #两种保温材料界面处温度，不得高于0.9倍外层材料最高使用温度，建议0.8倍
Q_max = 204 #GB50264-2013附录B要求400摄氏度管道的允许热损（W/m2）
T_s_max = 45 #最高表面温度；有时就是防烫伤温度
epsilon = 0.25 #镀锌钢板的黑度，其他外护层材料参考GB50264-2013的5.8.9
W = 0 #风速
k1_0 = 0.03 #内层材料导热系数方程多项式0次项系数
k1_1 = 5*10**-5 #内层材料导热系数方程多项式1次项系数
k1_2 = 2*10**-7 #内层材料导热系数方程多项式2次项系数
k2_0 = 0.035 #外层材料导热系数方程多项式0次项系数
k2_1 = 5*10**-5 #外层材料导热系数方程多项式1次项系数
k2_2 = 3*10**-7 #外层材料导热系数方程多项式2次项系数

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
                     (T_1-T_s)/np.log(D_2/D_1)*2*math.pi*(k2_0 + (T_1/2+T_s/2)*k2_1 + (T_1/2+T_s/2)**2*k2_2)-q, #外层材料传热方程
                     (T_0-T_1)/np.log(D_1/D_0)*2*math.pi*(k1_0 + (T_0/2+T_1/2)*k1_1 + (T_0/2+T_1/2)**2*k1_2)-q, #内层材料传热方程
                     math.pi*D_2*Q-q #线热损与面热损换算
                     ])

#传热方程组设立，由表面温度T_s求解 q、D_1、D_2，外表面换热系数的选取方法参照标准《GB50264-2013》
def GB50264_T_s_to_q(x):
    q,D_1,D_2 = x[0],x[1],x[2]
    alpha_r = 5.669*epsilon/(T_s-T_a)*(((273+T_s)/100)**4-((273+T_a)/100)**4)
    if W == 0:
        alpha_c = 26.4/(297+0.5*(T_s+T_a))**0.5*((T_s-T_a)/D_2)**0.25
    elif W*x[2] > 0.8:
        alpha_c = 4.53*W**0.805/D_2**0.195 
    else:
        alpha_c = 0.008/D_2+4.2*W**0.618/D_2**0.382 
    alpha_s = alpha_r + alpha_c
    return np.array([(T_s-T_a)*math.pi*D_2*alpha_s-q,
                     (T_1-T_s)/np.log(D_2/D_1)*2*math.pi*(k2_0 + (T_1/2+T_s/2)*k2_1 + (T_1/2+T_s/2)**2*k2_2)-q,
                     (T_0-T_1)/np.log(D_1/D_0)*2*math.pi*(k1_0 + (T_0/2+T_1/2)*k1_1 + (T_0/2+T_1/2)**2*k1_2)-q
                     ])

#传热方程组设立，由D_1、D_2求解 q、T_1、T_s，外表面换热系数的选取方法参照标准《GB50264-2013》
def GB50264_D_to_T_s(x):
    q,T_1,T_s = x[0],x[1],x[2]
    alpha_r = 5.669*epsilon/(T_s-T_a)*(((273+T_s)/100)**4-((273+T_a)/100)**4)
    if W == 0:
        alpha_c = 26.4/(297+0.5*(T_s+T_a))**0.5*((T_s-T_a)/D_2)**0.25
    elif W*D_2 > 0.8:
        alpha_c = 4.53*W**0.805/D_2**0.195 
    else:
        alpha_c = 0.008/D_2+4.2*W**0.618/D_2**0.382 
    alpha_s = alpha_r + alpha_c
    return np.array([(T_s-T_a)*math.pi*D_2*alpha_s-q,
                     (T_1-T_s)/np.log(D_2/D_1)*2*math.pi*(k2_0 + (T_1/2+T_s/2)*k2_1 + (T_1/2+T_s/2)**2*k2_2)-q,
                     (T_0-T_1)/np.log(D_1/D_0)*2*math.pi*(k1_0 + (T_0/2+T_1/2)*k1_1 + (T_0/2+T_1/2)**2*k1_2)-q
                     ])

def out_alpha_s():
    alpha_r = 5.669*epsilon/(T_s-T_a)*(((273+T_s)/100)**4-((273+T_a)/100)**4) #辐射换热系数
    if W == 0:
        alpha_c = 26.4/(297+0.5*(T_s+T_a))**0.5*((T_s-T_a)/D_2)**0.25 #无风时的对流换热系数
    elif W*D_2 > 0.8:
        alpha_c = 4.53*W**0.805/D_2**0.195 
    else:
        alpha_c = 0.008/D_2+4.2*W**0.618/D_2**0.382 
    alpha_s = alpha_r + alpha_c #外表面换热系数应为辐射换热系数与对流换热系数之和
    return alpha_s

#传热方程组求解，由最大热损q_max求解 T_s、D_1、D_2
Q = Q_max  #这里与以线热损q为要求的求解程序不同
sol_root1 = root(GB50264_q_to_T_s,[50,D_0+0.04,D_0+0.20,math.pi*(D_0+0.20)*Q]) 
sol_fsolve1 = fsolve(GB50264_q_to_T_s,[50,D_0+0.04,D_0+0.20,math.pi*(D_0+0.20)*Q])
T_s = sol_fsolve1[0]
D_1 = sol_fsolve1[1]
D_2 = sol_fsolve1[2]
q = sol_fsolve1[3]
print(" ")
print("solution :",sol_root1)
print(">>>>>>>>>>> When the Q is ",Q," :")
print("  q   = ",q)
print("  T_s = ",T_s)
print("  D_1 = ",D_1)
print("  D_2 = ",D_2)
print("  alpha_s = ",out_alpha_s())

#当由热损q求解 出的T_s>T_s_max时，重新由表面温度T_s_max求解 q、D_1、D_2
if T_s > T_s_max:
    T_s=T_s_max
    sol_root2 = root(GB50264_T_s_to_q,[q,D_1,D_2]) 
    sol_fsolve2 = fsolve(GB50264_T_s_to_q,[q,D_1,D_2])
    q = sol_fsolve2[0]
    D_1 = sol_fsolve2[1]
    D_2 = sol_fsolve2[2]
    Q=q/math.pi/D_2
    print(" ")
    print("solution :",sol_root2)
    print(">>>>>>>>>>> When the T_s is ",T_s_max," :")
    print("  q = ",q)
    print("  Q = ",Q)
    print("  D_1 = ",D_1)
    print("  D_2 = ",D_2)
    print("  alpha_s = ",out_alpha_s())

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
print("  alpha_s = ",out_alpha_s())


#以列表的形式展现相近结构方案的热损、界面温度、表面温度计算结果
#传热方程组求解，由D_1、D_2-0.04求解 q、T_1、T_s
D_1 = D_0 + 2*delta_1
D_2 = D_0 + 2*(delta-0.02)
sol_root3 = root(GB50264_D_to_T_s,[q,T_1,T_s]) 
sol_fsolve3 = fsolve(GB50264_D_to_T_s,[q,T_1,T_s])
q = sol_fsolve3[0]
T_1 = sol_fsolve3[1]
T_s = sol_fsolve3[2]
Q = q/math.pi/D_2
print(" ")
print("structure delta_1 delta_2 q T_1 T_s Q alpha_s")
#print("solution :",sol_root3)
print(1,round(delta_1,2),round(delta_2-0.02,2),round(q,2),round(T_1,2),round(T_s,2),round(Q,2),round(out_alpha_s(),2))

#传热方程组求解，由D_1、D_2-0.02求解 q、T_1、T_s
D_1 = D_0 + 2*delta_1
D_2 = D_0 + 2*(delta-0.01)
sol_root3 = root(GB50264_D_to_T_s,[q,T_1,T_s]) 
sol_fsolve3 = fsolve(GB50264_D_to_T_s,[q,T_1,T_s])
q = sol_fsolve3[0]
T_1 = sol_fsolve3[1]
T_s = sol_fsolve3[2]
Q = q/math.pi/D_2
#print("solution :",sol_root3)
print(2,round(delta_1,2),round(delta_2-0.01,2),round(q,2),round(T_1,2),round(T_s,2),round(Q,2),round(out_alpha_s(),2))

#传热方程组求解，由D_1、D_2求解 q、T_1、T_s
D_1 = D_0 + 2*delta_1
D_2 = D_0 + 2*delta
sol_root3 = root(GB50264_D_to_T_s,[q,T_1,T_s]) 
sol_fsolve3 = fsolve(GB50264_D_to_T_s,[q,T_1,T_s])
q = sol_fsolve3[0]
T_1 = sol_fsolve3[1]
T_s = sol_fsolve3[2]
Q = q/math.pi/D_2
#print("solution :",sol_root3)
print(3,round(delta_1,2),round(delta_2,2),round(q,2),round(T_1,2),round(T_s,2),round(Q,2),round(out_alpha_s(),2))

#传热方程组求解，由D_1、D_2+0.02求解 q、T_1、T_s
D_1 = D_0 + 2*delta_1
D_2 = D_0 + 2*(delta +0.01)
sol_root3 = root(GB50264_D_to_T_s,[q,T_1,T_s]) 
sol_fsolve3 = fsolve(GB50264_D_to_T_s,[q,T_1,T_s])
q = sol_fsolve3[0]
T_1 = sol_fsolve3[1]
T_s = sol_fsolve3[2]
Q = q/math.pi/D_2
#print("solution :",sol_root3)
print(4,round(delta_1,2),round(delta_2+0.01,2),round(q,2),round(T_1,2),round(T_s,2),round(Q,2),round(out_alpha_s(),2))

#传热方程组求解，由D_1、D_2+0.04求解 q、T_1、T_s
D_1 = D_0 + 2*delta_1
D_2 = D_0 + 2*(delta +0.02)
sol_root3 = root(GB50264_D_to_T_s,[q,T_1,T_s]) 
sol_fsolve3 = fsolve(GB50264_D_to_T_s,[q,T_1,T_s])
q = sol_fsolve3[0]
T_1 = sol_fsolve3[1]
T_s = sol_fsolve3[2]
Q = q/math.pi/D_2
#print("solution :",sol_root3)
print(5,round(delta_1,2),round(delta_2+0.02,2),round(q,2),round(T_1,2),round(T_s,2),round(Q,2),round(out_alpha_s(),2))
