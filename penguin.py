#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created: May, 20, 2022 14:27:54 by Hiroto Akira
# $Author: $
# $Date: $
# $URL: $
__giturl__ = "$URL: $"

def main():
  import numpy
  # 海水
  T_s = 0  # 温度 ℃
  rho_s = 1028  # 密度 kg/m^3
  mu_s = 0.0185  # 粘度 Pa s

  # 空気
  T_p = 38  # 体温 ℃
  rho_a = 1.0985  # 文献値

  # ペンギンの肺
  r_l_0 = 0.12  # 半径 m
  V_l_0 = 4/3*numpy.pi*r_l_0**3  # 肺の体積kg/m^3
  print('肺の体積',V_l_0)
  M_l = rho_a*V_l_0  # 肺の空気の重さ（一定）
  eta = 0.0025

  # ペンギン
  h_p = 1.1  # 身長 m
  d_p = 0.35  # 直径 m
  rho_p = rho_s*1.01  # 密度 kg/m^3
  V_p_0 = numpy.pi/6*h_p*d_p**2  # 外見の体積 m^3
  tV_p = V_p_0-V_l_0  # ペンギンの身の体積（一定） m^3
  tM_p = rho_p*tV_p  # ペンギンの身の重さ（一定） kg
  print('ペンギンの体重 =',tM_p)
  u_max = 35*1000/3600  # ペンギンの最大速度 m/s
  Re = rho_s*u_max*h_p/mu_s
  print('(2) Re =',Re)
  M_p = tM_p+M_l
  
  b = 0.2  # 羽の長さ m
  c = 0.08  # 羽の幅

  print(V_l_0/tV_p)

  # その他
  P_0 = 101325  # 大気圧 Pa
  g = 9.81  # 重力加速度 m/s^2

  # 条件を満たしているかの確認
  print('='*20)
  print(rho_s*tV_p)
  print(M_p)
  print(rho_s*(tV_p+V_l_0))
  print('='*20)

  # 問3
  # 釣り合いの位置z
  ## z_balance = (rho_s*V_l_0-M_l)/((rho_p-rho_s)*tV_p+M_l)*P_0/(rho_s*g)
  z_balance = (rho_s*V_l_0/(M_p-rho_s*tV_p)-1)*P_0/(rho_s*g)
  
  z_balance2 = (rho_s*(tV_p+V_l_0)-M_p)/(M_p-rho_s*tV_p)*P_0/(rho_s*g)
  print('(3) z =',z_balance)
  print('(3) z =',z_balance2)
  ## print(z)
  # 問4
  # 巡航速度を決める
  u_p = u_max*0.9
  # まずは上下方向の釣り合いを考える
  ## V_p_4 = get_V_p(z_balance,rho_s,g,P_0,tV_p,V_l_0)
  V_p_4,V_l_4 = get_V(0.99*z_balance,rho_s,g,P_0,tV_p,V_l_0)
  print('確認',rho_s*V_p_4-M_p)
  print('確認',rho_s*V_p_4/M_p)
  ## print(V_p_4)
  f_wing = (rho_s*V_p_4-M_p)*g  # 上下方向に釣り合うために羽から得なければならない下向きの力（二枚分） N
  ## f_wing = (tM_p-rho_s*V_p_4)*g  # 釣り合うために羽から得なければならない上向きの力 N
  print('f_wing =',f_wing)

  C_L = f_wing/(0.5*rho_s*2*c*b*u_p**2)
  alpha = C_L/2/numpy.pi
  print('alpha =',alpha)

  # 次に，モーメントを考える．
  # その準備として，式(19)より，
  zeta = (rho_p*V_l_4-M_l)*eta/M_p
  print('zeta =',zeta)
  
  # ここで，外形の図心まわりのモーメントを考えると，
  gamma = zeta*h_p*M_p*g/(h_p*f_wing)
  print('gamma =',gamma)


def get_V(z,rho_s,g,P_0,tV_p,V_l_0):
  P_l = rho_s*z*g+P_0
  V_l = V_l_0*P_0/P_l
  ## V_p = V_p_0-(V_l_0-V_l)   
  V_p = tV_p + V_l
  return V_p,V_l

if(__name__ == '__main__'):
  main()
