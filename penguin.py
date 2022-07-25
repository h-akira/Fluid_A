#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created: May, 20, 2022 14:27:54 by Hiroto Akira
# $Author: $
# $Date: $
# $URL: $
__giturl__ = "$URL: $"

def get_V(z,rho_s,g,P_0,tV_p,V_l_0):
  P_l = rho_s*z*g+P_0
  V_l = V_l_0*P_0/P_l
  ## V_p = V_p_0-(V_l_0-V_l)   
  V_p = tV_p + V_l
  return V_p,V_l

def main():
  import numpy
  # 海水
  T_s = 0  # 温度 ℃
  rho_s = 1028  # 密度 kg/m^3
  mu_s = 0.00185  # 粘度 Pa s

  # 空気
  T_p = 38  # 体温 ℃
  ## rho_a = 1.0985  # 密度 kg/s^3
  rho_a = 0.8*1.091+0.2*1.166  # 密度 kg/s^3
  # 以下のサイトにおける30度と40度の密度を内分した
  # https://www.hakko.co.jp/qa/qakit/html/h01040.htm

  print('(1) 海水の温度（文献値） =',T_s,'℃')
  print('(1) 海水の密度 （文献値） =',rho_s,'kg/m^3')
  print('(1) 海水の粘度 （文献値） =',mu_s,'Pa・s')

  # ペンギンの肺
  r_l_0 = 0.12  # 半径 m
  V_l_0 = 4/3*numpy.pi*r_l_0**3  # 肺の体積kg/m^3
  M_l = rho_a*V_l_0  # 肺の空気の重さ（一定）
  eta = 0.0025

  # ペンギン
  h_p = 1.1  # 身長 m
  d_p = 0.35  # 直径 m
  rho_p = rho_s*1.01  # 密度 kg/m^3
  V_p_0 = numpy.pi/6*h_p*d_p**2  # 外見の体積 m^3
  tV_p = V_p_0-V_l_0  # ペンギンの身の体積（一定） m^3
  tM_p = rho_p*tV_p  # ペンギンの身の重さ（一定） kg
  M_p = tM_p+M_l
  print('(1) ペンギンの身の密度（一定） =',rho_p,'kg')
  print('(1) ペンギンの身の体積（一定） =',tV_p,'kg')
  print('(1) ペンギンの身の体重（一定） =',tM_p,'kg')
  print('(補) ペンギンの身長h_p =',h_p,'m')
  print('(補) ペンギンの腹の直径d_p =',d_p,'m')
  print('(補) 大気圧，体温における空気の密度=',rho_s,'kg/m^3')
  print('(補) 地上における肺の直径 =',r_l_0,'m')
  print('(補) 地上における肺の体積 =',V_l_0,'m^3')
  print('(補) 肺の空気の重さ（深度によらず一定） =',M_l,'kg')
  u_max = 35*1000/3600  # ペンギンの最大速度 m/s
  Re = rho_s*u_max*h_p/mu_s
  print('(2) ペンギンの最大速度u_max（文献値） =',u_max,'m/s')
  print('(2) Re =',Re)
  
  b = 0.2  # 羽の長さ m
  c = 0.08  # 羽の幅

  ## print(V_l_0/tV_p)

  # その他
  P_0 = 101325  # 大気圧 Pa
  g = 9.81  # 重力加速度 m/s^2

  # 条件を満たしているかの確認
  ## print('='*10,'zが正になる条件を','='*10)
  ## print(rho_s*tV_p)
  ## print(M_p)
  ## print(rho_s*(tV_p+V_l_0))
  ## print('='*20)

  # 問3
  # 釣り合いの位置z
  ## z_balance = (rho_s*V_l_0-M_l)/((rho_p-rho_s)*tV_p+M_l)*P_0/(rho_s*g)  # 古い式(15)=間違い
  ## z_balance = (rho_s*V_l_0/(M_p-rho_s*tV_p)-1)*P_0/(rho_s*g)  # 古い式(14)
  z_balance = (rho_s*(tV_p+V_l_0)-M_p)/(M_p-rho_s*tV_p)*P_0/(rho_s*g)  # 新しい式(16)
  print('(3) z =',z_balance)
  
  # 問4
  # 巡航速度を決める
  u_p = u_max*0.8
  print('='*25+'(4)'+'='*25)
  print('ペンギンの巡航速度を0.8*u_max =',u_p,'m/s と仮定すると，')
  # まずは上下方向の釣り合いを考える
  V_p_4,V_l_4 = get_V(0.99*z_balance,rho_s,g,P_0,tV_p,V_l_0)
  print('0.99*z_balanceにおける肺の体積はボイルの法則より',V_l_4,'m^3 となり，')
  print('ペンギンの身の体積は一定より，外形の体積はそれとの足し算で',V_p_4,'m^3 となる．')
  f_wing = (rho_s*V_p_4-M_p)*g  # 上下方向に釣り合うために羽から得なければならない下向きの力（二枚分） N
  print('よって，上下方向の力のつりあいより，')
  print('羽根から得る必要がある下向きの力は',f_wing,'N となる．')
  print('一枚の羽根の長さb =',b,'m，')
  print('一枚の羽根の幅c =',c,'m と仮定し，羽根が二枚あることに注意して，')
  C_L = f_wing/(0.5*rho_s*2*c*b*u_p**2)
  print('揚力係数C_L =',C_L)
  alpha = C_L/2/numpy.pi
  print('-----[','羽根の角度alpha =',alpha,'rad?',']-----')

  print('次に，モーメントを考える．新しい式(21)=古い式(19)より，')
  zeta = (rho_p*V_l_4-M_l)*eta/M_p
  print('zeta =',zeta)
  print('これを用いて，外形の図心まわりのモーメントを考えることで')
  gamma = zeta*h_p*M_p*g/(h_p*f_wing)
  print('-----[','gamma =',gamma,']-----')

  # 問5
  print('='*25+'(5)'+'='*25)
  print('ペンギンの身長と身の密度が一定である場合，肺の体積の減少分，外形は痩せる．')
  print('これに基づくと，この深度におけるペンギンの腹の直径は')
  d = (V_p_4/h_p/numpy.pi*6)**0.5
  print('d =',d,'m となり，これを用いることで')
  A_ref = numpy.pi*d**2/4
  print('A_ref =',A_ref)
  print('また，')
  C_D = 0.1
  C_Di = C_L**2/(numpy.pi*b/c)
  print('C_D =',C_D)
  print('C_Di =',C_Di)
  print('これらを用いて，')
  F_D = 0.5*rho_s*C_D*A_ref*u_p**2
  F_Di = 0.5*rho_s*C_Di*2*c*b*u_p**2
  print('F_D =',F_D,'N')
  print('F_Di =',F_Di,'N')
  print('以上より，')
  print('-----[','推進に必要なpower =',(F_D+F_Di)*u_p,'W',']-----')
  
  # 問6
  print('='*25+'(6)'+'='*25)
  print('ペンギンが受ける上向きの力は，重力と浮力の合力より計算される．')
  print('浮力を計算するために各深度におけるペンギンの体積を計算すると，')
  V_p_6_075,V_l_6_075 = get_V(0.75*z_balance,rho_s,g,P_0,tV_p,V_l_0)
  V_p_6_025,V_l_6_025 = get_V(0.25*z_balance,rho_s,g,P_0,tV_p,V_l_0)
  print('z=0.75zのとき',V_p_6_075,'m^3')
  print('z=0.25zのとき',V_p_6_025,'m^3')
  print('これを用いて，ペンギンが受ける上向きの力は')
  F1 = rho_s*V_p_6_075*g - M_p*g
  F2 = rho_s*V_p_6_025*g - M_p*g
  print('z=0.75zのとき',F1,'N')
  print('z=0.25zのとき',F2,'N となる．')
  print('仕事とエネルギーの関係を用いて考えると，水面に出るまでにペンギンが受ける仕事は')
  W = F1*(0.5*z_balance) + F2*(0.5*z_balance)  # J
  print('W =',W,'J')
  print('水面に出た時のペンギンの運動エネルギーはWに等しい．')
  print('ゆえに，45°に飛び出た後の最大の高さは')
  h = 0.5*W/(M_p*g)
  print('-----[','h =',h,'m',']-----')

  # 先生の式の確認
  s_p = rho_p/rho_s
  s_z = rho_a/rho_s
  l_z = V_l_0/tV_p
  beta = 0.99
  ## gamma_check = (s_p-(1-beta)*s_z)*(((s_p-1)+l_z*s_z)-beta*s_z*(2+l_z))/((1-beta)*((s_p-1)+l_z*s_z)*(l_z*(1-s_z)-(s_p-1)))*l_z*eta
  gamma_check = (s_p-(1-beta)*s_z)*(((s_p-1)+l_z*s_z)-beta*s_z*(0+l_z))/((1-beta)*((s_p-1)+l_z*s_z)*(l_z*(1-s_z)-(s_p-1)))*l_z*eta
  ## gamma_check = (s_p-(1-beta)*s_z)*(((s_p-1)+l_z*s_z)-beta*s_z*(2+l_z))/((1-beta)*((s_p-1)+l_z*s_z)*(l_z*(1-s_z)-(s_p-1)))*l_z*0.003
  print('gamma_ =',gamma)
  print('gamma_check =',gamma_check)

  beta = 0.75
  Fbg_child = (1-beta)*(M_p-rho_s*tV_p)*(rho_s*(tV_p+V_l_0)-M_p)
  Fbg_mother = (1-beta)*M_p+beta*rho_s*(tV_p+V_l_0)-rho_s*tV_p
  Fbg = Fbg_child/Fbg_mother*g
  print(Fbg,F1)
  child = (1-beta)*((s_p-1)+l_z*s_z)*(l_z*(1-s_z)-(s_p-1))
  ## mother = (1-beta)*((s_p-1)+l_z*s_z)+beta*(2+l_z)
  mother = (1-beta)*((s_p-1)+l_z*s_z)+beta*(l_z)
  ## mother = (1-beta)*(s_p+l_z*s_z)+beta*(1+l_z)-1
  Fbg = child/mother/(s_p+l_z*s_z)*M_p*g
  print(Fbg,F1)
  ## Fbg=Fbg*(tV_p*rho_s)**0.5
  ## print(Fbg,F1)

## def get_Fbg(beta,M_p,):




if(__name__ == '__main__'):
  main()
