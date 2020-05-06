# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 20:48:49 2020

@author: Daisuke
"""

#課題1 等加速度運動のシミュレーション　
# ポイント　リストの使い方、ループの書き方

v = [ 0 for i in range(100)]
x = [ 0 for i in range(100)]

g = 9.80665
v[0] = 10
i = 1
 
while i < 100:
    v[i] = v[i-1] - g/100
    x[i] = x[i-1] + v[i]/100
    i = i + 1

for i in range(100):
    print("x位置 : %4f v速度 : %4f" % (x[i],v[i]))
    
'''
配列をnp.array等で実装すれば、matplotlib等で簡単に運動を可視化することが可能です。
実際は常微分方程式のシミュレートはオイラー法ではなく、ルンゲクッタ法を用いるべきです。
'''

#課題2
#(1)　
# ポイント　外部ライブラリの使い方
    
import math
print(112 + math.sqrt(math.pi)*abs(math.exp(2)))

#(2) 解は0.3413程度
# X~N(0,1) のとき、P(0<X<1)なる確率を求める積分計算を行う
# f(x)が標準正規分布であることに気づけば、np.random.randnを生成し、0<x<1の割合を数えてもよい

import numpy as np

def fun(x):
    return(np.exp(-pow(x,2)/2)/(np.sqrt(2*np.pi)))

# モンテカルロ法
from numpy import random

fx = fun(random.rand(1000000))
S = sum(fx)/1000000
print(S) #0.3412979230495587

# 台形公式
S = (fun(0) + fun(1))*1/2
print(S) #0.320456502460288

# 100分割した台形公式
S = 0
for i in range(100):
    S += (fun(i/100) + fun((i+1)/100))/200 #区間を100等分
print(S) #0.34134272963911727
 
# シンプソンの公式
S = (1-0)*(fun(0) + 4*fun((0+1)/2) + fun(1)) / 6
print(S) #0.3415290519962957

#課題3
#(1) 解答35430270439
# Chinese Reminder Theoremである。高校の整数分野でありがちな問題
# 計算量を減らすための工夫?を考えよう

print(1584891 * 3438478) #5449612835898間隔で解が存在

a = [i for i in range(3438478) if (i*1584891+32134)%3438478 == 193127]
print(a[0]*1584891+32134)

#(2) 
# 線型合同法に関する設問 疑似乱数のイメージをつかんでいただく。
# seedにはepochからの起動時間を使うことが多い。

def LCG(A = 48271,B = 0,M = pow(2,31)-1,n = 100,seed = 1111):
    x = [0 for i in range(n)]
    x[0] = seed
    i = 0
    for i in range(n):
        if (i == 99):
            break
        x[i+1] = (A*x[i]+B) % M
    return(x)

print(LCG())
        
#(3) おまけ
def saikoro(seed=111):
    x = LCG(seed = seed)
    [ print(i%6 + 1) for i in x[90:] ]


# 課題4
# (1)　幾何ブラウン運動のシミュレーション
'''
確率微分方程式は解析的に解くには、伊藤の公式などの確率解析の知識が必要だが、
シミュレートするだけなら常微分方程式の時に、正規分布の乱数の要素を追加するだけ!
更にシミュレーションなら解析的に難しいエキゾチックオプションのプライシングも楽にできる~
ただし、同値マルチンゲール測度でドライブすること、muではなくリスクフリーレートr、ボラティリティは同じ。
'''

from numpy.random import randn
import numpy as np

S = [0 for i in range(1000)]
S[0] = 100
mu = 0.05
sigma = 0.4

i = 1
while  i <= 999:
    S[i] = S[i-1] + mu*S[i-1]/1000 + sigma*S[i-1]*randn()*np.sqrt(1/1000)
    i = i + 1

print(S[999])

# (2)
print(S[0]*np.exp(mu))
# 計算すればわかるが、S(T)の期待値は、S(0)e^(rT)である。これは実質正規分布の積率母関数を求める問題であり容易
# あるいは同値マルチンゲールでは取引可能資産の期待収益率は全て共通にリスクフリーレートに等しいことを思い出せば自明である
