import numpy as np
import pandas as pd
import pylab as plt


def Zr(i,k):
    if community[i][k] in [7,5,3]:
        return 500  #cm
    elif community[i][k] in [2,1]:
        return 100  #cm
    else:
        return 0.1

def I(S,i,k):  # Precipitation
    return min(n*Zr(i,k)*(1-S),sheet['P'][day+1]*2.54/500)

def L(S):  # Leakage
    return (Ksat/500)*(np.exp(beta*(S-Sfc))-1)/(np.exp(beta*(1-Sfc))-1)

def Sw(i,k):  # S_w
    if community[i][k]==7:
        return 0.15
    elif community[i][k]==5:
        return 0.195
    elif community[i][k]==3:
        return 0.195
    elif community[i][k]==2:
        return 0.21
    elif community[i][k]==1:
        return 0.135

def Tmax(i,k):  # T_max
    if community[i][k]==7:
        return 8.128e-4
    elif community[i][k]==5:
        return 5.08e-4
    elif community[i][k]==3:
        return 5.08e-4
    elif community[i][k]==2:
        return 4.064e-4
    elif community[i][k]==1:
        return 9.144e-4


def alpha(i,k):  # Shielding Coefficient
    if 1<=i<=48 and 1<=k<=48:
        if community[i+1][k-1]-community[i][k]>=2 or\
            community[i+1][k]-community[i][k]>=2 or \
            community[i + 1][k+1] - community[i][k] >= 2 or\
            community[i][k-1]-community[i][k]>=2 or \
            community[i][k+1]-community[i][k]>=2 or\
            community[i-1][k-1]-community[i][k]>=2 or\
            community[i-1][k]-community[i][k]>=2 or\
            community[i-1][k+1]-community[i][k]>=2:
            return 0.5
        else:
            return 1
    else:
        return 1

def water_endure(i,k):
    if community[i][k]==7:
        return 180  # days
    if community[i][k] == 5:
        return 65
    if community[i][k] == 3:
        return 35
    if community[i][k] == 2:
        return 35
    if community[i][k] == 1:
        return 90

def T(S,i,k,dead):  # Transpiration
    if community[i][k] in [3,5,7]:
        if S<=Sw(i,k):
            water_deficit[i][k]+=1
            if water_deficit[i][k]>=water_endure(i,k):
                dead[community[i][k]] += 1
                community[i][k]=0

            return 0
        elif Sw(i,k)<=S<=Sstar:
            water_deficit[i][k]=0
            return (S-Sw(i,k))/(Sstar-Sw(i,k))*Sw(i,k)*Tmax(i,k)*alpha(i,k)
        elif S>=Sstar:
            water_deficit[i][k] = 0
            return Tmax(i,k)*alpha(i,k)
    elif community[i][k] in [1,2]:
        if sheet['P'][day+1]>0:
            if S <= Sw(i, k):
                water_deficit[i][k] += 1
                if water_deficit[i][k] >= water_endure(i, k):
                    dead[community[i][k]] += 1
                    community[i][k] = 0

                return 0
            elif Sw(i, k) <= S <= Sstar:
                water_deficit[i][k] = 0
                return (S - Sw(i, k)) / (Sstar - Sw(i, k)) * Sw(i, k) * Tmax(i, k) * alpha(i, k)
            elif S >= Sstar:
                water_deficit[i][k] = 0
                return Tmax(i, k) * alpha(i, k)
        else:
            if (S-0.8)/0.2<=Sw(i,k):
                water_deficit[i][k] += 1
                if water_deficit[i][k]>=water_endure(i,k):
                    dead[community[i][k]] += 1
                    community[i][k]=0

                return 0
            elif Sw(i, k) <=(S-0.8)/0.2 <= Sstar:
                water_deficit[i][k] = 0
                return ((S-0.8)/0.2 - Sw(i, k)) / (Sstar - Sw(i, k)) * Sw(i, k) * Tmax(i, k) * alpha(i, k)
            elif (S-0.8)/0.2>=Sstar:
                water_deficit[i][k] = 0
                return Tmax(i, k) * alpha(i, k)
    else:
        return 0


def E(S,i,k):  # Evaporation
    if community[i][k] in [3,5,7,0]:
        if S<=Sh:
            return 0
        elif Sh<S<Sstar:
            return (S-Sh)/(Sstar-Sh)*Emax*cita()/500
        elif S>=Sstar:
            return Emax * cita()/500
    elif community[i][k] in [1,2]:
        if (S-0.8)/0.2<=Sh:
            return 0
        elif Sh<(S-0.8)/0.2<Sstar:
            return ((S-0.8)/0.2-Sh)/(Sstar-Sh)*Emax*cita()/500
        elif (S-0.8)/0.2>=Sstar:
            return Emax * cita()/500

def cita():  # Temperature Coefficient
    return (sheet['T'][day+1]-32)/(95-32)


def water_discharge(S):
    Snew=[[np.nan for _ in range(50)] for _ in range(50)]
    for i in range(50):
        for k in range(50):
            if 1<=i<=48 and 1<=k<=48:
                Snew[i][k]=np.mean([S[i][k],S[i][k+1],S[i][k-1],S[i+1][k-1],S[i+1][k],S[i+1][k+1],S[i-1][k-1],S[i-1][k],S[i-1][k+1]])
            else:
                Snew[i][k]=S[i][k]
    return Snew.copy()


sheet=pd.read_excel(r'GZ precipitation.xlsx')  # weather data
sheet.drop([0],inplace=True)
sheet['P']=sheet['P'].astype(float)
dead=[0,0,0,0,0,0,0,0]

# constant
n=0.42
Ksat=109.8  # cm/day
beta=9
Sfc=0.29
Sstar=0.105
Sh=0.02
Emax=0.15  # cm/day
S0=0.2


community=[[0 for _ in range(50)] for _ in range(50)]  # Note down plants distribution
water_deficit=[[0 for _ in range(50)] for _ in range(50)]

# Randomly place the plants
np.random.seed(1)
for i in range(len(community)):
    for k in range(len(community[0])):
        community[i][k]=np.random.choice([*[0 for _ in range(18)],7,5,5,5,3,3,3,3,2,2,2,2,2,2,2,2,2,2])


S=[[S0 for _ in range(50)] for _ in range(50)]  # Note down Saturation
temp=[]
# Main Function
for day in range(365):
    for i in range(len(S)):
        for k in range(len(S[0])):
            S[i][k]=S[i][k]+I(S[i][k],i,k)/(n*Zr(i,k))
            S[i][k]=S[i][k]+(-L(S[i][k])-T(S[i][k],i,k,dead)-E(S[i][k],i,k))/(n*Zr(i,k))
            if S[i][k]<=0:
                S[i][k]=0
    S=water_discharge(S)
    print(S[30][30],day)
    temp.append(np.mean(S[25]))

# Death rate output
print('dead=',dead)
print(np.mean(temp))


# Figure Output

plt.figure()
for i in range(len(S)):
    for k in range(len(S[0])):
        if S[i][k]>=0.4:
            plt.plot(i,k,'o',c='#341CD4',markersize=3)  # dark blue
        elif S[i][k]>=0.35:
            plt.plot(i, k,'o', c='#2181FF',markersize=3)  # light blue
        elif S[i][k]>=0.3:
            plt.plot(i, k, 'o',c='#1AEEFF',markersize=3)  # more light blue
        elif S[i][k]>=0.25:
            plt.plot(i, k, 'o',c='#C7FFF1',markersize=3)  # light green
        elif S[i][k]>=0.2:
            plt.plot(i, k, 'o',c='#EBFFCB',markersize=3)  # more light yellow
        elif S[i][k]>=0.15:
            plt.plot(i, k, 'o',c='#FEFFBD',markersize=3)  # light yellow
        elif S[i][k]>=0.1:
            plt.plot(i, k, 'o',c='#FFF286',markersize=3)  # yellow
        else:
            plt.plot(i, k,'o', c='#FFB260',markersize=3)  # orange

plt.figure()
for i in range(len(S)):
    for k in range(len(S[0])):
        if community[i][k]==7:
            plt.plot(i,k,'o',c='#1A8C14',markersize=3)  # dark green
        elif community[i][k]==5:
            plt.plot(i, k,'o', c='#615C3A',markersize=3)  # brown
        elif community[i][k]==3:
            plt.plot(i, k, 'o',c='#F1F233',markersize=3)  # yellow
        elif community[i][k]==2:
            plt.plot(i, k, 'o',c='#99FF80',markersize=3)  # light green
        elif community[i][k] == 1:
            plt.plot(i, k,'o', c='#FF4F4F',markersize=3)  # red


plt.figure()
plt.plot(np.arange(1,len(temp)+1),temp,'-',color='skyblue',label='Water Saturation')
plt.xlabel('Days of a Year')
plt.ylabel('Saturation')
plt.legend()
plt.show()
