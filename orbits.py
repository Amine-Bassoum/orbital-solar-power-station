Nouveau ! Raccourcis clavier … 
Les raccourcis clavier de Drive ont été mis à jour pour vous permettre de naviguer à l'aide des premières lettres

#V2 all orbits one satelite + shadow positions (takes a long execution time) with menu =)
#lien utile : http://orbitsimulator.com/formulas/OrbitalElements.html
#imports

import numpy as np
import math as m
#import matplotlib.pyplot as plt
#useful functions

def plan(a,p01,p02):
    R = m.sqrt(p01**2+p02**2)
    if p01 == 0 :
        a1 = 0;   a2 = a
    elif p02 == 0 :
        a1 = a;   a2 = 0
    else :
        alpha = abs(m.asin(p02/R))
        a1 = abs(a * m.cos(alpha))
        a2 = abs(a * m.sin(alpha))
        if p01 > 0 :
            a1 = -a1
        if p02 > 0 :
            a2 = -a2
    return(a1,a2)

def in_shadow (A,B,C,D,E,P):
    n1 = np.cross(B - A, C - A)
    n2 = np.cross(C - A, D - A)
    n3 = np.cross(D - A, E - A)
    n4 = np.cross(E - A, B - A)
    n5 = np.cross(D - B, C - B)
    
    r1 = np.dot(P - A, n1)
    r2 = np.dot(P - A, n2)
    r3 = np.dot(P - A, n3)
    r4 = np.dot(P - A, n4)
    r5 = np.dot(P - B, n5)
    if r1 < 0 or r2 < 0 or r3 < 0 or r4 < 0 or r5 < 0  :
        return(False)
    else :
        return(True)

#conditions initiales des sattelites
#      p0x     p0y p0z v0x     v0y               v0z             T
ni = [[8529600 ,0 ,0 ,0 ,6020.842885982914 ,3269.0509615719548 ,133],#Low MEO
      [8529600 ,0 ,0 ,0 ,-1876.915445942415 ,6588.962866478181 ,133],#Sun Sync Repeat 2
      [7334000 ,0 ,0 ,0 ,-3779.886734454954 ,7548.259388099452 ,181],#HEO Elliptic
      [42157000 ,0 ,0 ,0 ,3081.6852028624658 ,0                ,1437],#Geo
      [6871000 ,0 ,0 ,0 ,6708.2872786680455 ,3642.30281276483  ,102],#LEO
      [8047500,0,0,0,-1574.6510333028455,6875.287453742872     ,120],#sun sunc Repeat 1
      [9090900,0,0,0,-2280.5955967035825,6232.020367705032     ,144],#sun Sunc Repeat 3
      [9754600,0,0,0,-2808.4089657039244,5758.091691108525     ,160],#Sun Sunc Repeat 4
      [13371000,0,0,0,3138.574767284475,4482.34929790841       ,257],#Low Van allen Gap
      [18371000,0,0,0,2677.6143368376547,3824.0295779787984    ,414],#High Van Allen Gap
      ]

#configuration
Dt = 10 #TimeStep

G = 6.67 * 10 ** -11
Mt = 6e24
Ms = 2e30
rs = 696340e3
rt = 6371e3

p0xt = 147493710000
p0yt = 0
p0zt = 0

v0xt = 0
v0yt = 30318.643648737947
v0zt = 835.5913145610228

pt = [np.array([p0xt,p0yt,p0zt])]

pst = 3.8651*10**26
energie = 0

orbit = int(input('''Orbit Selection : 
        0/ Low MEO
        1/ Sun Sync Repeat 2
        2/ HEO Elliptic
        3/ GEO
        4/ LEO
        5/ Sun Sync Repeat 1
        6/ Sun Sync Repeat 3
        7/ Sun Sync Repeat 4
        8/ Low Van Allen Gap
        9/ High Van Allen Gap
        '''))
        
print("Warning Reseting is necessary if the TimeStep has changed\n Default DT = 15")
reset = int(input("Reset Orbital Data ? \nYes : 1    No : 0\n"))

#primary sun secondary earth
e1 = 0
e2 = 0
P1 = False
P2 = False
recreateE = False
recreateO = False
if reset == 0 :
    
    while(e1 == 0 and P1 == False): 
        try :
            print("loading Earth's orbit\npress Ctrl+C to interupt command")
            pt = np.loadtxt('earth.csv', delimiter=',')
            P1 = True
        except :
            e1 = int(input('Error : File "earth.csv" not found / Corrupted ...\n 1 to create a new one \n 0 to try again\n t to terminate\n'))
            if e1 == 1 :
                recreateE = True
    
    while(e2 == 0 and P2 == False):
        try :
            print("loading Satellite's orbit\npress Ctrl+C to interupt command")
            p = np.loadtxt(f'orbit{orbit}.csv',delimiter=',')
            P2 = True
        except :
            e2 = int(input(f'Error : File "orbit{orbit}.csv" not found / Corrupted ...\n 1 to create a new one\n 0 to try again\n t to terminate\n'))
            if e2 == 1 :
                recreateO = True

if recreateE == True or reset == 1 :
    print('making Earth Orbit')
    for n in range(365*24*60*60//Dt):
    
        R = m.sqrt(p0xt**2+p0yt**2+p0zt**2)
        a = Ms * G / (R**2)
        
        if p0xt == 0 : #probleme plan sur Y, Z
            ax = 0
            ay , az = plan(a,p0yt,p0zt)
        elif p0yt == 0 : # probleme plan sur Z, X
            ay = 0
            az , ax = plan(a,p0zt,p0xt)
        elif p0zt == 0 : #probeleme plan sur X,Y
            az = 0
            ax , ay = plan(a,p0xt,p0yt)
        else :
            theta = m.atan2(p0yt,p0xt)
            alpha = m.atan2(m.sqrt(p0xt**2+p0yt**2),p0zt)
            ax = abs(a * m.cos(theta) * m.sin(alpha))
            ay = abs(a * m.sin(theta) * m.sin(alpha))
            az = abs(a * m.cos(alpha))
            
            if p0xt > 0 :
                ax = -ax
            
            if p0yt > 0 :
                ay = -ay
            
            if p0zt > 0 :
                az = -az 
        
        p0xt = p0xt + v0xt * Dt + (ax * Dt**2)/2
        p0yt = p0yt + v0yt * Dt + (ay * Dt**2)/2
        p0zt = p0zt + v0zt * Dt + (az * Dt**2)/2
                
        v0xt = v0xt + ax * Dt
        v0yt = v0yt + ay * Dt
        v0zt = v0zt + az * Dt
        
        pt.append(np.array([p0xt,p0yt,p0zt]))
        print(f"{(n/(364*24*60*60/Dt))*100}% complete")
        
    np.savetxt('earth.csv', pt, delimiter=',')
        
#satellite : Sat / primary : earth

if recreateO == True or reset == 1 :
    print("making Satellite orbit")
    p = [np.array([ni[orbit][0],ni[orbit][1],ni[orbit][2]],dtype = float)]
        
    for n in range(ni[orbit][6]*60//Dt):
        
        R = m.sqrt(ni[orbit][0]**2+ni[orbit][1]**2+ni[orbit][2]**2)
        a = Mt * G / (R**2)
        
        if ni[orbit][0] == 0 : #probleme plan sur Y, Z
            ax = 0
            ay , az = plan(a,ni[orbit][1],ni[orbit][2])
        elif ni[orbit][1] == 0 : # probleme plan sur Z, X
            ay = 0
            az , ax = plan(a,ni[orbit][2],ni[orbit][0])
        elif ni[orbit][3] == 0 : #probeleme plan sur X,Y    
            az = 0
            ax , ay = plan(a,ni[orbit][0],ni[orbit][1])
        else :
            theta = m.atan2(ni[orbit][1],ni[orbit][0])
            alpha = m.atan2(m.sqrt(ni[orbit][0]**2+ni[orbit][1]**2),ni[orbit][2])
            ax = abs(a * m.cos(theta) * m.sin(alpha))
            ay = abs(a * m.sin(theta) * m.sin(alpha))
            az = abs(a * m.cos(alpha))
            
            if ni[orbit][0] > 0 :
                ax = -ax
            
            if ni[orbit][1] > 0 :
                ay = -ay
            
            if ni[orbit][2] > 0 :
                az = -az 
        
        ni[orbit][0] = ni[orbit][0] + ni[orbit][3] * Dt + (ax * Dt**2)/2
        ni[orbit][1] = ni[orbit][1] + ni[orbit][4] * Dt + (ay * Dt**2)/2
        ni[orbit][2] = ni[orbit][2] + ni[orbit][5] * Dt + (az * Dt**2)/2
        
        ni[orbit][3] = ni[orbit][3] + ax * Dt
        ni[orbit][4] = ni[orbit][4] + ay * Dt
        ni[orbit][5] = ni[orbit][5] + az * Dt
        
        p.append(np.array([ni[orbit][0],ni[orbit][1],ni[orbit][2]],dtype = float))
        print(f"{n/(ni[orbit][6]*60/Dt)*100}% complete")
    #ajustement temporelle & spatial
    
    multiple = round(len(pt)/len(p))+1
    p *= multiple
    del p[len(pt):]
    
    for s in range(len(pt)):
        p[s][0] += pt[s][0]
        p[s][1] += pt[s][1]
        p[s][2] += pt[s][2]
        
    np.savetxt(f'orbit{orbit}.csv', pt, delimiter=',')

# #plotting
# x=[]
# y=[]
# z=[]
# for i in range(len(pt)):
#     x.append(p[i][0])
#     y.append(p[i][1])
#     z.append(p[i][2])
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter3D(x,y,z,color = 'hotpink',s = 1)

#power increment & shadows part2
#initial shadow points :

        
B1 = np.array([0,rs,0])
C1 = np.array([0,-rs,0])
B = np.array([p0xt,rt,0])
C = np.array([p0xt,-rt,0])
D = np.array([p0xt,0,rt])
E = np.array([p0xt,0,-rt])
A = np.array([rs*p0xt/(rs-rt),0,0]) #intersection of B1 B and C1 C
DST = m.sqrt(p[0][0]**2+p[0][1]**2+p[0][2]**2)

if in_shadow(A,B,C,D,E,p[0]) == False :
    Pss = pst/(4 * 3.14 * DST**2) #puissance soleil-satellite /m²
    energie += Pss * Dt               #/metre²
        
#mouvement of shadow position (rotate points):
s = 1
while(s<len(pt)):
    DST = m.sqrt(p[s][0]**2+p[s][1]**2+p[s][2]**2)
    v = pt[s] - pt[s-1]
    A += v
    B += v
    C += v
    D += v
    E += v
    if in_shadow(A,B,C,D,E,p[s]) == False :
        Pss = pst/(4 * 3.14 * DST**2) #puissance soleil-satellite /m²
        energie += Pss * Dt #/metre² it's also the fitness score for the neat algorithm... maybe
    s += 1
    print(f"{s/len(pt)*100}% complete")
print(f"Energie annuelle = {energie}J/m²")
print(f"Puissance annuelle moyenne = {energie/(365*24*60*60)}W/m²")
