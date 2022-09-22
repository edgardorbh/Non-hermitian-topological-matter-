Nobo
import matplotlib.pyplot as plt
import cmath
from numpy import*
import time 

def ssh(u,v,w,g,N):
    A=zeros(2*N, dtype=complex)
    B=zeros(2*N, dtype=complex)
    A[-1],A[1],A[3]=v,u+g,w
    B[-2],B[0],B[2]=w,u-g,v
    H=zeros((2*N,2*N),dtype=complex)
    for k in range(N):
        H[2*k]=roll(A,2*k)
        H[2*k+1]=roll(B,2*k)
    H[0,2*N-1]=0.
    H[1,2*N-2]=0.
    H[2*N-1,0]=0.
    H[2*N-2,1]=0.
    return H

u,v,w,g,N=1.1,1,2/5,2/3,70
H=ssh(u,v,w,g,N)
E,EV=linalg.eig(H)
Ros=[]
for e in E:
    R=roots([v*w,(u+g)*v+(u-g)*w,v**2+w**2+(u+g)*(u-g)-e**2,(u-g)*v+(u+g)*w,v*w])
    for i in arange(1,4):
        z=R[i]
        Ros.append(z)
Ros=list(sort(Ros)[::-1])
plt.scatter(real(Ros),imag(Ros),color="red",s=2)
plt.show()

def intersection(z0,z1,z2):
    z0_x,z1_x,z2_x=real([z0,z1,z2])
    z0_y,z1_y,z2_y=imag([z0,z1,z2])
    D=(z2_x-z1_x)*z0_y-(z2_y-z1_y)*z0_x
    if D==0:
        b=False
    else:
        t=(z0_x*z1_y-z0_y*z1_x)/D
        zz=z1+t*(z2-z1)
        if 0<=t<=1 and sign(real(zz))==sign(real(z0)):
            b=True
        else:
            b=False
    return b
                
def curva1(x):
    xx0=x.copy()
    #T=[]
    z0,l,curve=x[0],len(x),[]
    curve.append(z0)
    aux=x[1:]
    dis=list(abs(z0*ones(l-1)-aux))
    j=dis.index(min(dis))
    z1=aux[j]
    aux.remove(z1)
    aux=[z1]+aux
    curve.append(z1)
    t=angle(z1/z0)
    s0=sign(t)
    
    while True:
        xx,l=aux,len(aux) 
        aux=xx[1:]
        dis=list(abs(z1*ones(l-1)-aux))
        j=dis.index(min(dis))
        z2=aux[j]
        aux.remove(z2)
        t=angle(z2/z1)
        s=sign(t)
        if s==s0:
            curve.append(z2)
            aux=[z2]+aux
            break
        else:
            aux=[z1]+aux
            
    for z in [z0,z1,z2]:
        xx0.remove(z)
    xx=[z2]+xx0
    
    while(len(xx)>1):    
        z1,l=xx[0],len(xx)
        aux=xx[1:] 
        dis=list(abs(z1*ones(l-1)-aux))
        j=dis.index(min(dis))
        z2=aux[j]
        aux.remove(z2)
        t=angle(z2/z1)
        s=sign(t)
        if s==s0:
            b=intersection(z0,z1,z2)
            if b==True:
                curve.append(z0) 
                break
            else:
                curve.append(z2)
                aux=[z2]+aux
        else:
            aux=[z1]+aux
        xx=aux
    if curve[-1]!=z0:
        curve=curve+[z0]    
    return(curve,s0)

def curva2(x,s0):
    curve,xx,xx0,z0=[],x.copy(),x.copy(),x[0]
    curve.append(z0) 
    while len(curve)<3:
        aux=xx[1:]
        dis=list(abs(z0*ones(len(aux))-aux))
        j=dis.index(min(dis))
        z1=aux[j]
        aux.remove(z1)
        t=angle(z1/z0)
        s=sign(t)
        if s==-s0:
            curve.append(z1)
            z0=z1
            aux=[z1]+aux
            xx=aux
        else:
            xx.remove(z1)
    for z in curve:
        xx0.remove(z)
    xx=[curve[2]]+xx0
    z0=x[0]
    
    while(len(xx)>1):    
        z1,l=xx[0],len(xx)
        aux=xx[1:] 
        dis=list(abs(z1*ones(l-1)-aux))
        j=dis.index(min(dis))
        z2=aux[j]
        aux.remove(z2)
        t=angle(z2/z1)
        s=sign(t)
        if s==-s0:
            if intersection(z0,z1,z2)==True:
                curve.append(z0) 
                break
            else:
                curve.append(z2)
                aux=[z2]+aux
        else:
            aux=[z1]+aux
            #R.append(z2)
        xx=aux
    if curve[-1]!=z0:
        curve=curve+[z0]  
    return(curve)
    
def anglesort1(x,y):
    if angle(x)<=angle(y):
        return [x,y]
    else:
        return [y,x]
    
def angle_sort(Z):
    ZZ=[]
    for z in Z:
        ZZ.append(z)
    ZZ=array(ZZ)
    I=angle(ZZ).argsort()
    return ZZ[I]
            
def anglesort(U):
    UU=[]
    for uu in U:
        UU.append(uu)
    UU=array(UU)
    V=angle(list(UU))
    I=V[:,0].argsort()
    UU=UU[I]
    return UU

def palos(X):
    P=[]
    l=len(X)
    for i in range(l):
        P.append([X[i],X[(i+1)%l]])
    return P

def to_set(A):
    A_set=set()
    for a in A:
        A_set.add(tuple(a))
    return A_set

def separa(Z):
    DD=[]
    PP=anglesort(Z)
    l=len(PP)
    i=0
    while i<l:
        D=[]
        D.append(PP[i])
        while PP[i%l][1]==PP[(i+1)%l][0]:
            D.append(PP[(i+1)%l])
            i=i+1
        DD.append(D)
        i=i+1   
    if to_set(DD[0]).issubset(to_set(DD[-1]))==True:
        DD=DD[1:]
        
    return DD

def ordena(A,B):
    while A[0][0][0]!=B[0][0][0]:
        A=A[1:]+[A[0]]
    return A,B

def selecciona(A,B):
    l=min(len(B),len(A))
    C=[]
    for i in range(l):
        a1=A[i][0][1]-A[i][0][0]
        b1=B[i][0][1]-B[i][0][0]
        t1=angle(a1/b1)
        if t1>=0:
            C.append(B[i])
        else:
            C.append(A[i])
    return C

def limpia(A,B):
    U=A.union(B)
    I=A.intersection(B)
    R=U.difference(I)
    A1=A.intersection(R)
    B1=B.intersection(R)
    if (A1!=set()) and (B1!=set()):
        AA1=separa(A1)
        BB1=separa(B1)
        AA1,BB1=ordena(AA1,BB1)
        C=selecciona(AA1,BB1)
        for c in C:
            I=I.union(to_set(c))   
    return I

def envolvente(Z):
    c1,s0=curva1(Z)
    c2=curva2(Z,s0)
    c1=angle_sort(c1[0:-1])
    c2=angle_sort(c2[0:-1])
    P1=to_set(palos(c1))
    P2=to_set(palos(c2))
    return limpia(P1,P2)


c1,s0=curva1(Ros)
c2=curva2(Ros,s0)
C=envolvente(Ros)

plt.scatter(real(Ros),imag(Ros),color="black")
plt.scatter([0],[0],color="black")
plt.plot(real(c1),imag(c1),color="blue")
plt.plot(real(c2),imag(c2),color="orange")
for c in C:
    plt.plot(real(c),imag(c),color="red",linewidth=3.)
plt.show()

def enlaces(Z):
    l=len(Z)-1
    bonds=[]
    for z in Z:
        W=Z.remove(z)
        j=min(abs(z*ones(l)-W)).index()
        z_c=W[j]
        bonds.append([z,z_c])

def filtrar(Z):
    QQ=[]
    Q=Z.copy()
    QQ.append(Q[0])
    for z in Q:
        W=Q.remove(z)
        l=len(W)
        dis=lis(abs(z*onles(l)-W))
        j=dis.index(min(dis))
        z_c=W[j]
        if abs(z-z_c)<1e-14:  
            Q.remove(z_c)
            Q.remove(z)
            QQ.append([z,z_c])
    return QQ

u,v,w,g,N=1.1,1,2/5,2/3,70
H=ssh(u,v,w,g,N)
E,EV=linalg.eig(H)
Ros=[]
for e in E:
    R=roots([v*w,(u+g)*v+(u-g)*w,v**2+w**2+(u+g)*(u-g)-e**2,(u-g)*v+(u+g)*w,v*w])
    for i in arange(1,4):
        z=R[i]
        Ros.append(z)
Ros=list(sort(Ros)[::-1])
filtrar(Ros,1e-13)

def filtrar(Z,d):
    Q=set()
    QQ=set()
    for z in Q:
        W=Z.copy()
        W.remove(z)
        l=len(W)
        dis=abs(z*ones(l)-array(W))
        j=dis.index(min(dis))
        x_c=W[j]
        Q.add(tuple([z,x_c]))
    for q in Q:
        q=list(q)
        if abs(q[1]-q[0])<d:
            QQ.add(q[0])
        else:
            QQ.add(q[0])
            QQ.add(q[1])
    return list(QQ)

def enlaces(Z):
    enla=set()
    for z in Z:
        W=Z.copy()
        W.remove(z)
        l=len(W)
        dis=list(abs(z*ones(l)-W))
        j=dis.index(min(dis))
        z_c=W[j]
        enla.add(tuple(sort([x,x_c])))
    enla=array(list(enla))
    return enla 

x=random.rand(3)-.5
y=random.rand(3)-.5
Z=x+1j*y
x=list(x)+[x[0]]
y=list(y)+[y[0]]
z1,z2,z3=Z
z0=random.rand()+1j*random.rand()-.5-.5j
plt.xlim(-.5,.5)
plt.ylim(-.5,.5)
plt.plot(x,y)
if inside(z0,z1,z2,z3)==True:
    plt.scatter(real(z0),imag(z0),color="red")
if inside(z0,z1,z2,z3)==False:
    plt.scatter(real(z0),imag(z0),color="blue")
plt.show() 

def descarta(D):
    z1=0j
    W=D.copy()
    Q=D.copy()
    while len(W)>0:
        d=W[0]
        W.remove(d)
        Q.remove(d)
        for i in [1,2]:
            b=False
            z0=d[i]
            for q in Q:
                if inside(z0,0j,q[0],q[1])==True:
                    W.remove(d)
                    Q.remove(d)
                    b=True
                    break
            if b==True:
                break
    return W

def intersection(z1,z2,w1,w2):
    D=imag(z2-z1)*real(w2-w1)-real(z2-z1)*imag(w2-w1)
    D1=imag(w1-z1)*real(w2-w1)-real(w1-z1)*imag(w2-w1)
    D2=imag(w1-z1)*real(z2-z1)-real(w1-z1)*imag(z2-z1)
    if D==0:
        return False 
    else:
        t,s=D1/D,D2/D
        if (0<t<1) and (0<s<1):
            #z0=z1+t*(z2-z1)
            #w0=w1+t*(w2-w1)
            return True 
        else:
            return False 

def descarta(D):
    W=D.copy()
    Q=D.copy()
    w1=0+0j
    while(len(W)>0):
        d0=W[0]
        z1,z2=d0[0],d0[1]
        W.remove(d0)
        for d in Q:
            for w2 in d:
                if intersection(z1,z2,w1,w2)==True:
                    Q.remove(d)
                    break 
    return Q
def neighbor(z0,Z):
    l=len(Z)
    dis=list(abs(z0*ones(l)-array(Z)))
    j=dis.index(min(dis))
    z_c=Z[j]
    return z_c 

def encuentra(Z,N):
    for z0 in Z:
        L1,L2=[],[]
        Q=Z.copy()
        Q.remove(z0)
        z1=neighbor(z0,Q)
        s0=sign(angle(z1/z0))
        Q.remove(z1)  
        c1,c2=1,0
        while True:
            z_1=neighbor(z0,Q)
            s=sign(angle(z_1/z0))
            Q.remove(z_1)
            if s==s0:
                if c1<N: 
                    c1=c1+1      
                    L1.append(z_1)
            elif s==-s0:
                if c2<N:
                    c2=c2+1
                    L2.append(z_1)
            if (c1>=N) and (c2>=N):
                break
        L=L1+L2
        I=angle(L).argsort()
        L=L[I]
        z_min,z_max=L[0],L[-1]
        b=True
        for q in Q: 
            if trap(q,z_min,z_max)==True:
                b=False
                break 
    if b==True:
        break
    return z

def trap(z0,z1,z2):
    t1=angle(z1/z0)
    t2=angle(z2/z0)
    if sign(t1)==-sign(t2):
        return True 
    elif (t1==0) or (t2==0):
        return True 
    else:
        return False

A=array([1,4,5,3,1,2])
A.argsort()

from numpy import*
import matplotlib.pyplot as plt 
import pandas as pd 

M=random.randint(3,4,(3,3))
M












    


