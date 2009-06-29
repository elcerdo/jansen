#!/usr/bin/env python

import scipy as s
import scipy.linalg as l
from pylab import *
from matplotlib.widgets import Slider

def build_params(r,a,b,ap,bp,d):
    return [r,a,b,ap,bp,d]

def resolve_tri(A,B,a,b,up=True):
    AB=A-B
    c=l.norm(AB)
    aa=s.arctan2(AB[1],AB[0])
    bb=s.arccos((b**2+c**2-a**2)/(2*b*c))
    if up: return B+b*s.array((s.cos(aa+bb),s.sin(aa+bb)))
    else: return B+b*s.array((s.cos(aa-bb),s.sin(aa-bb)))

def resolve_leg(B,C,D,d):
    uCB=C-B
    uCB/=l.norm(uCB)
    return D-d*uCB

def plot_leg(theta,params):
    r,a,b,ap,bp,d=params
    A=s.array([1+r*s.cos(theta),r*s.sin(theta)]).transpose()
    B=s.array([0,0])
    C=resolve_tri(A,B,a,b)
    CC=C-B
    CC=s.array((-CC[1],CC[0]))
    CC/=l.norm(CC)
    CC*=.4
    BB=B+CC
    D=resolve_tri(A,B,ap,bp,up=False)
    DD=D+CC
    E=resolve_leg(B,C,D,d)
    points=s.array([A,C,B,D,E,DD,BB,C,B,BB,DD,D,A])
    return plot(points[:,0],points[:,1])

def update_leg(leg,theta,params):
    r,a,b,ap,bp,d=params
    A=s.array([1+r*s.cos(theta),r*s.sin(theta)]).transpose()
    B=s.array([0,0])
    C=resolve_tri(A,B,a,b)
    CC=C-B
    CC=s.array((-CC[1],CC[0]))
    CC/=l.norm(CC)
    CC*=.4
    BB=B+CC
    D=resolve_tri(A,B,ap,bp,up=False)
    DD=D+CC
    E=resolve_leg(B,C,D,d)
    points=s.array([A,C,B,D,E,DD,BB,C,B,BB,DD,D,A])
    leg.set_xdata(points[:,0])
    leg.set_ydata(points[:,1])

def compute_trajectory(params):
    r,a,b,ap,bp,d=params
    thetas=s.linspace(0,2*s.pi,256,endpoint=False)
    A=s.array((1+r*s.cos(thetas),r*s.sin(thetas))).transpose()
    C=s.array([resolve_tri(Ak,[0,0],a,b) for Ak in A])
    D=s.array([resolve_tri(Ak,[0,0],ap,bp,up=False) for Ak in A])
    E=s.array([resolve_leg([0,0],Ck,Dk,d) for Ck,Dk in zip(C,D)])
    return A,C,D,E

def analyse_trajectory(params):
    r,a,b,ap,bp,d=params
    foo,foo,foo,traj=compute_trajectory(params)

    #total length of tubes
    total_length=a+b+ap+2*bp+d

    #distance walked for each step
    step=traj[:,0].max()-traj[:,0].min()

    #mode quality
    speed=[]
    speed.append(traj[-1,:]-traj[1,:])
    for k in xrange(1,traj.shape[0]-1):
        speed.append(traj[k-1,:]-traj[k+1,:])
    speed.append(traj[-2,:]-traj[0,:])
    speed=s.array(speed)

    signe_change=0
    mode={}
    mode[-1]=None #contact mode
    mode[1]=None  #come back
    k=0
    while k<speed.shape[0]:
        ss=sign(speed[k,0])
        current=[traj[k,:]]
        k+=1

        while k<speed.shape[0] and ss==sign(speed[k,0]):
            current.append(traj[k,:])
            k+=1

        if signe_change<2:
            mode[ss]=s.array(current)
        elif signe_change==2:
            mode[ss]=s.vstack((current,mode[ss]))
        else:
            assert(False)

        signe_change+=1

    #crop transition
    mode[-1]=mode[-1][10:-10,:]

    #mode quality
    flatness=mode[-1][:,1].std()
    contact_height=mode[-1][:,1].mean()
    bad_comeback=0.
    for pos in mode[1]:
        if pos[1]<contact_height:
            bad_comeback+=20
        elif pos[1]<contact_height+.1:
            bad_comeback+=1.5
    bad_comeback/=mode[1].shape[0]
    
    score=-(total_length-6.)/2+step/1.-flatness/1e-2-bad_comeback
    print "length: %f step: %f flatness: %f height: %f comeback: %f score: %f" % (total_length,step,flatness,contact_height,bad_comeback,score)
    return score


params=build_params(r=.2,a=1.6,b=.8,ap=1.3,bp=.6,d=1.)
analyse_trajectory(params)

figure()
subplot(111)
subplots_adjust(bottom=0.2)

axtheta=axes([0.1, 0.05, 0.8, 0.1])
stheta=Slider(axtheta,'t',0,4*s.pi,valinit=s.pi/3)

subplot(111)
A,C,D,E=compute_trajectory(params)
plot(A[:,0],A[:,1])
plot(C[:,0],C[:,1])
plot(D[:,0],D[:,1])
plot(E[:,0],E[:,1])
leg,=plot_leg(s.pi/3,params)
axis("image")
xlim(-1,1.5)
ylim(-2,1)

def update(val):
    update_leg(leg,val,params)
    draw()

stheta.on_changed(update)
show()


