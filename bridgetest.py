#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 13:03:26 2018

@author: nikoogivehchian
"""

import bridgesss

#### kk time for my actual bridge!! ####
def main():
    obj = Bridges(15+1060+190)
    #obj.setI(1,"all")
    obj.setPins(15,15+1060+15)
    """ i set the external forces, but the reaction forces are calculated by the program. """
    #"""
    obj.setforce(-66.7,15)
    obj.setforce(-66.7,15+52+176)
    obj.setforce(-66.7,15+52+176+164)
    obj.setforce(-66.7,15+52+176*2+164)
    obj.setforce(-66.7,15+52+176*2+164*2)
    obj.setforce(-66.7,15+52+176*3+2*164)
    #"""
    #obj.setforce(-66.7,1060)
    #obj.setforce(-66.7,1060-176)
    #obj.setforce(-66.7,1060-176-164)
    ###
    #obj.calcPinsF()
    #obj.printSFD()
    #obj.printBMD()
    
    #Ls = [(-400/6) for i in range(0,6)]
    #xs = [52, 52+176, 52+176+164, 52+176*2+164, 52+176*2+164*2, 52+176*3+164*2, 52+176*3+164*2+52]
    ## WAIIIT WTF IS THIS SPACING
    #obj.simTrain(Ls,xs)
    #obj.Baldwin(100)
    
    #lol shoulda accumed..
    ## first xs is dist from edge of object, last is dist from front to end of obj. ##
    #obj.setMovingLoads(Ls,xs)
    #obj.genMLsim()
    
    
    mm=52+15
    """
    try:
        test = obj.setforce(-66.7,1060+190-52+mm)
        test = obj.setforce(-66.7,1060+190-176-52+mm)
        test = obj.setforce(-66.7, 1060+190-176-164-52+mm)
        test = obj.setforce(-66.7,1060+190-176*2-164-52+mm)
        test = obj.setforce(-66.7,1060+190-176*2-164*2-52+mm)
        test = obj.setforce(-66.7,1060+190-176*3-164*2-52+mm)
    except:
        print('oof')
    #"""
    
    #obj.calcPinsF()
    #obj.printSFD()
    #obj.printBMD()
  
    a = 110 # rip limited matboard...
    x = 20
    b = 75+2*x
    c = 1.27
    d = 1.69
    dd = 313 # 8? 
    
    #"""
    out = obj.createPseudoPi(a,b,c,d,x)
    Itot = out[0]
    
    ytop = out[1]
    ytop = abs(ytop)
    ybot = out[2]
    ybot = abs(ybot)
    ybar = out[3] 
    start = 0
    end = 300
    step = 100
    #"""
    
    #test = obj.printSFD()
    #test = obj.printBMD()
    #test = obj.BaldwinSetF(200)
    
    #test = obj.CrossSectionTester(a,b,c,d,x,dd)
    tetest = obj.FOStrain(a,b,c,d,ytop,ybot,Itot,ybar,start,end,step)
    
    #print(obj.getforces())
    ##obj.setDistributedLoad() <-- cool but wasnt this sorta useless here?? nikoo.... why... meh.
    
    ##first time setting forces
    #print("F", obj.getforces())
    #print("V", obj.getSFD())
    #print("M",obj.getBMD())
    #obj.printSFD()
    #obj.printBMD()
    #obj.printFI()
    return 1


main()