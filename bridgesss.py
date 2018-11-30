#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 18:44:23 2018

A program to help a couple of engsci kids build a strong matboard bridge in CIV102!

@author: nikoogivehchian
"""

import numpy as np
import matplotlib.pyplot as plt

class Bridges:
    def __init__(self, L):
        """
        inputted L is the length of the bridge, in mm. Should be integer- no decimal values!
        intervals of 1 mm will be considered along the bridge (ie. no decimals in between, 
        as the program currently uses python's indexing to count intervals along the bridge and i too lazy
        to go beyond these anyways small intervals.)
        
        """
        if L != int(L) or L < 0:
            print("err: L must be positive integer!")
        self.len = L
        self.x = [i for i in range(self.len+1)]
        self.f = [0 for i in range(self.len+1)] ##ALL FORCES (EXT + RXN)
        self.fext = [0 for i in range(self.len+1)] ##NON-RXN aka EXT FORCES
        self.pins = 0 #PIN FORCES - saved from last use.
        ### PLZ NOTE: NO OTHER FORCES ARE SAVED IN CLASS. that is, anything for the sim fcns vanishes after use.
        self.I = [0 for i in range(self.len+1)] #will be used? ...  ...  ...
        self.o = [0 for i in range(self.len+1)] #lol wtf is this
        ##add comments to the user about the different avaiable fcns, recomended procedure, how to... etc.
        ##print("PLEASE NOTE: all units are in terms of mm and N (and MPa when it goes through some calcs).")
        ##print("have fun not having to calculate these yourself!")

    """def setI(self,ii,rng):
        print("please note I has units mm^4")
        #split = input("is I constant? Y/N: ")
        if rng == "all" or rng == "ALL":
            self.I = [int(ii) for i in range(self.len+1)]
        else: #then rng is startpos, endpos list.
            first=rng[0]
            end=rng[1]
            for i in range(first,end+1,1):
                self.I[i] = ii
            print("i too lazy to check but make sure u initialize all ur 'I' values k?")
        #print(self.I)
        return 1"""
    
    """to start... this stuffs, which user can always call to get some infos."""
    def getx(self):
        return self.x
    def getforces(self):
        return self.f   
    def getPOIs(self, n):
        """takes in a list of generated values, outputs points of interest alongside these nums"""
        outVals = []
        outInds = []
        for i in range(len(n)):
            if (n[i] - n[i-1] != n[i-1] - n[i-2]) and n[i] != n[i-1]:
                outVals += [n[i]]
                outInds += [i]
        outVals += [n[len(n)-1]]
        outInds += [len(n)]
        for i in range(len(outInds)):
            print("ind: ",outInds[i], " val: ",outVals[i])
        return 1
    
    """manually input loads and applied forces, to be stored to the object"""
    def setforce(self, force, pos): #in newtons!!! NOT kN!
        if pos < 0:
            print("error: pos < 0")
            return 0
        if pos != int(pos):
            print("err: position invalid")
            return 0
        try:
            self.f[pos] = force
            self.fext[pos] = force
            #self.poi += [pos]
            return 1
        except:
            print("error: invalid position on the bridge")
            return 0
    def setDistributedLoad(self):
        print("(again, no error checking, so dont break me plz)")
        w = input("so.. what's the load? ")
        w=int(w)
        wtot=self.len * w
        halfw = wtot/2
        ws = [0 for i in range(0,self.len+1)]
        ws[0] = halfw
        #print("ws ",ws)
        for i in range(1,self.len+1):
            ws[i] = float(ws[i-1]) - float(w)
        #print("ws ",ws)
        for i in range(0,self.len+1):
            self.f[i] = ws[i] + self.f[i]
            self.fext[i] = ws[i] + self.fext[i]
        return 1
        
    """def setLoads(self):
        print("(again, no error checking, so dont break me plz)")
        n = input("time to add forces to our bridge! type 'w' for self-weight, or 'l' for other loads.")
        if n == "w":
            w = input("distributed load, in N/mm : ")
            whalf = w/2
            for i in range(0,self.len):
                self.f[i] = [self.f[i] + self.x[]]""" ##fin me later!

    """Dealing with the pins"""
    def setPins(self,p1,p2):
        """sets pin locations"""
        self.pins=[p1,p2]
        return 1
        
    def calcPinsF(self):
        """calculates reaction forces in the pins (which must have been set earlier)"""
        print("---------------- CalcPinsF ----------------")
        print("current pin locations: ", self.pins)
        #print("u better have inputted all the forces dammit.")
        #print("otherwise i totally wrong and u screwed.")
        #print("ps, once again, i very ez to break. be nice plz!!")
        ## NOTE: clockwise and up are positive.
        sumFy = 0
        #print(self.fext)
        for i in self.fext:
            sumFy += i
        ##we calculate moments about first pin.
        ##we also assume the system is very simple with two pins, otherwise i break.       
        sumMs = 0
        cntr=-1
        for i in self.fext:
            cntr+=1
            if i != 0:
                dist = cntr - self.pins[0]
                if dist > 0: #to the right
                    sumMs -= i*abs(dist) ##since dn force is pos moment
                elif dist < 0:
                    sumMs += i*abs(dist)
        ##yay for we have sum of forces and forces! since i know this is a two-pin system, let's make life ez...
        dist2 = self.pins[1] - self.pins[0] ##asume u aint stupid and gave me pins in order...
        F2 = sumMs/dist2
        F1 = -sumFy - F2
        print("Force at pin 1: ",F1)
        print("Force at pin 2: ",F2)
        self.f[self.pins[0]]= F1
        self.f[self.pins[1]]= F2 
        self.pinsf = [F1, F2]
        #print("finished saving pin forces")
        #print(self.f)
        return 1

    def calcPinsFint(self, L):
        print("---------------- CalcPinsFint ----------------")
        """
        input L is list of non-rxn forces
        """
        #print("u better have inputted all the forces dammit.")
        #print("otherwise i totally wrong and u screwed.")
        #print("ps, once again, i very ez to break. be nice plz!!")
        ##clockwise and up are positive.
        sumFy = 0
        pins = list(self.pins)
        #print(self.fext)
        for i in L:
            sumFy += i
        print("sumFY :", sumFy)
        ##we do moments about first pin.
        ##we also assume the system is very simple with two pins, otherwise i break.       
        sumMs = 0
        cntr=-1
        for i in L:
            #print(i)
            cntr+=1
            if i != 0:
                dist = cntr - pins[0]
                #print("dist ", dist, "force ",i)
                if dist > 0: #to the right, should always be this!
                    sumMs -= i*dist ##since dn force is pos moment, and i is negative.
                elif dist < 0:
                    print("WHY I NEGATIVEEE???????")
                    sumMs += (i)*abs(dist)
        print("sumMS: ",sumMs)
        ##yay for we have sum of forces and forces! since i know this is a two-pin system, let's make life ez...
        dist2 = pins[1] - pins[0] ##asume u aint stupid and gave me pins in order...
        F2 = sumMs/abs(dist2)
        F1 = -sumFy - F2
        print("Force at pin 1: ",F1)
        print("Force at pin 2: ",F2)
        pinsF = [F1,F2]
        #L[pins[0]]= F1
        #L[pins[1]]= F2 
        #self.pinsf = [F1, F2]
        #print("finished saving pin forces")
        #print(self.f)
        return [pins,pinsF]

    def simTrain(self, Ls, xs):
        print("----------- Train Simulator -----------")
        tempf = [0 for i in range(self.len+1)]
        loads = list(Ls)
        print(xs)
        ctr=-1
        Mpos=[0,0,0] ##  index assignments: 0 is location, 1 is Mval, 3 is tracker for where i was (where train located)
        Mneg=[0,0,0]
        for i in range(0,self.len,500): ### each sim!! ###
            print("ith turn: ",i)
            ctr = -1
            print(ctr)
            tempf = list(self.f)
            while ctr < len(xs):
                ctr+=1
                print(ctr)
                for n in range(0,self.len-1):
                    if ctr == len(xs):
                        print("broken")
                        break
                    if (n+i) == xs[ctr]:
                        print("made it to point of trying...")
                        print(n+i)
                        print("ctr :", ctr)
                        try:
                            tempf[n+i] = loads[ctr]
                            ctr+=1
                            print("tried and could")
                        #tempf[n+i] = loads[n]
                        except:
                            print("except")
                            ctr += 1
            print("tempf ",tempf)
            #solve rxn forces
            pstuffs = self.calcPinsFint(tempf)
            pinpos = pstuffs[0]
            pinf = pstuffs[1]
            #print("pinpos ", pinpos)
            #print("pinf ",pinf)
            tempf[pinpos[0]] = pinf[0]
            tempf[pinpos[1]] = pinf[1]
            #print("tempf :", tempf)
            #now we have list of complete forces! proceed to generate SFD, BMD, and POIs...
            V = self.getSFDint(tempf)
            #print("V ",V)
            self.printSFDint(V)
            M = self.getBMDint(V)
            self.printBMDint(M)
            ## let's quickly get max positive and negative moments...
            for mmmm in range(len(M)):
                if M[mmmm] > Mpos[1]:
                    Mpos[0] = mmmm
                    Mpos[1] = M[mmmm]
                    Mpos[2] = i #where art thou train...
                if M[mmmm] < Mneg[1]:
                    Mneg[0] = mmmm
                    Mneg[1] = M[mmmm]
                    Mneg[2] = i
            print("okie next simulation!")## end each sim! ##
            print("")
        out = [Mpos,Mneg]
        print("M max pos: ",Mpos)
        print("M max neg: ",Mneg)
        return 1
    
    def simTrainint(self, Ls, xs,start,end,step):
        print("----------- Train Simulator Int -----------")
        """
        Ls is list of point loads the train (or any moving obj) has been approximated to (ie. -66.7 N)
        xs is the spacing between each of these point loads. 
        """
        tempf = [0 for i in range(self.len+1)]
        loads = list(Ls)
        #print(loads)
        ctr=-1
        Mpos=[0,0,0] ## INDEX 0 is location, 1 IS Mval, 3 is tracker for what i was (where train located)
        Mneg=[0,0,0]
        Vmax = [0,0,0]
        for i in range(start,end+1,step): ### each sim!! ###
            print("+++++ the train has shifted ",i, "mm from the first pin. +++++")
            ctr = -1
            tempf = [0 for i in range(self.len+1)]
            # NIKOO U DUMBASS FIX THIS !!
            while ctr < len(loads)-1: ## keep from going over
                ctr+=1
                #print(ctr)
                print('load ',loads[ctr])
                for n in range(0,self.len-1): ## iterate over positions, stop when position + shift is in xs list
                    if ctr == len(loads):
                        print("broken")
                        break
                    if (n+i) == xs[ctr]:
                        print("made it to point of trying...")
                        #print(n+i)
                        #print("ctr :", ctr)
                        try:
                            tempf[n+i] = loads[ctr]
                            ctr+=1
                            print("tried and could")
                        except:
                            print("except")
                            ctr += 1
            #print("tempf ",tempf)
            #solve rxn forces
            pstuffs = self.calcPinsFint(tempf)
            pinpos = pstuffs[0]
            pinf = pstuffs[1]
            print("pinpos ", pinpos)
            print("pinf ",pinf)
            tempf[pinpos[0]] = pinf[0]
            tempf[pinpos[1]] = pinf[1]
            #print("tempf :", tempf)
            #now we have list of complete forces! proceed to generate SFD and BMD for this train position...
            V = self.getSFDint(tempf)
            self.printSFDint(V)
            M = self.getBMDint(V)
            self.printBMDint(M)
            ## let's also get max positive and negative moments... and V... 
            for vvv in range(len(V)):
                if V[vvv] > Vmax[1]:
                    Vmax[0] = vvv
                    Vmax[1] = V[vvv]
                    Vmax[2] = i
            for mmmm in range(len(M)):
                if M[mmmm] > Mpos[1]:
                    Mpos[0] = mmmm
                    Mpos[1] = M[mmmm]
                    Mpos[2] = i #where art thou train...
                if M[mmmm] < Mneg[1]:
                    Mneg[0] = mmmm
                    Mneg[1] = M[mmmm]
                    Mneg[2] = i
            print("okie next simulation!")## end each sim! ##
            print("")
        out = [Mpos,Mneg,Vmax]
        return out
    
    def FOStrain(self,a,b,c,d,ytop,ybot,Itot,ybar,start,end,step):
        Ls = [(-400/6) for i in range(0,6)]
        xs = [52, 52+176, 52+176+164, 52+176*2+164, 52+176*2+164*2, 52+176*3+164*2, 52+176*3+164*2+52]
        MMM = self.simTrainint(Ls,xs,start,end,step)
        Mpos = MMM[0]
        Mneg = MMM[1]
        Vs = MMM[2]
        print(" == Max M +ve Infos == ")
        print(" location: ",Mpos[0])
        print(" M: ",Mpos[1])
        print(" how much train shifted: ",Mpos[2]) ## this is relative to the leftmost pin
        
        MPbot = Mpos[1]*ybot/Itot #moment positive bottom -- tension
        MPtop = Mpos[1]*ytop/Itot # comp
        
        print("")
        print(" == Max M -ve Infos == ")
        print(" location: ",Mneg[0])
        print(" M: ",Mneg[1])
        print(" how much train shifted: ",Mneg[2]) ## this is relative to the leftmost pin
        
        MNbot = Mneg[1]*ybot/Itot # comp
        MNtop = Mneg[1]*ytop/Itot # tens
        
        Mtens = max([MPbot,MNtop])
        Mcomp = max([MPtop,MNbot])
        
        FOSt = 30/Mtens
        FOSc = 6/Mcomp
        
        print("")
        print(" == FOS stuffs == ")
        print(" FOS tens: ", FOSt)
        print(" FOS comp: ", FOSc)
        
        """ THIS USES TRAIN LOADING, BUT THE SHEAR PART THINGY """
        A = ybar*d
        y = 0.5*ybar
        Q = A*y
        
        print("")
        print(" == Max shear infos == \n","location: ", Vs[0], "\n magnitude: ",Vs[1]," \n shift: ",Vs[2])
        taoo = (Vs[1])*Q/(2*d*Itot)
        FOStaoo = 4/taoo
        
        print("")
        print(" calculated shear: ", taoo, "\n ... vs. max shear: 4")
        print(" FOS shear: ", FOStaoo)
        
        return 1

    def BaldwinGetF(self, F): ##F IS EACH POINT LOAD SO TEST 400 AND UPPPP
        holes = [550,(550+510+190)]
        tempf = [0 for i in range(self.len+1)]
        tempf[holes[0]] = -F
        tempf[holes[1]] = -F       
        ## set reaction forces
        pstuffs = self.calcPinsFint(tempf)
        pinpos = pstuffs[0]
        pinf = pstuffs[1]
        tempf[pinpos[0]] = pinf[0]
        tempf[pinpos[1]] = pinf[1]        
        return tempf
        
    def BaldwinSetF(self, F): #F EACH POINT force applied by the baldwin
        holes = [550,(550+510+190)]
        tempf = [0 for i in range(self.len+1)]
        tempf[holes[0]] = -F
        tempf[holes[1]] = -F
        #print(tempf)
        
        ## set reaction forces
        pstuffs = self.calcPinsFint(tempf)
        pinpos = pstuffs[0]
        pinf = pstuffs[1]
        tempf[pinpos[0]] = pinf[0]
        tempf[pinpos[1]] = pinf[1]
        #print(tempf)
        
        ## make the diagrams!
        V = self.getSFDint(tempf)
        self.printSFDint(V)
        M = self.getBMDint(V)
        self.printBMDint(M)
        #self.calcStrength() ##ADD THIS WHEN I HAVE Y VALS AND SUCH
        return 1        

    def BaldwinResults(self,F): #F EACH POINT force applied by the baldwin
        print("------------ baldwin results ------------")
        holes = [550,(550+510+190)]
        tempf = [0 for i in range(self.len+1)]
        tempf[holes[0]] = -F
        tempf[holes[1]] = -F
        #print(tempf)
        
        ## set reaction forces
        pstuffs = self.calcPinsFint(tempf)
        pinpos = pstuffs[0]
        pinf = pstuffs[1]
        tempf[pinpos[0]] = pinf[0]
        tempf[pinpos[1]] = pinf[1]
        #print(tempf)
        
        ## make the diagrams!
        V = self.getSFDint(tempf)
        self.printSFDint(V)
        M = self.getBMDint(V)
        self.printBMDint(M)
        #self.calcStrength() ##ADD THIS WHEN I HAVE Y VALS AND SUCH
        return 1        

    
    """SFD stuffs"""
    """INT VERSION"""
    def getSFDint(self, F): 
        print("------------ SFD int ------------")
        V = [0 for i in range(self.len+1)]
        V[0] = F[0]
        #self.calcPinsF() ##lol just do it here to be sure XD
        for i in range(1,self.len+1):
            V[i]=V[i-1] + F[i]
        #print("points of interest (from getSFD): ", self.poi)
        #print(V)
        return V
    def printSFDint(self,V):
        #the following is for INTERNAL USE, where I choose V myself.
        plt.plot(self.x,V,"-b",self.x,self.o,"--k")
        plt.grid(b=True, which='both', axis='both', linestyle=':')
        plt.xlabel("x [mm]")
        plt.ylabel("V [N]")
        plt.show()
        print("SFD POIs: ") ##not pretty rn but whatever
        self.getPOIs(V)
        return 1

    """USER VERSION"""
    def getSFD(self): 
        V = [0 for i in range(self.len+1)]
        self.calcPinsF() ## here ONLY for user version.
        V[0] = 0
        for i in range(1,self.len+1):
            V[i]=V[i-1]+self.f[i-1] ##issue in diragram: needs be shifted 1 over... rip
        #print("points of interest (from getSFD): ", self.poi)
        #print(V)
        return V
    def printSFD(self):
        #the following is bad code. please revise later!
        #print("fs",self.f)
        V = self.getSFD()
        #print("V",V)
        plt.plot(self.x,V,"-b",self.x,self.o,"--k")
        plt.grid(b=True, which='both', axis='both', linestyle=':')
        plt.xlabel("x [mm]")
        plt.ylabel("V [N]")
        plt.show()
        print("SFD POIs: ") ##not pretty rn but whatever
        self.getPOIs(V)
        return 1
    
    """BMD stuffs"""
    """INT VERSION"""
    def getBMDint(self,V):
        print("------------ BMD int ------------")
        #print("VVV",V)
        M = [0 for i in range(self.len+1)]
        for i in range(1,self.len+1):
            M[i] = M[i-1] + (V[i-1])*1 #incrementing distance by 1
        return M
    def printBMDint(self,M):
        #the following is for INTERNAL USE ONLY.
        plt.plot(self.x,M,"-b",self.x,self.o,"--k")
        plt.grid(b=True, which='both', axis='both', linestyle=':')
        plt.xlabel("x [mm]")
        plt.ylabel("M [Nmm]")
        plt.show()
        print("BMD POIs: ")
        self.getPOIs(M)
        return 1
    
    """USER VERSION"""
    def getBMD(self):
        V = self.getSFD()
        M = [0 for i in range(self.len+1)]
        for i in range(1,self.len+1):
            M[i] = M[i-1]+ (V[i-1])*1 #incrementing distance by 1
        return M
    def getMmax(self, M):
        return max(M)
    def getMmin(self,M):
        return min(M)
    def printBMD(self):
        #the following is bad code. please revise later!
        #print("fs",self.f)
        M = self.getBMD()
        #print("M",M)
        plt.plot(self.x,M,"-b",self.x,self.o,"--k")
        plt.grid(b=True, which='both', axis='both', linestyle=':')
        plt.xlabel("x [mm]")
        plt.ylabel("M [Nmm]")
        plt.show()
        print("BMD POIs: ")
        self.getPOIs(M)
        return 1
      
    """curvature diagram stuffs"""
    def getFI(self, I): 
        """
        I is a list of length self.len, with I values at each position given.
        """
        E = 4000 #[MPa]
        M = self.getBMD()
        FI = list(M)
        for i in range(0,self.len+1):
            FI[i] = FI[i]/(E*self.I[i-1])
        return FI
    def printFI(self):
        FI = self.getFI()
        plt.plot(self.x,FI,"-b",self.x,self.o,"--k")
        plt.grid(b=True, which='both', axis='both', linestyle=':')
        plt.xlabel("x [mm]")
        plt.ylabel("Curvature [rad/mm]")
        plt.show()
        print("FI POIs: ",self.getPOIs(FI))
        return 1
    
    
    def calcStrength(self, M, ytop, ybot): 
        """
        GIVEN: a list of M values from some BMD, ytop, ybot, and type of loading.
        """
        if typef != "b" and typef != "t":
            print("invalid type of loading defined.")
            return -1
        
        ## max        
        mmax = self.getMmax(M)
        maxAt = M.index(mmax)
        Imax = self.I[maxAt]
        
        print("for M max...")
        smaxt = abs(mmax*ytop/Imax)#shear max top
        smaxb = abs(mmax*ybot/Imax)
        print("shear top, compression: ", smaxt)
        print("shear bot, tension: ", smaxb)
        
        ## min
        mmin = self.getMmin(M)
        minAt = M.index(mmin)
        Imin = self.I[minAt]

        print("for M min...")
        smint = abs(mmin*ytop/Imin)
        sminb = abs(mmin*ybot/Imin)
        print("shear top, compression: ", smint)
        print("shear bot, tension: ", sminb)
        
        compShears = [smaxt, smint]
        tensShears = [smaxb, sminb]
        
        compMax = max(compShears)
        if compMax == smaxt: #max compression at the top always, but using max M
            Isc = Imax
        elif compMax == smint: #max compression, min M
            Isc = Imin
        tensMax = max(tensShears)
        if tensMax == smaxb: 
            Ist = Imax
        else:
            Ist = Imin
        
        FOScomp = 6/compMax
        FOStens = 30/tensMax
        
        print("FOS comp: ",FOScomp)
        print("FOS tens: ", FOStens)
               
        return 1
    
    def predictFailureLoad(self, Itot, ybot, ytop):
        print("--------- predicting failure load round 1... ---------")
        
        print("Mmax +ve is 165P...")
        PMt = (30*Itot)/(165*ybot)
        PMc = (6*Itot)/(165*ytop)
        print("Max compressive load (P): ",PMc)
        print("Max tensile load (P): ",PMt)
        
        print("Mmax -ve is 190P...")
        pmt = (30*Itot)/(190*ytop)            
        pmc = (6*Itot)/(190*ybot)
        print("Max compressive load (P): ",pmc)
        print("Max tensile load (P): ",pmt)
        
        Ps = [PMt,PMc,pmt,pmc]
        P = min(Ps)
        print("Predicted Failure Loads (P) from Geometry: \n",Ps)
        #print("max load our poor little bridge can take according to predictFailureLoad: ", P)
        print("")
        return Ps
            
    """MAKE THE GEOMETRY"""
    def createPseudoPi(self,a,b,c,d,x): ## Returns Itot, Ytop, Ybot
        print("--------- creating pseudo pi... ---------")
        td4 =  d #0.75*d ## 5 is sweet
        A1 = c*b
        A2 = a*d
        A3 = a*d
        #A4 = 2*d*(((x*x)+(a*a))**(0.5))## this is the diaphragms
        A4 = (b-2*x-2*d)*(td4)
        #A4 = (td4)*((b-2*x-2*d)
        h = a+c
        Y1 = a + 0.5*c
        Y2 = 0.5*a
        Y3 = 0.5*a
        #Y4 = 2*a/
        Y4 = a-0.5*td4 ###THIS IS THE ORIGINAL CHANGE ME BACKKKK
        #Y4 = 0.5*d
        #ybar = (A1*Y1+A2*Y2+A3*Y3)/(A1+A2+A3)
        ybar = (A1*Y1+A2*Y2+A3*Y3+A4*Y4)/(A1+A2+A3+A4)
        
        y1 = Y1 - ybar
        y2 = Y2 - ybar
        y3 = Y3 - ybar
        y4 = Y4 - ybar
        
        I1 = A1*(y1*y1) + (b*c*c*c)/12
        I2 = A2*(y2*y2) + (d*a*a*a)/12
        I3 = A3*(y3*y3) + (d*a*a*a)/12
        I4 = A4*(y4*y4) + ((b-2*x-2*d)*(td4*td4*td4))/12
        Itot = I1 + I2 + I3 + I4
        
        ytop = a+c-ybar
        ybot = ybar
        out = [Itot, ytop, ybot, ybar]
        print("Itot: ",Itot)
        print("ytop: ",ytop)
        print("ybot: ",ybot)
        print("")
        return out
    
    
    """TEST THE GEOMETRY"""
    def ShearFailureMat(self,a,b,c,d,ybar,Itot): ## d
        print("---------- Shear Failure ----------")
        A = ybar*d*2
        y = 0.5*ybar
        Q = A*y
        print("Q: ",Q)
        #tao = 2*((Vmax)*Q/(Itot*d)) ## 2 bc two legs of the pi
        #FOStao = 4/tao
        #print("tao: ",tao)
        #print("FOS for tao: ",FOStao)
        pmaxpot = 4*Itot*(2*d)/(Q) ## ps. 2*d = b in the formula but dont wanna replace b here... yep.
        print("Shear Failure load (P): ", pmaxpot)
        print("")
        return pmaxpot
        
    def ShearGlueFailure(self,a,b,c,d,ybar,Itot):
        print("---------- Shear Glue Failure ----------")
        A = c*b
        y = (a+0.5*c) - ybar
        Q = A*y
        pmaxpot = 4*Itot*(2*d)/(Q)
        print("Shear Failure Load (P): ", pmaxpot)
        print("")
        return pmaxpot

    def ThinPlateBuckling(self,a,b,c,d,x,dd,Itot,ybar):
        print("---------- Thin Plate Buckling ----------")
        
        sigmaCR1 = 4*(3.14*3.14)*4000*(d*d)/((b-2*x)*(b-2*x)*12*(1-(0.2*0.2)))
        ## this is stress on the centre of the beam (4 sides fixed)
        sigmaCR2 = 0.425*(3.14*3.14)*4000*(d*d)/((x*x)*12*(1-(0.2*0.2)))
        ## this is stress on the edges of the beam (3 sides fixed)
        sigmaCR3 = 6*(3.14*3.14)*4000*(d*d)/((ybar)*(ybar)*12*(1-(0.2*0.2)))
        ## this is stresses on the webs (triangular stress, 4 sides fixed)
        taoCR = (5*(3.14*3.14)*4000/(12*(1-(0.2*0.2))))*((d/a)*(d/a)+(d/dd)*(d/dd))
        ## this is shear buckling
        
        pCR1 = sigmaCR1*Itot/(165*(ybar*0.5))
        pCR2 = sigmaCR2*Itot/(165*(ybar*0.5))
        pCR3 = sigmaCR3*Itot/(165*(ybar*0.5))
        
        A = c*b
        y = (a+0.5*c) - ybar
        Q = A*y
        
        ptaoCR = taoCR*Itot*(2*d)/Q
        
        print("for positive moment... ")
        print("crushing stress #1: ",sigmaCR1)
        print("Corresponding P: ", pCR1)
        print("crushing stress #2: ",sigmaCR2)
        print("Corresponding P: ", pCR2)
        print("crushing stress #3: ",sigmaCR3)
        print("Corresponding P: ", pCR3)
        print("Shear buckling: ",taoCR)
        print("Corresponding P: ",ptaoCR)
        print("")
        
        y = (a+c - ybar) - c

        sigmaCR3n = 6*(3.14*3.14)*4000*(d*d)/((y)*(y)*12*(1-(0.2*0.2)))
        ## this is stresses on the webs (triangular stress, 4 sides fixed)

        pCR3n = sigmaCR3n*Itot/(165*(0.5*ybar))
        
        print("for negative moment...")
        print("crushing stress #3: ",sigmaCR3n)
        print("Corresponding P: ", pCR3n)
        print("")
        
        out = [pCR1,pCR2,pCR3,ptaoCR]
        
        return out



    
    def CrossSectionTester(self,a,b,c,d,x,dd): 
        """ 
        a,b,c,d as per written on the paper... u better have the paper XD
        x: width of flange
        d: thickness of the matboard
        dd: spacing of diaphragms
        """
        out = self.createPseudoPi(a,b,c,d,x) 
        Itot = out[0]
        ytop = out[1]
        ytop = abs(ytop)
        ybot = out[2]
        ybot = abs(ybot)
        ybar = out[3]
        
        print("note: all predicted loads are P, NOT 2P (which is the total baldwin load on the bridge).")
        print("thus, max loading is truly double these predictions. \n")
        #I = [Itot for i in range(self.len+1)] ##this would exist to vary the cross section... but i dont wanna.
        Pfailss = self.predictFailureLoad(Itot, ybot, ytop)
        Pfail8 = Pfailss[0]
        Pfail9 = Pfailss[1]
        Pfail10 = Pfailss[2]
        Pfail11 = Pfailss[3]

        #self.BaldwinResults(Pfail)
        Pfail2 = self.ShearFailureMat(a,b,c,d,ybar,Itot)
        Pfail3 = self.ShearGlueFailure(a,b,c,d,ybar,Itot)
        Pfailssss = self.ThinPlateBuckling(a,b,c,d,x,dd,Itot,ybar)
        Pfail4 = Pfailssss[0]
        Pfail5 = Pfailssss[1]
        Pfail6 = Pfailssss[2]
        Pfail7 = Pfailssss[3]
        PPPPPP = [Pfail8,Pfail9,Pfail10,Pfail11,Pfail2,Pfail3,Pfail4,Pfail5,Pfail6,Pfail7]
        #print("Potential Failures summ: ",PPPPPP)
        PFAILURE = min(PPPPPP)
        print("if this cross-section was for the whole thing... final Pfail: ", PFAILURE, "N")
        #self.BaldwinResults(Pfail)
        return 1
    
        
    def MidspanDeflection(self):
        """barely started then abandoned... finish later!"""
        BaldwinLoading = 200
        forces = self.BaldwinGetF(BaldwinLoading) ## these are all forces in the bridge!
        V = self.printSFDint(forces)
        M = self.printBMDint(V)
        
        #print(forces)
        return 1
        


