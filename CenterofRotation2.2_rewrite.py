import scipy as sp
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import os
from decimal import Decimal as D

Gconsts = {'DEBUG' : False,
           'PixConv' : 1/1280, # To convert pixels to metres
           'FontSize' : 12, # Size of font on the graph generated
           'SubjectNum' : 15,
           'LoadNumber' : 2, # 1 for 8kg added, 2 for 16kg added
           'Save' : False, # Save the graph?
           'cut' : True, # Adjusts the limits of the graph to show the full section of spine
           'Labels' : ['L5','L4','L3','L2','L1','T12','T11','T10'] # Vertebral names
           }


if Gconsts['DEBUG'] == True:
    print("Current working directory: ", format(os.getcwd()))
    print("----------")

os.chdir('I:\MPhys Project\Data')
FileNs = np.genfromtxt('List_of_filenames.txt', dtype=str)

if Gconsts['DEBUG'] == True:
    print("File Names")
    print(FileNs)
    print("----------")
    
os.chdir('I:\MPhys Project\Data\LoadedSpineData')

if Gconsts['DEBUG'] == True:
    print("Current working directory: ", format(os.getcwd()))
    print("----------")

def ReadIn(Data, Gconsts):
    Horr=[]
    Vert=[]
    if Gconsts['DEBUG'] == True:
        print(Data)
        print("----------")
        
    # Convertes the values from pixels to metres and inverts the y axis. Vert is the y componant of the points Horr is the x componant. 
    # These are plotted with Vert on the x axis and Horr on the y axis due to the monge parameterisation.
    for i in range(len(Data)):
        Vert.append(Data[i][0]*Gconsts['PixConv'])
        Horr.append(Data[i][1]*Gconsts['PixConv'])
        
    for i in range (len(Horr)):
        Horr[i]=-1*Horr[i]
        
    vertCentreBot=str(Vert[14]) # 14th point in each set is the centre of the bottom endplate of the L5 disc.  
    horrCentreBot=str(Horr[14]) # This is moved to (0,0) for every dataset.
    
    if Gconsts['DEBUG'] == True:
        print(horrCentreBot,vertCentreBot)
        print(Horr,Vert)
        print("----------")
        
    for i in range (len(Vert)):
        # Makes the centre of the sacrum endplate (Point 14) the zero point and deals with the floating point error.
        Vert[i]=float(D(str(Vert[i]))-D(vertCentreBot))
        Horr[i]=float(D(str(Horr[i]))-D(horrCentreBot))
        
    if Gconsts['DEBUG'] == True:
        print(Horr)
        print(Vert)
        print("----------")
    return Horr,Vert

def Midpoint(A,A_):
    # A and A_ are the two points the midpoint is found between.
    Midx = (A[0]+A_[0])/2
    Midy = (A[1]+A_[1])/2
    Midway = [Midx,Midy]
    return Midway

def Bisector(A,A_,x):
    # Finds the perpendicular bisector of two corrosponding points on the loaded and unloaded top endplates (A and A_).
    #Returns the cofishents of a straight line (The perp Bisector).
    dydx = VctrDiff(A, A_)
    dx = dydx[0]
    dy = dydx[1]
    Fit = dy/dx
    if Gconsts['DEBUG'] == True:
        print('Bisector Run')
        print("----------")
    
    # M is the gradient of the perpendicular bisector (y = Mx + c)
    M = -1/Fit
    Mp =  Midpoint(A, A_)
    y = M*(x-Mp[0])+Mp[1]
    
    c = M*(-1*Mp[0])+Mp[1]
    if Gconsts['DEBUG'] == True:
        print([M,c])
        plt.plot(x,y,linewidth = 0.1)
        print("----------")
    return [M, c]

def Graph(Data,Data2,Gconsts):
    # This produces the plot of the unloaded and loaded datasets
        Horr,Vert = ReadIn(Data, Gconsts)
        Horr2,Vert2 = ReadIn(Data2, Gconsts)
        plt.scatter(Horr,Vert,marker='x',color='Red',s=15,linewidths=0.5,label='Raw Points Unloaded')
        plt.scatter(Horr2,Vert2,marker='x',color='Blue',s=15,linewidths=0.5,label='Raw Points Loaded')
        if Gconsts['cut'] == True:
            plt.xlim([-50*Gconsts['PixConv'],375*Gconsts['PixConv']])
            plt.ylim([-80*Gconsts['PixConv'],150*Gconsts['PixConv']])
        
        plt.title(Dataset[:5]+' '+'Participant '+str(Gconsts['SubjectNum']), fontsize = Gconsts['FontSize']+2)
        plt.xlabel('Vertical Position (m)', fontsize = Gconsts['FontSize'])
        plt.ylabel('Horrizontal Position (m)', fontsize = Gconsts['FontSize'])
        plt.grid()
        return

'''
Frame of refference works  by converting into the staionary endplate FoR (FoRin).
Then feeds into Bisector and returns the perpendicular bisector of of the endplate past into it.
'''
def FoRin (enplttop,enpltbot,Horr,Vert,Horr2,Vert2):
    # Changing FoR to that of the lower endplate.
    if Gconsts['DEBUG'] == True:
        print('FoRin run')
        print("----------")

    #The middle point of Sacrum endplate in the loaded and unloaded data (1 and 2).
    CentreBottom1 = [Horr[enpltbot[4]],Vert[enpltbot[4]]]
    CentreBottom2 = [Horr2[enpltbot[4]],Vert2[enpltbot[4]]]
    
    '''
    TRANSLATIONAL TRANSFORM
    '''
    # This loop is where the translational transform is preformed by 
    # minusing the centrebottom point from each coordinate on the endplate.
    endplt1TTFbot = []
    endplt2TTFbot = []
    endplt1TTFtop = []
    endplt2TTFtop = []
    
    endplt1botTTFX = []
    endplt1botTTFY = []
    endplt2botTTFX = []
    endplt2botTTFY = []
    
    endplt1topTTFX = []
    endplt1topTTFY = []
    endplt2topTTFX = []
    endplt2topTTFY = []  
    
    
    for i in range(len(enpltbot)):
        
        HorrTTFbot = Horr[enpltbot[i]] - CentreBottom1[0]
        VertTTFbot = Vert[enpltbot[i]] - CentreBottom1[1]
        Horr2TTFbot = Horr2[enpltbot[i]] - CentreBottom2[0]
        Vert2TTFbot = Vert2[enpltbot[i]] - CentreBottom2[1]
        
        HorrTTFtop = Horr[enplttop[i]] - CentreBottom1[0]
        VertTTFtop = Vert[enplttop[i]] - CentreBottom1[1]
        Horr2TTFtop = Horr2[enplttop[i]] - CentreBottom2[0]
        Vert2TTFtop = Vert2[enplttop[i]] - CentreBottom2[1] 
        
        
        endplt1TTFbot.append([HorrTTFbot,VertTTFbot])
        endplt2TTFbot.append([Horr2TTFbot,Vert2TTFbot])
        
        endplt1TTFtop.append([HorrTTFtop,VertTTFtop])
        endplt2TTFtop.append([Horr2TTFtop,Vert2TTFtop])       
        
        
        endplt1botTTFX.append(HorrTTFbot)
        endplt1botTTFY.append(VertTTFbot)
        endplt2botTTFX.append(Horr2TTFbot)
        endplt2botTTFY.append(Vert2TTFbot)
        
        endplt1topTTFX.append(HorrTTFtop)
        endplt1topTTFY.append(VertTTFtop)
        endplt2topTTFX.append(Horr2TTFtop)
        endplt2topTTFY.append(Vert2TTFtop)
    
    if Gconsts['DEBUG'] == True:
        print ('VVVVVVVVVVVV')
        print ('Bottom')
        print (endplt1TTFbot)
        print ('-------------')
        print (endplt2TTFbot)
        print ('Top')
        print (endplt1TTFtop)
        print ('-------------')
        print (endplt2TTFtop)
        print("----------")
    
    '''
    ROTATIONAL TRANSFORM
    '''
    #Preformed on the already translated lists endplt1TTF (Unloaded) and endplt2TTF (Loaded)
    #(TTF = Translationally Transformed)
    
    #Finding the difference in angle of the lower endplates before applying it to the upper.
    Plate1 = np.polyfit(endplt1botTTFX,endplt1botTTFY, 1)
    Plate2 = np.polyfit(endplt2botTTFX,endplt2botTTFY, 1)
    
    Theta1 = np.arctan(Plate1[0])
    Theta2 = np.arctan(Plate2[0])
    
    # dTheta is the amount the loaded data must be turned anti-clockwise by in order to be in the same FOR as the unloaded.
    dTheta = Theta2 - Theta1 # NOTE: Theta is in RADIANS!!
    
    # This is the rotation transform matrix.
    RTMatrix = np.array([[np.cos(dTheta),-1*np.sin(dTheta)],[np.sin(dTheta),np.cos(dTheta)]])
    # The inverse for Transforming back. (Not needed as stay in the unloaded Frame of reference)
    InvRTMatrix = np.linalg.inv(RTMatrix)
    
    # Total transform of loaded data set into the same reference frame as the first.
    endplt1TFbot = np.array(endplt1TTFbot)
    endplt1TFtop = np.array(endplt1TTFtop)
    
    endplt2TFbot = np.matmul(np.array(endplt2TTFbot),RTMatrix)
    endplt2TFtop = np.matmul(np.array(endplt2TTFtop),RTMatrix)
    
    return endplt1TFbot, endplt1TFtop, endplt2TFbot, endplt2TFtop, CentreBottom1, CentreBottom2

def VctrDiff(nu,mu):
    # Finds the difference between 2 vectors. Returns vector.
    DeltaX = nu[0] - mu[0]
    DeltaY = nu[1] - mu[1]
    V = [DeltaX,DeltaY]
    return V

def VctrAdd(nu,mu):
    # Adds two vectors. Returns Vector
    DeltaX = nu[0] + mu[0]
    DeltaY = nu[1] + mu[1]
    V = [DeltaX,DeltaY]
    return V

def Interceptor(line1,line2):
    # Finds the point of intersectoin of 2 straight lines given their cofishents.
    if Gconsts['DEBUG'] == True:
        print('Interceptor run')
        
    X = (line2[1]-line1[1])/(line1[0]-line2[0])
    Y = line1[0]*X + line1[1]
    return [X,Y]

def Main(Dataset,Dataset2,Gconsts):
    Labels = Gconsts['Labels']
    Data = np.genfromtxt(Dataset,skip_header=3,skip_footer=1)
    Data2 = np.genfromtxt(Dataset2,skip_header=3,skip_footer=1)
    Graph(Data, Data2, Gconsts)
    Horr,Vert = ReadIn(Data, Gconsts)
    Horr2,Vert2 = ReadIn(Data2, Gconsts)
    # l cycles through the vertebra. 0 = bottom
    for l in range (8):
        print(Labels[l])
        #l is the vertebra number. ZeroPoint is the zero point on each vertebra.
        ZeroPoint = 28*l
        #top and bottom are the numbers of the points along each endplate of THE DISC.
        bottom = [10,11,12,13,14,15,16,17,18]
        top = [32,31,30,29,28,55,54,53,52]
        #j loops over the points on the endplates
        for j in range(len(top)):
            #print('j ',j)
            temptop = []
            tempbot = []
            for i in range (len(top)):
                # These are the points on the endplates being adjusted for the apropreate vertebra.
                # temptop and tempbot contain the positions of points on the endplate required in the lists Horr, Vert and so Horr2, Vert2 as well.
                temptop.append (top[i] + ZeroPoint)
                tempbot.append (bottom[i] + ZeroPoint)
                
        if Gconsts['DEBUG'] == True:
            print(tempbot)
            print(temptop)
            
        endpltTFbot_UL, endpltTFtop_UL, endpltTFbot_L, endpltTFtop_L, CentreBottom1, CentreBottom2 = FoRin(temptop,tempbot,Horr,Vert,Horr2,Vert2)
    
        # Finding perpendicular bisectors of the loaded points.
        x = np.arange(-50, 375, 10)
        # x is a sample rate. Ish. Gives something to plot against.
        Cofish = []
        Cors = [] # CoR = Centre of rotation
        Corsx = []
        Corsy = []

        for i in range (len(endpltTFtop_UL)):
             #print(endpltTFtop_UL[i],endpltTFtop_L[i])
             PBCofish = Bisector(endpltTFtop_UL[i],endpltTFtop_L[i],x)
             Cofish.append(PBCofish)
             
        for i in range(len(Cofish)):
            for k in range(i + 1, len(Cofish)):
                IPs = Interceptor(Cofish[i], Cofish[k])
                IPs = np.round(IPs)
                #removing outliers
                if int(IPs[0]) in range(-50,375,1) and int(IPs[1]) in range(-80,150,1):
                    Cors.append(IPs)
                    Corsx.append(IPs[0])
                    Corsy.append(IPs[1])
                    v = IPs[0] + CentreBottom1[0]
                    h = IPs[1] + CentreBottom1[1]
                    if Gconsts['DEBUG'] == True:
                        plt.scatter(v,h,color = 'Green', marker = 'x')
                    
                # Cors should be the Centres of Rotation when they are all moved to 0 
                # in the FoR of the unloaded data. There should be 36 of them.
        
        CoR = [np.average(Corsx),np.average(Corsy)]
        CoR = VctrAdd(CoR, CentreBottom1) # Moving them back in to the unloaded FoR of the sacrum centred at (0,0)
        CoR = np.round(CoR,4)
        print(CoR)
        cut = Gconsts['cut']
        Save = Gconsts['Save']
        if cut == True:
            if CoR[0] < -50*Gconsts['PixConv'] or CoR[0] > 375*Gconsts['PixConv'] or CoR[1] < -80*Gconsts['PixConv'] or CoR[1] > 150*Gconsts['PixConv']:
                if Save == True:
                    name = Dataset[:4]+'_'+'Participant_'+str(Gconsts['SubjectNum'])+'_'+str(Gconsts['LoadNumber']*8)+'kg_CUT.png'
                    os.chdir('I:\MPhys Project\Research\Clean Dis References\Results\CoR_Graphs')
                    plt.savefig(name, dpi = 200)
                    os.chdir('I:\MPhys Project\Data\LoadedSpineData')
                    print('cut')
            else:
                plt.scatter(CoR[0],CoR[1],color = 'Green', label = 'CoR in unloaded FoR')
                txt = Labels[l]
                plt.annotate(txt, (CoR[0],CoR[1]), fontstyle = 'oblique', fontsize = Gconsts['FontSize']-2)
                if l == 0:
                    plt.legend(fontsize = Gconsts['FontSize'])
        
                if Save == True:
                    name = Dataset[:4]+'_'+'Participant_'+str(Gconsts['SubjectNum'])+'_'+str(Gconsts['LoadNumber']*8)+'kg_CUT.png'
                    os.chdir('I:\MPhys Project\Research\Clean Dis References\Results\CoR_Graphs')
                    plt.savefig(name, dpi = 200)
                    os.chdir('I:\MPhys Project\Data\LoadedSpineData')
        
        else:
            plt.scatter(CoR[0],CoR[1],color = 'Green', label = 'CoR in unloaded FoR')
    
            txt = Labels[l]
            plt.annotate(txt, (CoR[0],CoR[1]), fontstyle = 'oblique', fontsize = Gconsts['FontSize']-2)
            if l == 0:
                plt.legend(fontsize = Gconsts['FontSize'])
        
            if Save == True:
                name = Dataset[:4]+'_'+'Participant_'+str(Gconsts['SubjectNum'])+'_'+str(Gconsts['LoadNumber']*8)+'kg.png'
                os.chdir('I:\MPhys Project\Research\Clean Dis References\Results\CoR_Graphs\WideAngle')
                plt.savefig(name, dpi = 200)
                os.chdir('I:\MPhys Project\Data\LoadedSpineData')
        
    return 'Done'


Participant = 3*Gconsts['SubjectNum']
print(FileNs[Participant])
Dataset = FileNs[Participant]
Dataset2 = FileNs[Participant + Gconsts['LoadNumber']]


Runner = Main(Dataset,Dataset2,Gconsts)
print(Runner)
