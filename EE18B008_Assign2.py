"""
        EE2703 Applied Programming Lab - 2019
            Assignment 1: Sample solution
"""
from sys import argv, exit
import numpy as np
import cmath
import math

"""
It's recommended to use constant variables than hard-coding them everywhere.
For example, if you decide to change the command from '.circuit' to '.start' later,
    you only need to change the constant
"""
START = '.circuit'
END = '.end'
"""
It's a good practice to check if the user has given required and only the required inputs
Otherwise, show them the expected usage.
"""
if len(argv) != 2:
    print('\nUsage: %s <inputfile>' % argv[0])          # print('\nUsage: %s <inputfile>' % argv[0])
    exit()
if len(argv) == 2:
    if(argv[1][::-1].split('.')[0] !='tsilten'):            ## if your file type is netlist then give the correct commandline
        print("ERROR ERROR")
        print("give the correct format of netlist file eg: for file name circuit1 it is circuit1.netlist")
        exit()
"""
The user might input a wrong file name by mistake.
In this case, the open function will throw an IOError.
Make sure you have taken care of it using try-catch
"""
try:
    with open(argv[1]) as file:
        lines = file.readlines()
        start = -1; end = -2
        for line in lines:              # extracting circuit definition start and end lines
            if START == line[:len(START)]:
                start = lines.index(line)
            elif END == line[:len(END)]:
                end = lines.index(line)
                break
        if start >= end:                # validating circuit block
            print('Invalid circuit definition')
            exit(0)
except IOError:
    print('Invalid file')
    exit()
for i in range(len(lines)):
    lines[i]=lines[i].split('#')[0].split('\n')[0]
w=0
try:
    if(lines[end+1].split('#')[0].split()[0]=='.ac'):
        ac=lines[end+1].split('#')[0].split()[len(lines[end+1].split('#')[0].split())-1]
        w=float(ac)
except Exception:
    pass

class Resistor:
    name=''
    n1=''
    n2=''
    val=''
    def __init__(self,name,n1,n2,val):
        self.name=name
        if(n1=='GND'):             # Making GND nnode as ZERO
            self.n1='0'                                   
        else:
            self.n1=n1
        if(n2=='GND'):
            self.n2='0'
        else:                                            # Class for Impedance
            self.n2=n2                                   # name,node 1,node 2 and value of impedance
        if(name[0]=='R'):
            self.value=float(val)+(0j)
        elif(name[0]=='L'):
            self.value=w*float(val)*(1j)
        elif(name[0]=='C'):
            self.value=1/(w*float(val)*(1j))
class Voltage:
    name=''
    n1=''
    n2=''
    val=''
    phase=''
    def __init__(self,name,n1,n2,val,phase):                         # Class for Voltage Sources
        self.name=name                                               # name,node 1,node 2 and peak to peak value of voltage source and phase
        if(n1=='GND'):
            self.n1='0'
        else:
            self.n1=n1
        if(n2=='GND'):
            self.n2='0'
        else:
            self.n2=n2  
        self.value=(float(val))*(cmath.exp(float(phase)*(1j)))
class Current:
    name=''
    n1=''
    n2=''
    val=''
    phase=''
    def __init__(self,name,n1,n2,val,phase):                         # Class for Current Sources
        self.name=name                                               # name,node 1,node 2 and peak to peak value of current source and phase
        if(n1=='GND'):
            self.n1='0'
        else:
            self.n1=n1
        if(n2=='GND'):
            self.n2='0'
        else:
            self.n2=n2
        self.value=(float(val))*(cmath.exp(float(phase)*(1j)))

resistor=[]                     #list which stores objects of Resistors
voltage=[]                      #list which stores objects of voltage sources                     
current=[]                      #list which stores objects of current sources
totalnodes=[]

for line in lines[start+1:end]:
            reqline=line.split('#')[0]
            if(reqline[0]=='R' or reqline[0]=='L' or reqline[0]=='C'):
                name,n1,n2,val=reqline.split()
                resistor.append(Resistor(name,n1,n2,val))  
            elif(reqline[0]=='V' and w!=0):                                                # Check for ac voltage
                req=reqline.split('ac')[0].split()+reqline.split('ac')[1].split()
                name,n1,n2,val,phase=req
                val=str(float(val)/2)
                voltage.append(Voltage(name,n1,n2,val,phase))
            elif(reqline[0]=='V' and w==0):                                                 # Check for dc voltage
                name,n1,n2,val=reqline.split()                                             
                voltage.append(Voltage(name,n1,n2,val,'0'))
            elif(reqline[0]=='I' and w!=0):                                                # Check for ac current
                req=reqline.split('ac')[0].split()+reqline.split('ac')[1].split()            
                name,n1,n2,val,phase=req
                val=str(float(val)/2)
                current.append(Current(name,n1,n2,val,phase))
            elif(reqline[0]=='I'):                                                           # Check for dc current
                name,n1,n2,val=reqline.split()
                current.append(Current(name,n1,n2,val,'0'))

for obj in resistor+voltage+current:
    totalnodes=[obj.n1,obj.n2]+totalnodes

totalnodes=set(totalnodes)                 #which gives total list of nodes
nodes=list(totalnodes)
nodes.sort()
for p in range(len(voltage)):
    nodes.append('i'+str(p))                    #adds current variables passing through voltage sources
dictionary=dict()
for i in range(len(nodes)-len(totalnodes)):                         # To make sure which current is passing through which voltage source
    dictionary.update({voltage[i]:nodes[len(totalnodes)+i]})
a=list(dictionary.values())
b=list(dictionary.keys())
f1=[]
B=[[0] for i in range(len(nodes))]         # Matrix B
for i in range(len(nodes)):
    if(i==0):                                       # To fill Matrix A which contain current variables 
        f1.append([1])                                              
    else:
        f1.append([0])
    for j in range(len(nodes)-1):                                
        f1[i].append(0)
for i in range(1,len(nodes)):
    for obj1 in voltage:
        if(obj1.n1==nodes[i]):                                          # Filling Matrix A which contain currents through voltage sources as variables
            f1[i][nodes.index(dictionary[obj1])]=-1
        elif(obj1.n2==nodes[i]):
            f1[i][nodes.index(dictionary[obj1])]=1
    for obj2 in current:
        if(obj2.n1==nodes[i]):                                           
            B[i][0]=-(obj2.value)                                      # Filling Matrix B if there are any independent current sources
        elif(obj2.n2==nodes[i]):
            B[i][0]=obj2.value
    for obj in resistor:                                                
        if(obj.n1==nodes[i]):
            f1[i][nodes.index(obj.n1)]=(1/(obj.value))+f1[i][nodes.index(obj.n1)]
            f1[i][nodes.index(obj.n2)]=f1[i][nodes.index(obj.n2)]-(1/(obj.value))           # Filling Matrix A 
        elif(obj.n2==nodes[i]):                                                         
            f1[i][nodes.index(obj.n2)]=(1/(obj.value))+f1[i][nodes.index(obj.n2)]
            f1[i][nodes.index(obj.n1)]=f1[i][nodes.index(obj.n1)]-(1/(obj.value))
for obj in voltage:
    f1[nodes.index(dictionary[obj])][nodes.index(obj.n1)]=1
    f1[nodes.index(dictionary[obj])][nodes.index(obj.n2)]=-1                  # Filling Matrix A and B corresponding EMF equations 
    B[nodes.index(dictionary[obj])][0]=obj.value

A=np.array(f1)
B=np.array(B)
for i in range(len(totalnodes)):
    nodes[i]='v'+nodes[i]
for i in range(len(totalnodes),len(nodes)):
    nodes[i]=nodes[i]+' through '+b[a.index(nodes[i])].name

X=np.linalg.solve(A,B)                   # Solving  X which has unknown variables           
if(w==0):
    X=np.real(X)

for i in range(len(nodes)):
    print(nodes[i] ,"=", X[i][0])










