#!/usr/bin/env python

# Copyright 2016 - 2017 Bas van Meerten and Wouter Franssen

# This file is part of ssNake.
#
# ssNake is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ssNake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ssNake. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import spectrum_classes as sc
import re

def convertDIFDUB(dat):
    def checkWrite(dup,currentNum,step,numberList):
        if dup != '':
            for dupli in range(int(dup)):
                numberList.append(numberList[-1] + int(step))
            dup = ''
            step = ''
        elif currentNum != '': #Write if available
            numberList.append(int(currentNum)) #write old num
            currentNum = ''
        elif step != '':
            numberList.append(numberList[-1] + int(step))
            step = ''
            
        return dup,currentNum,step,numberList
        
        
    SQZ = {'@': 0, 'A': 1, 'B': 2,'C': 3,'D': 4,'E': 5,'F': 6,'G': 7,'H': 8,'I': 9,
           'a': -1,'b': -2,'c': -3,'d': -4,'e': -5,'f': -6,'g': -7,'h': -8,'i': -9}
    DIF = {'%': 0,'J': 1,'K': 2,'L': 3,'M': 4,'N': 5,'O': 6,'P': 7,'Q': 8,'R': 9,
           'j': -1,'k': -2,'l': -3,'m': -4,'n': -5,'o': -6,'p': -7,'q': -8,'r': -9}
    DUP = {'S': 1,'T': 2,'U': 3,'V': 4,'W': 5,'X': 6,'Y': 7,'Z': 8,'s': 9}
    
    
    currentNum = ''
    step = ''
    dup = ''
    numberList = []
    
    last = False
    for char in dat:
        if char in '0123456789':
                if dup != "":
                    dup = dup + char
                elif currentNum != '':
                    currentNum = currentNum + char
                elif step != '':
                    step = step + char
                else:
                    continue
    
        elif char in SQZ.keys():
            dup,currentNum,step,numberList = checkWrite(dup,currentNum,step,numberList)
            
            currentNum = currentNum + str(SQZ[char])
            
        elif char in DIF.keys():
            dup,currentNum,step,numberList = checkWrite(dup,currentNum,step,numberList)
                
            step = step + str(DIF[char])     
        elif char in DUP.keys():
            dup = dup + str(DUP[char]) #For now, assume no SQZ defore DUP
        elif char == ' ':
            last = True
            break
    
    dup,currentNum,step,numberList = checkWrite(dup,currentNum,step,numberList)
    if last:
       return np.array(numberList)
    else:
        return np.array(numberList)[:-1] 
        
        
def loadJCAMP(filePath,name):
    with open(filePath, 'r') as f:
        data = f.read().split('\n')


    realDataPos = []   
    imagDataPos = []
    spectDataPos = []
    currentPos = 0    
    for line in data:
        testline = re.sub('[\t ]*','', line)

        if '#.OBSERVEFREQUENCY=' in testline:
            freq = float(line[line.index('=')+1:]) * 1e6
        elif '##DATATYPE=' in testline:   
            dataType = line[line.index('=')+1:]
        elif '#VAR_DIM=' in testline:
            nPoints = line[line.index('=')+1:]
            nPoints = re.sub(',[\t ][\t ]*',' ', nPoints)
            nPoints = re.sub('[\t\r]*','', nPoints)
            nPoints = int(nPoints.split()[0])
        elif '#VAR_FORM=' in testline:
            varForm = line[line.index('=')+1:]
            varForm = re.sub(',[\t ][\t ]*',' ', varForm)
            varForm = re.sub('[\t\r]*','', varForm)
            varForm = varForm.split()
        elif '#UNITS=' in testline:
            units = line[line.index('=')+1:]
            units = re.sub(',[\t ][\t ]*',' ', units)
            units = re.sub('[\t\r]*','', units)
            units = units.split()[0]
        elif '#FIRST=' in testline:
            first = line[line.index('=')+1:]
            first = re.sub(',[\t ][\t ]*',' ', first)
            first = re.sub('[\t\r]*','', first)
            first = float(first.split()[0].replace(',', '.'))
        elif '#LAST=' in testline:
            last = line[line.index('=')+1:]
            last = re.sub(',[\t ][\t ]*',' ', last)
            last = re.sub('[\t\r]*','', last)
            last = float(last.split()[0].replace(',', '.'))
        elif '#FACTOR=' in testline:
            factor = line[line.index('=')+1:]
            factor = re.sub(',[\t ][\t ]*',' ', factor)
            factor = re.sub('[\t\r]*','', factor)
            factor = factor.split()
            for elem in range(len(factor)):
                factor[elem] = float(factor[elem].replace(',', '.'))
        elif '(X++(R..R))' in testline:
            realDataPos.append(currentPos + 1)
        elif '#PAGE=' in testline and len(realDataPos) == 1:
            realDataPos.append(currentPos - 1)
        elif '(X++(I..I))' in testline:
            imagDataPos.append(currentPos + 1)
        elif '#ENDNTUPLES=' in testline and len(imagDataPos) == 1:
            imagDataPos.append(currentPos - 1)  
        #Spectrum specific
        elif   '(X++(Y..Y))' in testline:  
            spectDataPos.append(currentPos + 1)   
        elif   '##END' in testline and len(spectDataPos) == 1:  
            spectDataPos.append(currentPos - 1) 
        elif   '##XUNITS=' in testline: 
            Xunit = line[line.index('=')+1:]
            Xunit = re.sub('[ \t]*','', Xunit)
        elif '##FIRSTX=' in testline:
            FirstX = float(line[line.index('=')+1:])
        elif '##LASTX=' in testline:
            LastX = float(line[line.index('=')+1:])
        elif '##YFACTOR=' in testline:
            YFactor= float(line[line.index('=')+1:])
        elif '##NPOINTS=' in testline:
            NPoints= int(line[line.index('=')+1:])    
        currentPos += 1

    #Convert the data
    if 'NMR FID' in dataType:
        realDat = np.array([])
        if varForm[1] == 'ASDF': #If DIFDUB form
            for line in data[realDataPos[0]:realDataPos[1]+1]:
                realDat = np.append(realDat,convertDIFDUB(line))
        elif varForm[1] == 'AFFN': #If regular list form
            for line in data[realDataPos[0]:realDataPos[1]+1]:
                realDat = np.append(realDat,np.fromstring(line,sep=' ')[1:])
        realDat = realDat * factor[1]
            
        imagDat = np.array([])  
        if varForm[2] == 'ASDF':
            for line in data[imagDataPos[0]:imagDataPos[1]+1]:
                imagDat = np.append(imagDat,convertDIFDUB(line))
        elif varForm[1] == 'AFFN':
            for line in data[imagDataPos[0]:imagDataPos[1]+1]:
                imagDat = np.append(imagDat,np.fromstring(line,sep=' ')[1:])       
        imagDat = imagDat * factor[2]  
        fullData = realDat - 1j *  imagDat
        sw = 1.0/((last - first)/(nPoints-1))
        masterData = sc.Spectrum(name, fullData, (10, filePath), [freq], [sw], [False])
    elif 'NMRSPECTRUM' in dataType:
        spectDat = np.array([])
        for line in data[spectDataPos[0]:spectDataPos[1]+1]:
            spectDat = np.append(spectDat,convertDIFDUB(line))
        spectDat = np.flipud(spectDat) * YFactor
        if Xunit == 'HZ':
            sw = abs(FirstX - LastX) 
            sw = sw + sw / NPoints
        elif Xunit == 'PPM':
            sw = abs(FirstX - LastX) * freq
            sw = sw + sw / NPoints
        masterData = sc.Spectrum(name, spectDat, (10, filePath), [freq], [sw], [True], ref = [None])
    
    return masterData
