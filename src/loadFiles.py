#!/usr/bin/env python

# Copyright 2016 - 2018 Bas van Meerten and Wouter Franssen

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
import re
import os
import spectrum_classes as sc
import hypercomplex as hc

def loading(num, filePath, name=None, realpath=False, dialog=None):
    if num == 0:
        masterData = loadVarianFile(filePath, name)
    elif num == 1:
        masterData = loadBrukerTopspin(filePath, name)
    elif num == 2:
        masterData = loadChemFile(filePath, name)
    elif num == 3:
        masterData = loadMagritek(filePath, name, realpath)
    elif num == 4:
        masterData = loadSimpsonFile(filePath, name)
    elif num == 5:
        masterData = loadJSONFile(filePath, name)
    elif num == 6:
        masterData = loadMatlabFile(filePath, name)
    elif num == 7:
        masterData = loadBrukerSpectrum(filePath, name)
    elif num == 8:
        masterData = loadPipe(filePath, name)
    elif num == 9:
        masterData = loadJEOLDelta(filePath, name)
    elif num == 10:
        masterData = loadJCAMP(filePath, name)
    elif num == 11:
        if dialog is None:
            return None
        masterData = loadAscii(filePath, name, dialog.dataDimension, dialog.dataSpec, dialog.dataOrder, dialog.delim, dialog.sw)
    elif num == 12:
        masterData = loadMinispec(filePath, name)
    elif num == 13:
        masterData = loadBrukerEPR(filePath, name)
    return masterData

def fileTypeCheck(filePath):
    returnVal = 0
    fileBase = ''
    direc = filePath
    if os.path.isfile(filePath):
        filename = os.path.basename(filePath)
        fileBase = os.path.splitext(filename)[0]
        direc = os.path.dirname(filePath)
        if filename.endswith('.fid') or filename.endswith('.spe'):
            with open(filePath, 'r') as f:
                check = int(np.fromfile(f, np.float32, 1))
            if check == 0:
                return (8, filePath, returnVal)  # Suspected NMRpipe format
            else:  # SIMPSON
                return (4, filePath, returnVal)
        if filename.endswith('.ft') or filename.endswith('.ft1') or filename.endswith('.ft2') or filename.endswith('.ft3') or filename.endswith('.ft4'):
            with open(filePath, 'r') as f:
                check = int(np.fromfile(f, np.float32, 1))
            if check == 0:
                return (8, filePath, returnVal)  # Suspected NMRpipe format
        elif filename.endswith('.json') or filename.endswith('.JSON'):
            return (5, filePath, returnVal)
        elif filename.endswith('.mat') or filename.endswith('.MAT'):
            return (6, filePath, returnVal)
        elif filename.endswith('.jdf'):  # JEOL delta format
            return (9, filePath, returnVal)
        elif filename.endswith('.dx') or filename.endswith('.jdx') or filename.endswith('.jcamp'):  # JCAMP format
            return (10, filePath, returnVal)
        elif filename.endswith('.sig'):  # Bruker minispec
            return (12, filePath, returnVal)
        returnVal = 1
        direc = os.path.dirname(filePath)
    if os.path.exists(direc + os.path.sep + 'procpar') and os.path.exists(direc + os.path.sep + 'fid'):
        return (0, direc, returnVal)
        # And for varian processed data
    if (os.path.exists(direc + os.path.sep + '..' + os.path.sep + 'procpar') or os.path.exists(direc + os.path.sep + 'procpar')) and os.path.exists(direc + os.path.sep + 'data'):
        return (0, direc, returnVal)
    elif os.path.exists(direc + os.path.sep + 'acqus') and (os.path.exists(direc + os.path.sep + 'fid') or os.path.exists(direc + os.path.sep + 'ser')):
        return (1, direc, returnVal)
    elif os.path.exists(direc + os.path.sep + 'procs') and (os.path.exists(direc + os.path.sep + '1r') or os.path.exists(direc + os.path.sep + '2rr') or os.path.exists(direc + os.path.sep + '3rrr')):
        return (7, direc, returnVal)
    elif os.path.exists(direc + os.path.sep + 'acq') and os.path.exists(direc + os.path.sep + 'data'):
        return (2, direc, returnVal)
    elif os.path.exists(direc + os.path.sep + 'acqu.par'):
        dirFiles = os.listdir(direc)
        files2D = [x for x in dirFiles if '.2d' in x]
        files1D = [x for x in dirFiles if '.1d' in x]
        if len(files2D) != 0 or len(files1D) != 0:
            return (3, direc, returnVal)
    elif os.path.exists(direc + os.path.sep + fileBase + '.spc') and os.path.exists(direc + os.path.sep + fileBase + '.par'):
        return (13, direc + os.path.sep + fileBase, returnVal)
    elif os.path.isfile(filePath):  # If not recognised, load as ascii
        return (11, filePath, returnVal)
    return (None, filePath, 2)

def loadVarianFile(filePath, name=''):
    from struct import unpack
    if os.path.isfile(filePath):
        Dir = os.path.dirname(filePath)
    else:
        Dir = filePath
    if os.path.exists(Dir + os.path.sep + 'procpar'):
        file = Dir + os.path.sep + 'procpar'
    elif os.path.exists(Dir + os.path.sep + '..' + os.path.sep + 'procpar'):
        file = Dir + os.path.sep + '..' + os.path.sep + 'procpar'
    else:
        file = None
    Out = [0,1,1,0,0,0,0,0] 
    indirectRef = 'dfrq'
    if file is not None:
        with open(file, 'r') as f:
            data = f.read().split('\n')
        for s in range(0, len(data)): #First check indirect ref name
            if data[s].startswith('refsource1' + " "):
                indirectRef = data[s + 1].split()[1][1:-1]
        f, fM, S = lambda x: float(x), lambda x: float(x) * 1e6, lambda x: x
        Elem = [['sfrq', fM],['sw', f],['sw1', f],['reffrq',fM],['reffrq1',fM],['rp',f],['phfid',f],[indirectRef,fM]]
        for s in range(0, len(data)):
            for index in range(len(Elem)):
                if data[s].startswith(Elem[index][0] + ' '):
                   Out[index] = Elem[index][1]((data[s + 1].split()[1]))
            #elif data[s].startswith('phfid '):
            #    if int(data[s].split()[-2]): #if on
            #        phfid = float(data[s + 1].split()[1]) 
    freq, sw, sw1, reffreq, reffreq1, rp, phfid, freq1 = Out 
    if os.path.exists(Dir + os.path.sep + 'fid'):
        filePath = Dir + os.path.sep + 'fid'
    elif os.path.exists(Dir + os.path.sep + 'data'):
        filePath = Dir + os.path.sep + 'data'
    with open(filePath, "rb") as f:
        nblocks, ntraces, npoints, ebytes, tbytes, bbytes  = np.fromfile(f, np.int32, 6).newbyteorder('>l')
        status = np.fromfile(f, np.int16, 2).newbyteorder('>h')[1]
        status = '{0:016b}'.format(status)[::-1] #send to zeropadded string, and put order correct
        spec, fid32, fidfloat, hypercomplex, flipped = np.array([bool(int(x)) for x in status])[[1,2,3,5,9]]
        nbheaders = np.fromfile(f, np.int32, 1).newbyteorder('>l')[0]
        SizeTD2 = npoints
        SizeTD1 = nblocks * ntraces
        if fidfloat:
            bitType = ['>f', np.float32, 7] # [bitorder, read as, number of elements in nbheader]
        else:
            if fid32:
                bitType = ['>l', np.int32, 7]
            else:
                bitType = ['>h', np.int16, 14]
        totalpoints = (ntraces * npoints + nbheaders**2 * bitType[2])*nblocks
        fid = np.fromfile(f, bitType[1], totalpoints).newbyteorder(bitType[0]).astype(np.complex128)
        if not spec or (spec and not hypercomplex):
            fid = fid.reshape(nblocks, int(totalpoints / nblocks))
            fid = fid[:, bitType[2]::] #Cut off block headers
            fid = fid[:, ::2] - 1j * fid[:, 1::2]
        else:
            fid = fid[nbheaders * bitType[2]::4] - 1j * fid[nbheaders * bitType[2] + 1::4]
            fid = fid.reshape(int(SizeTD1 / 4), SizeTD2)
            if flipped:
                fid = np.fliplr(fid)
    if spec == 0:
        fid = fid * np.exp((rp + phfid) / 180 * np.pi * 1j)  # apply zero order phase
    if SizeTD1 == 1:
        fid = fid[0][:]
        if spec:  # flip if spectrum
            fid = np.flipud(fid)
        masterData = sc.Spectrum(name, fid, (0, filePath), [freq], [sw], [bool(int(spec))], ref=[reffreq])
    else:
        masterData = sc.Spectrum(name, fid, (0, filePath), [freq1, freq], [sw1, sw], [bool(int(spec))] * 2, ref=[reffreq1, reffreq])
    masterData.addHistory("Varian data loaded from " + filePath)
    return masterData

def loadPipe(filePath, name=''):
    with open(filePath, 'r') as f:
        header = np.fromfile(f, np.float32, 512)
    NDIM = int(header[9])
    SIZE = [int(header[32]),int(header[15]),int(header[219]),int(header[99])]
    quadFlag = [int(header[54]),int(header[51]),int(header[55]),int(header[56])] #0 complex, 1 real        
    spec = [int(header[31]),int(header[13]),int(header[222]),int(header[220])]  # 1 if ft, 0 if time
    freq = np.array([header[28],header[10],header[218],header[119]]) * 1e6
    sw = [header[29],header[11],header[229],header[100]]
    ref = [header[30],header[12],header[249],header[101]]  # frequency of last point in Hz
    numFiles = int(header[442]) #Number of files to be loaded
    pipeFlag = int(header[57]) #Indicates stream, or separate files
    cubeFlag = int(header[447]) #1 if 4D, and saved as sets of 3D
    for i in range(len(spec)): #get reference frequencies
        sidefreq = -np.floor(SIZE[i] / 2) / SIZE[i] * sw[i]  # frequency of last point on axis
        ref[i] = sidefreq + freq[i] - ref[i]
    TotP = SIZE[3] * SIZE[2] #Max file size
    if quadFlag[3] == 0:  # if complex direct axis
        TotP = TotP * 2
    if NDIM == 4 and cubeFlag == 1 and pipeFlag == 0:
        TotP *= SIZE[1] #If file 3D format, load more data points
    elif NDIM == 4 and pipeFlag > 0:
        TotP *= SIZE[0] * SIZE[1]
    if NDIM == 3 and pipeFlag > 0:
        TotP *= SIZE[1]
    if numFiles  == 1:
        files = [filePath]
    else: #Get the names of the files, if more than 1
        dir,file = os.path.split(filePath)
        base, ext = os.path.splitext(file)
        start = re.search('[0-9]+$',base).start()
        basename = base[:start]
        numbers = len(base[start:])
        files = [dir + os.path.sep + basename + str(x + 1).zfill(numbers) + ext for x in range(numFiles)]
    data = []
    for file in files: #Load all the data from the files
        with open(file, 'r') as f:
            tmp = np.fromfile(f, np.float32, 512)
            data.append( np.fromfile(f, np.float32, TotP))
    for i in range(len(data)): #Reshape all the data
        if NDIM > 1 and cubeFlag == 0 and pipeFlag == 0: #Reshape 2D sets if needed
            data[i] = np.reshape(data[i], (SIZE[2],int(TotP/SIZE[2])))
        elif NDIM == 4 and cubeFlag == 1 and pipeFlag == 0: #For 4D, in sets of 3D
            data[i] = np.reshape(data[i], (SIZE[1],SIZE[2],int(TotP/SIZE[2]/SIZE[1])))
        elif NDIM == 4 and pipeFlag > 0: #For stream
            data[i] = np.reshape(data[i], (SIZE[0],SIZE[1],SIZE[2],int(TotP/SIZE[2]/SIZE[1]/SIZE[0])))
        elif NDIM ==3 and pipeFlag > 0: #For stream
            data[i] = np.reshape(data[i], (SIZE[1],SIZE[2],int(TotP/SIZE[2]/SIZE[1])))
    #For 3D or 4D data, merge the datasets
    if NDIM > 2 and pipeFlag == 0:
        data = [np.array(data)]
    eS = (slice(None),) #empty slice
    if quadFlag[3] == 0: #If complex along last dim
        useSlice = eS * (NDIM - 1)
        data[0] = data[0][useSlice + (slice(None,SIZE[3],None),)] + 1j * data[0][useSlice + (slice(SIZE[3],None,None),)]
    hyper = np.array([0])
    if NDIM > 1: #Reorder data, if hypercomplex along an axis
        for dim in range(NDIM - 1):
            newdata = []
            if quadFlag[4 - NDIM + dim] == 0:
                useSlice1 = eS * dim + (slice(None,None,2),) + eS * (NDIM - dim - 1)
                useSlice2 = eS * dim + (slice(1,None,2),) + eS * (NDIM - dim - 1)
                for dat in data:
                    newdata.append(dat[useSlice1])
                    newdata.append(dat[useSlice2])
                data = newdata
                hyper = np.append(hyper, hyper + 2**dim)
    for k in range(len(data)): #Flip LR if spectrum axis
        for i in range(NDIM):
            if spec[-1 - i] == 1: 
                data[k] = np.flip(data[k], NDIM -1 - i)
    masterData = sc.Spectrum(name, hc.HComplexData(data, hyper), (8, filePath), freq[4 - NDIM:4], sw[4 - NDIM:4], spec[4 - NDIM:4], ref=ref[4 - NDIM:4])
    masterData.addHistory("NMRpipe data loaded from " + filePath)
    return masterData

def loadJEOLDelta(filePath, name=''):
    from struct import unpack
    multiUP = lambda typ, bit, num, start: np.array([unpack(typ,header[start + x:start + bit + x])[0] for x in range(0,num * bit,bit)])
    with open(filePath, "rb") as f:
        header = f.read(1296)
    endian =['>d','<d'][multiUP('>B', 1, 1, 8)[0]]
    NDIM =  multiUP('>B', 1, 1, 12)[0]
    #data_dimension_exist = multiUP('>B', 1, 1, 13)[0]
    #data_type = multiUP('>B', 1, 1, 14)[0]
    #translate = multiUP('>B', 1, 8, 16)
    dataType = multiUP('>B', 1, 8, 24)
    dataUnits = multiUP('>B', 1, 16, 32).reshape(8, 2)
    NP = multiUP('>I', 4, 8, 176)
    #dataStart = multiUP('>I', 4, 8, 208)
    dataStop = multiUP('>I', 4, 8, 240)
    axisStart = multiUP('>d', 8, 8, 272)
    axisStop = multiUP('>d', 8, 8, 336)
    baseFreq = multiUP('>d', 8, 8, 1064)
    #zero_point = multiUP('>d', 8, 8, 1128)
    reverse = multiUP('>B', 1, 8, 1192)
    readStart = multiUP('>I', 4, 1, 1284)[0]
    #data_length = multiUP('>Q', 8, 1, 1288)[0]
    loadSize = np.cumprod(NP[:NDIM])[-1]
    if NDIM == 1 and dataType[0] == 3: #Complex 1D
        loadSize *= 2
    elif NDIM == 2 and dataType[0] == 4: #2D Real-Complex (non-Hypercomplex) 
        loadSize *= 2
    elif NDIM == 2 and dataType[0] == 3: #2D Complex-Complex (Hypercomplex) 
        loadSize *= 4
    with open(filePath, "rb") as f:
        f.seek(readStart) #Set read start to position of data
        data = np.fromfile(f, endian, loadSize)
    hyper = np.array([0])
    if NDIM == 1 and dataType[0] == 1: #Real 1D
        data = [data]
    elif NDIM == 1 and dataType[0] == 3: #Complex 1D
        data = data[:NP[0]] - 1j * data[NP[0]:]
        data = [data[0:data_offset_stop[0] + 1]]
    elif NDIM == 2 and dataType[0] == 4: #2D Real-Complex (non-Hypercomplex) 
        Step = 4
        data = data[:int(loadSize/2)] - 1j * data[int(loadSize/2):]
        data = np.reshape(data, [int(NP[1] / Step), int(NP[0] / Step), Step, Step])
        data = [np.concatenate(np.concatenate(data, 1), 1)]
    elif NDIM == 2 and dataType[0] == 3: #2D Complex-Complex (Hypercomplex) 
        hyper = np.array([0, 1])
        Step = 32  # Step size of block
        tmp = np.split(data,4)
        data = [tmp[0] - 1j * tmp[1], tmp[2] - 1j * tmp[3]]
        del tmp
        for i in range(len(data)):
            data[i] = np.reshape(data[i], [int(NP[1] / Step), int(NP[0] / Step), Step, Step])
            data[i] = np.concatenate(np.concatenate(data[i], 1), 1)
    eS = (slice(None),) #empty slice
    for dim in range(NDIM): #Cut data for every dim
        useSlice = eS * (NDIM - dim - 1) +(slice(0,dataStop[dim],None),) + eS * dim 
        for i in range(len(data)):
            data[i] = data[i][useSlice]
    freq = baseFreq[0:NDIM][::-1] * 1e6
    spec = dataUnits[0:NDIM,1][::-1] != 28 #If not 28 (Hz), then spec = true
    sw = []
    ref = []
    for axisNum in reversed(range(NDIM)):
        axisType = dataUnits[axisNum][1]  # Sec = 28, Hz = 13, PPM = 26
        axisScale = dataUnits[axisNum][0]
        if axisType == 28:  # Sec
            scale = (axisScale >> 4) & 15
            if scale > 7:
                scale = scale - 16
            dw = (axisStop[axisNum] - axisStart[axisNum]) / (dataStop[axisNum] + 1 - 1) * 10.0**(-scale * 3)  # minus one to give same axis as spectrum???
            # scale for SI prefix
            sw.append(1.0 / dw)
            sidefreq = -np.floor((dataStop[axisNum] + 1) / 2) / dataStop[axisNum] + 1 * sw[-1]  # frequency of last point on axis
            ref.append(baseFreq[axisNum] * 1e6)
        if axisType == 13:  # Hz
            sw.append(np.abs(axisStart[axisNum] - axisStop[axisNum]))
            sidefreq = -np.floor((dataStop[axisNum] + 1) / 2) / (dataStop[axisNum] + 1) * sw[-1]  # frequency of last point on axis
            ref.append(sidefreq + baseFreq[axisNum] * 1e6 - axisStop[axisNum])
        if axisType == 26:  # ppm
            sw.append(np.abs(axisStart[axisNum] - axisStop[axisNum]) * baseFreq[axisNum])
            sidefreq = -np.floor((dataStop[axisNum] + 1) / 2) / (dataStop[axisNum] + 1) * sw[-1]  # frequency of last point on axis
            ref.append(sidefreq + baseFreq[axisNum] * 1e6 - axisStop[axisNum] * baseFreq[axisNum])
    for k in range(len(data)): #Flip LR if spectrum axis
        for i in range(NDIM):
            if spec[-1 - i] == 1: 
                data[k] = np.flip(data[k], NDIM -1 - i)
    masterData = sc.Spectrum(name, hc.HComplexData(np.array(data), hyper), (9, filePath), freq, sw, spec, ref=ref)
    masterData.addHistory("JEOL Delta data loaded from " + filePath)
    return masterData

def saveJSONFile(filePath, spectrum):
    import json
    struct = {}
    item = spectrum.data
    struct['dataReal'] = np.real(item.data).tolist()
    struct['dataImag'] = np.imag(item.data).tolist()
    struct['hyper'] = item.hyper.tolist()
    struct['freq'] = spectrum.freq.tolist()
    struct['sw'] = list(spectrum.sw)
    struct['spec'] = list(1.0 * np.array(spectrum.spec))
    struct['wholeEcho'] = list(1.0 * np.array(spectrum.wholeEcho))
    struct['ref'] = np.array(spectrum.ref, dtype=np.float).tolist()
    struct['history'] = spectrum.history
    tmpXax = []
    for i in spectrum.xaxArray:
        tmpXax.append(i.tolist())
    struct['xaxArray'] = tmpXax
    with open(filePath, 'w') as outfile:
        json.dump(struct, outfile)

def loadJSONFile(filePath, name=''):
    import json
    with open(filePath, 'r') as inputfile:
        struct = json.load(inputfile)
    if 'hyper' in struct.keys():
        hyper = list(struct['hyper'])
        tmpReal = struct['dataReal']
        tmpImag = struct['dataImag']
        data = np.array(tmpReal) + 1j * np.array(tmpImag)
    else:
        hyper = [0]
        data = np.array([np.array(struct['dataReal']) + 1j * np.array(struct['dataImag'])])
    ref = np.where(np.isnan(struct['ref']), None, struct['ref'])
    if 'history' in struct.keys():
        history = struct['history']
    else:
        history = None
    xaxA = []
    for i in struct['xaxArray']:
        xaxA.append(np.array(i))
    masterData = sc.Spectrum(name,
                             hc.HComplexData(data, hyper),
                             (5, filePath),
                             list(struct['freq']),
                             list(struct['sw']),
                             list(struct['spec']),
                             list(np.array(struct['wholeEcho'], dtype=bool)),
                             list(ref),
                             xaxA,
                             history=history)
    masterData.addHistory("JSON data loaded from " + filePath)
    return masterData

def saveMatlabFile(filePath, spectrum, name='spectrum'):
    import scipy.io
    struct = {}
    struct['dim'] = spectrum.ndim()
    struct['data'] = spectrum.data.data
    struct['hyper'] = spectrum.data.hyper
    struct['freq'] = spectrum.freq
    struct['sw'] = spectrum.sw
    struct['spec'] = spectrum.spec
    struct['wholeEcho'] = spectrum.wholeEcho
    struct['ref'] = np.array(spectrum.ref, dtype=np.float)
    struct['history'] = spectrum.history
    struct['xaxArray'] = spectrum.xaxArray
    matlabStruct = {name: struct}
    scipy.io.savemat(filePath, matlabStruct)

def loadMatlabFile(filePath, name=''):
    import scipy.io
    with open(filePath, 'rb') as inputfile:  # read first several bytes the check .mat version
        teststring = inputfile.read(13)
    version = float(teststring.decode("utf-8")[7:10])  # extract version from the binary array
    if version < 7.3:  # all versions below 7.3 are supported
        matlabStruct = scipy.io.loadmat(filePath)
        var = [k for k in matlabStruct.keys() if not k.startswith('__')][0]
        mat = matlabStruct[var]
        if 'hyper' in mat.dtype.names:
            if len(mat['hyper'][0,0]) == 0:
                hyper = None
            else:
                hyper = mat['hyper'][0,0][0]
        else:
            hyper = None
        data = []
        if mat['dim'] == 1:
            if hyper is None:
                data = np.array(mat['data'][0][0][0])
            else:
                data = np.array(mat['data'][0][0])
            xaxA = [k[0] for k in (mat['xaxArray'][0])]
        else:
            if hyper is None: #If old format
                data = np.array(mat['data'][0, 0])
            else: #If new format
                data = np.array(mat['data'][0][0])
            if all(x == data[0].shape[0] for x in data[0].shape):
                xaxA = [k for k in (mat['xaxArray'][0, 0])]
            else:
                xaxA = [k[0] for k in (mat['xaxArray'][0, 0][0])]
        # insert some checks for data type
        ref = mat['ref'][0, 0][0]
        ref = np.where(np.isnan(ref), None, ref)
        if 'history' in mat.dtype.names:
            history = list(np.array(mat['history'][0, 0], dtype=str))
        else:
            history = None
        masterData = sc.Spectrum(name,
                                 hc.HComplexData(data, hyper),
                                 (6, filePath),
                                 list(mat['freq'][0, 0][0]),
                                 list(mat['sw'][0, 0][0]),
                                 list(mat['spec'][0, 0][0]),
                                 list(np.array(mat['wholeEcho'][0, 0][0]) > 0),
                                 list(ref),
                                 xaxA,
                                 history=history)
        masterData.addHistory("Matlab data loaded from " + filePath)
        return masterData
    else:  # If the version is 7.3, use HDF5 type loading
        import h5py  # For .mat v 7.3 support
        f = h5py.File(filePath, 'r')
        Groups = []
        for name in f:
            if name != '#refs#':
                Groups.append(name)
        DataGroup = Groups[0]  # get the group name
        mat = f[DataGroup]
        if 'hyper' in mat.dtype.names:
            hyper = np.array(mat['hyper'])
        else:
            hyper = None
        if np.array(mat['dim'])[0][0] == 1:
            xaxA = list([np.array(mat['xaxArray'])[:, 0]])
            if hyper is None: #Old data format
                data = np.array(mat['data'])
                data = [(data['real'] + data['imag'] * 1j)[:, 0]]  # split and use real and imag part
            else: #New format (hypercomplex)
                data = np.transpose(np.array(mat['data']))
                data = list(data['real'] + data['imag'] * 1j)
        else:
            if hyper is None:
                data = np.transpose(np.array(mat['data']))
                data = [data['real'] + data['imag'] * 1j]
            else:
                tmp = np.transpose(np.array(mat['data']))
                data = []
                for index in range(tmp.shape[0]):
                    data.append(tmp['real'][index] + tmp['imag'][index] * 1j)
            if all(x == data[0].shape[0] for x in data[0].shape):
                xaxA = [np.array(mat[k]) for k in (mat['xaxArray'])]
            else:
                xaxA = [np.array(mat[k[0]]) for k in (mat['xaxArray'])]
        ref = np.array(mat['ref'])[:, 0]
        ref = np.where(np.isnan(ref), None, ref)
        if 'history' in mat.keys():
            history = list()
            history.append([item.astype(np.int8).tostring().decode("ascii") for item in np.array(mat['history']).transpose()])
            history = history[0]
        else:
            history = None
        masterData = sc.Spectrum(name,
                                 hc.HComplexData(data, hyper),
                                 (6, filePath),
                                 list(np.array(mat['freq'])[:, 0]),
                                 list(np.array(mat['sw'])[:, 0]),
                                 list(np.array(mat['spec'])[:, 0]),
                                 list(np.array(mat['wholeEcho'])[:, 0] > 0),
                                 list(ref),
                                 xaxA,
                                 history=history)
        masterData.addHistory("Matlab data loaded from " + filePath)
        return masterData

def loadBrukerTopspin(filePath, name=''):
    if os.path.isfile(filePath):
        Dir = os.path.dirname(filePath)
    else:
        Dir = filePath
    f,fM, i = lambda x: float(x), lambda x: float(x) * 1e6, lambda x: int(x) #Conversion functions
    Elem = [['TD', i, []],['SFO1', fM ,[]],['SW_h', f, []],['O1',f, []], ['BYTORDA',i,[]]] #The elements to be found [Name, conversion, list with hits]
    for File in ['acqu','acqu2','acqu3']:
        if os.path.exists(Dir + os.path.sep + File):
            with open(Dir + os.path.sep + File, 'r') as f:
                data = f.read().split('\n')
            for s in range(0, len(data)):
                for var in Elem:
                    if data[s].startswith('##$' + var[0] + '='):
                        var[2].append( var[1](re.sub('##\$' + var[0] + '=', '', data[s])))
    SIZE, FREQ, SW, REF, BYTE = [x[2] for x in Elem] #Unpack results
    ByteOrder = ['l','b'][BYTE[0]] #The byte orders that is used 
    REF = list(- np.array(REF) + np.array(FREQ))
    totsize = np.cumprod(SIZE)[-1]
    dim = len(SIZE)
    directSize = int(np.ceil(float(SIZE[0]) / 256)) * 256 #Size of direct dimension including
    #blocking size of 256 data points
    for file in ['fid','ser']:
        if os.path.exists(Dir + os.path.sep + file):
            if file == 'ser':
                totsize = int(totsize / SIZE[0]) * directSize #Always load full 1024 byte blocks (256 data points) for >1D
            with open(Dir + os.path.sep + file, "rb") as f:
                raw = np.fromfile(f, np.int32, totsize)
            raw = raw.newbyteorder(ByteOrder) #Load with right byte order
    ComplexData = np.array(raw[0:len(raw):2]) + 1j * np.array(raw[1:len(raw):2])
    if dim >= 2:
        newSize = list(SIZE)
        newSize[0] = int(directSize / 2)
        ComplexData = ComplexData.reshape(*newSize[-1::-1])
    if dim == 2:
        ComplexData = ComplexData[:,0:int(SIZE[0]/2)] #Cut off placeholder data
    elif dim == 3:
        ComplexData = ComplexData[:,:,0:int(SIZE[0]/2)] #Cut off placeholder data
    masterData = sc.Spectrum(name, ComplexData, (1, filePath), FREQ[-1::-1], SW[-1::-1], [False] * dim, ref = REF[-1::-1])
    masterData.addHistory("Bruker data loaded from " + filePath)
    return masterData

def loadBrukerSpectrum(filePath, name=''):
    if os.path.isfile(filePath):
        Dir = os.path.dirname(filePath)
    else:
        Dir = filePath
    f,fM, i = lambda x: float(x), lambda x: float(x) * 1e6, lambda x: int(x) #Conversion functions
    Elem = [['SI', i, []],['XDIM', i ,[]],['SW_p', f, []],['SF',fM, []],['OFFSET', f, []], ['BYTORDP',i,[]]] #The elements to be found [Name, conversion, list with hits]
    #OFFSET: The highest ppm values of the axis
    #XDIM: The blocking size along an axis
    for File in ['procs','proc2s','proc3s']:
         if os.path.exists(Dir + os.path.sep + File):
            with open(Dir + os.path.sep + File, 'r') as f:
                data = f.read().split('\n')
            for s in range(0, len(data)):
                for var in Elem:
                    if data[s].startswith('##$' + var[0] + '='):
                        var[2].append( var[1](re.sub('##\$' + var[0] + '=', '', data[s])))
    SIZE, XDIM, SW, FREQ, OFFSET, BYTE = [x[2] for x in Elem] #Unpack results
    ByteOrder = ['l','b'][BYTE[0]] #The byte orders that is used 
    REF = []
    for index in range(len(SIZE)): #For each axis
        pos = np.fft.fftshift(np.fft.fftfreq(SIZE[index], 1.0 / SW[index]))[-1] #Get last point of axis
        pos2 = OFFSET[index] * 1e-6 * FREQ[index] #offset in Hz
        REF.append(FREQ[index] + pos - pos2)
    totsize =  np.cumprod(SIZE)[-1] 
    dim = len(SIZE)
    DATA = []
    files = [['1r','1i'],['2rr','2ir','2ri','2ii'],['3rrr','3irr','3rir','3iir','3rri','3iri','3rii','3iii']]
    counter = 0
    for file in files[dim - 1]: # For all the files
        if os.path.exists(Dir + os.path.sep + file): 
            with open(Dir + os.path.sep + file, "rb") as f:
                raw = np.fromfile(f, np.int32, totsize)
                raw = raw.newbyteorder(ByteOrder) # Set right byteorder
                if counter % 2 == 0: # If even, data is real part
                    DATA.append(np.flipud(raw))
                else: # If odd, data is imag, and needs to be add to the previous
                    DATA[-1] = DATA[-1] - 1j * np.flipud(raw)
                counter += 1 # only advance counter when file is found
    del raw
    hyper = np.array([0])
    if dim == 2: # If 2D data has more than 1 part: hypercomplex along the first axis
        if len(DATA) != 1:
            hyper = np.array([0, 1])
    if len(SIZE) == 2:
        for index in range(len(DATA)): # For each data set
            # Reshape DATA to 4D data using the block information
            # Twice concat along axis 1 constructs the regular x-y data
            DATA[index] = np.reshape(DATA[index],[int(SIZE[1]/XDIM[1]),int(SIZE[0]/XDIM[0]),XDIM[1],XDIM[0]])
            DATA[index] = np.concatenate(np.concatenate(DATA[index],1),1)
    elif len(SIZE) == 3:
        for index in range(len(DATA)):
            # The same as 2D, but now split to 6D data, and concat along 2
            DATA[index] = np.reshape(DATA[index],[int(SIZE[2]/XDIM[2]),int(SIZE[1]/XDIM[1]),int(SIZE[0]/XDIM[0]),XDIM[2],XDIM[1],XDIM[0]])
            DATA[index] = np.concatenate(np.concatenate(np.concatenate(DATA[index],2),2),2)
    spec = [True]
    masterData = sc.Spectrum(name, hc.HComplexData(DATA, hyper), (7, filePath), FREQ[-1::-1], SW[-1::-1], spec * dim, ref=REF[-1::-1])
    masterData.addHistory("Bruker spectrum data loaded from " + filePath)
    return masterData

def loadChemFile(filePath, name=''):
    if os.path.isfile(filePath):
        Dir = os.path.dirname(filePath)
    else:
        Dir = filePath
    sizeTD1 = 1
    sw1 = 1
    H = dict(line.strip().split('=') for line in open(Dir + os.path.sep + 'acq', 'r'))
    sizeTD2 = int(float(H['al']))
    freq = float(H['sf' + str(int(float(H['ch1'])))])
    sw = 1 / float(H['dw'][:-1])
    if any('array_num_values_' in s for s in H.keys()):
        if 'use_array=1' in open(Dir + '/acq_2').read():
            for s in H.keys():
                if ('array_num_values_' in s):
                    sizeTD1 = sizeTD1 * int(H[s])
        else:
            if 'al2' in H:
                sizeTD1 = int(float(H['al2']))
                if 'dw2' in H:
                    sw1 = 1 / float(H['dw2'][:-1])
    else:
        if 'al2' in H:
            sizeTD1 = int(float(H['al2']))
            if 'dw2' in H:
                sw1 = 1 / float(H['dw2'][:-1])
    with open(Dir + os.path.sep + 'data', 'rb') as f:
        raw = np.fromfile(f, np.int32)
        b = np.complex128(raw.byteswap())
    filePath = Dir + os.path.sep + 'data'
    fid = b[:int(len(b) / 2)] + 1j * b[int(len(b) / 2):]
    fid = np.reshape(fid, (sizeTD1, sizeTD2))
    data = np.array(fid)
    spec = [False]
    if sizeTD1 is 1:
        data = data[0][:]
        masterData = sc.Spectrum(name, data, (2, filePath), [freq * 1e6], [sw], spec)
    else:
        data = data.reshape((sizeTD1, sizeTD2))
        masterData = sc.Spectrum(name, data, (2, filePath), [freq * 1e6] * 2, [sw1, sw], spec * 2)
    masterData.addHistory("Chemagnetics data loaded from " + filePath)
    return masterData

def loadMagritek(filePath, name='', realPath=''):
    # Magritek load script based on some Matlab files by Ole Brauckman
    if os.path.isfile(filePath):
        Dir = os.path.dirname(filePath)
    else:
        Dir = filePath
    if realPath:
        rememberPath = realPath
    else:
        rememberPath = filePath
    DirFiles = os.listdir(Dir)
    Files2D = [x for x in DirFiles if '.2d' in x]
    Files1D = [x for x in DirFiles if '.1d' in x]
    # initialize 2D values to some dummy value
    sizeTD1 = 0
    sw1 = 50e3
    lastfreq1 = None
    ref1 = None
    # Start pars extraction
    H = dict(line.strip().split('=') for line in open(Dir + os.path.sep + 'acqu.par', 'r'))
    for key in H.keys():
        if key.startswith("bandwidth "):
            sw = float(H[key]) * 1000
        elif key.startswith("nrPnts "):
            sizeTD2 = int(H[key])
        elif key.startswith("b1Freq "):
            freq = float(H[key]) * 1e6
        elif key.startswith("lowestFrequency "):
            lastfreq = float(H[key])
        elif key.startswith("nrSteps "):
            sizeTD1 = int(H[key])
        elif key.startswith("bandwidth2 "):
            sw1 = float(H[key]) * 1000
        elif key.startswith("lowestFrequency2 "):
            lastfreq1 = float(H[key])
    sidefreq = -np.floor(sizeTD2 / 2) / sizeTD2 * sw  # freqeuency of last point on axis
    ref = sidefreq + freq - lastfreq
    if len(Files2D) == 1:
        File = Files2D[0]
        if lastfreq1 is not None:
            sidefreq1 = -np.floor(sizeTD1 / 2) / sizeTD1 * sw1  # freqeuency of last point on axis
            ref1 = sidefreq1 + freq - lastfreq1
        with open(Dir + os.path.sep + File, 'rb') as f:
            raw = np.fromfile(f, np.float32)
        Data = raw[-2 * sizeTD2 * sizeTD1::]
        ComplexData = Data[0:Data.shape[0]:2] - 1j * Data[1:Data.shape[0]:2]
        ComplexData = ComplexData.reshape((sizeTD1, sizeTD2))
        masterData = sc.Spectrum(name, ComplexData, (3, rememberPath), [freq] * 2, [sw1, sw], [False] * 2, ref=[ref1, ref])
    elif len(Files1D) != 0:
        File = 'data.1d'
        with open(Dir + os.path.sep + File, 'rb') as f:
            raw = np.fromfile(f, np.float32)
        Data = raw[-2 * sizeTD2::]
        ComplexData = Data[0:Data.shape[0]:2] - 1j * Data[1:Data.shape[0]:2]
        masterData = sc.Spectrum(name, ComplexData, (3, rememberPath), [freq], [sw], [False], ref=[ref])
    masterData.addHistory("Magritek data loaded from " + rememberPath)
    return masterData

def saveSimpsonFile(filePath, spectrum):
    data = spectrum.getHyperData(0) # SIMPSON does not support hypercomplex
    with open(filePath, 'w') as f:
        f.write('SIMP\n')
        if data.ndim == 2:
            f.write('NP=' + str(data.shape[1]) + '\n')
            f.write('NI=' + str(data.shape[0]) + '\n')
            f.write('SW=' + str(spectrum.sw[1]) + '\n')
            f.write('SW1=' + str(spectrum.sw[0]) + '\n')
        else:
            f.write('NP=' + str(data.shape[0]) + '\n')
            f.write('SW=' + str(spectrum.sw[0]) + '\n')
        if spectrum.spec[0]:
            f.write('TYPE=SPE' + '\n')
        else:
            f.write('TYPE=FID' + '\n')
        f.write('DATA' + '\n')
        if data.ndim == 1:
            for line in data:
                f.write(str(line.real) + ' ' + str(line.imag) + '\n')
        if data.ndim == 2:
            points = data.shape
            for i in range(0, points[0]):
                for j in range(0, points[1]):
                    f.write(str(data[i][j].real) + ' ' + str(data[i][j].imag) + '\n')
        f.write('END')

def loadSimpsonFile(filePath, name=''):
    with open(filePath, 'r') as f:
        Lines = f.read().split('\n')
    NP, NI, SW, SW1, TYPE, FORMAT = 0, 1, 0, 0, '', 'Normal'
    DataStart = Lines.index('DATA')
    DataEnd = Lines.index('END')
    for s in range(0, DataStart):
        if Lines[s].startswith('NP='):
            NP = int(re.sub('NP=', '', Lines[s]))
        elif Lines[s].startswith('NI='):
            NI = int(re.sub('NI=', '', Lines[s]))
        elif Lines[s].startswith('SW='):
            SW = float(re.sub('SW=', '', Lines[s]))
        elif Lines[s].startswith('SW1='):
            SW1 = float(re.sub('SW1=', '', Lines[s]))
        elif Lines[s].startswith('TYPE='):
            TYPE = re.sub('TYPE=', '', Lines[s])
        elif Lines[s].startswith('FORMAT='):
            FORMAT = re.sub('FORMAT=', '', Lines[s])
    if 'Normal' in FORMAT:
        length = DataEnd - DataStart - 1
        data = np.zeros(length, dtype=complex)
        for i in range(length):
            temp = Lines[DataStart + 1 + i].split()
            data[i] = float(temp[0]) + 1j * float(temp[1])
    elif 'BINARY' in FORMAT:
        # Binary code based on:
        # pysimpson: Python module for reading SIMPSON files
        # By: Jonathan J. Helmus (jjhelmus@gmail.com)
        # Version: 0.1 (2012-04-13)
        # License: GPL
        chardata = ''.join(Lines[DataStart + 1:DataEnd])
        nquads, mod = divmod(len(chardata), 4)
        assert mod == 0     # character should be in blocks of 4
        BASE = 33
        charst = np.fromstring(chardata, dtype=np.uint8)
        charst = charst.reshape(nquads, 4) - BASE

        def FIRST(f, x): return ((x) & ~(~0 << f))

        def LAST(f, x): return ((x) & (~0 << (8 - f)))
        first = FIRST(6, charst[:, 0]) | LAST(2, charst[:, 1] << 2)
        second = FIRST(4, charst[:, 1]) | LAST(4, charst[:, 2] << 2)
        third = FIRST(2, charst[:, 2]) | LAST(6, charst[:, 3] << 2)
        Bytes = np.ravel(np.transpose(np.array([first, second, third]))).astype('int64')
        # convert every 4 'bytes' to a float
        num_points, num_pad = divmod(len(Bytes), 4)
        Bytes = np.array(Bytes)
        Bytes = Bytes[:-num_pad]
        Bytes = Bytes.reshape(num_points, 4)
        mantissa = ((Bytes[:, 2] % 128) << 16) + (Bytes[:, 1] << 8) + Bytes[:, 0]
        exponent = (Bytes[:, 3] % 128) * 2 + (Bytes[:, 2] >= 128) * 1
        negative = Bytes[:, 3] >= 128
        e = exponent - 127
        m = np.abs(mantissa) / np.float64(1 << 23)
        data = np.float32((-1)**negative * np.ldexp(m, e))
        data = data.view('complex64')
    if NI != 1:  # 2D data, reshape to NI, NP
        data = data.reshape(int(NI), -1)
    if 'FID' in TYPE:
        spec = [False]
    elif 'SPE' in TYPE:
        spec = [True]
    if NI is 1:
        masterData = sc.Spectrum(name, data, (4, filePath), [0], [SW], spec)
    else:
        masterData = sc.Spectrum(name, data, (4, filePath), [0, 0], [SW1, SW], spec * 2)
    masterData.addHistory("SIMPSON data loaded from " + filePath)
    return masterData

def convertDIFDUB(dat):
    def checkWrite(dup, currentNum, step, numberList):
        if dup != '':
            for dupli in range(int(dup)):
                numberList.append(numberList[-1] + int(step))
            dup = ''
            step = ''
        elif currentNum != '':  # Write if available
            numberList.append(int(currentNum))  # write old num
            currentNum = ''
        elif step != '':
            numberList.append(numberList[-1] + int(step))
            step = ''
        return dup, currentNum, step, numberList
    SQZ = {'@': 0, 'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7, 'H': 8, 'I': 9,
           'a': -1, 'b': -2, 'c': -3, 'd': -4, 'e': -5, 'f': -6, 'g': -7, 'h': -8, 'i': -9}
    DIF = {'%': 0, 'J': 1, 'K': 2, 'L': 3, 'M': 4, 'N': 5, 'O': 6, 'P': 7, 'Q': 8, 'R': 9,
           'j': -1, 'k': -2, 'l': -3, 'm': -4, 'n': -5, 'o': -6, 'p': -7, 'q': -8, 'r': -9}
    DUP = {'S': 1, 'T': 2, 'U': 3, 'V': 4, 'W': 5, 'X': 6, 'Y': 7, 'Z': 8, 's': 9}
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
            dup, currentNum, step, numberList = checkWrite(dup, currentNum, step, numberList)
            currentNum = currentNum + str(SQZ[char])
        elif char in DIF.keys():
            dup, currentNum, step, numberList = checkWrite(dup, currentNum, step, numberList)
            step = step + str(DIF[char])
        elif char in DUP.keys():
            dup = dup + str(DUP[char])  # For now, assume no SQZ defore DUP
        elif char == ' ':
            last = True
            break
    dup, currentNum, step, numberList = checkWrite(dup, currentNum, step, numberList)
    if last:
        return np.array(numberList)
    else:
        return np.array(numberList)[:-1]

def loadJCAMP(filePath, name):
    with open(filePath, 'r') as f:
        data = f.read().split('\n')
    realDataPos = []
    imagDataPos = []
    spectDataPos = []
    currentPos = 0
    for line in data:
        testline = re.sub('[\t ]*', '', line)
        if '#.OBSERVEFREQUENCY=' in testline:
            freq = float(line[line.index('=') + 1:]) * 1e6
        elif '##DATATYPE=' in testline:
            dataType = line[line.index('=') + 1:]
        elif '#VAR_DIM=' in testline:
            nPoints = line[line.index('=') + 1:]
            nPoints = re.sub(',[\t ][\t ]*', ' ', nPoints)
            nPoints = re.sub('[\t\r]*', '', nPoints)
            nPoints = int(nPoints.split()[0])
        elif '#VAR_FORM=' in testline:
            varForm = line[line.index('=') + 1:]
            varForm = re.sub(',[\t ][\t ]*', ' ', varForm)
            varForm = re.sub('[\t\r]*', '', varForm)
            varForm = varForm.split()
        elif '#UNITS=' in testline:
            units = line[line.index('=') + 1:]
            units = re.sub(',[\t ][\t ]*', ' ', units)
            units = re.sub('[\t\r]*', '', units)
            units = units.split()[0]
        elif '#FIRST=' in testline:
            first = line[line.index('=') + 1:]
            first = re.sub(',[\t ][\t ]*', ' ', first)
            first = re.sub('[\t\r]*', '', first)
            first = float(first.split()[0].replace(',', '.'))
        elif '#LAST=' in testline:
            last = line[line.index('=') + 1:]
            last = re.sub(',[\t ][\t ]*', ' ', last)
            last = re.sub('[\t\r]*', '', last)
            last = float(last.split()[0].replace(',', '.'))
        elif '#FACTOR=' in testline:
            factor = line[line.index('=') + 1:]
            factor = re.sub(',[\t ][\t ]*', ' ', factor)
            factor = re.sub('[\t\r]*', '', factor)
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
        # Spectrum specific
        elif '(X++(Y..Y))' in testline:
            spectDataPos.append(currentPos + 1)
        elif '##END' in testline and len(spectDataPos) == 1:
            spectDataPos.append(currentPos - 1)
        elif '##XUNITS=' in testline:
            Xunit = line[line.index('=') + 1:]
            Xunit = re.sub('[ \t]*', '', Xunit)
        elif '##FIRSTX=' in testline:
            FirstX = float(line[line.index('=') + 1:])
        elif '##LASTX=' in testline:
            LastX = float(line[line.index('=') + 1:])
        elif '##YFACTOR=' in testline:
            YFactor = float(line[line.index('=') + 1:])
        elif '##NPOINTS=' in testline:
            NPoints = int(line[line.index('=') + 1:])
        currentPos += 1
    # Convert the data
    if 'NMR FID' in dataType:
        realDat = np.array([])
        if varForm[1] == 'ASDF':  # If DIFDUB form
            for line in data[realDataPos[0]:realDataPos[1] + 1]:
                realDat = np.append(realDat, convertDIFDUB(line))
        elif varForm[1] == 'AFFN':  # If regular list form
            for line in data[realDataPos[0]:realDataPos[1] + 1]:
                realDat = np.append(realDat, np.fromstring(line, sep=' ')[1:])
        realDat = realDat * factor[1]
        imagDat = np.array([])
        if varForm[2] == 'ASDF':
            for line in data[imagDataPos[0]:imagDataPos[1] + 1]:
                imagDat = np.append(imagDat, convertDIFDUB(line))
        elif varForm[1] == 'AFFN':
            for line in data[imagDataPos[0]:imagDataPos[1] + 1]:
                imagDat = np.append(imagDat, np.fromstring(line, sep=' ')[1:])
        imagDat = imagDat * factor[2]
        fullData = realDat - 1j * imagDat
        sw = 1.0 / ((last - first) / (nPoints - 1))
        masterData = sc.Spectrum(name, fullData, (10, filePath), [freq], [sw], [False])
    elif 'NMRSPECTRUM' in dataType:
        spectDat = np.array([])
        for line in data[spectDataPos[0]:spectDataPos[1] + 1]:
            spectDat = np.append(spectDat, convertDIFDUB(line))
        spectDat = np.flipud(spectDat) * YFactor
        if Xunit == 'HZ':
            sw = abs(FirstX - LastX)
            sw = sw + sw / NPoints
        elif Xunit == 'PPM':
            sw = abs(FirstX - LastX) * freq
            sw = sw + sw / NPoints
        masterData = sc.Spectrum(name, spectDat, (10, filePath), [freq], [sw], [True], ref=[None])
    return masterData

def saveASCIIFile(filePath, spectrum, axMult=1):
    axis = np.array([spectrum.xaxArray[-1] * axMult]).transpose()
    tmpData = spectrum.data.getHyperData(0)
    if tmpData.ndim == 1:  # create nx1 matrix if it is a 1d data set
        data = np.array([tmpData]).transpose()
    else:
        data = tmpData.transpose()
    splitdata = np.zeros([data.shape[0], data.shape[1] * 2])
    for line in np.arange(data.shape[1]):
        splitdata[:, line * 2] = np.real(data[:, line])
        splitdata[:, line * 2 + 1] = np.imag(data[:, line])
    data = np.concatenate((axis, splitdata), axis=1)
    np.savetxt(filePath, data, delimiter='\t')

def loadAscii(filePath, name, dataDimension, dataSpec, dataOrder, delimitor, swInp=0.0):
    freq = 0.0
    delimChar = ''
    if delimitor == 'Tab':
        delimChar = '\t'
    elif delimitor == 'Space':
        delimChar = ' '
    elif delimitor == 'Comma':
        delimChar = ','
    else:
        return
    matrix = np.genfromtxt(filePath, dtype=None, delimiter=delimChar)
    if dataOrder == 'XRI' or dataOrder == 'XR' or dataOrder == 'XI':
        if not dataSpec:
            sw = 1.0 / (matrix[1, 0] - matrix[0, 0])
        else:
            sw = abs(matrix[0, 0] - matrix[-1, 0]) / (matrix.shape[0] - 1) * matrix.shape[0]
    else:
        sw = swInp * 1000
    if dataDimension == 1:
        if dataOrder == 'XRI':
            data = matrix[:, 1] + 1j * matrix[:, 2]
        elif dataOrder == 'XR':
            data = matrix[:, 1]
        elif dataOrder == 'XI':
            data = 1j * matrix[:, 1]
        elif dataOrder == 'RI':
            data = matrix[:, 0] + 1j * matrix[:, 1]
        elif dataOrder == 'R':
            data = matrix
        masterData = sc.Spectrum(name, data, (11, filePath), [freq], [sw], [dataSpec], ref=[None])
    elif dataDimension == 2:
        if dataOrder == 'XRI':
            data = np.transpose(matrix[:, 1::2] + 1j * matrix[:, 2::2])
        elif dataOrder == 'XR':
            data = np.transpose(matrix[:, 1:])
        elif dataOrder == 'XI':
            data = 1j * np.transpose(matrix[:, 1:])
        elif dataOrder == 'RI':
            data = np.transpose(matrix[:, 0::2] + 1j * matrix[:, 1::2])
        elif dataOrder == 'RI':
            data = np.transpose(matrix)
        masterData = sc.Spectrum(name, data, (11, filePath), [freq, freq], [1, sw], [False, dataSpec], ref=[None, None])
    else:
        return
    masterData.addHistory("ASCII data loaded from " + filePath)
    return masterData

def loadMinispec(filePath, name):
    with open(filePath, 'r') as f:
        data = f.read().split('\n')
    dataType = int(data[1][data[1].index('=') + 1:])
    dataLimits = np.fromstring(data[2][data[2].index('=') + 1:], sep=',')
    dw = (dataLimits[1] - dataLimits[0]) / (dataLimits[2] - 1)
    if 'Time/ms' in data[3]:
        sw = 1.0 / dw * 1000
    elif 'Time/s' in data[3]:
        sw = 1.0 / dw
    totaldata = np.array([])
    if dataType == 1:  # real data?
        for line in data[7:]:
            if len(line) > 0:
                totaldata = np.append(totaldata, float(line))
    if dataType == 2:  # Complex data
        for line in data[7:]:
            if len(line) > 0:
                temp = np.fromstring(line, sep='\t')
                totaldata = np.append(totaldata, temp[0] + 1j * temp[1])
    masterData = sc.Spectrum(name, totaldata, (12, filePath), [0], [sw], [False])
    masterData.addHistory("Minispec data loaded from " + filePath)
    return masterData

def loadBrukerEPR(filePath, name=''):
    with open(filePath + '.par', mode='r') as f:
        textdata = [row.split() for row in f.read().replace('\r', '\n').split('\n')]
    for row in textdata:
        if len(row) < 2:
            continue
        if row[0] == 'ANZ':
            numOfPoints = int(row[1])
        elif row[0] == 'GSI':
            sweepWidth = float(row[1])
        elif row[0] == 'GST':
            leftX = float(row[1])
    with open(filePath + '.spc', mode='rb') as f:
        data = np.fromfile(f, np.float32, numOfPoints)
    masterData = sc.Spectrum(name, data, (13, filePath), [(sweepWidth + 2 * leftX) / 2], [sweepWidth], [True], ref=[0])
    masterData.addHistory("Bruker EPR data loaded from " + filePath)
    return masterData
