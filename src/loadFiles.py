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
import os

def LoadVarianFile(filePath, name=''):
    from struct import unpack
    if os.path.isfile(filePath):
        Dir = os.path.dirname(filePath)
    else:
        Dir = filePath
    freq = 0
    sw = 1
    sw1 = 1
    freq1 = 0  # intialize second dimension freqency as 0

    if os.path.exists(Dir + os.path.sep + 'procpar'):
        file = Dir + os.path.sep + 'procpar'
    elif os.path.exists(Dir + os.path.sep + '..' + os.path.sep + 'procpar'):
        file = Dir + os.path.sep + '..' + os.path.sep + 'procpar'
    else:
        file = None
    if file is not None:
        with open(file, 'r') as f:
            data = f.read().split('\n')
        
        indirectRef = ''
        for s in range(0, len(data)):
            if data[s].startswith('sfrq '):
                freq = float(data[s + 1].split()[1]) * 1e6
            elif data[s].startswith('refsource1 '):
                indirectRef = data[s + 1].split()[1][1:-1]
            elif data[s].startswith('sw '):
                sw = float(data[s + 1].split()[1])
            elif data[s].startswith('sw1 '):
                sw1 = float(data[s + 1].split()[1])
            elif data[s].startswith('reffrq '):
                reffreq = float(data[s + 1].split()[1]) * 1e6
            elif data[s].startswith('reffrq1 '):
                reffreq1 = float(data[s + 1].split()[1]) * 1e6
            elif data[s].startswith('phfid '):
                if int(data[s].split()[-2]): #if on
                    phfid = float(data[s + 1].split()[1]) 
                else:
                    phfid = 0
            elif data[s].startswith('rp '):
                rp = float(data[s + 1].split()[1])  
        if indirectRef:
            for s in range(0, len(data)): #Extra loop to get freq in indirect dimension
                if data[s].startswith(indirectRef + " "):
                    freq1 = float(data[s + 1].split()[1]) * 1e6
                
            

    if os.path.exists(Dir + os.path.sep + 'fid'):
        datafile = Dir + os.path.sep + 'fid'
        filePath = datafile
    elif os.path.exists(Dir + os.path.sep + 'data'):
        datafile = Dir + os.path.sep + 'data'
        filePath = datafile

    with open(datafile, "rb") as f:
        raw = np.fromfile(f, np.int32, 6)
        nblocks = unpack('>l', raw[0])[0]
        ntraces = unpack('>l', raw[1])[0]
        npoints = unpack('>l', raw[2])[0]
        ebytes = unpack('>l', raw[3])[0]
        tbytes = unpack('>l', raw[4])[0]
        bbytes = unpack('>l', raw[5])[0]
        raw = np.fromfile(f, np.int16, 2)
        vers_id = unpack('>h', raw[0])[0]
        status = unpack('>h', raw[1])[0]
        spec = bool(int(bin(status)[-2]))
        raw = np.fromfile(f, np.int32, 1)
        nbheaders = unpack('>l', raw[0])[0]
        SizeTD2 = npoints
        SizeTD1 = nblocks * ntraces
        a = []
        fid32 = int(bin(status)[-3])
        fidfloat = int(bin(status)[-4])
        hypercomplex = bool(bin(status)[-5])
        

        if not fid32 and fidfloat:  # only for `newest' format, use fast routine
            flipped = bool(bin(status)[-10])
            totalpoints = (ntraces * npoints + nbheaders**2 * 7)*nblocks
            raw = np.fromfile(f, np.float32, totalpoints)
            a = raw.newbyteorder('>f')
#                print(bin(status)[-10])
            if not spec or (spec and not hypercomplex):
                a = a.reshape(nblocks, int(totalpoints / nblocks))
                a = a[:, 7::]
                fid = a[:, ::2] - 1j * a[:, 1::2]
            else:
                fid = a[nbheaders*7::4] - 1j * a[nbheaders*7+1::4]
                fid = fid.reshape(int(SizeTD1/4),SizeTD2)
                if flipped:
                    fid = np.fliplr(fid)
 
        elif fid32 and not fidfloat:  # for VNMRJ 2 data
            totalpoints = (ntraces * npoints + nbheaders**2 * 7)*nblocks
            raw = np.fromfile(f, np.int32, totalpoints)
            a = raw.newbyteorder('>l')
            a = a.reshape(nblocks, int(totalpoints / nblocks))
            a = a[:, 7::]
            fid = a[:, ::2] - 1j * a[:, 1::2]
        else:  # use slow, but robust routine
            for iter1 in range(0, nblocks):
                b = []
                for iter2 in range(0, nbheaders):
                    raw = np.fromfile(f, np.int16, nbheaders * 14)
                if not fid32 and not fidfloat:
                    raw = np.fromfile(f, np.int16, ntraces * npoints)
                    for iter3 in raw:
                        b.append(unpack('>h', iter3)[0])
                elif fid32 and not fidfloat:
                    raw = np.fromfile(f, np.int32, ntraces * npoints)
                    for iter3 in raw:
                        b.append(unpack('>l', iter3)[0])
                else:
                    raw = np.fromfile(f, np.float32, ntraces * npoints)
                    for iter3 in raw:
                        b.append(unpack('>f', iter3)[0])
                b = np.array(b)
                if(len(b) != ntraces * npoints):
                    b.append(np.zeros(ntraces * npoints - len(b)))
                a.append(b)
            a = np.complex128(a)
            fid = a[:, ::2] - 1j * a[:, 1::2]
    
    
    fid = fid * np.exp((rp + phfid) / 180 * np.pi * 1j) #apply zero order phase
    if SizeTD1 is 1:
        fid = fid[0][:]
        if spec:  # flip if spectrum
            fid = np.flipud(fid)
        masterData = sc.Spectrum(name, fid, (0, filePath), [freq], [sw], [bool(int(spec))],ref = [reffreq])
    else:
        masterData = sc.Spectrum(name, fid, (0, filePath), [freq1, freq], [sw1, sw], [bool(int(spec))] * 2,ref = [reffreq1,reffreq])
    masterData.addHistory("Varian data loaded from " + filePath)
    return masterData
    
        
def LoadPipe(filePath, name=''):
        with open(filePath, 'r') as f:
            header = np.fromfile(f, np.float32, 512)

            NumberofPoints = int(header[99])
            data = np.fromfile(f, np.float32, NumberofPoints)
            if int(header[106]) == 0:  # if complex
                data = data + 1j * np.fromfile(f, np.float32, NumberofPoints)

            spec = int(header[220])  # 1 if ft, 0 if time
            freq = header[119] * 1e6
            sw = header[100]
            reference = header[101]  # frequency of last point in Hz

        sidefreq = -np.floor(NumberofPoints / 2) / NumberofPoints * sw  # freqeuency of last point on axis
        ref = sidefreq + freq - reference
        if spec == 1:
            data = np.flipud(data)

        masterData = sc.Spectrum(name, data, (8, filePath), [freq], [sw], [spec], ref=[ref])
        masterData.addHistory("NMRpipe data loaded from " + filePath)
        return masterData        


def LoadJEOLDelta(filePath, name=''):
    from struct import unpack
    with open(filePath, "rb") as f:
        file_identifier = f.read(8)
        endian = unpack('>B',f.read(1))[0]
        f.read(3) #placeholder to get rid of unused data
        data_dimension_number = unpack('>B',f.read(1))[0]
        data_dimension_exist = unpack('>B',f.read(1))[0]
        data_type = unpack('>B',f.read(1))[0]
        f.read(1) #placeholder to get rid of unused data
        translate = np.fromfile(f, '>B', 8)
        data_axis_type = np.fromfile(f, '>B', 8)
        data_units = np.fromfile(f, '>B', 16).reshape(8,2)#Reshape for later unit extraction
        title = f.read(76)
        f.read(52)
        #176
        data_points = np.fromfile(f, '>I', 8)
        data_offset_start = np.fromfile(f, '>I', 8)
        data_offset_stop = np.fromfile(f, '>I', 8)
        data_axis_start = np.fromfile(f, '>d', 8)
        data_axis_stop = np.fromfile(f, '>d', 8)
        #400
        f.read(664)
        base_freq = np.fromfile(f, '>d', 8)
        #1128
        zero_point = np.fromfile(f, '>d', 8)
        reverse = np.fromfile(f, '>B', 8)
        #1200
        f.read(12)
        param_start = np.fromfile(f, '>I', 1)[0]
        param_length = np.fromfile(f, '>I', 1)[0]
        #1220
        f.read(64)
        data_start = np.fromfile(f, '>I', 1)[0]
        data_length = np.fromfile(f, '>Q', 1)[0]
        #1296
        f.read(data_start-1296) #skip to data_start
        #start reading the data
        if endian:
            dataendian = '<d'
        else:
            dataendian = '>d'
        if data_dimension_number == 1:
            if data_axis_type[0] == 1: #if real
                data = np.fromfile(f, dataendian, data_points[0])
            elif data_axis_type[0] == 3: #if complex
                data = np.fromfile(f, dataendian, data_points[0]) - 1j*np.fromfile(f, dataendian, data_points[0])
                data = data[0:data_offset_stop[0]+1]
        elif data_dimension_number == 2:
            if data_axis_type[0] == 4: #if real-complex (no hypercomplex)
                Step = 4 #Step size of block
                pointsD2 = data_points[0]
                pointsD1 = data_points[1]
                datalength = pointsD2 * pointsD1
                datareal = np.fromfile(f, dataendian, datalength)
                dataimag = np.fromfile(f, dataendian, datalength)
                data = datareal - 1j*dataimag
                data = np.reshape(data,[pointsD1/Step,datalength/pointsD1/Step,Step,Step])
                data = np.concatenate(np.concatenate(data,1),1)
                data = data[0:data_offset_stop[1]+1,0:data_offset_stop[0]+1] #cut back to real size
            if data_axis_type[0] == 3: #if complex (i.e. hypercomplex)
                Step = 32 #Step size of block
                pointsD2 = data_points[0]
                pointsD1 = data_points[1]
                datalength = pointsD2 * pointsD1
                datareal1 = np.fromfile(f, dataendian, datalength)
                dataimag1 = np.fromfile(f, dataendian, datalength)
                datareal2 = np.fromfile(f, dataendian, datalength)
                dataimag2 = np.fromfile(f, dataendian, datalength)
                data1 = datareal1 - 1j*dataimag1
                data1 = np.reshape(data1,[pointsD1/Step,datalength/pointsD1/Step,Step,Step])
                data1 = np.concatenate(np.concatenate(data1,1),1)
                data2 = datareal2 - 1j*dataimag2
                data2 = np.reshape(data2,[pointsD1/Step,datalength/pointsD1/Step,Step,Step])
                data2 = np.concatenate(np.concatenate(data2,1),1)
                data = np.zeros([data2.shape[0]*2,data2.shape[1]],dtype=complex)
                if reverse[0] == 0:
                    data[::2,:] = data1 #Interleave both types in D1
                    data[1::2,:] = data2
                else:
                    data[::2,:] = data2 #Interleave both types in D1
                    data[1::2,:] = data1           
                data = data[0:(data_offset_stop[1] + 1)*2,0:data_offset_stop[0]+1] #Cut back to real size
                if reverse[1] == 1:
                    data = np.conjugate(data)
        else:
                return 'ND error'
	#unitExp = data_units[0][0] &15#unit_exp: if (unitExp > 7) >unitExp -= 16;
	#scaleType = (data_units[0][1]>>4) &15 #scaleType: if (scaleType > 7) scaleType -= 16;
    sw = []
    spec = []
    ref = []
    freq = []
    for axisNum in reversed(range(data_dimension_number)):
        freq.append(base_freq[0]*1e6)
        axisType = data_units[axisNum][1] #Sec = 28, Hz = 13, PPM = 26
        axisScale = data_units[axisNum][0]
        if axisType == 28: #Sec
            scale = (axisScale >> 4) & 15
            if scale > 7:
                scale = scale -16
            dw = (data_axis_stop[axisNum]-data_axis_start[axisNum])/(data_offset_stop[axisNum]+1 -1)*10**(-scale*3)#minus one to give same axis as spectrum???
            #scale for SI prefix
            sw.append(1.0 / dw)
            spec.append(False)
            sidefreq = -np.floor((data_offset_stop[axisNum]+1) / 2) / data_offset_stop[axisNum]+1 * sw[-1]  # frequency of last point on axis
            ref.append(base_freq[axisNum]*1e6)
        if axisType == 13: #Hz
            sw.append(np.abs(data_axis_start[axisNum]-data_axis_stop[0]))
            spec.append(True)
            if data_dimension_number == 1:
                data = np.flipud(data)
            sidefreq = -np.floor(data_points[axisNum] / 2) / (data_offset_stop[axisNum]+1) * sw[-1]  # frequency of last point on axis
            ref.append(sidefreq + base_freq[axisNum]*1e6 - data_axis_stop[axisNum])
        if axisType == 26: #ppm
            sw.append(np.abs(data_axis_start[axisNum]-data_axis_stop[axisNum])*base_freq[axisNum])
            spec.append(True)
            if data_dimension_number == 1:
                data = np.flipud(data)
            sidefreq = -np.floor((data_offset_stop[axisNum]+1) / 2) / (data_offset_stop[axisNum]+1) * sw[-1]  # frequency of last point on axis
            ref.append(sidefreq + base_freq[axisNum]*1e6 - data_axis_stop[axisNum]*base_freq[axisNum])
    masterData = sc.Spectrum(name, data, (9, filePath), freq, sw, spec,ref=ref)
    masterData.addHistory("JEOL Delta data loaded from " + filePath)
    return masterData

def loadJSONFile(filePath, name=''):
    import json
    with open(filePath, 'r') as inputfile:
        struct = json.load(inputfile)
    data = np.array(struct['dataReal']) + 1j * np.array(struct['dataImag'])
    ref = np.where(np.isnan(struct['ref']), None, struct['ref'])
    if 'history' in struct.keys():
        history = struct['history']
    else:
        history = None
    xaxA = []
    for i in struct['xaxArray']:
        xaxA.append(np.array(i))
    masterData = sc.Spectrum(name,
                             data,
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
    
def loadMatlabFile(filePath, name=''):
    import scipy.io
    with open(filePath, 'rb') as inputfile:  # read first several bytes the check .mat version
        teststring = inputfile.read(13)
    version = float(teststring.decode("utf-8")[7:10])  # extract version from the binary array
    if version < 7.3:  # all versions below 7.3 are supported
        matlabStruct = scipy.io.loadmat(filePath)
        var = [k for k in matlabStruct.keys() if not k.startswith('__')][0]
        mat = matlabStruct[var]
        if mat['dim'] == 1:
            data = mat['data'][0, 0][0]
            xaxA = [k[0] for k in (mat['xaxArray'][0])]
        else:
            data = mat['data'][0, 0]
            if all(x == data.shape[0] for x in data.shape):
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
                                 data,
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
        DataGroup = Groups[0]  # get the groupo name
        mat = f[DataGroup]
        if np.array(mat['dim'])[0][0] == 1:
            xaxA = list([np.array(mat['xaxArray'])[:, 0]])
            data = np.array(mat['data'])
            data = (data['real'] + data['imag'] * 1j)[:, 0]  # split and use real and imag part
        else:
            data = np.transpose(np.array(mat['data']))
            data = data['real'] + data['imag'] * 1j
            if all(x == data.shape[0] for x in data.shape):
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
                                 data,
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


def LoadBrukerTopspin(filePath, name=''):
    if os.path.isfile(filePath):
        Dir = os.path.dirname(filePath)
    else:
        Dir = filePath
    if os.path.exists(Dir + os.path.sep + 'acqus'):
        with open(Dir + os.path.sep + 'acqus', 'r') as f:
            data = f.read().split('\n')
        for s in range(0, len(data)):
            if data[s].startswith('##$TD='):
                sizeTD2 = int(data[s][6:])
            if data[s].startswith('##$SFO1='):
                freq2 = float(data[s][8:]) * 1e6
            if data[s].startswith('##$SW_h='):
                SW2 = float(data[s][8:])
            if data[s].startswith('##$BYTORDA='):
                ByteOrder = int(data[s][11:])
    sizeTD1 = 1
    if os.path.exists(Dir + os.path.sep + 'acqu2s'):
        with open(Dir + os.path.sep + 'acqu2s', 'r') as f:
            data2 = f.read().split('\n')
        SW1 = 10e3 #pre initialize
        for s in range(0, len(data2)):
            if data2[s].startswith('##$TD='):
                sizeTD1 = int(data2[s][6:])
            if data2[s].startswith('##$SFO1='):
                freq1 = float(data2[s][8:]) * 1e6
            if data2[s].startswith('##$SW_h='):
                SW1 = float(data2[s][8:])
    if os.path.exists(Dir + os.path.sep + 'fid'):
        filePath = Dir + os.path.sep + 'fid'
        with open(Dir + os.path.sep + 'fid', "rb") as f:
            raw = np.fromfile(f, np.int32, sizeTD1* sizeTD2)
    elif os.path.exists(Dir + os.path.sep + 'ser'):
        filePath = Dir + os.path.sep + 'ser'
        with open(Dir + os.path.sep + 'ser', "rb") as f:
            raw = np.fromfile(f, np.int32, sizeTD1 * int(np.ceil(sizeTD2 / 256))*256) #Always load full 1024 byte blocks (256 data points)
    if ByteOrder:
        RawInt = raw.newbyteorder('b')
    else:
        RawInt = raw.newbyteorder('l')
    ComplexData = np.array(RawInt[0:len(RawInt):2]) + 1j * np.array(RawInt[1:len(RawInt):2])
    spec = [False]
    if sizeTD1 is 1:
        masterData = sc.Spectrum(name, ComplexData, (1, filePath), [freq2], [SW2], spec)
    else:
        ComplexData = ComplexData.reshape(sizeTD1, int(np.ceil(sizeTD2 / 256) * 256 / 2))
        ComplexData = ComplexData[:,0:int(sizeTD2/2)] #Cut off placeholder data
        masterData = sc.Spectrum(name, ComplexData, (1, filePath), [freq1, freq2], [SW1, SW2], spec * 2)
    masterData.addHistory("Bruker data loaded from " + filePath)
    return masterData

def LoadBrukerSpectrum(filePath, name=''):
    if os.path.isfile(filePath):
        Dir = os.path.dirname(filePath)
    else:
        Dir = filePath
    if os.path.exists(Dir + os.path.sep + 'procs'):  # Get D2 parameters
        with open(Dir + os.path.sep + 'procs', 'r') as f:
            data = f.read().split('\n')
        for s in range(0, len(data)):
            if data[s].startswith('##$SI='):
                sizeTD2 = int(re.findall("\#\#\$SI= (.*.)", data[s])[0])
#                if data[s].startswith('##$XDIM='):
#                    blockingD2 = int(data[s][8:])
            if data[s].startswith('##$BYTORDP='):
                ByteOrder = int(data[s][11:])
            if data[s].startswith('##$SW_p='):
                SW2 = float(data[s][8:])
            if data[s].startswith('##$SF='):
                Ref2 = float(data[s][6:])*1e6 
                
    freq2 = 0
    if os.path.exists(Dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'acqus'):  # Get D2 parameters from fid directory, if available
        with open(Dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'acqus', 'r') as f:
            data = f.read().split('\n')
        for s in range(0, len(data)):
            if data[s].startswith('##$SFO1='):
                freq2 = float(data[s][8:]) * 1e6
    sizeTD1 = 1
    if os.path.exists(Dir + os.path.sep + 'proc2s'):  # Get D1 parameters
        with open(Dir + os.path.sep + 'proc2s', 'r') as f:
            data2 = f.read().split('\n')
        for s in range(0, len(data2)):
            if data2[s].startswith('##$SI='):
                sizeTD1 = int(data2[s][6:])
#                if data2[s].startswith('##$XDIM='):
#                    blockingD1 = int(data[s][8:])
            if data2[s].startswith('##$SW_p='):
                SW1 = float(data2[s][8:])
            if data2[s].startswith('##$SF='):
                Ref1 = float(data2[s][6:])*1e6
    freq1 = 0
    if os.path.exists(Dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'acqu2s'):  # Get D1 parameters from fid directory, if available
        with open(Dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'acqu2s', 'r') as f:
            data = f.read().split('\n')
        for s in range(0, len(data)):
            if data[s].startswith('##$SFO1='):
                freq1 = float(data[s][8:]) * 1e6
    if os.path.exists(Dir + os.path.sep + '1r'):  # Get D2 data
        filePath = Dir + os.path.sep + '1r'
        with open(Dir + os.path.sep + '1r', "rb") as f:
            RawReal = np.fromfile(f, np.int32, sizeTD1 * sizeTD2)
        RawImag = np.zeros([sizeTD1 * sizeTD2])
        if os.path.exists(Dir + os.path.sep + '1i'):
            with open(Dir + os.path.sep + '1i', "rb") as f:
                RawImag = np.fromfile(f, np.int32, sizeTD1 * sizeTD2)
    elif os.path.exists(Dir + os.path.sep + '2rr'):  # Get D1 data
        filePath = Dir + os.path.sep + '2rr'
        with open(Dir + os.path.sep + '2rr', "rb") as f:
            RawReal = np.fromfile(f, np.int32, sizeTD1 * sizeTD2)
        RawImag = np.zeros([sizeTD1 * sizeTD2])
        if os.path.exists(Dir + os.path.sep + '2ir'):  # If hypercomplex
            with open(Dir + os.path.sep + '2ir', "rb") as f:
                RawImag = np.fromfile(f, np.int32, sizeTD1 * sizeTD2)
        elif os.path.exists(Dir + os.path.sep + '2ii'):
            with open(Dir + os.path.sep + '2ii', "rb") as f:
                RawImag = np.fromfile(f, np.int32, sizeTD1 * sizeTD2)
    if ByteOrder:
        RawReal = RawReal.newbyteorder('b')
        RawImag = RawImag.newbyteorder('b')
    else:
        RawReal = RawReal.newbyteorder('l')
        RawImag = RawImag.newbyteorder('l')
    Data = np.flipud(RawReal) - 1j * np.flipud(RawImag)
    spec = [True]
    if sizeTD1 is 1:
        masterData = sc.Spectrum(name, Data, (7, filePath), [freq2], [SW2], spec,ref=[Ref2])
    else:
        Data = Data.reshape(sizeTD1, sizeTD2)
        masterData = sc.Spectrum(name, Data, (7, filePath), [freq1, freq2], [SW1, SW2], spec * 2,ref=[Ref1,Ref2])
    masterData.addHistory("Bruker spectrum data loaded from " + filePath)
    return masterData


def LoadChemFile(filePath, name=''):
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

def LoadMagritek(filePath, name='',realPath=''):
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
    
    #initialize 2D values to some dummy value
    sizeTD1 = 0
    sw1 = 50e3
    lastfreq1 = None
    ref1 = None
    
    #Start pars extraction
    H = dict(line.strip().split('=') for line in open(Dir + os.path.sep + 'acqu.par', 'r'))
    for key in H.keys():
        if key.startswith("bandwidth "):
            sw = float(H[key]) * 1000
        elif key.startswith("nrPnts "):
            sizeTD2 = int(H[key])
        elif key.startswith("b1Freq "):
            freq = float(H[key])*1e6
        elif key.startswith("lowestFrequency "):
            lastfreq = float(H[key])
        elif key.startswith("nrSteps "):   
            sizeTD1 = int(H[key])
        elif key.startswith("bandwidth2 "):
            sw1 = float(H[key]) *1000
        elif key.startswith("lowestFrequency2 "):
            lastfreq1 = float(H[key])
    
    sidefreq = -np.floor(sizeTD2 / 2) / sizeTD2 * sw  # freqeuency of last point on axis
    ref = sidefreq + freq - lastfreq
    if len(Files2D) == 1:
        File = Files2D[0]
        
        if lastfreq1 != None:   
            sidefreq1 = -np.floor(sizeTD1 / 2) / sizeTD1 * sw1  # freqeuency of last point on axis
            ref1 = sidefreq1 + freq - lastfreq1

        
        with open(Dir + os.path.sep + File, 'rb') as f:
            raw = np.fromfile(f, np.float32)
        Data = raw[-2 * sizeTD2 * sizeTD1::]
        ComplexData = Data[0:Data.shape[0]:2] - 1j * Data[1:Data.shape[0]:2]
        ComplexData = ComplexData.reshape((sizeTD1, sizeTD2))
        masterData = sc.Spectrum(name, ComplexData, (3, rememberPath), [freq] * 2, [sw1,sw], [False] * 2,ref=[ref1,ref])
    elif len(Files1D) != 0:
        File = 'data.1d'
        with open(Dir + os.path.sep + File, 'rb') as f:
            raw = np.fromfile(f, np.float32)
        Data = raw[-2 * sizeTD2::]
        ComplexData = Data[0:Data.shape[0]:2] - 1j * Data[1:Data.shape[0]:2]
        masterData = sc.Spectrum(name, ComplexData, (3, rememberPath), [freq], [sw], [False],ref=[ref])
    masterData.addHistory("Magritek data loaded from " + rememberPath)
    return masterData
 


def LoadSimpsonFile(filePath, name=''):
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
        charst =  np.fromstring(chardata, dtype=np.uint8)
        charst = charst.reshape(nquads,4) - BASE
        FIRST = lambda f, x: ((x) & ~(~0 << f))
        LAST = lambda f, x: ((x) & (~0 << (8 - f)))
        
        first = FIRST(6, charst[:,0]) | LAST(2, charst[:,1] << 2)
        second  = FIRST(4, charst[:,1]) | LAST(4, charst[:,2] << 2)
        third = FIRST(2, charst[:,2]) | LAST(6, charst[:,3] << 2)
        
        Bytes = np.ravel(np.transpose(np.array([first,second,third]))).astype('int64')
        
        # convert every 4 'bytes' to a float
        num_points, num_pad = divmod(len(Bytes), 4)
        Bytes = np.array(Bytes)
        Bytes=Bytes[:-num_pad]
        Bytes=Bytes.reshape(num_points,4)
        mantissa = ((Bytes[:,2] % 128) << 16) + (Bytes[:,1] << 8) + Bytes[:,0]
        exponent = (Bytes[:,3] % 128) * 2 + (Bytes[:,2] >= 128) * 1
        negative = Bytes[:,3] >= 128
        e = exponent - 127
        m = np.abs(mantissa) / np.float64(1 << 23)
        data = np.float32((-1)**negative*np.ldexp(m,e))
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


def LoadAscii(filePath, name, dataDimension, dataSpec, dataOrder, delimitor, swInp = 0.0):
    
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
    
    matrix = np.genfromtxt(filePath,dtype=None, delimiter = delimChar)
    if dataOrder == 'XRI' or dataOrder == 'XR' or dataOrder == 'XI':
        if dataSpec == False:
            sw = 1.0 / (matrix[1,0] - matrix[0,0])
        else:
            sw = abs(matrix[0,0] - matrix[-1,0])/(matrix.shape[0] - 1) * matrix.shape[0]
    else:
        sw = swInp * 1000
    
    if dataDimension == 1:
        if dataOrder == 'XRI':
            data = matrix[:,1] + 1j * matrix[:,2]
        elif dataOrder == 'XR':
            data = matrix[:,1]
        elif dataOrder == 'XI':
            data = 1j * matrix[:,1]
        elif dataOrder == 'RI':
            data = matrix[:,0] + 1j * matrix[:,1]
        elif dataOrder == 'R':
            data = matrix
        masterData = sc.Spectrum(name, data, (11, filePath), [freq], [sw], [dataSpec], ref = [None])
    elif dataDimension == 2:
        if dataOrder == 'XRI':
            data = np.transpose(matrix[:,1::2] + 1j * matrix[:,2::2])
        elif dataOrder == 'XR':
           data = np.transpose(matrix[:,1:]) 
        elif dataOrder == 'XI':
           data = 1j * np.transpose(matrix[:,1:]) 
        elif dataOrder == 'RI':  
            data = np.transpose(matrix[:,0::2] + 1j * matrix[:,1::2])
        elif dataOrder == 'RI':  
            data = np.transpose(matrix)
        masterData = sc.Spectrum(name, data, (11, filePath), [freq,freq], [1,sw], [False,dataSpec], ref = [None,None])
    else:
        return
            
    masterData.addHistory("ASCII data loaded from " + filePath)    
    return masterData
        
        
        
def LoadMinispec(filePath,name):
    with open(filePath, 'r') as f:
        data = f.read().split('\n')    
       
    dataType = int(data[1][data[1].index('=')+1:])
    dataLimits = np.fromstring(data[2][data[2].index('=')+1:],sep = ',')
    dw = (dataLimits[1] - dataLimits[0]) / ( dataLimits[2] - 1)
    if 'Time/ms' in data[3]:    
        sw = 1.0/dw * 1000
    elif 'Time/s' in data[3]: 
        sw = 1.0/dw
    totaldata = np.array([])
    if dataType == 1: #real data?
        for line in data[7:]:
            if len(line) > 0:
                totaldata = np.append(totaldata,float(line))
    if dataType == 2: #Complex data
        for line in data[7:]:
            if len(line) > 0:
                temp = np.fromstring(line,sep = '\t')
                totaldata = np.append(totaldata,temp[0] + 1j * temp[1])
    
    
    
    masterData = sc.Spectrum(name, totaldata, (12, filePath), [0], [sw], [False])
    return masterData       
        




def LoadBrukerEPR(filePath,name=''):
    from struct import unpack, calcsize
    if os.path.isfile(filePath):
        Dir = os.path.dirname(filePath)
    else:
        Dir = filePath

    with open(filePath + '.spc',mode='rb') as f:
        content = f.read()

    x = 'f'
    xsize = calcsize(x)
    #print(xsize)
    data = []
    for i in range(int(len(content)/xsize)):
        data.append(unpack(x,content[i*xsize:i*xsize+xsize])[0])

    data =np.array(data) 

    with open(filePath + '.par',mode='r')     as f:
        textdata = [row.split() for row in f.readlines()]

    for row in textdata:
        if row[0]=='ANZ':
            numOfPoints = int(row[1])
        elif row[0]=='GSI':
            sweepWidth = float(row[1])
        elif row[0]=='GST':
            leftX = float(row[1])


    xdata = np.arange(leftX,leftX+sweepWidth,sweepWidth/numOfPoints)        


    masterData = sc.Spectrum(name, data, (0, filePath), [(sweepWidth + 2 * leftX)/2], [sweepWidth], [True], ref = [0])
    masterData.addHistory("Bruker EPR data loaded from " + filePath)
    return masterData
    

