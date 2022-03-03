# -*- coding: cp1252 -*-
from pylab import *
from matplotlib.pyplot import *
from numpy import *
import os.path as osp

DIR_NAME = osp.dirname(osp.abspath(__file__))

def get_freqs_data(filename, Trace=1, unit='MHz', dirname=None):
    if unit == 'Hz':
        fact = 1
    elif unit == 'kHz':
        fact = 1e-3
    elif unit == 'MHz':
        fact = 1e-6
    else:
        raise ValueError('Wrong unit specified : '+str(unit))

    if dirname is None:
        dirname = DIR_NAME

    if not filename[-4:] == '.csv':
        filename += '.csv'
        
    with open(dirname + '\\' + filename, 'r') as ff:
        lines = ff.readlines()

    num_trace = 0
    II = 0
    while not num_trace == Trace:
        for ii in range(II, len(lines)):
            line = lines[ii]
            if line[-1] == '\n':
                line = line[:-1]
            """
            if not '[Trace]' in line.split(','):
                pass
            else:
                if '[Trace]' in line.split(','):
                    #print("Trace found, line", ii)
                    break
            """
            if '[Trace]' in line.split(','):
                if not 'Spectrogram' in lines[ii+1]:
                    print("Trace found, line", ii)
                    break
        else:
            raise ValueError('Trace not found')
        num_trace = int(lines[ii+1].split(',')[0].split('Trace ')[1])
        II = ii+1
        #print('Num trace :', num_trace, '- II :', II)
    numpoints = int(lines[ii+2].split(',')[1])
    freq_start = float(lines[ii+3].split(',')[1])
    freq_stop = float(lines[ii+4].split(',')[1])
    freqs = linspace(freq_start, freq_stop, numpoints)

    data = []
    for nn in range(numpoints):
        data.append(float(lines[ii+5+nn].split('\n')[0].split(',')[0]))

    #print('len(data) :', len(data))
    return array(freqs), array(data)

def get_times_data(filename, Trace=1, unit='s', dirname=None, linear=True):
    if unit == 's':
        fact = 1
    elif unit == 'ms':
        fact = 1e-3
    else:
        raise ValueError('Wrong unit specified : '+str(unit))

    if dirname is None:
        dirname = DIR_NAME

    if not filename[-4:] == '.csv':
        filename += '.csv'
        
    with open(dirname + '\\' + filename, 'r') as ff:
        lines = ff.readlines()

    num_trace = 0
    II = 0
    while not num_trace == Trace:
        for ii in range(II, len(lines)):
            line = lines[ii]
            if not line == '[Trace]\n':
                pass
            else:
                if line == '[Trace]\n':
                    #print("Trace found, line", ii)
                    break
        else:
            raise ValueError('Trace not found')
        num_trace = int(lines[ii+1].split(',')[0].split('Trace ')[1])
        II = ii+1
        #print('Num trace :', num_trace, '- II :', II)
    numpoints = int(lines[ii+2].split(',')[1])
    times_start = float(lines[ii+3].split(',')[1])
    times_stop = float(lines[ii+4].split(',')[1])
    times = linspace(times_start, times_stop, numpoints)

    data = []
    for nn in range(numpoints):
        data.append(float(lines[ii+5+nn].split('\n')[0]))

    if linear:
        data = [10**(x/20) for x in data]
    #print('len(data) :', len(data))
    return array(times), array(data)

def get_freqs_data_2(filename, unit='MHz', dirname=None):
    if unit == 'Hz':
        fact = 1
    elif unit == 'kHz':
        fact = 1e-3
    elif unit == 'MHz':
        fact = 1e-6
    else:
        raise ValueError('Wrong unit specified : '+str(unit))

    if dirname is None:
        dirname = DIR_NAME

    if not filename[-4:] == '.csv':
        filename += '.csv'
        
    with open(dirname + '\\' + filename, 'r') as ff:
        lines = ff.readlines()

    for ii, line in enumerate(lines):
        if line.split(',')[0].split('\n')[0] == '[Traces]':
            break
    else:
        raise ValueError('Traces not found')

    print('Traces found, starting line', ii+1)
    numpoints = []
    line_numpoints = lines[ii+3]
    for nn, x in enumerate(line_numpoints.split(',')):
        if x == 'NumberPoints':
            numpoints.append(int(line_numpoints.split(',')[nn+1]))
    print('Numpoints :', numpoints)

    xstarts = []
    line_xstarts = lines[ii+4]
    for nn, x in enumerate(line_xstarts.split(',')):
        if x == 'XStart':
            xstarts.append(float(line_xstarts.split(',')[nn+1]))
    print('xstarts :', xstarts)

    xstops = []
    line_xstops = lines[ii+5]
    for nn, x in enumerate(line_xstops.split(',')):
        if x == 'XStop':
            xstops.append(float(line_xstops.split(',')[nn+1]))
    print('xstops :', xstops)
    
    start = ii+6
    freqs = []
    traces = [[], [], []]
    for ii, line in enumerate(lines[start:]):
        tu = line.split(',')
        count = int(tu[0])
        freqs.append(float(tu[1]))
        jj = 0
        for x in tu[2:]:
            if not x in ('', '\n'):
                traces[jj].append(float(x))
                jj += 1

    while [] in traces:
        traces.remove([])

    if not product([x == count for x in numpoints]) or \
       not product([len(ll) == count for ll in traces]):
        print("""\n!!! Warning !!!
Num Points not corresponding
Count : """+str(count)+"""
Len(traces) : """+str([len(ll) for ll in traces])+"""
Numpoints : """+str(numpoints)+"""\n""")

    if not product([x == freqs[0] for x in xstarts]) or \
       not product([x == freqs[-1] for x in xstops]):
        print("""\n!!!Warning !!!
Freqs not corresponding\n""")

    return array(freqs), array(traces)

def get_IQ(filename, dirname):
    if not '.' in filename or not filename.split('.')[-1] == 'csv':
        filename += '.csv'
    with open(dirname+'/'+filename, 'r') as ff:
        lines = ff.readlines()

    I_found = False
    Q_found = False
    Is = []
    Qs = []
    for ii, line in enumerate(lines):
        if line.split(',')[0] == 'I Trace':
            I_found = True
            idx_I = ii
            print('I found line', ii)
        elif line.split(',')[0] == 'Q Trace':
            Q_found = True
            idx_Q = ii
            print('Q found line', ii)

        if I_found and Q_found and ii > idx_Q+3:
            Qs.append(float(line.split('\n')[0]))
        elif I_found and not Q_found and ii > idx_I+3 and not '[Trace]' in line:
            Is.append(float(line.split('\n')[0]))

    numpoints_I = int(lines[idx_I+1].split(',')[1])
    xstart_I = float(lines[idx_I+2].split(',')[1])
    xstop_I = float(lines[idx_I+3].split(',')[1])
    time_unit_I = lines[idx_I+2].split(',')[2].split('\n')[0]
    xI = linspace(xstart_I, xstop_I, numpoints_I)
    print('numpoints_I:', numpoints_I)
    print('xstart_I =', xstart_I, time_unit_I)
    print('xstop_I =', xstop_I, time_unit_I)

    numpoints_Q = int(lines[idx_Q+1].split(',')[1])
    xstart_Q = float(lines[idx_Q+2].split(',')[1])
    xstop_Q = float(lines[idx_Q+3].split(',')[1])
    time_unit_Q = lines[idx_Q+2].split(',')[2].split('\n')[0]
    xQ = linspace(xstart_Q, xstop_Q, numpoints_Q)
    print('numpoints_Q:', numpoints_I)
    print('xstart_Q =', xstart_I, time_unit_Q)
    print('xstop_Q =', xstop_I, time_unit_Q)

    return array(xI), array(Is), array(xQ), array(Qs)

def get_IQ_2(filename, dirname, start_idx=10):
    if not '.' in filename or not filename.split('.')[-1] == 'csv':
        filename += '.csv'
    with open(dirname+'/'+filename, 'r') as ff:
        lines = ff.readlines()

    for ii in range(start_idx):
        line = lines[ii]
        if 'SamplingFrequency' in line:
            freq = int(line.split(',')[1].split('\n')[0])
            print('SamplingFrequency,', freq)
        elif 'NumberSamples' in line:
            nb = int(line.split(',')[1].split('\n')[0])
            print('NumberSamples', nb)

    Is = []
    Qs = []
    for ii in range(start_idx, len(lines)):
        line = lines[ii]
        ll = line.split(',')

        if len(ll) == 4:
            I = float(ll[0]+'.'+ll[1])
            Q = float(ll[2]+'.'+ll[3].split('\n')[0])
        elif len(ll) == 2:
            I = float(ll[0])
            Q = float(ll[1].split('\n')[0])
        elif 'E' in ll[0]:
            I = float(ll[0])
            Q = float(ll[1]+'.'+ll[2].split('\n')[0])
        elif 'E' in ll[1] or ll[1][:2] == '00':
            I = float(ll[0]+'.'+ll[1])
            Q = float(ll[2].split('\n')[0])
        elif 'E' in ll[2]:
            I = float(ll[0])
            Q = float(ll[1]+'.'+ll[2].split('\n')[0])

        Is.append(I)
        Qs.append(Q)

    return array(Is), array(Qs), freq, nb

def get_spectrogram(filename, dirname):
    if not '.' in filename or not filename.split('.')[-1] == 'csv':
        filename += '.csv'
    with open(dirname+'/'+filename, 'r') as ff:
        lines = ff.readlines()

    spectrogram = []
    start = lines.index('[Traces]\n') + 1
    while start < len(lines):
        if not lines[start] == '[SpectrogramTrace]\n':
            print('Warning line '+str(start)+' : not Spectrogram Trace !')
            break
        else:
            ll = []
            freqs = []

            numline = int(lines[start+1].split(',')[0].split('SpectrogramLine')[1])
            numpoints = int(lines[start+2].split(',')[1].split('\n')[0])
            xstart = int(lines[start+3].split(',')[1])
            xstop = int(lines[start+4].split(',')[1])
            timestamp = float(lines[start+5].split(',')[1]+'.'+lines[start+5].split(',')[2])
            print('SpectrogramLine %i found line %i : numpoints = %i ; xstart = %f ; xstop = %f ; timestamp = %f'%(numline, start, numpoints, xstart, xstop, timestamp))

            for ii in range(numpoints):
                line = lines[start + 6 + ii]
                if len(line.split(',')) == 3:
                    val = float(line.split(',')[0]+'.'+line.split(',')[1])
                    freq = float(line.split(',')[2].split('\n')[0])
                else:
                    val = float(line.split(',')[0])
                    freq = float(line.split(',')[1].split('\n')[0])

                ll.append(val)
                freqs.append(freq)

            start += 6 + numpoints
            spectrogram.append(ll)

    return array(spectrogram), array(freqs)
