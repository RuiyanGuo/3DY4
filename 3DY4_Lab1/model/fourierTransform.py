#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
import numpy as np
import cmath, math
import random
import sys
from scipy import signal

def plotSpectrum(x, Fs, type = 'FFT'):

    n = len(x)             # length of the signal
    df = Fs/n              # frequency increment (width of freq bin)

    # compute Fourier transform, its magnitude and normalize it before plotting
    if type == 'FFT':
        Xfreq = np.fft.fft(x)
    elif type == 'DFT':
        Xfreq = myDFT(x)
    XMag = abs(Xfreq)/n

    # Note: because x is real, we keep only the positive half of the spectrum
    # Note also: half of the energy is in the negative half (not plotted)
    XMag = XMag[0:int(n/2)]

    # freq vector up to Nyquist freq (half of the sample rate)
    freq = np.arange(0, Fs/2, df)

    fig, ax = plt.subplots()
    ax.plot(freq, XMag)
    ax.set(xlabel='Frequency (Hz)', ylabel='Magnitude',
        title='Frequency domain plot')
    # fig.savefig("freq.png")
    plt.show()

def plotTime(x, time):

    fig, ax = plt.subplots()
    ax.plot(time, x)
    ax.set(xlabel='Time (sec)', ylabel='Amplitude',
            title='Time domain plot')
    # fig.savefig("time.png")
    plt.show()

def plotTime(x, time, type='IDFT'):
    x = myInDFT(myDFT(x))
    fig, ax = plt.subplots()
    ax.plot(time, x)
    ax.set(xlabel='Time (sec)', ylabel='Amplitude',
            title='Time domain plot')
    # fig.savefig("time.png")
    plt.show()

def generateSin(Fs, interval, frequency = 7.0, amplitude = 5.0, phase = 0.0):

    dt = 1.0/Fs                          # sampling period (increment in time)
    time = np.arange(0, interval, dt)    # time vector over interval

    # generate the sin signal
    x = amplitude*np.sin(2*math.pi*frequency*time+phase)

    return time, x

def cli_error_msg():

    # error message to provide the correct command line interface (CLI) arguments
    print('Valid arguments:')
    print('\trc:  reference code')
    print('\til1: in-lab 1')
    print('\til2: in-lab 2')
    print('\til3: in-lab 3')
    print('\tth:  take-home')
    sys.exit()

def myDFT(x):
    #Code for DFT
    N = len(x)
    n = np.arange(N)
    k = n.reshape((N, 1))
    e = np.exp(-2j * cmath.pi * k*n / N)
    x_dft = np.dot(e, x)

    return x_dft

def myInDFT(x):
    #Code for DFT
    Ni = len(x)
    ki = np.arange(Ni)
    ni = ki.reshape((Ni, 1))
    ei = 1/Ni * np.exp(2j * cmath.pi * ki*ni / Ni)

    x_indft = np.dot(ei, x)

    return x_indft

if __name__ == "__main__":

    if len(sys.argv[0:]) != 2:
        cli_error_msg()

    Fs = 100.0          # sampling rate
    interval = 1.0      # set up to one full second

    if (sys.argv[1] == 'rc'): # runs the reference code (rc)

        print('Reference code for the Fourier transform')

        # generate the user-defined sin function
        time, x = generateSin(Fs, interval)
        # plot the signal in time domain
        plotTime(x, time)
        # plot the signal in frequency domain
        plotSpectrum(x, Fs, type = 'FFT')

    elif (sys.argv[1] == 'il1'):

        print('In-lab experiment 1 for the Fourier transform')

        # compute the spectrum with your own DFT
        # you can use cmath.exp() for complex exponentials
        # plotSpectrum(x, Fs, type = 'your DFT name')
        time, x = generateSin(Fs, interval)
        # plot the signal in time domain
        #plotTime(x, time)
        # plot the signal in frequency domain
        plotSpectrum(x, Fs, type = 'DFT')
        plotTime(x, time, type = 'IDFT')
        # confirm DFT/IDFT correctness by checking if x == IDFT(DFT(x))
        print(np.allclose(x, myInDFT(myDFT(x))))
        print(np.allclose(np.fft.fft(x), myDFT(x)))
        # for further details, if any, check the lab document

    elif (sys.argv[1] == 'il2'):

        print('In-lab experiment 2 for the Fourier transform')

        print('In-lab experiment 2 for the Fourier transform')
        # use np.random.randn() for randomization
        # we can owverwrie the default values
        # frequency =  8.0                     # frequency of the signal
        # amplitude =  3.0                     # amplitude of the signal
        # phase = 1.0                          # phase of the signal
        # You should also numerically check if the signal energy
        # in time and frequency domains is identical
        # for further details, if any, check the lab document
        #Generate signal
        frequency = np.random.randn()
        amplitude = np.random.randn()
        phase = np.random.randn()
        time, x = generateSin(1000, 1, frequency, amplitude, phase)
        DFT = []
        IDFT = []
        DFT = myDFT(x)
        IDFT = myInDFT(DFT)
        sumDFT = 0
        sumIDFT = 0
        for i in range(1000):
            sumDFT = sumDFT + math.pow(abs(DFT[i]),2)
        sumDFT = sumDFT/1000
        for j in range(1000):
            sumIDFT = sumIDFT + math.pow(abs(IDFT[j]),2)

        # You should also numerically check if the signal energy
        # in time and frequency domains is identical
        print(np.allclose(sumDFT,sumIDFT))

        # for further details, if any, check the lab document

    elif (sys.argv[1] == 'il3'):

        print('In-lab experiment 3 for the Fourier transform')

        time, tone0 = generateSin(Fs, interval, 9, 17, 0.5)
        time, tone1 = generateSin(Fs, interval, 3, 2, 0.01)
        time, tone2 = generateSin(Fs, interval, 7, 30, 2.1)
        x = tone0 + tone1 + tone2
        # plot the signal in time domain
        plotTime(tone1, time)
        plotTime(x, time)
        # plot the signal in frequency domain
        plotSpectrum(x, Fs, type = 'FFT')
        # generate randomized multi-tone signals
        # plot them in both time and frequency domain

        # for further details, if any, check the lab document

    elif (sys.argv[1] == 'th'):

        print('Take-home exercise for the Fourier transform')
        x = 0
        # for specific details check the lab document
        interval = 20*interval

        duty_cyc = random.uniform(0,1)
        #print(duty_cyc)
        dt = 1.0/Fs
        time = np.arange(0, interval, dt)
        x = signal.square(time, duty_cyc)
        plotTime(x, time)
        plotSpectrum(x, Fs, type = 'DFT')

    else:

        cli_error_msg()

    plt.show()
