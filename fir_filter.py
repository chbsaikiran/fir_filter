from pylab import *
import scipy.signal as signal
import numpy as np
from matplotlib import pyplot as plt

#Plot frequency and phase response
def mfreqz(b,a=1):
    w,h = signal.freqz(b,a)
    h_dB = 20 * log10 (abs(h))
    subplot(211)
    plot(w/max(w),h_dB)
    ylim(-150, 5)
    ylabel('Magnitude (db)')
    xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    title(r'Frequency response')
    subplot(212)
    h_Phase = unwrap(arctan2(imag(h),real(h)))
    plot(w/max(w),h_Phase)
    ylabel('Phase (radians)')
    xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    title(r'Phase response')
    subplots_adjust(hspace=0.5)

#Plot step and impulse response
def impz(b,a=1):
    l = len(b)
    impulse = repeat(0.,l); impulse[0] =1.
    x = arange(0,l)
    response = signal.lfilter(b,a,impulse)
    subplot(211)
    stem(x, response)
    ylabel('Amplitude')
    xlabel(r'n (samples)')
    title(r'Impulse response')
    subplot(212)
    step = cumsum(response)
    stem(x, step)
    ylabel('Amplitude')
    xlabel(r'n (samples)')
    title(r'Step response')
    subplots_adjust(hspace=0.5)
    
def overlap_and_add_filter(b,x,y):
    y_temp1 = np.convolve(x[0:4000],b);
    y_temp2 = np.convolve(x[4000:8000],b);
    y_temp3 = np.convolve(x[8000:12000],b);
    y_temp4 = np.convolve(x[12000:16000],b);
    
    y[0:4000,0] = y_temp1[0:4000];
    y[4000:(4000+510),0] = y_temp1[4000:(4000+510)] + y_temp2[0:510];
    y[(4000+510):(4000+510+3490),0] = y_temp2[510:(510+3490)];
    i1 = (4000+510+3490);
    i2 = (4000+510+3490+510);
    y[i1:i2,0] = y_temp2[(510+3490):(510+3490+510)] + y_temp3[0:510];
    i1 = (4000+510+3490+510);
    i2 = (4000+510+3490+510+3490);
    y[i1:i2,0] = y_temp3[510:(510+3490)];
    i1 = (4000+510+3490+3490+510);
    i2 = (4000+510+3490+3490+510+510);
    y[i1:i2,0] = y_temp3[(510+3490):(510+3490+510)] + y_temp4[0:510];
    i1 = (4000+510+3490+510+3490+510);
    i2 = (4000+510+3490+510+3490+510+3490);
    y[i1:i2,0] = y_temp4[510:(510+3490)];
    

np.random.seed(42)  # for reproducibility
fs = 8000  # sampling rate, Hz
ts = np.arange(0, 2, 1.0 / fs)  # time vector - 5 seconds
ys = np.sin(2*np.pi * 200.0 * ts)  # signal @ 1.0 Hz, without noise
yerr = 0.5 * np.random.normal(size=len(ts))  # Gaussian noise
yraw = float32(ys + yerr)
output_file = open('input_pygen.bin', 'wb')
yraw.tofile(output_file)
output_file.close()
n = 511
b = float32(signal.firwin(n, cutoff = 0.2, window = "hamming"))
output_file = open('fir_ceoffs_pygen.bin', 'wb')
b.tofile(output_file)
output_file.close()

y_lfilter = np.zeros((len(ts),1));
zf = y_lfilter[0:510,0];

y_lfilter[0:4000,0], zf = signal.lfilter(b, 1,yraw[0:4000],zi=zf );
y_lfilter[4000:8000,0], zf = signal.lfilter(b, 1, yraw[4000:8000],zi=zf); 
y_lfilter[8000:12000,0], zf = signal.lfilter(b, 1,yraw[8000:12000] ,zi=zf);
y_lfilter[12000:16000,0], zf = signal.lfilter(b, 1,yraw[12000:16000] ,zi=zf);
    
tot_outs = 16000;
y_lfilter_full = np.zeros((len(ts),1));
#y_lfilter_full[0:tot_outs,0] = signal.lfilter(b, 1, yraw)
y_lfilter_full[0:tot_outs,0] = np.fromfile("out_msvc.bin", dtype=np.float32)
#overlap_and_add_filter(b,yraw,y_lfilter_full);
error = np.square(float32(y_lfilter_full[0:tot_outs,0])-float32(y_lfilter[0:tot_outs,0]))
abs_error = 10*log10(np.sum(error)) - 10*log10(len(ts))

plt.figure(figsize=[6.4, 2.4])

plt.plot(ts[0:tot_outs], y_lfilter[0:tot_outs,0], label="Raw signal")
plt.plot(ts[0:tot_outs], error[0:tot_outs], alpha=0.8, lw=3, label="SciPy lfilter")

plt.xlabel("Time / s")
plt.ylabel("Amplitude")
plt.legend(loc="lower center", bbox_to_anchor=[0.5, 1],
           ncol=2, fontsize="smaller")
#Frequency and phase response
#mfreqz(b)
#show()
#Impulse and step response
#figure(2)
#impz(a)
#show()