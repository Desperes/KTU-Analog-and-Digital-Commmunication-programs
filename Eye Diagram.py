import numpy as np
import matplotlib.pyplot as plt
import sympy
def rc(t,Ts,alpha):
    rc_signal=(1/Ts)*(np.sin(np.pi*t/Ts)/(np.pi*t/Ts))*((np.cos(np.pi*alpha*t/Ts))/(1-(2*alpha*t/Ts)**2))
    x=sympy.symbols('x')
    expr=(1/Ts)*(sympy.sin(np.pi*x/Ts)/(np.pi*x/Ts))*((sympy.cos(np.pi*alpha*x/Ts))/(1-(2*alpha*x/Ts)**2))
    u=np.isnan(rc_signal)
    v=np.isinf(rc_signal)
    rc_signal[u]=(np.pi/(4*Ts))*(np.sin(np.pi/(2*alpha))/(np.pi/(2*alpha)))
    rc_signal[v]=(np.pi/(4*Ts))*(np.sin(np.pi/(2*alpha))/(np.pi/(2*alpha)))
    z_index=np.where(t==0)
    rc_signal[z_index]=sympy.limit(expr,x,0)
    return rc_signal
Ts=10
span=8*Ts
L=8
t=np.arange(-span/2,span/2+1/L,1/L)
alpha=0.4
rc_pulse=rc(t,Ts,alpha)
plt.figure(1)
plt.plot(t,rc_pulse)
plt.show()
N=50
binary_seq=np.random.randint(0,2,N)
polarized_binary_seq=2*binary_seq-1
upscale_seq=np.zeros(N*len(rc_pulse))
upscale_seq[::len(rc_pulse)]=polarized_binary_seq
plt.figure(2)
plt.title("Bianry sequence")
plt.stem(polarized_binary_seq)
plt.show()
print(polarized_binary_seq)
plt.figure(3)
plt.title("Upscaled binary sequence")
plt.plot(upscale_seq)
plt.show()
trans_sig=np.convolve(upscale_seq,rc_pulse)
plt.figure(4)
plt.title("Transmitted signal")
plt.plot(trans_sig)
plt.show()
Snr_db=10
snr=10**(Snr_db/10)
sig_pow=(trans_sig@trans_sig)/len(trans_sig)
noise_pow=sig_pow/snr
noise_sig=np.sqrt(noise_pow/2)*np.random.standard_normal(len(trans_sig))
channel_sig=noise_sig+trans_sig
rec_sig=np.convolve(channel_sig,rc_pulse)
rec_sig=rec_sig[int(len(rc_pulse)/2):-int(len(rc_pulse)/2)]
plt.figure(5)
plt.title("Matched filter output")
plt.plot(rec_sig)
plt.show()
plt.figure(6)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Eye diagram")
for i in range(N):
    plt.plot(channel_sig[i*len(rc_pulse):(i+1)*len(rc_pulse)])
plt.show()
