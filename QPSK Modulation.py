import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import random
N=10**5
bit_stream=np.random.randint(0,2,N)
L=16;Eb=2;Ts=2
t=np.arange(0,Ts,1/L)
cos_sig=np.sqrt(2*Eb/Ts)*np.cos(2*np.pi*t/Ts)
sin_sig=np.sqrt(2*Eb/Ts)*np.sin(2*np.pi*t/Ts)
mod_energy=cos_sig@cos_sig
qpsk_out=[]
for i in range(0,N,2):
    qpsk_out.extend((np.sqrt(2*Eb/Ts)*np.cos(2*np.pi*t/Ts+(bit_stream[i]*np.pi)))-(np.sqrt(2*Eb/Ts)*np.sin(2*np.pi*t/Ts+bit_stream[i+1]*np.pi)))
qpsk_out=np.array(qpsk_out)

Snr_db=np.linspace(-30,-10,10)
snr=10**(Snr_db/10)
ber_calc=0.5*erfc(np.sqrt(snr))
ber_values=[]
for val in snr:
    noise_pow=mod_energy/val
    noise_signal=np.sqrt(noise_pow/2)*(np.random.standard_normal(len(qpsk_out)))
    channel_sig=noise_signal+qpsk_out
    demodulated_out=[]
    for j in range(int(N/2)):
        correlated_value1=channel_sig[j*L*Ts:(j+1)*(L*Ts)]@cos_sig
        correlated_value2=channel_sig[j*L*Ts:(j+1)*(L*Ts)]@sin_sig
        demodulated_out.append(correlated_value1<0)
        demodulated_out.append(correlated_value2>0)
    ber_values.append(np.sum(demodulated_out!=bit_stream)/N)
plt.semilogy(Snr_db,ber_calc)
plt.semilogy(Snr_db,ber_values,linestyle='',marker='*')
plt.show()
    