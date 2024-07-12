#BPSK MODULATION
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
N=10**5
bit_stream=np.random.randint(0,2,N)
Eb=2
Ts=2
L=16
t=np.arange(0,Ts,1/L)
mod_sig=np.sqrt(2*Eb/Ts)*(np.sin(2*np.pi*t/Ts))
trans_sig=[]
for i in bit_stream:
    trans_sig.extend(np.sqrt((2*Eb)/Ts)*(np.sin(2*np.pi*t/Ts+(i*np.pi))))
trans_sig=np.array(trans_sig)
plt.figure(2)
plt.title("transmitted signal")
plt.plot(trans_sig)
plt.show()

Snr_db=np.linspace(-30,-10,10)
snr=10**(Snr_db/10)
ber_calc=0.5*erfc(np.sqrt(snr))
ber_values=[]
for i in snr:
    noise_pow=(mod_sig@mod_sig)/i
    noise_sig=np.sqrt(noise_pow/2)*(np.random.standard_normal(len(trans_sig)))
    channel_sig=trans_sig+noise_sig
    demodulated_out=[]
    for j in range(N):
        correlated_value=mod_sig@channel_sig[j*(L*Ts):(j+1)*(L*Ts)]
        demodulated_out.append((int(correlated_value<1)))
        
    ber_values.append(np.sum(demodulated_out!=bit_stream)/N)
plt.figure(3)
plt.title("bit error rate vs snr")
plt.semilogy(Snr_db,ber_calc)
plt.semilogy(Snr_db,ber_values,linestyle="",marker='*')

plt.show()
    
        
        
