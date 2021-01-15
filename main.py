import numpy as np
import matplotlib.pyplot as plt
from sympy.combinatorics.graycode import GrayCode
import statistics as st
import math
from scipy import special
import gc
from scipy.io import wavfile
from scipy.io.wavfile import write
import scipy.io

AM1 = 7
fm1 = 7
AM2 = 9
fm2 = 3

#Question 1

#part a
def erwthma_1_ai(AM,fm,title):
  plt.clf()
  
  T=1/(2*fm) 
  fs1=20*fm
  t=np.arange(0,4*T,(1/fs1)) #ftiaxnei ton xroniko axona gia 4 periodous symfwna me th syxnothta deigmatolhpsias
  y=np.cos(2*np.pi*fm*t)*np.cos(2*np.pi*(AM+2)*fm*t)
  
  #dhmiourgia diagrammatos
  plt.figure(1)
  plt.xlabel("time(ms)")
  plt.ylabel("voltage(V)")
  plt.title(title)
  plt.plot(t,y,'-ob')
  plt.grid()
  plt.savefig(title)
  return
erwthma_1_ai(AM1,fm1,"1_a_i_Sample_el18007")
erwthma_1_ai(AM2,fm2,"1_a_i_Sample_el18129")

def erwthma_1_aii(AM,fm,title):
  plt.clf()
  T=1/(2*fm)
  fs2=100*fm
  t=np.arange(0,4*T,(1/fs2))
  y=np.cos(2*np.pi*fm*t)*np.cos(2*np.pi*(AM+2)*fm*t)

  plt.figure(2)
  plt.xlabel("time(ms)")
  plt.ylabel("voltage(V)")
  plt.title(title)
  plt.plot(t,y,'-ob')
  plt.grid()
  plt.savefig(title)
erwthma_1_aii(AM1,fm1,"1_a_ii_Sample_el18007")
erwthma_1_aii(AM2,fm2,"1_a_ii_Sample_el18129")
def erwthma_1_aiii(AM,fm,title):
  plt.clf()
  T=1/(2*fm)
  fs2=100*fm
  t=np.arange(0,4*T,(1/fs2))
  y=np.cos(2*np.pi*fm*t)*np.cos(2*np.pi*(AM+2)*fm*t)

  #deigmatolipsia 100fm
  plt.figure(1)
  plt.xlabel("time(ms)")
  plt.ylabel("voltage(V)")
  plt.title(title)
  plt.plot(t,y,'-or', label='fs2')
  plt.legend()
  plt.savefig(title)
  T=1/(2*fm) 
  fs1=20*fm
  t=np.arange(0,4*T,(1/fs1))
  y=np.cos(2*np.pi*fm*t)*np.cos(2*np.pi*(AM+2)*fm*t)

  #deigmatolipsia 20fm
  plt.figure(1)
  plt.xlabel("time(ms)")
  plt.ylabel("voltage(V)")
  plt.title(title)
  plt.plot(t,y,'-ob',label='fs1')
  plt.legend()
  plt.grid()
  plt.savefig(title)
erwthma_1_aiii(AM1,fm1,"1_a_iii_Sample_el18007")
erwthma_1_aiii(AM2,fm2,"1_a_iii_Sample_el18129")

#part b
def erwthma_1_b(AM,fm,title):
  plt.clf()
  T=1/(2*fm) 
  fs1=5*fm
  t=np.arange(0,4*T,(1/fs1))
  y=np.cos(2*np.pi*fm*t)*np.cos(2*np.pi*(AM+2)*fm*t)
  #5fm
  plt.figure(1)
  plt.xlabel("time(ms)")
  plt.ylabel("voltage(V)")
  plt.title(title)
  plt.plot(t,y,'-ob')
  plt.grid()
  plt.savefig(title)
  return
erwthma_1_b(AM1,fm1,"1_b_Sample_el18007")
erwthma_1_b(AM2,fm2,"1_b_Sample_el18129")

#Question 2

#part a

def erwthma2_a(AM,fm,title):
    T=1/(2*fm)
    fs1=20*fm
    t=np.arange(0,4*T,(1/fs1))
    y=np.cos(2*np.pi*fm*t)*np.cos(2*np.pi*(AM+2)*fm*t)

    L=2**5 
    D=2/L #megethos bhmatos kbantisti

    def midriser(x):
        return D*(math.floor(x/D)+1/2)
    i=midriser(-2)
    
    #dhmiourgia epipedwn kbantishs
    ql=[]
    for i in range(L):
        ql.append(round(i,6))
        i+=D
    y=((y+1)*16-1)
    
		#antistoixisi simatos sta epipeda kbantisis
    yq=list(map(midriser,y))
    for j in range(len(yq)):
        yq[j]=round(yq[j],0)
    gv=GrayCode(5)
    gv=list(gv.generate_gray())

    plt.clf()
    plt.figure()
    plt.xlabel("time(ms)")
    plt.ylabel("midriser output")
    plt.title(title)
    plt.yticks(ql,gv)
    plt.step(t,yq)
    plt.grid()
    plt.savefig(title)
    return(yq, y)
(yq1, y1)=erwthma2_a(AM1,fm1,"2_a_Sample_el18007")
(yq2, y2)=erwthma2_a(AM2,fm2,"2_a_Sample_el18129")

#part b
def erwthma2_b(yq, y, b):
    print("Erwthma 2.b.")
    for j in range(len(yq)):
        yq[j] = ((yq[j] + 1) / 16 - 1)
        y[j] = ((y[j] + 1) / 16 - 1)
    D=2/(2**5)
    print(b, "theoretical sigma:", np.sqrt(D**2/12))
    err = y - yq

    #erwthma 2 b (i)
    sn1 = st.stdev(err[:10])
    print(b, "sigma10=", sn1)

    #erwthma 2 b (ii)
    sn2 = st.stdev(err[:20])
    print(b, "sigma20=", sn2)
    
    #erwthma 2 b (iii)
    ss1 = st.stdev(y[:10])
    ss2 = st.stdev(y[:20])


    def db(x):
        return 10 * math.log10(x)


    SNR1 = ss1 / sn1
    SNR1 = SNR1**2
    print(b, "SNR10=", SNR1)
    SNR1 = db(SNR1)
    print(b, "SNR10(dB) = ", SNR1)
		
    SNR2= ss2 / sn2
    SNR2=SNR2**2
    print(b, "SNR20=", SNR2)
    SNR2=db(SNR2)
    print(b, "SNR20(dB)=", SNR2)

    return
erwthma2_b(yq1, y1,"el18007")
erwthma2_b(yq2, y2,"el18129")

#part c
def erwthma2_c(yq, fm, title):
  def polar_NRZ(x):
    b=[]
    for i in range(5):
        if x[i]=="0":
            b.append((-1)*fm)
        else:
            b.append(fm)
    return b

  #polar NRS bitstream
  def polar_NRZ_bs(qsig):
    strm=[]
    for i in range(len(qsig)):
        strm+=polar_NRZ(qsig[i])
		    
    return strm
  
	#synarthsh quantized to signal
  def q2g(qsig,gv,ql):
    a=[]
    for i in range(len(qsig)):
        x=ql.index(qsig[i])
        a.append(gv[x])
    return a

  
  gv=GrayCode(5)
  gv=list(gv.generate_gray())
  ql=[]
  
  for i in range(32): #Dhmiourgia Epipedwn Kbantismou
    z=-1+(i+1)*(1/16)
    ql.append(z)
  a=q2g(yq, gv, ql)
  b=polar_NRZ_bs(a)
  b=b[:50] #ta stoixeia tou gia mia periodo
  n=np.arange(0, 0.05, 0.001)

  plt.clf()
  plt.figure(2)
  plt.step(n, b)
  plt.xlabel("time(ms)")
  plt.ylabel("voltage(V)")
  plt.title(title)
  plt.grid()
  plt.savefig(title)
  return
erwthma2_c(yq1,fm1,"2_c_Sample_el18007")
erwthma2_c(yq2,fm2,"2_c_Sample_el18129")
yq1.clear()
yq2.clear()

#Question 3
A1 = 0 + 0 + 7
A2 = 3 #1 + 2 + 9 = 12, 1 + 2 = 3
Tb = 0.5  # bit - pulse duration
Eb1 = (A1**2) * Tb # mean of symbol energy Eb=Integral(si^2dt)from 0 to Tb
Eb2 = (A2**2) * Tb # mean of symbol energy Eb=Integral(si^2dt)from 0 to Tb
R = 1/Tb  # transmission
samples = 50  # samples of each pulse
fs = int(samples/Tb)  # sampling frequency


#function to simulate adding "Additive White Gaussian Noise" to the signal
def add_noise_to_signal(signal, SNR, Eb):
   No = Eb/(10**(SNR/10)) #SNR is given in db
   noise=np.random.normal(0,np.sqrt(No/2), size=len(signal))+1j*np.random.normal(0,np.sqrt(No/2), size=len(signal))
   #Αdd the generated noise to the signal 
   r = noise + signal

   return r

#part a



def part3_a(A,Title):
   plt.clf()
   bits = 46
   bitstream = np.random.randint(2, size=46)
   # BPAM, symbol assignment
   BPAM = []
   BPAM.append(0)
   for i in range(bits):
       bitstream[i] = -A if bitstream[i] == 0 else A
       BPAM.append(bitstream[i])
   t1=np.arange( 0,23.5, 0.5)
   plt.xlabel("time")
   plt.ylabel("voltage(V)")
   plt.step(t1,BPAM,"r")
   plt.grid()
   plt.title(Title)
   plt.savefig(Title)
   return(BPAM) 

BPAM1 = part3_a(A1,'3_a_Sample_el18007')
BPAM2 = part3_a(A2,'3_a_Sample_el18129')   

#part b


def part3_b(A,Eb,Title):
   plt.clf()
   #constellation of the BPAM
   s0 = -np.sqrt(Eb)
   s1 = np.sqrt(Eb)
   x = [s0.real, s1.real]
   y = [s0.imag, s1.imag]
   plt.figure(figsize=(6,6))
   plt.grid(True)
   plt.plot(x, y, 'o')
   plt.title('Constellation Diagram for B-PAM')
   plt.xlabel('Re')
   plt.ylabel('Im')
   plt.axhline(y=0, linewidth=0.5, color='black')
   plt.axvline(x=0, linewidth=0.5, color='black')
   bits = ['0', '1']
   for i, bit in enumerate(bits):
       plt.annotate(bit, (x[i] + 0.01, y[i] + 0.003))
   plt.axvline(color="blue",label="Decision Threshold")
   plt.legend(loc=0,fontsize="small")
   plt.savefig(Title)

part3_b(A1,Eb1,'3_b_Conste03118007')
part3_b(A2,Eb2,'3_b_Conste03118129')

#part c


def part3_c(BPAM,SNR,Eb,Title,timi):
   plt.clf()
   #labame 50 digmata gia na petyxoume orthi apikonisi me th synarthsh plot
   BPAM_50x=[]
   for i in range(46):
     i = i+1
     for j in range(50):
         BPAM_50x.append(BPAM[i])
   t1=np.arange(0,23,0.5/50)
   noisy = add_noise_to_signal(BPAM_50x, SNR, Eb)
   plt.figure(figsize = (20, 8))
   plt.clf()
   plt.xlabel("time")
   plt.ylabel("voltage(V)")
   plt.grid()
   plt.plot( t1,noisy.real, label=timi)
   plt.title(Title)
   plt.savefig(Title)
   return noisy

noisy_5_1 = part3_c(BPAM1,5,Eb1,'3_c_awgn+signal5_el18007','SNR1 = Eb/N0 = 5 dB')
noisy_15_1 = part3_c(BPAM1,15,Eb1,'3_c_awgn+signal15_el18007','SNR2 = Eb/N0 = 15 dB')
noisy_5_2 = part3_c(BPAM2,5,Eb2,'3_c_awgn+signal5_el18129','SNR1 = Eb/N0 = 5 dB')
noisy_15_2 = part3_c(BPAM2,15,Eb2,'3_c_awgn+signal15_el18129','SNR2 = Eb/N0 = 15 dB')

#part d

def part3_d(noisy,Eb,Title,AM):
   plt.clf()
   s0 = -np.sqrt(Eb)
   s1 = np.sqrt(Eb)
	 #diagramma asterismou
   plt.figure(figsize = (10, 10))
   plt.grid()
   plt.scatter(np.real(noisy)*np.sqrt(Tb),np.imag(noisy)*np.sqrt(Tb),label="Received signal")
   plt.scatter(s1,0,label="Trasmitted signal")
   plt.scatter(s0,0,label="Trasmitted signal")
   plt.axvline(color="black",label="Decision Threshold")
   plt.legend()
   plt.xlabel("Re")
   plt.ylabel("Im")
   plt.xlim((-1)*np.sqrt(Eb)-np.sqrt(AM)-3,np.sqrt(Eb)+np.sqrt(AM)+3)
   plt.ylim((-1)*np.sqrt(Eb)-np.sqrt(AM)-3,np.sqrt(Eb)+np.sqrt(AM)+3)
   plt.title(Title)
   plt.savefig(Title)
   
part3_d(noisy_5_1,Eb1,'3_d_Sample_el18007_5',AM1)
part3_d(noisy_15_1,Eb1,'3_d_Sample_el18007_15',AM1)
part3_d(noisy_5_2,Eb2,'3_d_Sample_el18129_5',AM2)
part3_d(noisy_15_2,Eb2,'3_d_Sample_el18129_15',AM2)
#part e
noisy_5_1=[]
noisy_15_1=[]
noisy_5_2=[]
noisy_15_2=[]
BPAM1=[]
BPAM2=[]

def part3_e(Eb,A,Title): 
   gc.collect()
   N1=10**6
   s0 = -np.sqrt(Eb)
   s1 = np.sqrt(Eb)
   b=np.random.randint(2,size=N1)
   const=[]
   for i in range(len(b)):
    const.append(s0) if b[i] == 0 else const.append(s1)
   #dhmiourgia pinaka gia ta pososta
   percentage=[]
   for d in range(0,16):
     noisy = add_noise_to_signal(const, d, Eb)
     counter = 0
     for i in range(len(noisy)):
      if const[i]*(np.real(noisy[i]))<0: 
          counter+=1  
     percentage.append(counter/N1)
   #dhmiourgia pinaka gia ta thewritika pososta
   th=[]
   for d in range(0,16):
	   th.append(0.5*special.erfc(np.sqrt(10**(d/10))))
   plt.clf()
   dB=np.arange(0,16,1)
   plt.yscale("log")
   plt.plot(dB, th, '-or')
   plt.plot(dB, percentage, 'ob')
   plt.grid()
   plt.title(Title)
   plt.savefig(Title)

part3_e(Eb1,A1,'3_e_Sample_el18007')
part3_e(Eb2,A2,'3_e_Sample_el18129')

#Question 4

A1 = 7
A2 = 3
bits = 46
Es1 = (A1**2)
Es2 = (A2**2)
Eb1 = Es1/2
Eb2 = Es2/2

#part a
def part1(Es,Title):
   currsym=[]
   sym2plot=[]
   bs = np.random.randint(2, size=46)
   #anagnwsh tou bitstream kai dhmiourgia akoloythias symbolwn
   for i in range(0, 45, 2):
     if (bs[i]==0 and bs[i+1]==0):
       currsym.append("00")
       sym2plot.append((np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es))
     elif (bs[i]==0 and bs[i+1]==1):
       currsym.append("01")
       sym2plot.append((np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es))
     elif (bs[i]==1 and bs[i+1]==0):
       currsym.append("10")
       sym2plot.append(-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es))
     else:
       currsym.append("11")  
       sym2plot.append(-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es))
   #ta 4 simeia tou asterismou
   z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
   x=[]
   y=[]
   for i in range(23):
     x.append(sym2plot[i].real)
     y.append(sym2plot[i].imag)
   plt.clf()
   plt.figure(1)
   plt.grid(True)
   plt.plot(x, y, 'o')
   plt.title('Constellation Diagram for pi/4 QPSK')
   plt.xlabel('Re')
   plt.ylabel('Im')
   plt.axhline(y=0, linewidth=0.5, color='blue')
   plt.axvline(x=0, linewidth=0.5, color='blue')
   plt.annotate('00', (z00.real+0.1, z00.imag+0.1))
   plt.annotate('01', (z01.real+0.1, z01.imag+0.1))
   plt.annotate('11', (z11.real+0.1, z11.imag+0.1))
   plt.annotate('10', (z10.real+0.1, z10.imag+0.1))
   plt.title(Title)
   plt.savefig(Title)
   return currsym
currsym1 = part1(Es1,'4_a_Conste_el18007')
currsym2 = part1(Es2,'4_a_Conste_el18129')

Tb=0.5


def signal (A,Es,Title,currsym):
   t=np.arange(0,2*Tb,0.05)
   s00=[]
   s01=[]
   s10=[]
   s11=[]
	 #xrhomopoieitai gia th diamorfwsh qpsk:
   for i in range(len(t)):
     r1=np.sqrt(2)*A*np.cos(2*np.pi*2*t[i]+3*np.pi/4)
     s00.append(r1)
     r2=np.sqrt(2)*A*np.cos(2*np.pi*2*t[i]+np.pi/4)
     s10.append(r2)
     r3=np.sqrt(2)*A*np.cos(2*np.pi*2*t[i]-3*np.pi/4)
     s01.append(r3)
     r4=np.sqrt(2)*A*np.cos(2*np.pi*2*t[i]-np.pi/4)
     s11.append(r4)
   z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)

   QPSK=[] #   QPSK Modulation
   qpsk_conste=[]
   for i in range(23):
     if (currsym[i]=="00"):
       for j in range(len(t)):
         QPSK.append(s00[j])
         qpsk_conste.append(z00)
     elif (currsym[i]=="01"):
       for j in range(len(t)):
         QPSK.append(s01[j])
         qpsk_conste.append(z01)
     elif (currsym[i]=="10"):

       for j in range(len(t)):
         QPSK.append(s10[j])
         qpsk_conste.append(z10)
     else:
       for j in range(len(t)):
         QPSK.append(s11[j])
         qpsk_conste.append(z11)

   plt.clf()
   plt.figure(figsize=(20,8))
   t1=np.arange(0,23,0.05)     
   plt.plot(t1,QPSK)
   plt.savefig(Title)
   return (QPSK,qpsk_conste)

(QPSK1,qpsk_conste_1) = signal (A1,Es1,'QPSK_el18007',currsym1)
(QPSK2,qpsk_conste_2) = signal (A2,Es2,'QPSK_el18129',currsym2)

#part b

#function to simulate adding "Additive White Gaussian Noise" to the signal
def add_noise_to_signal(signal, SNR, Es):
    No = Es/(10**(SNR/10)) #SNR is given in db
    noise=np.random.normal(0,np.sqrt(No/2), size=len(signal))+1j*np.random.normal(0,np.sqrt(No/2), size=len(signal))
    #Αdd the generated noise to the signal 
    r = noise + signal
    return r

def noisy_signal(QPSK,SNR,Es,Title):
   noisy = add_noise_to_signal(QPSK, SNR, Es)
   t2=np.arange(0,23,0.05)
   plt.clf()
   plt.plot(t2,noisy.real)
   plt.savefig(Title)

noisy_signal(QPSK1,5,Es2,'Signal+Noise_5dB_el18007')
noisy_signal(QPSK1,15,Es2,'Signal+Noise_15dB_el18007')
noisy_signal(QPSK2,5,Es2,'Signal+Noise_5dB_el18129')
noisy_signal(QPSK2,15,Es2,'Signal+Noise_15dB_el18129')

def noisy_conste(qpsk_conste_1,Es,SNR,Title):
   plt.figure(figsize = (8, 8))

   z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)

   qpsk_conste=add_noise_to_signal(qpsk_conste_1, SNR , Es)
   plt.grid()
   plt.scatter(np.real(qpsk_conste),np.imag(qpsk_conste),label="Received signal")
   plt.scatter(z00.real,z00.imag,label="Trasmitted signal")
   plt.scatter(z01.real,z01.imag,label="Trasmitted signal")
   plt.scatter(z10.real,z10.imag,label="Trasmitted signal")
   plt.scatter(z11.real,z11.imag,label="Trasmitted signal")
   plt.axvline(color="black",label="Decision Threshold")
   plt.axhline(color="black",label="Decision Threshold")
   plt.legend()
   plt.xlabel("I")
   plt.ylabel("Q")
   plt.title(Title)
   plt.savefig(Title)

noisy_conste(qpsk_conste_1,Es1,5,'4_b_Constellation+Noise_5dB_el18007')
noisy_conste(qpsk_conste_1,Es1,15,'4_b_Constellation+Noise_15dB_el18007')
noisy_conste(qpsk_conste_2,Es2,5,'4_b_Constellation+Noise_5dB_el18129')
noisy_conste(qpsk_conste_2,Es2,15,'4_b_Constellation+Noise_15dB_el18129')

#part c

def part_c(Eb,Es,Title1,Title2):
   N1=10**6
   N2=int(np.floor(N1/2))
   bs=np.random.randint(2,size=N1)
   QPSK = []
   QPSK.append(0)
   currsym=[]
   z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
   for i in range(0, N1, 2):
     if (bs[i]==0 and bs[i+1]==0):
       currsym.append("00")
     elif (bs[i]==0 and bs[i+1]==1):
       currsym.append("01")
     elif (bs[i]==1 and bs[i+1]==0):
       currsym.append("10")
     else:
       currsym.append("11")  
   #pinakas asterismou
   QPSKc=[]
   for i in range(N2):
     if (currsym[i]=="00"):
         QPSKc.append(z00)
     elif (currsym[i]=="01"):
         QPSKc.append(z01)
     elif (currsym[i]=="10"):
         QPSKc.append(z10)
     else:
         QPSKc.append(z11)
   percentage=[]
   for d in range(0,16):
     noisy = add_noise_to_signal(QPSKc, d, Eb)
     counter = 0
     for i in range(len(noisy)):
        if (((currsym[i]=="00") and ((noisy[i].real<0) or (noisy[i].imag>0))) or ((currsym[i]=="01") and ((noisy[i].real<0) or (noisy[i].imag<0))) or ((currsym[i]=="11") and ((noisy[i].real>0) or (noisy[i].imag<0))) or ((currsym[i]=="10") and ((noisy[i].real>0) or (noisy[i].imag>0)))):
          counter+=1 #auxisi metriti kata 1 gia diaforetika symbola
        if (((currsym[i]=="00") and ((noisy[i].real<0) and (noisy[i].imag>0))) or ((currsym[i]=="01") and ((noisy[i].real<0) and (noisy[i].imag<0))) or ((currsym[i]=="11") and ((noisy[i].real>0) and (noisy[i].imag<0))) or ((currsym[i]=="10") and ((noisy[i].real>0) and (noisy[i].imag>0)))):
         counter+=1 #peretero auxisi gia dio lathos bits
          
     percentage.append(counter/N1)
   th=[]
   for d in range(0,16):
	   th.append(0.5*special.erfc(np.sqrt(10**(d/10))))

   b=np.random.randint(2,size=N1)
   BPSK = []
   for i in range(N1):
     if(b[i]==1):
       BPSK.append(np.sqrt(Eb))
     else:
       BPSK.append(-np.sqrt(Eb))
     
   percentagebpsk=[]
   for d in range(0,16):
     noisy = add_noise_to_signal(BPSK, d, Eb)
     counter = 0
     for i in range(len(noisy)):
        if (((bool(b[i]==0)) & bool(noisy[i] > 0)) | ((bool(b[i]==1) & bool(noisy[i] < 0)))):
          counter+=1
     percentagebpsk.append(counter/N1)
    
   thbpsk=[]
   for d in range(0,16):
   	thbpsk.append(0.5*special.erfc(np.sqrt(10**(d/10))))

   plt.clf()
   dB=np.arange(0,16,1)
   plt.yscale("log")
   plt.plot(dB, th, '-or', label="Theoretical QPSK")
   plt.plot(dB, percentage, 'ob', label="Percentage QPSK")
   plt.legend()
   plt.grid()
   plt.savefig(Title1)

   plt.clf()
   plt.yscale("log")
   plt.plot(dB, th, '-or', label="Theoretical QPSK")
   plt.plot(dB, thbpsk, '-oy', label="Theoretical BPSK")
   plt.plot(dB, percentage, 'ob', label="Percentage QPSK")
   plt.plot(dB, percentagebpsk, 'og', label="Percentage BPSK")
   plt.legend()
   plt.grid()
   plt.title(Title2)
   plt.savefig(Title2)

part_c(Eb1,Es1,'4_c_BEP_el18007','4_c_BEP_Compare_el18007')
part_c(Eb2,Es2,'4_c_BEP_el18129','4_c_BEP_Compare_el18129')

#part d

#(i)

def text_to_bits(text, encoding='ascii'):
   bits=bin(int.from_bytes(text.encode(encoding), 'big'))[2:]
   return bits.zfill(8*((len(bits)+7)//8))

#el18007, 0+0+7=7 odd
#el18129, 1+2+9=12 even
#ascii to bits

text1=open('shannon_odd.txt', 'r')
tx1=text1.read()
text1.close()
bnr1=text_to_bits(tx1)
bnr1=bnr1.replace("0b","")

text2=open('shannon_even.txt', 'r')
tx2=text2.read()
text2.close()
tx2=tx2.replace("—","")
bnr2=text_to_bits(tx2)
bnr2=bnr2.replace("0b","")

#(ii)

def kvantismos(bnr, Title):
   t=np.arange(0,len(bnr),8)
   dec=[]
   dec2=[]
   for i in range(0, len(bnr), 8):
     k=0
     for j in range(8):
        k = k + int(bnr[i+j])*(2**(7-j))
     dec.append(k-128)
     dec2.append(k)

   L=2**8
   D=2/L

   def midriser(x):
     return D*(np.floor(x/D)+1/2)

   i=midriser(-2)
   ql=[]
   for i in range(L):
     ql.append(round(i,6))
   yq=list(map(midriser,dec))
   gv=GrayCode(8)
   gv=list(gv.generate_gray())

   plt.clf()
   plt.figure(figsize=(30, 5))
   plt.xlabel("time(s)")
   plt.ylabel("midriser output")
   plt.yticks(ql,gv) #antistoixizw epipeda Gray apo 0-256
   plt.step(t,yq)
   plt.grid()
   plt.title(Title)
   plt.savefig(Title)
   return dec2,gv

(dec2_1,gv) = kvantismos(bnr1, '4_d_ii_Binary_el18007')
(dec2_2,gv) = kvantismos(bnr2, '4_d_ii_Binary_el18129')

#(iii)

Es = 1**2

def qpsk(dec2,gv):
   signalwithGray=[]
   for i in range(len(dec2)):
     a = dec2[i]
     signalwithGray.append(gv[a])

   z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)

   #diamorfwsh qpsk
   t=np.arange(0,2*Tb,0.05)
   s00=[]
   s01=[]
   s10=[]
   s11=[]
   for i in range(len(t)):
     r1=np.sqrt(2)*np.cos(2*np.pi*2*t[i]+3*np.pi/4)
     s00.append(r1)
     r2=np.sqrt(2)*np.cos(2*np.pi*2*t[i]+np.pi/4)
     s10.append(r2)
     r3=np.sqrt(2)*np.cos(2*np.pi*2*t[i]-3*np.pi/4)
     s01.append(r3)
     r4=np.sqrt(2)*np.cos(2*np.pi*2*t[i]-np.pi/4)
     s11.append(r4)
   
	 #symbol sequence to bit list
   bs=''.join(signalwithGray)
   currsym=[]
   for i in range(0, len(bs), 2):
     if (bs[i]=='0' and bs[i+1]=='0'):
       currsym.append("00")
     elif (bs[i]=='0' and bs[i+1]=='1'):
       currsym.append("01")
     elif (bs[i]=='1' and bs[i+1]=='0'):
       currsym.append("10")
     else:
       currsym.append("11")  

   QPSK=[] #   QPSK Modulation
   qpsk_conste_1=[]
   for i in range(len(currsym)):
     if (currsym[i]=="00"):
       qpsk_conste_1.append(z00)
       for j in range(len(t)):
         QPSK.append(s00[j])

     elif (currsym[i]=="01"):
       qpsk_conste_1.append(z01)
       for j in range(len(t)):
         QPSK.append(s01[j])
      
     elif (currsym[i]=="10"):
       qpsk_conste_1.append(z10)   
       for j in range(len(t)):
         QPSK.append(s10[j])
      
     else:
       qpsk_conste_1.append(z11)
       for j in range(len(t)):
         QPSK.append(s11[j])
   return QPSK,qpsk_conste_1,currsym

(QPSK1,qpsk_conste_1_1,currsym1)=qpsk(dec2_1,gv)
(QPSK2,qpsk_conste_1_2,currsym2)=qpsk(dec2_2,gv)

# (iv), (v)

def noisy_qpsk(QPSK,qpsk_conste_1,SNR,Title): 
   noisy = add_noise_to_signal(QPSK, SNR, Es)

   plt.figure(figsize = (10, 10))
   qpsk_conste_AWGN=add_noise_to_signal(qpsk_conste_1,SNR,Es)
   z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
   z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
   
   plt.grid()
   plt.scatter(np.real(qpsk_conste_AWGN),np.imag(qpsk_conste_AWGN),label="Received signal")
   plt.scatter(z00.real,z00.imag,label="Trasmitted signal")
   plt.scatter(z01.real,z01.imag,label="Trasmitted signal")
   plt.scatter(z10.real,z10.imag,label="Trasmitted signal")
   plt.scatter(z11.real,z11.imag,label="Trasmitted signal")
   plt.axvline(color="black",label="Decision Threshold")
   plt.axhline(color="black",label="Decision Threshold")
   plt.legend()
   plt.xlabel("I")
   plt.ylabel("Q")
   plt.title(Title)
   plt.savefig(Title)
   return qpsk_conste_AWGN

qpsk_conste_AWGN_5_1=noisy_qpsk(QPSK1,qpsk_conste_1_1,5,'4_d_v_Conste+Noise_5dB_Partd_el18007')
qpsk_conste_AWGN_15_1=noisy_qpsk(QPSK1,qpsk_conste_1_1,15,'4_d_v_Conste+Noise_15dB_Partd_el18007')
qpsk_conste_AWGN_5_2=noisy_qpsk(QPSK2,qpsk_conste_1_2,5,'4_d_v_Conste+Noise_5dB_Partd_el18129')
qpsk_conste_AWGN_15_2=noisy_qpsk(QPSK2,qpsk_conste_1_2,15,'4_d_v_Conste+Noise_15dB_Partd_el18129')

def decoding_conste(qpsk_conste):
	 #apodiamorfwsh
   decurrsym=[]
   for i in range(0, len(qpsk_conste), 1):
     if (qpsk_conste[i].real>0 and qpsk_conste[i].imag<0):
       decurrsym.append('00')
    
     elif (qpsk_conste[i].real>0 and qpsk_conste[i].imag>0):
       decurrsym.append('01')
    
     elif (qpsk_conste[i].real<0 and qpsk_conste[i].imag<0):
       decurrsym.append('10')
    
     else:
       decurrsym.append('11')  
   return decurrsym

decurrsym_5_1=decoding_conste(qpsk_conste_AWGN_5_1)
decurrsym_15_1=decoding_conste(qpsk_conste_AWGN_15_1)
decurrsym_5_2=decoding_conste(qpsk_conste_AWGN_5_2)
decurrsym_15_2=decoding_conste(qpsk_conste_AWGN_15_2)

#(vi)

def possibility(qpsk_conste,currsym,SNR,a,b):
   percentage=[]
   counter = 0
   for i in range(len(currsym)):
     if (((currsym[i]=="00") and ((qpsk_conste[i].real<0) or (qpsk_conste[i].imag>0))) or ((currsym[i]=="01") and ((qpsk_conste[i].real<0) or (qpsk_conste[i].imag<0))) or ((currsym[i]=="11") and ((qpsk_conste[i].real>0) or (qpsk_conste[i].imag<0))) or ((currsym[i]=="10") and ((qpsk_conste[i].real>0) or (qpsk_conste[i].imag>0)))):
       counter+=1
     if (((currsym[i]=="00") and ((qpsk_conste[i].real<0) and (qpsk_conste[i].imag>0))) or ((currsym[i]=="01") and ((qpsk_conste[i].real<0) and (qpsk_conste[i].imag<0))) or ((currsym[i]=="11") and ((qpsk_conste[i].real>0) and (qpsk_conste[i].imag<0))) or ((currsym[i]=="10") and ((qpsk_conste[i].real>0) and (qpsk_conste[i].imag>0)))):
       counter+=1     
   percentage.append(counter/len(currsym))
   print(a, 0.5*counter*100/len(currsym),"%")

   print(b, 50*special.erfc(np.sqrt(0.5*10**((SNR/10)))), "%")

possibility(qpsk_conste_AWGN_5_1,currsym1,5,'4_d_vi_Percentage(5dB,el18007)=','Tpercentage(5dB,el18007)=')
possibility(qpsk_conste_AWGN_15_1,currsym1,15,'4_d_vi_Percentage(15dB,el18007)=','Tpercentage(15dB,el18007)=')
possibility(qpsk_conste_AWGN_5_2,currsym2,5,'4_d_vi_Percentage(5dB,el18129)=','Tpercentage(5dB,el18129)=')
possibility(qpsk_conste_AWGN_15_2,currsym2,15,'4_d_vi_Percentage(15dB,el18129)=','Tpercentage(15dB,el18129)=')   

#(vii)

def reconstruction(decurrsym,Title):
   string=[]
   for i in range(0,len(decurrsym),4):
     a=''
     for j in range(4):
       a+=decurrsym[i+j]
     string.append(a)
    
   #diadikasia metatrophs gray se dyadiko
   degray2=[]
   for i in range(len(string)):
     degray=gv.index(string[i])#vriskei th thesi tou stoixeiou ston pinaka gray
     inter=bin(degray)

     interlendiff=10-len(str(inter))
     for j in range(interlendiff):
       degray2.append('0')
     degray2.append(bin(degray))

   degray3=''.join(degray2)

   decode=degray3.replace("0b","")

   def bits_to_text(bits):    
       debytes=[bits[8*i:8*i+8] for i in range(len(bits)//8)]
       return "".join([chr(int(b,2)) for b in debytes])

   detext=bits_to_text(decode)

   print(Title)
   print(detext)

reconstruction(decurrsym_5_1,'4_d_vii_Reconstructed shannon_odd, SNR=5dB, el18007')
reconstruction(decurrsym_15_1,'4_d_vii_Reconstructed shannon_odd, SNR=15dB, el18007')
reconstruction(decurrsym_5_2,'4_d_vii_Reconstructed shannon_even, SNR=5dB, el18129')
reconstruction(decurrsym_15_2,'4_d_vii_Reconstructed shannon_even, SNR=15dB, el18129')

#Question 5

#el18007, 0+0+7=7 odd

#part a

#diabazoume to arxeio
audio=scipy.io.wavfile.read('soundfile1_lab2.wav')
audio2=[]
audio2=audio[1] #pinakas me tis lhftheises times
t1=np.arange(0,len(audio2),1)
plt.figure(figsize=(20,5))
plt.plot(t1,audio2)
plt.grid()
plt.xlabel("time")
plt.ylabel("voltage")
plt.title('Sounfile1 el18007')
plt.savefig('5_a_Soundfile1 (el18007)')

#part b

L=2**8
max1=max(audio2)-min(audio2)
D=max1/L

#synarthsh omoiomorfou kbantisth
def midriser(x):
  if(x>9985):
    return D*(np.floor(x/D)-1/2)
  else:
    return D*(np.floor(x/D)+1/2)
j=midriser(min(audio2))

#orismos epipedwn kbantishs
ql=[]
for i in range(L):
  ql.append(np.round(j,6))
  j+=D

#antistoixhsh shmatos sta epipeda kbantishs
audio_quantized=list(map(midriser,audio2))
for i in range(len(audio_quantized)):
  audio_quantized[i]=np.round(audio_quantized[i],6)


plt.clf()
plt.figure(figsize=(20, 5))
plt.title("midrised")
plt.xlabel("time(s)")
plt.ylabel("midriser output")
plt.step(t1,audio_quantized)
plt.grid()
plt.title('Quantized Signal el18007')
plt.savefig('5_b_YQ1 el18007')

#part c

#synarthsh metatrophs kbantismenou shmatos se kwdika gray
def q2g(qsig,gv,ql):
    a=[]
    for i in range(len(qsig)):
     x=ql.index(qsig[i])
     a.append(gv[x])
    return a

#energeies symbolwn kai bits
Es = (1**2)
Eb = Es/2
Tb=0.5

#shmeia asterismou
z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)

gv=GrayCode(8)
gv=list(gv.generate_gray())

#metatroph kbantismenou shmatos se kwdika gray
audiogray=q2g(audio_quantized,gv,ql)
bs=''.join(audiogray)

#antistoixisi bits se symbola
currsym=[]
for i in range(0, len(bs), 2):
  if (bs[i]=='0' and bs[i+1]=='0'):
    currsym.append("00")
  elif (bs[i]=='0' and bs[i+1]=='1'):
    currsym.append("01")
  elif (bs[i]=='1' and bs[i+1]=='0'):
    currsym.append("10")
  else:
    currsym.append("11")  

#antistoixisi symbolwn se asterismo
qpsk_conste=[]
for i in range(len(currsym)):
  if (currsym[i]=="00"):
    qpsk_conste.append(z00)

  elif (currsym[i]=="01"):
    qpsk_conste.append(z01)
      
  elif (currsym[i]=="10"):
    qpsk_conste.append(z10)   

  else:
    qpsk_conste.append(z11)

#part d

#synarthsh prosthikis thorybou sto shma
def add_noise_to_signal(signal, SNR, Es):
    No = Es/(10**(SNR/10))
    noise=np.random.normal(0,np.sqrt(No/2), size=len(signal))+1j*np.random.normal(0,np.sqrt(No/2), size=len(signal))
    r = noise + signal
    
    return r

#4 dB
qpsk_conste1=add_noise_to_signal(qpsk_conste,4,Es)
qpsk_consteAWGN4=[]
#bhma 2 giati allios den to trexei
for i in range(0, len(qpsk_conste1), 2):
	qpsk_consteAWGN4.append(qpsk_conste1[i])

#14 dB
qpsk_conste2=add_noise_to_signal(qpsk_conste,14,Es)
qpsk_consteAWGN14=[]
#bhma 2 giati allios den to trexei
for i in range(0, len(qpsk_conste2), 2):
	qpsk_consteAWGN14.append(qpsk_conste2[i])

#part e
#4dB
plt.clf()
plt.grid()
plt.figure(figsize=(10,10))
plt.scatter(np.real(qpsk_consteAWGN4),np.imag(qpsk_consteAWGN4))
plt.scatter(z00.real,z00.imag,label="Trasmitted signal")
plt.scatter(z01.real,z01.imag,label="Trasmitted signal")
plt.scatter(z10.real,z10.imag,label="Trasmitted signal")
plt.scatter(z11.real,z11.imag,label="Trasmitted signal")
plt.axvline(color="black",label="Decision Threshold")
plt.axhline(color="black",label="Decision Threshold")
plt.xlabel("I")
plt.ylabel("Q")
plt.title('Constellation + Noise 4dB')
plt.savefig('5_d_Constellation+Noise 4dB')

#14dB
plt.clf()
plt.grid()
plt.figure(figsize=(10,10))
plt.scatter(np.real(qpsk_consteAWGN14),np.imag(qpsk_consteAWGN14))
plt.scatter(z00.real,z00.imag,label="Trasmitted signal")
plt.scatter(z01.real,z01.imag,label="Trasmitted signal")
plt.scatter(z10.real,z10.imag,label="Trasmitted signal")
plt.scatter(z11.real,z11.imag,label="Trasmitted signal")
plt.axvline(color="black",label="Decision Threshold")
plt.axhline(color="black",label="Decision Threshold")
plt.xlabel("I")
plt.ylabel("Q")
plt.title('Constellation + Noise ')
plt.savefig('5_d_Constellation+Noise 14dB')

#part st

#4 dB
counter = 0
for i in range(len(currsym)):
  if (((currsym[i]=="00") and ((qpsk_conste1[i].real<0) or (qpsk_conste1[i].imag>0))) or ((currsym[i]=="01") and ((qpsk_conste1[i].real<0) or (qpsk_conste1[i].imag<0))) or ((currsym[i]=="11") and ((qpsk_conste1[i].real>0) or (qpsk_conste1[i].imag<0))) or ((currsym[i]=="10") and ((qpsk_conste1[i].real>0) or (qpsk_conste1[i].imag>0)))):
    counter+=1
  if (((currsym[i]=="00") and ((qpsk_conste1[i].real<0) and (qpsk_conste1[i].imag>0))) or ((currsym[i]=="01") and ((qpsk_conste1[i].real<0) and (qpsk_conste1[i].imag<0))) or ((currsym[i]=="11") and ((qpsk_conste1[i].real>0) and (qpsk_conste1[i].imag<0))) or ((currsym[i]=="10") and ((qpsk_conste1[i].real>0) and (qpsk_conste1[i].imag>0)))):
    counter+=1
print("4 dB el18007")
print("percentage=", 0.5*counter*100/len(currsym),"%")


print("tpercentage=", 50*special.erfc(np.sqrt(0.5*10**((4/10)))), "%")

#14 dB

counter = 0
for i in range(len(currsym)):
  if (((currsym[i]=="00") and ((qpsk_conste2[i].real<0) or (qpsk_conste2[i].imag>0))) or ((currsym[i]=="01") and ((qpsk_conste2[i].real<0) or (qpsk_conste2[i].imag<0))) or ((currsym[i]=="11") and ((qpsk_conste2[i].real>0) or (qpsk_conste2[i].imag<0))) or ((currsym[i]=="10") and ((qpsk_conste2[i].real>0) or (qpsk_conste2[i].imag>0)))):
    counter+=1
  if (((currsym[i]=="00") and ((qpsk_conste2[i].real<0) and (qpsk_conste2[i].imag>0))) or ((currsym[i]=="01") and ((qpsk_conste2[i].real<0) and (qpsk_conste2[i].imag<0))) or ((currsym[i]=="11") and ((qpsk_conste2[i].real>0) and (qpsk_conste2[i].imag<0))) or ((currsym[i]=="10") and ((qpsk_conste2[i].real>0) and (qpsk_conste2[i].imag>0)))):
    counter+=1
     
print("14 dB el18007")
print("percentage=", 0.5*counter*100/len(currsym),"%")


print("tpercentage=", 50*special.erfc(np.sqrt(0.5*10**((14/10)))), "%")

#part z
#apodiamorfwsi

#antistoixisi simeion asterismou se symbola
decurrsym=[]
for i in range(0, len(qpsk_conste1), 1):
  if (qpsk_conste1[i].real>0 and qpsk_conste1[i].imag<0):
    decurrsym.append('00')
  elif (qpsk_conste1[i].real>0 and qpsk_conste1[i].imag>0):
    decurrsym.append('01')
  elif (qpsk_conste1[i].real<0 and qpsk_conste1[i].imag<0):
    decurrsym.append('10')
  else:
    decurrsym.append('11') 

string=[]
for i in range(0,len(decurrsym),4):
  a=''
  for j in range(4):
    a+=decurrsym[i+j]
  string.append(a)

degray2=[]
for i in range(len(string)):
  degray=gv.index(string[i])#vriskei th thesi tou stoixeiou ston pinaka 
  degray2.append(ql[degray])


data = np.array(degray2, dtype=np.int16)
scipy.io.wavfile.write("soundfile_result4_el18007.wav", 44100, data)

#apodiamorfwsi

#antistoixisi simeion asterismou se symbola
decurrsym=[]
for i in range(0, len(qpsk_conste2), 1):
  if (qpsk_conste2[i].real>0 and qpsk_conste2[i].imag<0):
    decurrsym.append('00')
  elif (qpsk_conste2[i].real>0 and qpsk_conste2[i].imag>0):
    decurrsym.append('01')
  elif (qpsk_conste2[i].real<0 and qpsk_conste2[i].imag<0):
    decurrsym.append('10')
  else:
    decurrsym.append('11') 

string=[]
for i in range(0,len(decurrsym),4):
  a=''
  for j in range(4):
    a+=decurrsym[i+j]
  string.append(a)

degray2=[]
for i in range(len(string)):
  degray=gv.index(string[i])#vriskei th thesi tou stoixeiou ston pinaka 
  degray2.append(ql[degray])


data = np.array(degray2, dtype=np.int16)
scipy.io.wavfile.write("soundfile_result14_el18007.wav", 44100, data)

#el18129, 1+2+9=12 even

#part a

#diabazoume to arxeio
audio=scipy.io.wavfile.read('soundfile2_lab2.wav')
audio2=[]
audio2=audio[1] #pinakas me tis lhftheises times
t1=np.arange(0,len(audio2),1)
plt.figure(figsize=(20,5))
plt.plot(t1,audio2)
plt.grid()
plt.xlabel("time")
plt.ylabel("voltage")
plt.title('Sounfile1 el18129')
plt.savefig('5_a_Soundfile2 (el18129)')

#part b

L=2**8
max1=max(audio2)-min(audio2)
D=max1/L

#synarthsh omoiomorfou kbantisth
def midriser(x):
  if(x>7925):
    return D*(np.floor(x/D)-1/2)
  else:
    return D*(np.floor(x/D)+1/2)
j=midriser(min(audio2))

#orismos epipedwn kbantishs
ql=[]
for i in range(L):
  ql.append(np.round(j,6))
  j+=D

#antistoixhsh shmatos sta epipeda kbantishs
audio_quantized=list(map(midriser,audio2))
for i in range(len(audio_quantized)):
  audio_quantized[i]=np.round(audio_quantized[i],6)


plt.clf()
plt.figure(figsize=(20, 5))
plt.title("midrised")
plt.xlabel("time(s)")
plt.ylabel("midriser output")
plt.step(t1,audio_quantized)
plt.grid()
plt.title('Quantized Signal el18129')
plt.savefig('5_b_YQ2 el18129')

#part c

#synarthsh metatrophs kbantismenou shmatos se kwdika gray
def q2g(qsig,gv,ql):
    a=[]
    for i in range(len(qsig)):
     x=ql.index(qsig[i])
     a.append(gv[x])
    return a

#energeies symbolwn kai bits
Es = (1**2)
Eb = Es/2
Tb=0.5

#shmeia asterismou
z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)

gv=GrayCode(8)
gv=list(gv.generate_gray())

#metatroph kbantismenou shmatos se kwdika gray
audiogray=q2g(audio_quantized,gv,ql)
bs=''.join(audiogray)

#antistoixisi bits se symbola
currsym=[]
for i in range(0, len(bs), 2):
  if (bs[i]=='0' and bs[i+1]=='0'):
    currsym.append("00")
  elif (bs[i]=='0' and bs[i+1]=='1'):
    currsym.append("01")
  elif (bs[i]=='1' and bs[i+1]=='0'):
    currsym.append("10")
  else:
    currsym.append("11")  

#antistoixisi symbolwn se asterismo
qpsk_conste=[]
for i in range(len(currsym)):
  if (currsym[i]=="00"):
    qpsk_conste.append(z00)

  elif (currsym[i]=="01"):
    qpsk_conste.append(z01)
      
  elif (currsym[i]=="10"):
    qpsk_conste.append(z10)   

  else:
    qpsk_conste.append(z11)

#part d

#synarthsh prosthikis thorybou sto shma
def add_noise_to_signal(signal, SNR, Es):
    No = Es/(10**(SNR/10))
    noise=np.random.normal(0,np.sqrt(No/2), size=len(signal))+1j*np.random.normal(0,np.sqrt(No/2), size=len(signal))
    r = noise + signal
    
    return r

#4dB
qpsk_conste1=add_noise_to_signal(qpsk_conste,4,Es)
qpsk_consteAWGN4=[]
#bhma 2 giati allios den to trexei
for i in range(0, len(qpsk_conste1), 2):
	qpsk_consteAWGN4.append(qpsk_conste1[i])

#14dB
qpsk_conste2=add_noise_to_signal(qpsk_conste,14,Es)
qpsk_consteAWGN14=[]
#bhma 2 giati allios den to trexei
for i in range(0, len(qpsk_conste2), 2):
	qpsk_consteAWGN14.append(qpsk_conste2[i])

#part e
#4dB
plt.clf()
plt.grid()
plt.figure(figsize=(10,10))
plt.scatter(np.real(qpsk_consteAWGN4),np.imag(qpsk_consteAWGN4))
plt.scatter(z00.real,z00.imag,label="Trasmitted signal")
plt.scatter(z01.real,z01.imag,label="Trasmitted signal")
plt.scatter(z10.real,z10.imag,label="Trasmitted signal")
plt.scatter(z11.real,z11.imag,label="Trasmitted signal")
plt.axvline(color="black",label="Decision Threshold")
plt.axhline(color="black",label="Decision Threshold")
plt.xlabel("I")
plt.ylabel("Q")
plt.title('Constellation + Noise 4dB')
plt.savefig('5_d_Constellation+Noise 4dB')

#14dB
plt.clf()
plt.grid()
plt.figure(figsize=(10,10))
plt.scatter(np.real(qpsk_consteAWGN14),np.imag(qpsk_consteAWGN14))
plt.scatter(z00.real,z00.imag,label="Trasmitted signal")
plt.scatter(z01.real,z01.imag,label="Trasmitted signal")
plt.scatter(z10.real,z10.imag,label="Trasmitted signal")
plt.scatter(z11.real,z11.imag,label="Trasmitted signal")
plt.axvline(color="black",label="Decision Threshold")
plt.axhline(color="black",label="Decision Threshold")
plt.xlabel("I")
plt.ylabel("Q")
plt.title('Constellation + Noise ')
plt.savefig('5_d_Constellation+Noise 14dB')

#part st

#4 dB
counter = 0
for i in range(len(currsym)):
  if (((currsym[i]=="00") and ((qpsk_conste1[i].real<0) or (qpsk_conste1[i].imag>0))) or ((currsym[i]=="01") and ((qpsk_conste1[i].real<0) or (qpsk_conste1[i].imag<0))) or ((currsym[i]=="11") and ((qpsk_conste1[i].real>0) or (qpsk_conste1[i].imag<0))) or ((currsym[i]=="10") and ((qpsk_conste1[i].real>0) or (qpsk_conste1[i].imag>0)))):
    counter+=1
  if (((currsym[i]=="00") and ((qpsk_conste1[i].real<0) and (qpsk_conste1[i].imag>0))) or ((currsym[i]=="01") and ((qpsk_conste1[i].real<0) and (qpsk_conste1[i].imag<0))) or ((currsym[i]=="11") and ((qpsk_conste1[i].real>0) and (qpsk_conste1[i].imag<0))) or ((currsym[i]=="10") and ((qpsk_conste1[i].real>0) and (qpsk_conste1[i].imag>0)))):
    counter+=1

print("4 dB el18129")
print("percentage=", 0.5*counter*100/len(currsym),"%")


print("tpercentage=", 50*special.erfc(np.sqrt(0.5*10**((4/10)))), "%")

#14 dB

counter = 0
for i in range(len(currsym)):
  if (((currsym[i]=="00") and ((qpsk_conste2[i].real<0) or (qpsk_conste2[i].imag>0))) or ((currsym[i]=="01") and ((qpsk_conste2[i].real<0) or (qpsk_conste2[i].imag<0))) or ((currsym[i]=="11") and ((qpsk_conste2[i].real>0) or (qpsk_conste2[i].imag<0))) or ((currsym[i]=="10") and ((qpsk_conste2[i].real>0) or (qpsk_conste2[i].imag>0)))):
    counter+=1
  if (((currsym[i]=="00") and ((qpsk_conste2[i].real<0) and (qpsk_conste2[i].imag>0))) or ((currsym[i]=="01") and ((qpsk_conste2[i].real<0) and (qpsk_conste2[i].imag<0))) or ((currsym[i]=="11") and ((qpsk_conste2[i].real>0) and (qpsk_conste2[i].imag<0))) or ((currsym[i]=="10") and ((qpsk_conste2[i].real>0) and (qpsk_conste2[i].imag>0)))):
    counter+=1
     
print("14 dB el18129")
print("percentage=", 0.5*counter*100/len(currsym),"%")


print("tpercentage=", 50*special.erfc(np.sqrt(0.5*10**((14/10)))), "%")

#part z
#apodiamorfwsi
#4dB

#antistoixisi simeion asterismou se symbola
decurrsym=[]
for i in range(0, len(qpsk_conste1), 1):
  if (qpsk_conste1[i].real>0 and qpsk_conste1[i].imag<0):
    decurrsym.append('00')
  elif (qpsk_conste1[i].real>0 and qpsk_conste1[i].imag>0):
    decurrsym.append('01')
  elif (qpsk_conste1[i].real<0 and qpsk_conste1[i].imag<0):
    decurrsym.append('10')
  else:
    decurrsym.append('11') 

string=[]
for i in range(0,len(decurrsym),4):
  a=''
  for j in range(4):
    a+=decurrsym[i+j]
  string.append(a)

degray2=[]
for i in range(len(string)):
  degray=gv.index(string[i])#vriskei th thesi tou stoixeiou ston pinaka 
  degray2.append(ql[degray])


data = np.array(degray2, dtype=np.int16)
scipy.io.wavfile.write("soundfile_result4_el18129.wav", 44100, data)

#apodiamorfwsi
#14dB
#antistoixisi simeion asterismou se symbola
decurrsym=[]
for i in range(0, len(qpsk_conste2), 1):
  if (qpsk_conste2[i].real>0 and qpsk_conste2[i].imag<0):
    decurrsym.append('00')
  elif (qpsk_conste2[i].real>0 and qpsk_conste2[i].imag>0):
    decurrsym.append('01')
  elif (qpsk_conste2[i].real<0 and qpsk_conste2[i].imag<0):
    decurrsym.append('10')
  else:
    decurrsym.append('11') 

string=[]
for i in range(0,len(decurrsym),4):
  a=''
  for j in range(4):
    a+=decurrsym[i+j]
  string.append(a)

degray2=[]
for i in range(len(string)):
  degray=gv.index(string[i])#vriskei th thesi tou stoixeiou ston pinaka 
  degray2.append(ql[degray])


data = np.array(degray2, dtype=np.int16)
scipy.io.wavfile.write("soundfile_result14.wav", 44100, data)