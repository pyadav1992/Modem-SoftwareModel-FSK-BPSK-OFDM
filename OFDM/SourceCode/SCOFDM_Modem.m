%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Name: Pratik Yadav
%%  Project: Software model for SCOFDM modem 
%%  Description: This project is presented as a final exam for the course  
%%               Modem Design (EE-655)
%%               This is Part 1 of the Final: Single Carrier OFDM Modem  
%%                 
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  
%%  If not, see <https://www.gnu.org/licenses/>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%% Part A
% Implement an OFDM Modem as a 128 point transform with 65 occupied
% frequency bins as shown above.In particular, bins –27-to-1 and +1 to +27,
% (bin 0 empty). Each bin is modulated 16-QAM. The time series formed at the 
% output of the transform is separated by a guard band of 32 zero valued 
% samples located at the beginning of the OFDM symbol. Form an OFDM packet 
% composed of 50 OFDM symbol.

symbols = 50;

SCOFDMcyc = [];
for k = 1:symbols
    cp0=(floor(4*rand(1,64))-1.5)/1.5+1j*(floor(4*rand(1,64))-1.5)/1.5;
    fcp0=(fft(cp0));
    fcp1=zeros(1,128);
    fcp1(65+(-32:31))=fcp0;

    cp1=2*ifft(fftshift(fcp1));
    SCOFDMcyc=[SCOFDMcyc zeros(1,32) cp1];
end

%% Part B
% Plot a power spectrum of the time series formed by a windowed 2048 point FFT.
% Comment on the side lobe structure of the spectrum.

window = kaiser(2048,8)';
window = window/sum(window);
window = 30*window;

figure(1)
subplot(2,1,1)
plot(0:999,real(SCOFDMcyc(1:1000)),'linewidth',2)
axis([0 999 -2 2])
title('Real part of the time series SC-OFDM')
xlabel('Samples')
ylabel('amplitude')
grid on 
    
subplot(2,1,2)
plot((-0.5:1/2048:0.50-1/2048),fftshift(20*log10(abs(fft(SCOFDMcyc(1:2048).*window)))),'linewidth',1.5)
grid on
axis([-0.5 0.5 -50 10])
title('Power spectrum of the time series SC-OFDM')
xlabel('Normalized frequency')
ylabel('20*log10(Mag)')

%% Part C

% Demodulate the OFDM packet and display the (time) overlaid real part of 
% the successive transforms as well as the (time and frequency) overlaid 
% constellation diagram when there is no channel distortion.

figure(2)
subplot(2,1,1)
plot(0,0)
hold on
for k=1:symbols
    v2=SCOFDMcyc((33:160)+(k-1)*160)/2;
    ff1 = fft(v2);
    ff2 = fftshift(ff1);
    ff3 = fftshift(ff2((-32:31)+65));
    ff4 = ifft(ff3);
    plot(-32:31,real(ff4),'o','linewidth',1.5) 
end
hold off
grid on
title('Overlay 50 SC-OFDM Symbols Real Part of demodulated signal without Channel, Without Cyclic Prefix')
% axis([0 130 -1.5 1.5])

subplot(2,3,4)
plot(0,0)
hold on
for k=1:symbols
    v2=SCOFDMcyc((33:160)+(k-1)*160)/2;
    ff1 = fft(v2);
    ff2 = fftshift(ff1);
    ff3 = fftshift(ff2((-32:31)+65));
    ff4 = ifft(ff3);
    plot(ff4,'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, 64 Bins')
axis('equal')
axis([-1.5 1.5 -1.5 1.5])

subplot(2,3,5)
plot(0,0)
hold on
for k=1:symbols
    v2=SCOFDMcyc((33:160)+(k-1)*160)/2;
    ff1 = fft(v2);
    ff2 = fftshift(ff1);
    ff3 = fftshift(ff2((-32:31)+65));
    ff4 = ifft(ff3);
    plot(ff4(20),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, -20 Bin')
axis('equal')
axis([-1.5 1.5 -1.5 1.5])

subplot(2,3,6)
plot(0,0)
hold on
for k=1:symbols
    v2=SCOFDMcyc((33:160)+(k-1)*160)/2;
    ff1 = fft(v2);
    ff2 = fftshift(ff1);
    ff3 = fftshift(ff2((-32:31)+65));
    ff4 = ifft(ff3);
    plot(ff4(44),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, +20 Bin')
axis('equal')
axis([-1.5 1.5 -1.5 1.5])

%% Part D
% Pass the OFDM packet through a channel with impulse response 
% [1   0   j 0.2   0   0   0   0   0   0   0.1   ] Now demodulate the 
% received OFDM packet and display the (time) overlaid real part of the 
% successive transforms as well as the (time and frequency) overlaid 
% constellation diagram.

ch = [1 0 1j*0.2 0 0 0 0 0 0 0.1];
SCOFDMch = filter(ch,1,SCOFDMcyc);

figure(3)
subplot(2,1,1)
plot(0,0)
hold on
for k=1:symbols
    v2=SCOFDMch((33:160)+(k-1)*160)/2;
    ff1 = fft(v2);
    ff2 = fftshift(ff1);
    ff3 = fftshift(ff2((-32:31)+65));
    ff4 = ifft(ff3);
    plot(-32:31,real(ff4),'o','linewidth',1.5) 
end
hold off
grid on
title('Overlay 50 SC-OFDM Symbols Real Part of demodulated signal with Channel, Without Cyclic Prefix')
% axis([0 130 -1.5 1.5])

subplot(2,1,2)
plot(0,0)
hold on
for k=1:symbols
    v2=SCOFDMch((33:160)+(k-1)*160)/2;
    ff1 = fft(v2);
    ff2 = fftshift(ff1);
    ff3 = fftshift(ff2((-32:31)+65));
    ff4 = ifft(ff3);
    plot(ff4,'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, 64 Bins')
axis('equal')
axis([-1.5 1.5 -1.5 1.5])

%% Part E
% Replace the 32-sample guard band with a 32 sample cyclic prefix taken 
% from the last 32 samples of each OFDM symbol. Pass this new packet through 
% the same channel. Now demodulate the received OFDM packet and display the 
% (time) overlaid real part of the successive transforms as well as the 
% (time and frequency) overlaid constellation diagram.

SCOFDMcyc=[];
symbols = 50;

for k=1:symbols
    cp0=(floor(4*rand(1,64))-1.5)/1.5+j*(floor(4*rand(1,64))-1.5)/1.5;
    fcp0=(fft(cp0));    
    fcp1=zeros(1,128);
    fcp1(65+(-32:31))=fcp0;

    cp1=2*ifft(fftshift(fcp1));
    SCOFDMcyc=[SCOFDMcyc cp1(97:128) cp1];
end
ch=[1 0 j*0.2 0 0 0 0 0 0 0.1];
SCOFDMcyc_ch=filter(ch,1,SCOFDMcyc);

figure(4)

subplot(2,1,1)
plot(0,0)
hold on
for k=1:160:8000-160
    v2=SCOFDMcyc_ch(k:k+159)/2;
    ff1=v2(33:160);
    ff2=fftshift(fft(ff1));
    ff3=(fftshift(ff2((-32:31)+65)));
    ff4 = ifft(ff3);
    plot(-32:31,real(ff4),'o')
end
hold off
grid on
title('Overlay 50 SC-OFDM Symbols Real Part of demodulated signal with Channel, With Cyclic Prefix')

subplot(2,1,2)
plot(0,0)
hold on
for k=1:160:8000-160
v2=SCOFDMcyc_ch(k:k+159)/2;
ff1=v2(33:160);
ff2=fftshift(fft(ff1));
ff3=(fftshift(ff2((-32:31)+65)));
ff4 = ifft(ff3);
plot((ff4),'o')
end
hold off
axis('equal')
title('Constellations, 64 Bins')
axis('equal')
axis([-1.5 1.5 -1.5 1.5])

%% Part F
% Insert 50 samples of delay and a single OFDM symbol length preamble at 
% the beginning of the OFDM packet formed with 4-repeated time segments. 
% Do this with a random qpsk time sequence of length 16, repeated 4-times. 
% Previously we obtained the periodicity by zero packing the spectrum but 
% here we simply repeat the segment. The modulation process interpolates 
% these segments from 16 samples per segment to 32 samples per segment. 
% This symbol has an appended cyclic prefix. Add white noise with standard
% deviation 0.01 to the received signal. Spin the data being delivered 
% through the channel at 1 degree per sample. Use a 32 sample delay line 
% and form the cross correlation between input and output of the delay line
% and form the autocorrelation at the output of the delay line. Verify that
% the sample by sample ratio of the cross to the auto correlation identifies
% the start of the OFDM frame. Present a figure showing magnitude of cross, auto,
% and ratio. Also show the angle obtained at the output of the 32 sample averager
% processing the cross correlation. Verify this is a reasonable estimate of the 
% input rotation rate. Present a figure showing output of 32-point averager
% and the estimated rate of rotation. Use this estimated rotation rate to 
% de-spin subsequent signal samples following this short preamble. 
% For this problem, we will not use the estimated rotation rate to 
% guide the de-spinning. This was included to illustrate the process. 
% For the remainder of this problem, remove the input spin but not the input noise.

QPSK = (floor(2*rand(1,16))-0.5)/0.5 + 1j*(floor(2*rand(1,16))-0.5)/0.5;

temp1 = zeros(1,64);
temp1(1:16) = QPSK;
temp1(17:32) = QPSK;
temp1(33:48) = QPSK;
temp1(49:64) = QPSK;

v1 = fft(temp1);
v2 = zeros(1,128);
v2(65+(-32:31))=v1;
v3 = fftshift(v2);
v4 = 2*ifft(v3);

ShPreamble = [v4(97:end) v4];

SCOFDMcyc = [];
for i = 1:symbols
    cp0=(floor(4*rand(1,64))-1.5)/0.5+j*(floor(4*rand(1,64))-1.5)/0.5;
    fcp0=fft(cp0);   
    fcp1=zeros(1,128);
    fcp1(65+(-32:31))=fftshift(fcp0);
  
    cp1=2*ifft(fftshift(fcp1));
    SCOFDMcyc =[SCOFDMcyc cp1(97:end) cp1];
end

sigma = 0.01;
SCOFDMrec = [zeros(1,50) ShPreamble SCOFDMcyc];
SCOFDMrec = sigma*(randn(1,length(SCOFDMrec))+1j*(randn(1,length(SCOFDMrec)))) + SCOFDMrec;
SCOFDMrec11 = SCOFDMrec.*exp(1j*2*pi*(1/360)*(1:8210));

DelayLine = zeros(1,33);
regCCr = zeros(1,32);
regACr = zeros(1,32);
ratio_sv = [];
avgCCR_sv = [];
avgACR_sv = [];

ShPreambleEx = SCOFDMrec11(51:210);

for n = 1 : length(ShPreambleEx)

    DelayLine = [ShPreambleEx(n) DelayLine(1:32)];
    CCr = DelayLine(1)*conj(DelayLine(33));
    ACr = DelayLine(33)*conj(DelayLine(33));
    
    regCCr = [CCr regCCr(1:31) ];
    regACr = [ACr regACr(1:31) ];
    
    avgCCR = sum(regCCr)/32;
    avgACR = sum(regACr)/32;
    
    avgCCR_sv = [ avgCCR_sv avgCCR];
    avgACR_sv = [ avgACR_sv avgACR];
    
    ratio = avgCCR/(avgACR+sigma);
    ratio_sv = [ratio_sv ratio];
end

Despin = (angle(avgCCR_sv)*360)/(2*pi*32);
SCOFDMrec1 = SCOFDMrec11(115:end);
SCOFDMdespin = SCOFDMrec1.*exp(j*2*pi*(-1/360)*(115:8210));

figure(5)
subplot(2,2,1)
plot(abs(avgACR_sv),'linewidth',1.5)
title('Auto correlation')
xlabel('Samples')
ylabel('Magnitude')
grid on

subplot(2,2,2)
plot(abs(avgCCR_sv),'linewidth',1.5)
title('Cross correlation')
xlabel('Samples')
ylabel('Magnitude')
grid on

subplot(2,2,3)
plot(abs(ratio_sv),'linewidth',1.5)
title('Ratio of Auto-correlation and Cross-correlation')
xlabel('Samples')
ylabel('Magnitude')
axis([0 160 0 1.1])
title('Ratio of autocorrelation and cross correlation')
grid on

subplot(2,2,4)
plot(abs(avgCCR_sv),'linewidth',1.5)
hold on
plot(abs(avgACR_sv),'ro')
hold off
grid on
title('Auto-correlation and Cross-correlation comparison')
xlabel('Samples')
ylabel('Magnitude')

figure(6)
plot(Despin,'linewidth',1.5)
title('Angle of cross correlation')
xlabel('Samples')
ylabel('Degrees')
grid on

%% Part G
% Form the long preamble used to estimate the channel frequency response.
% Since SC-OFDM forms the series in the time domain we have to form a signal 
% with the cascade transforms with approximately constant spectral amplitude. 
% We do this with a linear FM sweep known as CAZAC (constant amplitude, zero
% auto correlation) as shown below.  

NN=64;
rr=(1:NN).*(1:NN)/2;
prb0=exp(j*2*pi*rr/NN);
fprb0=fft(prb0);
fprb1=zeros(1,128);
fprb1(65+(-32:31))=fprb0;
prb1=ifft(fftshift(fprb1));
LnPreamble=[prb1(97:128) prb1];

figure(7);

subplot(3,1,1);
plot(real(LnPreamble),'linewidth',1.5);
title('Real part of long preamble time series');
xlabel('Samples')
ylabel('Amplitude')
grid on

subplot(3,1,2);
plot(imag(LnPreamble),'linewidth',1.5);
title('Imaginary part of long preamble time series');
xlabel('Samples')
ylabel('Amplitude')
grid on

subplot(3,1,3);
plot(-0.5:1/128:0.5-1/128,fftshift(abs(fft(LnPreamble(33:160),128))),'linewidth',1.5);
title('frequency response of long preamble');
xlabel('Normalized frequency')
ylabel('20*log10(Mag)  dB')
grid on

%% Part H
% Pass the long preamble through the channel as shown below
% xprb=filter(cc,1,prb1a);
% xprb=xprb+0.01*(randn(1,160)+j*randn(1,160));	
% On one of two subplots plot the FFT of the long preamble, 
% the FFT of the channel response to the long preamble, and the 
% FFT of the Channel. All responses are for series of length 128 without
% their cyclic prefix. On the second subplot plot the estimated channel
% frequency response and the actual frequency response. 

ch = [1 0 j*0.2 0 0 0 0 0 0 0.1];
SCOFDMcyc_Noise = filter(ch,1,LnPreamble);
SCOFDMcyc_Noise = SCOFDMcyc_Noise + sigma*((randn(1,160)+j*randn(1,160)));

LnPreambleRx=fftshift(((fft(SCOFDMcyc_Noise(33:160)))));
LnPreambleTx=fftshift(((fft(LnPreamble(33:160)))));
ChEst=zeros(1,128);
ChEst((-32:31)+65)=LnPreambleRx((-32:31)+65)./LnPreambleTx((-32:31)+65);

figure(8)

subplot(2,1,1);
plot((-0.5:1/128:0.5-1/128),fftshift((abs(fft(SCOFDMcyc_Noise(33:160),128)))),'k','linewidth',1.5);
hold on;
plot((-0.5:1/128:0.5-1/128),fftshift((abs(fft(LnPreamble(33:160),128)))),'linewidth',1.5);
plot((-0.5:1/128:0.5-1/128),fftshift((8*abs(fft(ch,128)))),'r','linewidth',1.5);
hold off;
title('Frequency response');
legend('FFT of Channel response to long preamble','FFT of the long preamble','FFT of the Channel');
xlabel('Normalized frequency');
ylabel('20*log10(Mag)  dB');
grid on;

subplot(2,1,2);
plot((-0.5:1/128:0.5-1/128),abs(ChEst),'linewidth',1.5);
hold on;
plot((-0.5:1/128:0.5-1/128),fftshift((abs(fft(ch,128)))),'r','linewidth',1.5);
hold off;
grid on;
title('Frequency response');
legend('FFT of Estimated Channel','FFT of Channel');
xlabel('Normalized frequency');
ylabel('20*log10(Mag)  dB');

%% Part I
% Pass the 50 SC-OFDM symbols through the channel and repeat demodulate 
% and channel equalize using the inverse channel of the channel estimate 
% formed in the previous section.

SCOFDMcyc = [];
for k=1:50
    cp0=(floor(4*rand(1,64))-1.5)/1.5+j*(floor(4*rand(1,64))-1.5)/1.5;
    fcp0=(fft(cp0));
    fcp1=zeros(1,128);
    fcp1(65+(-32:31))=fcp0;
    cp1=2*ifft(fftshift(fcp1));
    SCOFDMcyc=[SCOFDMcyc cp1(97:128) cp1];
end

SCOFDMrec = [LnPreamble SCOFDMcyc];
ch=[1 0 1j*0.2 0 0 0 0 0 0 0.1];
SCOFDMcyc_Noise=SCOFDMrec + sigma*(randn(1,length(SCOFDMrec))+1j*randn(1,length(SCOFDMrec)));
SCOFDMcyc_ChNoise=filter(ch,1,SCOFDMcyc_Noise);

figure(9)
subplot(2,1,1)
plot(0,0)
hold on

for k=161:160:8000-160   
    v2=SCOFDMcyc_ChNoise(k:k+159)/2;
    ff1 = v2(33:160);
    ff2 = fft(ff1);
    ff3 = fftshift(ff2)./ChEst;
    ff4 = ifft(fftshift(ff3((-32:31)+65)));
    plot(-32:31,real(ff4),'o','linewidth',1.5);
end

hold off
grid on
title('Real part overlaid,after channel estimation')

figure(9)
subplot(2,1,2)
plot(0,0)
hold on

for k=161:160:8000-160   
    v2=SCOFDMcyc_ChNoise(k:k+159)/2;
    ff1 = v2(33:160);
    ff2 = fft(ff1);
    ff3 = fftshift(ff2)./ChEst;
    ff4 = ifft(fftshift(ff3((-32:31)+65)));
    plot(ff4,'o','linewidth',1.5);
end
hold off
axis('equal')
axis([-1.5 1.5 -1.5 1.5])
grid on
title('Constellation Diagram for all bins,after channel estimation')
