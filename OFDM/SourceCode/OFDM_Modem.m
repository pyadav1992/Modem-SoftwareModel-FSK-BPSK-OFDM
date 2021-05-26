%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Name: Pratik Yadav
%%  Project: Software model for OFDM modem 
%%  Description: This project is presented as a final exam for the course  
%%               Modem Design (EE-655)
%%               This is Part 1 of the Final: OFDM Modem  
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

%% Part A and Part B
% Implement an OFDM Modem as a 128 point transform with 54 occupied frequency bins. 
% In particular, bins –27-to-1 and +1 to +27, (bin 0 empty). Each bin is modulated 16-QAM.
% The time series formed at the output of the transform is separated by a guard band of 32 
% zero valued samples located at the beginning of the OFDM symbol. Form an OFDM packet composed 
% of 50 OFDM symbol.

% Plot a power spectrum of the time series formed by a windowed 2048 point FFT.
% Comment on the side lobe structure of the spectrum.

% Random Signal Generator for 16 QAM
OFDM = [];
symbols = 50;

for k = 1:symbols
    QAM = (floor(4*rand(1,54))-1.5)/1.5 + 1j*(floor(4*rand(1,54))-1.5)/1.5;
    v1 = zeros(1,128);
    v1((-27:-1)+65) = QAM(1:27);
    v1((1:27)+65) = QAM(28:54);
    v1 = fftshift(v1);
    v2 = 10*ifft(v1);
    OFDM = [OFDM zeros(1,32) v2];
end    

window = kaiser(2048,8)';
window = window/sum(window);
window = 25*window;

figure(1)

subplot(2,1,1)
plot(0:999,real(OFDM(1:1000)),'linewidth',2)
title('Real part of the time series SC-OFDM')
xlabel('Samples')
ylabel('amplitude')
grid on 
    
subplot(2,1,2)
plot((-0.5:1/2048:0.50-1/2048)*128,fftshift(20*log10(abs(fft(OFDM(1:2048).*window)))),'linewidth',1.5)
grid on
axis([-64 64 -50 10])
title('Power spectrum of the time series SC-OFDM')
xlabel('Normalized frequency')
ylabel('20*log10(Mag)')

%% Part C
% Demodulate the OFDM packet and display the (time) overlaid real part of the successive 
% transforms as well as the (time and frequency) overlaid constellation diagram when there is 
% no channel distortion.

figure(2)

subplot(2,1,1)
plot(0,0)
hold on
for k=1:50
    v2=OFDM((33:160)+(k-1)*160)/10;
    ff1 = fft(v2);
    ff2=fftshift(ff1);
    ff3=ff2((-60:60)+65);
    plot(-60:60,real(ff3),'o','linewidth',1.5) 
end
hold off
grid on
title('Overlay 50 OFDM Symbols Real Part of demodulated signal without Channel, Without Cyclic Prefix')

subplot(2,3,4)
plot(0,0)
hold on
for k=1:50
    v2=OFDM((33:160)+(k-1)*160)/10;
    ff1 = fft(v2);
    ff2=fftshift(ff1);
    plot(ff2,'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, 54 Bins')
axis('equal')
axis([-1.5 1.5 -1.5 1.5])

subplot(2,3,5)
plot(0,0)
hold on
for k=1:50
    v2=OFDM((33:160)+(k-1)*160)/10;
    ff1 = fft(v2);
    ff2=fftshift(ff1);
    plot(ff2(45),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, -20 Bin')
axis('equal')
axis([-1.5 1.5 -1.5 1.5])

subplot(2,3,6)
plot(0,0)
hold on
for k=1:50
    v2=OFDM((33:160)+(k-1)*160)/10;
    ff1 = fft(v2);
    ff2=fftshift(ff1);
    plot(ff2(85),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, +20 Bin')
axis('equal')
axis([-1.5 1.5 -1.5 1.5])

%% Part D
% Pass the OFDM packet through a channel with impulse response 
% [1   0   j 0.2   0   0   0   0   0   0   0.1   ]
% Now demodulate the received OFDM packet and display the (time) overlaid 
% real part of the successive transforms as well as the (time and frequency) 
% overlaid constellation diagram.

ch = [1 0 1*j*0.2 0 0 0 0 0 0 0.1];
OFDMch = filter(ch,1,OFDM);

figure(3)

subplot(2,1,1)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMch((33:160)+(k-1)*160)/10;
    ff1 = fft(v2);
    ff2=fftshift(ff1);
    ff3=ff2((-60:60)+65);
    plot(-60:60,real(ff3),'o','linewidth',1.5) 
end
hold off
grid on
title('Overlay 50 OFDM Symbols Real Part of demodulated signal with Channel, Without Cyclic Prefix')

subplot(2,3,4)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMch((33:160)+(k-1)*160)/10;
    ff1=fft(v2);
    ff2=fftshift(ff1);  
    plot(ff2,'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, 54 Bins')
axis('equal')
axis([-2 2 -2 2])

subplot(2,3,5)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMch((33:160)+(k-1)*160)/10;
    ff1=fft(v2);
    ff2=fftshift(ff1);
    plot(ff2(45),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, -20 Bin')
axis('equal')
axis([-2 2 -2 2])

subplot(2,3,6)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMch((33:160)+(k-1)*160)/10;
    ff2=fftshift(fft(v2));
    plot(ff2(85),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, +20 Bin')
axis('equal')
axis([-2 2 -2 2])

%% Part E
% Replace the 32-sample guard band with a 32 sample cyclic prefix taken from
% the last 32 samples of each OFDM symbol. Pass this new packet through the same
% channel. Now demodulate the received OFDM packet and display the (time) overlaid
% real part of the successive transforms as well as the (time and frequency) 
% overlaid constellation diagram.

OFDMcyc = [];
for i = 1:symbols
    QAM = (floor(4*rand(1,54))-1.5)/1.5 + 1j*(floor(4*rand(1,54))-1.5)/1.5;
    v1 = zeros(1,128);
    v1((-27:-1)+65) = QAM(1:27);
    v1((1:27)+65) = QAM(28:54);
    v1 = fftshift(v1);
    v2 = 10*ifft(v1);
    OFDMcyc = [OFDMcyc v2(97:end) v2];
end    

OFDMcyc_ch = filter(ch,1,OFDMcyc);

figure(4)

subplot(2,1,1)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMcyc_ch((33:160)+(k-1)*160)/10;
    ff1 = fft(v2);
    ff2=fftshift(ff1);
    ff3=ff2((-60:60)+65);
    plot(-60:60,real(ff3),'o','linewidth',1.5) 
end
hold off
grid on
title('Overlay 50 OFDM Symbols Real Part of demodulated signal with Channel, With Cyclic Prefix')

subplot(2,3,4)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMcyc_ch((33:160)+(k-1)*160)/10;
    ff1=fft(v2);
    ff2=fftshift(ff1);  
    plot(ff2,'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, 54 Bins')
axis('equal')
axis([-2 2 -2 2])

subplot(2,3,5)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMcyc_ch((33:160)+(k-1)*160)/10;
    ff1=fft(v2);
    ff2=fftshift(ff1);
    plot(ff2(45),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, -20 Bin')
axis('equal')
axis([-2 2 -2 2])

subplot(2,3,6)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMcyc_ch((33:160)+(k-1)*160)/10;
    ff2=fftshift(fft(v2));
    plot(ff2(85),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, +20 Bin')
axis('equal')
axis([-2 2 -2 2])

%% Part F
% Insert 50 samples of delay and a single OFDM symbol length preamble at the beginning 
% of the OFDM packet formed by zero packing the spectrum 1-to-4. Only insert random QPSK 
% data in bins ?4,  ?8, ?12,…., ?24. This symbol has an appended cyclic prefix. 
% Add white noise with standard deviation 0.01 to the received signal. 
% Spin the data being delivered through the channel at 1 degree per sample. Use a 32 sample 
% delay line and form the cross correlation between input and output of the delay line and
% form the autocorrelation at the output of the delay line. Verify that the sample by sample
% ratio of the cross to the auto correlation identifies the start of the OFDM frame. 
% Present a figure showing magnitude of cross, auto, and ratio.
% Also show the angle obtained at the output of the 32 sample averager processing the
% cross correlation. Verify this is a reasonable estimate of the input rotation rate.
% Present a figure showing output of 32-point averager and the estimated rate of rotation.
% Use this estimated rotation rate to de-spin subsequent signal samples following this short
% preamble. For this problem, we will not use the estimated rotation rate to guide the de-spinning.
% This was included to illustrate the process. For the remainder of this problem, 
% remove the input spin but not the input noise.

QPSK = (floor(2*rand(1,12))-0.5)/0.5 + 1j*(floor(2*rand(1,12))-0.5)/0.5;
temp1 = zeros(1,128);
temp1((-24:4:-4)+65) = QPSK(1:6);
temp1((4:4:24)+65) = QPSK(7:12);
v1 = fftshift(temp1);
v2 = 20*ifft(v1);
ShPreamble = [v2(97:end) v2];

figure(5)

subplot(2,1,1)
plot(real(ShPreamble),'linewidth',1.5)
grid on
title('Short preamble time series')
xlabel('Samples')
ylabel('Amplitude')

subplot(2,1,2)
plot((-0.5:1/128:0.5-1/128),abs(temp1),'linewidth',1.5)
grid on
title('Short Preamble')
xlabel('Normalized frequency')
ylabel('Amplitude')

OFDMcyc = [];
for i = 1:symbols
    QAM = (floor(4*rand(1,54))-1.5)/1.5 + 1j*(floor(4*rand(1,54))-1.5)/1.5;
    v1 = zeros(1,128);
    v1((-27:-1)+65) = QAM(1:27);
    v1((1:27)+65) = QAM(28:54);
    v1 = fftshift(v1);
    v2 = 10*ifft(v1);
    OFDMcyc = [OFDMcyc v2(97:end) v2];
end  

sigma = 0.01;
OFDMrec = [zeros(1,50) ShPreamble OFDMcyc];
OFDMrec = sigma*(randn(1,length(OFDMrec))+1j*(randn(1,length(OFDMrec)))) + OFDMrec;
OFDMrec11 = OFDMrec.*exp(1j*2*pi*(1/360)*(1:8210));

DelayLine = zeros(1,33);
regCCr = zeros(1,32);
regACr = zeros(1,32);
ratio_sv = [];
avgCCR_sv = [];
avgACR_sv = [];

ShPreambleEx = OFDMrec11(51:210);

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
OFDMrec1 = OFDMrec11(211:end);
OFDMdespin = OFDMrec1.*exp(j*2*pi*(-1/360)*(211:8210));

figure(6)
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

figure(7)
plot(Despin,'linewidth',1.5)
title('Angle of cross correlation')
xlabel('Samples')
ylabel('Degrees')
grid on

figure(8)
subplot(2,1,1)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMdespin((33:160)+(k-1)*160)/10;
    ff1 = fft(v2);
    ff2=fftshift(ff1);
    ff3=ff2((-60:60)+65);
    plot(-60:60,real(ff3),'o','linewidth',1.5) 
end
hold off
grid on
title('Overlay 50 despined OFDM Symbols Real Part of demodulated signal with noise, With Cyclic Prefix')

subplot(2,3,4)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMdespin((33:160)+(k-1)*160)/10;
    ff1=fft(v2);
    ff2=fftshift(ff1);  
    plot(ff2,'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, 54 Bins')
axis('equal')
axis([-2 2 -2 2])

subplot(2,3,5)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMdespin((33:160)+(k-1)*160)/10;
    ff1=fft(v2);
    ff2=fftshift(ff1);
    plot(ff2(45),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, -20 Bin')
axis('equal')
axis([-2 2 -2 2])

subplot(2,3,6)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDMdespin((33:160)+(k-1)*160)/10;
    ff2=fftshift(fft(v2));
    plot(ff2(85),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, +20 Bin')
axis('equal')
axis([-2 2 -2 2])

%% Part G
% Insert a second preamble following the short preamble formed by QPSK 
% modulation in all the bins (-27 through +27 but not bin 0). This symbol has 
% an appended cyclic prefix as do all data symbols. Pass the ofdm packet with 
% two preambles through the channel and add the noise described earlier. When 
% the receiver receives this preamble symbol it demodulates it and uses the known
% spectral values to measure and correct the distortion caused by the channel. 
% The relationship between channel and received signal at each frequency is of 
% the form Frcvd(K)=Fxmtd(k)* Fchannel(k). From this relationship we can estimate 
% the channel and remove it in subsequent observations of the received signal. 
% From the channel estimate use the inverse of the channel at each frequency 
% to demodulate the data segment of the symbol frames and plot the (time) 
% overlaid real part of the successive transforms as well as the (time and frequency) 
% overlaid constellation diagram.

QPSK = (floor(2*rand(1,54))-0.5)/0.5 + 1j*(floor(2*rand(1,54))-0.5)/0.5;
temp2 = zeros(1,128);
temp2((-27:-1)+65) = QPSK(1:27);
temp2((1:27)+65) = QPSK(28:54);
v1 = fftshift(temp2);
v2 = 20*ifft(v1);
LnPreamble = [v2(97:end) v2];

figure(9)

plot(real(LnPreamble),'linewidth',1.5)
grid on
title('Long preamble time series')
xlabel('Samples')
ylabel('Amplitude')


OFDMcyc = [];
for i = 1:symbols
    QAM = (floor(4*rand(1,54))-1.5)/1.5 + 1j*(floor(4*rand(1,54))-1.5)/1.5;
    v1 = zeros(1,128);
    v1((-27:-1)+65) = QAM(1:27);
    v1((1:27)+65) = QAM(28:54);
    v1 = fftshift(v1);
    v2 = 10*ifft(v1);
    OFDMcyc = [OFDMcyc v2(97:end) v2];
end  

sigma = 0.01;
OFDMrec = [ShPreamble LnPreamble OFDMcyc];
OFDMcyc_Noise = sigma*(randn(1,length(OFDMrec))+1j*(randn(1,length(OFDMrec)))) + OFDMrec;
OFDMcyc_ChNoise = filter(ch,1,OFDMcyc_Noise);

LnPreambleTx = OFDMcyc_Noise(161:320);
LnPreambleRx = OFDMcyc_ChNoise(161:320);
ChEst= fftshift(fft(LnPreambleRx(33:160)))./fftshift(fft(LnPreambleTx(33:160)));

InvCh = 1./ChEst;
OFDM_Rx = OFDMcyc_ChNoise(321:end);

figure(10)
subplot(2,1,1)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDM_Rx((33:160)+(k-1)*160)/10;
    ff1 = fft(v2);
    ff2=fftshift(ff1).*InvCh;
    ff3=ff2((-60:60)+65);
    plot(-60:60,real(ff3),'o','linewidth',1.5) 
end
hold off
grid on
title('Overlay 50 OFDM Symbols(After channel estimation) Real Part of demodulated signal with noise, With Cyclic Prefix')

subplot(2,3,4)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDM_Rx((33:160)+(k-1)*160)/10;
    ff1=fft(v2);
    ff2=fftshift(ff1).*InvCh;  
    plot(ff2,'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, 54 Bins')
axis('equal')
axis([-2 2 -2 2])

subplot(2,3,5)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDM_Rx((33:160)+(k-1)*160)/10;
    ff1=fft(v2);
    ff2=fftshift(ff1).*InvCh;
    plot(ff2(45),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, -20 Bin')
axis('equal')
axis([-2 2 -2 2])

subplot(2,3,6)
plot(0,0)
hold on
for k=1:symbols
    v2=OFDM_Rx((33:160)+(k-1)*160)/10;
    ff2=fftshift(fft(v2)).*InvCh;
    plot(ff2(85),'o','linewidth',1.5) 
end
hold off
grid on
title('Constellations, +20 Bin')
axis('equal')
axis([-2 2 -2 2])