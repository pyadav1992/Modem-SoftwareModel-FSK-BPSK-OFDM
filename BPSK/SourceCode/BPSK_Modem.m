%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Name: Pratik Yadav
%%  Project: Software model for BPSK modem 
%%  Description: This project is presented as midterm exam for the course  
%%               Modem Design (EE-655)
%%               This is Part E of the Midterm: Complete BPSK Modem  
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

%% Part E 
%Modulator

% Basic Declarations
SymRate = 1;
Fs = 8;
TypeFlag = 'sqrt';
Alpha = 0.40;
Delay = 10; 
Window = kaiser(2048,8)';
Window = 10*(Window/sum(Window));

%Random Signal Generator (Bi-Polar)
X0 = (floor(2*rand(1,(2000)))-0.5)/0.5;

%Raised Cosine Filter
h = rcosine(SymRate, Fs, TypeFlag, Alpha, Delay);  
h1 = h/max(h);
h11 = reshape([h1 0 0 0 0 0 0 0],8,21);

reg = zeros(1,21);
m = 1;
for n = 1:2000
    reg = [X0(n), reg(1:20)];
    for k = 1 : 8
        BPSK(m) = reg*h11(k,:)';                       % BPSK Output 
        m = m + 1;
    end
end

b = [1 0 0 0 0 0 0.2 0 0 0 0 0 1j*0.005];              % channel

BPSK_Ch = filter(b,1,BPSK);
%Phase shift
BPSK_Ch_Offset = BPSK_Ch.*exp(1j*(2*pi/100)*(1:length(BPSK_Ch)));
%Downsample
BPSK_Ch_Offset_DownSample = BPSK_Ch_Offset(2:2:end);

%% Plots for Modulator
figure(11)
subplot(2,1,1)
plot(real(BPSK_Ch_Offset_DownSample(1:800)),'linewidth',1.5); 
hold on
plot([1:4:length(BPSK_Ch_Offset_DownSample(1:800))],real(BPSK_Ch_Offset_DownSample(1:4:800)),'ro','linewidth',1.5);
hold off
title('BPSK signal passed through channel')
xlabel('samples');
ylabel('Amplitude');
grid on;

subplot(2,1,2)
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(BPSK_Ch_Offset_DownSample(1:2048).*Window)))),'linewidth',1.5);
hold on
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(h1(1:2:end)/sum(h1(1:2:end)),2048)))),'r','linewidth',1.5);
hold off
axis([-0.5 0.5 -80 10]);
title('Frequency spectrum of BPSK signal passed through channel')
xlabel('Normalized frequency');
ylabel('20*log10(Mag)  dB');
grid on;
    
figure(12)
subplot(2,1,1)
plot(0,0)
% axis([-1 1 -1.2 1.2])
hold on;
for n = 161:8:8000-8
    plot(-1:1/4:1,real(BPSK_Ch_Offset_DownSample(n:n+8)))
end 
title('EYE Diagram for BPSK signal passed through channel')
hold off
grid on

subplot(2,1,2)
plot(0,0)
axis('equal')
axis([-2 2 -2 2])
hold on;
plot([-2 2],[0 0],'--k','linewidth',1.5)
plot([0 0],[-2 2],'--k','linewidth',1.5)
plot((BPSK_Ch_Offset_DownSample(4001:4:8000)),0,'rx')
hold off;
title('Constellation mapping for BPSK signal passed through channel')
xlabel('Real');
ylabel('Imaginary')
grid on

%% FLL
[ggPos,ggNeg] = band_edge_harris(4,Alpha,Delay);

theta=1/25;
eta = sqrt(2)/2;
eta=1*eta;
denom = (1+(2*eta*theta)+(theta*theta));
k_i = (4*theta*theta)/denom;
k_p = (4*eta*theta)/denom;

temp=zeros(1,81);
int=0;
loopAcc=0;
tempB = 0;
tempA = 0;
c = 0;
s = 0;
m1 = 1;
tempProd = zeros(1,length(BPSK_Ch_Offset_DownSample));
loopAcci = 0;

for i1 = 1:1:length(BPSK_Ch_Offset_DownSample)
   tempProd(i1) = BPSK_Ch_Offset_DownSample(i1).*exp(1j*2*pi*loopAcc);
   temp = [tempProd(i1) temp(1:80)];
   c = temp*ggPos';
   c_sv(m1) = c;
   s = temp*ggNeg';
   s_sv(m1) = s;
   
   err = abs(c)^2 - abs(s)^2;
   err_sv(m1) = err;
   loopAcci = loopAcci + err;
   loopAcci_sv(m1) = loopAcci;
   
   tempA=err*k_p;
   tempB=err*k_i+tempB;
   int = tempA+tempB;
   int_sv(m1) = int;
   
   loopAcc = loopAcc + int;
   loopAcc_sv(m1) = loopAcc;
   m1 = m1 + 1;
end

LoopOut = tempProd(2:2:end); %Input to filter

%Match Filter
h2 = rcosine(SymRate, 64, TypeFlag, Alpha, Delay);
h3 = h2/(h2(1:32:length(h2))*h1(1:4:length(h1))');
h33 = reshape([h3 zeros(1,(1312-length(h3)))],32,41);
%Best Arm
FLLout = filter(h33(18,:),1,LoopOut);

%% Plots for FLL
figure(21)
subplot(3,1,1)
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(tempProd(1:2048).*Window)))),'linewidth',1.5);
hold on
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(h1(1:2:end)/sum(h1(1:2:end)),2048)))),'r','linewidth',1.5);
hold off
axis([-0.5 0.5 -80 10]);
title([' Frequency spectrum of Complex BPSK time series (Loop enabled) for \alpha =' ,num2str(Alpha)])
xlabel('Normalized frequency');
ylabel('20*log10(Mag) dB');
grid on;

subplot(3,1,2)
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(c_sv(1:2048).*Window)))),'linewidth',1.5);
hold on
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(ggNeg,2048)))),'r','linewidth',1.5);
hold off
axis([-0.5 0.5 -80 10]);
title(' Frequency spectrum of negetive band edge filter')
xlabel('Normalized frequency');
ylabel('20*log10(Mag) dB');
grid on;

subplot(3,1,3)
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(s_sv(1:2048).*Window)))),'linewidth',1.5);
hold on
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(ggPos,2048)))),'r','linewidth',1.5);
hold off
axis([-0.5 0.5 -80 10]);
title(' Frequency spectrum of positive band edge filter')
xlabel('Normalized frequency');
ylabel('20*log10(Mag) dB');
grid on;

figure(22)
subplot(1,3,1)
plot(0,0)
axis('equal')
axis([-2 2 -2 2])
hold on;
plot([-2 2],[0 0],'--k','linewidth',1.5)
plot([0 0],[-2 2],'--k','linewidth',1.5)
plot((BPSK_Ch_Offset_DownSample(2000:4:8000)),'rx','linewidth',1.2)
hold off;
title('FLL input time series Constellation Mapping')
xlabel('Real');
ylabel('Imaginary')
hold off
grid on

subplot(1,3,2)
plot(0,0)
axis('equal')
axis([-2 2 -2 2])
hold on;
plot([-2 2],[0 0],'--k','linewidth',1.5)
plot([0 0],[-2 2],'--k','linewidth',1.5)
plot((LoopOut(1000:2:4000)),'rx','linewidth',1.2)
hold off;
title('Poly-phase match filter input time series Constellation Mapping')
xlabel('Real');
ylabel('Imaginary')
hold off
grid on

subplot(1,3,3)
plot(0,0)
axis('equal')
axis([-2 2 -2 2])
hold on;
plot([-2 2],[0 0],'--k','linewidth',1.5)
plot([0 0],[-2 2],'--k','linewidth',1.5)
plot((FLLout(1000:2:4000)),'rx','linewidth',1.2)
hold off;
title('Poly-phase match filter output(Best Arm) time series Constellation Mapping')
xlabel('Real');
ylabel('Imaginary')
hold off
grid on

figure(23)
subplot(2,1,1)
plot(20*log10(abs((err_sv))),'linewidth',1.5)
hold on
plot(20*log10(abs((int_sv))),'r','linewidth',1.5)
hold off
xlabel('samples')
ylabel('20*log10(Mag) dB')
title('Input and output to Carrier loop filter')
grid on;
axis([1 8000 -400 5]);
legend('Input to carrier loop filter','Output of carrier loop filter')

subplot(2,1,2)
plot(abs(loopAcci_sv),'linewidth',1.5)
hold on
plot(abs(loopAcc_sv),'r','linewidth',1.5)
hold off
grid on
axis([0 8000 -10 200])
title('Phase profile for FLL')
xlabel('Samples')
ylabel('Phase')
legend('Input phase','Output phase','location','NorthWest')

figure(24)
subplot(3,1,1)
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(BPSK_Ch_Offset_DownSample(1:2048).*Window)))),'linewidth',1.5);
hold on
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(h1(1:2:end)/sum(h1(1:2:end)),2048)))),'r','linewidth',1.5);
hold off
axis([-0.5 0.5 -80 10]);
title([' Frequency spectrum of Complex BPSK time series (Loop disabled) for \alpha =' ,num2str(Alpha)])
xlabel('Normalized frequency');
ylabel('20*log10(Mag) dB');
grid on;

subplot(3,1,2)
y2 = filter(ggNeg,1,BPSK_Ch_Offset_DownSample);
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(y2(1:2048).*Window)))),'linewidth',1.5);
hold on
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(ggNeg,2048)))),'r','linewidth',1.5);
hold off
axis([-0.5 0.5 -80 10]);
title(' Frequency spectrum of negetive band edge filter')
xlabel('Normalized frequency');
ylabel('20*log10(Mag) dB');
grid on;

subplot(3,1,3)
y1 = filter(ggPos,1,BPSK_Ch_Offset_DownSample);
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(y1(1:2048).*Window)))),'linewidth',1.5);
hold on
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(ggPos,2048)))),'r','linewidth',1.5);
hold off
axis([-0.5 0.5 -80 10]);
title(' Frequency spectrum of positive band edge filter')
xlabel('Normalized frequency');
ylabel('20*log10(Mag) dB');
grid on;

%% Timing Recovery

%derivative Filter
h4 = conv(h3,[1 0 -1]*32);
h4 = h4(2:end-1);
h4 = h4/max(h4);
h44 = reshape([h4 zeros(1,1312-length(h4))],32,41);

theta = 2*pi/200;
eta = sqrt(2)/2;
eta = 7*eta;
denom = (1+(2*eta*theta)+(theta*theta));
k_i = (4*theta*theta)/denom;
k_p = (4*eta*theta)/denom;

tempA = 0;
tempB = 0;
int = 0;   
accum = 10;      
regT = zeros(1,41);
path_sv = zeros(1,2000);
prod_sv = zeros(1,2000);
accum_sv = zeros(1,2000);
int_sv = zeros(1,2000);

prod = 0;

m = 1;
for n=1:2:4000
    regT=[LoopOut(n) regT(1:40)];  
    
    path=floor(real(accum));
    path_sv(m)=path;
    
    MF(n)=regT*h33(path,:)';
    DMF(n)=regT*h44(path,:)';
    
    regT=[LoopOut(n+1) regT(1:40)];
    MF(n+1)=regT*h33(path,:)';
    DMF(n+1)=regT*h44(path,:)';
    
    prod = real(MF(n+1))*real(DMF(n+1));
    prod_sv(m) = prod;
    
    tempA=prod*k_p;
    tempB=prod*k_i + tempB;
    int = tempA + tempB;
    int_sv(m) = int;
    
    accum=accum+int;
    accum_sv(m)=accum;
    
    if accum>=33
        accum=accum-1;
    elseif accum<=1
        accum=accum+1;
    end
    m=m+1;
end

%% Plots for timing recovery
figure(31)
subplot(2,2,1)
plot(h3,'linewidth',1.5);
grid on
axis([0 1281 -0.3 1.1])
ylabel('Amplitude')
xlabel('Samples')
title('Impulse response of prototype poly-phase match filter')

subplot(2,2,3)
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(h3/sum(h3),2048)))),'linewidth',1.5);
grid on
axis([-0.2 0.2 -80 10]);
title(' Frequency spectrum of prototype poly-phase match filter')
xlabel('Normalized frequency');
ylabel('20*log10(Mag) dB');

subplot(2,2,2)
plot(h4,'linewidth',1.5);
grid on
axis([0 1281 -1.1 1.1])
ylabel('Amplitude')
xlabel('Samples')
title('Impulse response of prototype poly-phase derivative filter')

subplot(2,2,4)
plot((-0.5:(1/2048):0.5-(1/2048)),fftshift(20*log10(abs(fft(h4/80,2048)))),'linewidth',1.5);
grid on
axis([-0.2 0.2 -80 10]);
title(' Frequency Spectrum of prototype poly-phase derivative filter')
xlabel('Normalized frequency');
ylabel('20*log10(Mag) dB');

figure(32)
for k=1:32   
    phs=unwrap(fftshift(angle(fft(h33(k,:),1000))))/(2*pi);
    dphs_x=-conv(phs,[1 0 -1]*1000/2);
    dphs=dphs_x(2:1001);    
    hold on;
    plot(-0.5:1/1000:0.5-1/1000,dphs,'color',rand(1,3),'linewidth',2);  
    grid on;
    axis([-0.4 0.4 16 23]);
    title('Group Delay Response of the poly-phase match filter arms ');
    xlabel('Normalized frequency');
    ylabel('Delay (Samples)');    
end
hold off;

figure(33)
plot(path_sv,'linewidth',1.1);
hold on;
plot(accum_sv,'r','linewidth',1.1)
hold off;
axis([0 2000 7 25])
grid on
xlabel('Samples')
ylabel('Filter arm')
title('Timing recovery loop poly-phase filter pointer')
legend('Integer part of phase accumulator','Phase accumulator')

figure(34)
subplot(2,1,1)
plot(prod_sv)
grid on
axis([0 2000 -2.5 2.5])
xlabel('Samples')
ylabel('Amplitude')
title('Input to Loop filter')

subplot(2,1,2)
plot(int_sv)
grid on
axis([0 2000 -2.6 2.6])
xlabel('Samples')
ylabel('Amplitude')
title('Output to Loop filter')

figure(35)
subplot(2,1,1)
plot(0,0)
axis('equal')
axis([-2 2 -2 2])
hold on;
plot([-2 2],[0 0],'--k','linewidth',1.5)
plot([0 0],[-2 2],'--k','linewidth',1.5)
plot((LoopOut(3500:4:4000)),'rx','linewidth',1.5)
hold off;
title('Timing recovery input time series Constellation Mapping')
xlabel('Real');
ylabel('Imaginary')
hold off
grid on

subplot(2,1,2)
plot(0,0)
axis('equal')
axis([-2 2 -2 2])
hold on;
plot([-2 2],[0 0],'--k','linewidth',1.5)
plot([0 0],[-2 2],'--k','linewidth',1.5)
plot((MF(3500:4:4000)),'rx','linewidth',1.5)
hold off;
title('Timing recovery output time series Constellation Mapping')
xlabel('Real');
ylabel('Imaginary')
hold off
grid on
%% Equalizer / Timing Aquisition

%CarrierLoop and DDS
theta=pi/5000;
eta = sqrt(2)/2;
eta=7*eta;
denom = (1+(2*eta*theta)+(theta*theta));
k_i = (4*theta*theta)/denom;
k_p = (4*eta*theta)/denom;
accum = 0;
tempA = 0;
tempB = 0;
int = 0;
m = 1;

% LMS and Equalizer
W = zeros(1,40);
W(11) = 1;
reg1 = zeros(1,40);
mu = 0.01;

%Downsample
Flag = 0;

BPSK_Match = MF;

for i = 1:length(BPSK_Match)

    prod1(i) = BPSK_Match(i)*exp(-1j*2*pi*accum);
    reg1 = [prod1(i) reg1(1:39)];
    Ucap = reg1*W';
    Ucap_sv(i) = Ucap;
    
    if(Flag == 0)
        DetReal = real(Ucap);
        DetOut = 1;
        if (DetReal < 0)
            DetOut = -1;
        end
        DetImg = 0;
        Detector = DetOut + 1j*DetImg;
        errNeg =  Detector - Ucap;
        errNeg_sv(m) = errNeg;     
       
        W = W + reg1*mu*conj(errNeg);
        prod2 = conj(Ucap)*Detector;

        phi = -angle(prod2)/2*pi;
        phi1_sv(m) = phi;
        
        tempA=phi*k_p;
        tempB=phi*k_i+tempB;
        int = tempA+tempB;
        int1_sv(m) = int;
        
        m = m + 1;
        Flag = Flag + 1;
    else
        Flag = Flag + 1;
        
        if Flag == 2
            Flag = 0;
        end
    end
    accum = accum + int;
    accum_sv1(i) = accum; 
end
%% Plots for timing aqisition

figure(51)
plot(20*log10(abs((errNeg_sv))),'linewidth',1.2)
grid on;
xlabel('samples')
ylabel('20*log10(Mag)  dB')
title('Learning curve of equalizer')
axis([1 2000 -80 5 ]);

figure(52)
subplot(2,1,1)
plot(phi1_sv,'linewidth',1.5)
grid on
axis([0 2000 -3 3])
xlabel('Samples')
ylabel('Amplitude')
title('Input time series of PLL loop')

subplot(2,1,2)
plot(int1_sv,'linewidth',1.5)
grid on
axis([0 2000 -0.04 0.04])
xlabel('Samples')
ylabel('Amplitude')
title('Output time series of PLL loop')

figure(53)
subplot(2,1,1)
plot(accum_sv1,'linewidth',1.5)
title('Input time series to equalizer constellation mapping')
xlabel('Samples');
ylabel('Amplitude')
grid on

subplot(2,2,3)
plot(0,0)
axis('equal')
axis([-2 2 -2 2])
hold on;
plot([-2 2],[0 0],'--k','linewidth',1.5)
plot([0 0],[-2 2],'--k','linewidth',1.5)
plot((BPSK_Match(3001:2:4000)),'rx','linewidth',1.5)
hold off;
title('Input time series to timing aquisition constellation mapping')
xlabel('Real');
ylabel('Imaginary')
grid on

subplot(2,2,4)
plot(0,0)
axis('equal')
axis([-2 2 -2 2])
hold on;
plot([-2 2],[0 0],'--k','linewidth',1.5)
plot([0 0],[-2 2],'--k','linewidth',1.5)
plot((Ucap_sv(3901:2:4000)),'rx','linewidth',1.5)
hold off;
title('Output time series of timing aquisition constellation mapping')
xlabel('Real');
ylabel('Imaginary')
grid on