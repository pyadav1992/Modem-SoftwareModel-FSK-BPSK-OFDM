%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Name: Pratik Yadav
%%  Project: Software model for BPSK modem 
%%  Description: This project is presented as midterm exam for the course  
%%               Modem Design (EE-655)
%%               This is Part D of the Midterm   
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

%% PartD
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

b = [1 0 0 0 0 0 0.2 0 0 0 0 0 1j*0.05];
a = [1 0 0 0 0 0 0 0 0 0 0 0 0];

BPSK_Ch = filter(b,1,BPSK);
% Downsample
BPSK_Ch_DownSample = BPSK_Ch(2:4:length(BPSK));
% Offset
BPSK_Ch_DownSample_Offset = BPSK_Ch_DownSample.*(exp(1j*pi/10));

%Match Filter
h2 = rcosine(SymRate, 2, TypeFlag, Alpha, Delay);  
h3 = h2/(h2*h1(1:4:end)');
BPSK_Match = filter(h3,1,BPSK_Ch_DownSample_Offset);

%CarrierLoop and DDS
theta=pi/5000;
eta = sqrt(2)/2;
eta=1*eta;
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
W(20) = 1;
reg1 = zeros(1,40);
mu = 0.01;

%Downsample
Flag = 0;

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
        phi_sv(m) = phi;
        
        tempA=phi*k_p;
        tempB=phi*k_i+tempB;
        int = tempA+tempB;
        int_sv(m) = int;
        
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


figure(51)
plot(20*log10(abs((errNeg_sv))),'linewidth',1.2)
grid on;
xlabel('samples')
ylabel('20*log10(Mag)  dB')
title('Learning curve of equalizer')
axis([1 2000 -80 5 ]);


figure(52)
subplot(2,1,1)
plot(phi_sv,'linewidth',1.5)
grid on
axis([0 2000 -2 2])
xlabel('Samples')
ylabel('Amplitude')
title('Input time series of PLL loop')

subplot(2,1,2)
plot(int_sv,'linewidth',1.5)
grid on
axis([0 2000 -4*10^-3 4*10^-3])
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
axis([-1.3 1.3 -1.3 1.3])
hold on;
plot([-2 2],[0 0],'--k','linewidth',1.5)
plot([0 0],[-2 2],'--k','linewidth',1.5)
plot((BPSK_Match(3801:2:4000)),'rx','linewidth',1.5,'MarkerSize',10)
hold off;
title('Input time series to timing aquisition constellation mapping')
xlabel('Real');
ylabel('Imaginary')
grid on

subplot(2,2,4)
plot(0,0)
axis('equal')
axis([-1.3 1.3 -1.3 1.3])
hold on;
plot([-2 2],[0 0],'--k','linewidth',1.5)
plot([0 0],[-2 2],'--k','linewidth',1.5)
plot((Ucap_sv(3901:2:4000)),'rx','linewidth',1.5,'MarkerSize',12)
hold off;
title('Output time series of timing aquisition constellation mapping')
xlabel('Real');
ylabel('Imaginary')
grid on
