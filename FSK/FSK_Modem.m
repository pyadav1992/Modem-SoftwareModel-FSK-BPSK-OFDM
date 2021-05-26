%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Name: Pratik Yadav
%%  Project: Software model for FSK modem 
%%  Description: This project is presented as midterm exam for the course  
%%               Digital Signal Processing (EE-556)
%%               Bianry data is modulated using FSK techniques and 
%%               detected using:
%%               a) Band pass filters
%%               b) Complex sinusoids
%%               c) Derivative functions   
%%               d) Hilbert transform
%%               e) Arc tangent   
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

clc
clear all
close all

%% PART A
f0=2000; %Frequency in FSK to represent bit '0'
f1=6000; %Frequency in FSK to represent bit '1'
fs=40000; % sampling frequency

%Signals with f0=2KHz frequency
h1a=cos(2*pi*(f0/fs)*(0:19));
h1b=sin(2*pi*(f0/fs)*(0:19));
h1c=exp(+1*j*2*pi*(f0/fs)*(0:19));

%Signals with f1=6KHz frequency
h2a=cos(2*pi*(f1/fs)*(0:19));
h2b=sin(2*pi*(f1/fs)*(0:19));
h2c=exp(+j*2*pi*(f1/fs)*(0:19));

figure(1)

subplot(3,2,1)
plot(h1a,'linewidth',2)
title('H1A signal')
xlabel('Time Index')
ylabel('Amplitude')
grid on;
subplot(3,2,2)
plot((-0.5:1/400:0.5-1/400)*40,(abs(fftshift(fft(h1a,400)))),'linewidth',2)
hold on
stem(2,max(abs(fftshift(fft(h1a,400)))),'r-o','linewidth',2)
stem(-2,max(abs(fftshift(fft(h1a,400)))),'r-o','linewidth',2)
hold off
title('H1A signal Frequency Spectrum')
xlabel('Frequency')
ylabel('Magnitude')
grid on;

%h1b=sin(2*pi*(f0/fs)*(0:19));
subplot(3,2,3)
plot(h1b,'linewidth',2)
title('H1B signal')
xlabel('Time Index')
ylabel('Amplitude')
grid on;
subplot(3,2,4)
plot((-0.5:1/400:0.5-1/400)*40,(abs(fftshift(fft(h1b,400)))),'linewidth',2)
hold on
stem(2,max(abs(fftshift(fft(h1a,400)))),'r-o','linewidth',2)
stem(-2,max(abs(fftshift(fft(h1a,400)))),'r-o','linewidth',2)
hold off
grid on
title('H1B signal Frequency Spectrum')
xlabel('Frequency')
ylabel('Magnitude')

%h1c=exp(-1*j*2*pi*(f0/fs)*(0:19))
subplot(3,2,5)
plot(0:19,h1c,'linewidth',1.5)
title('H1C signal')
xlabel('Time Index')
ylabel('Amplitude')
grid on;

subplot(3,2,6)
plot((-0.5:1/400:0.5-1/400)*40,(abs(fftshift(fft(h1c,400)))),'linewidth',2)
hold on
stem(2,max(abs(fftshift(fft(h1c,400)))),'r-o','linewidth',2)
hold off
axis([-20 20 0 25]);
grid on;
title('H1C signal Frequency Spectrum')
xlabel('Frequency')
ylabel('Magnitude')

% h2a=cos(2*pi*(f1/fs)*(1:20));
figure(2)

subplot(3,2,1)
plot(h2a,'linewidth',2)
title('H2A signal')
xlabel('Time Index')
ylabel('Amplitude')
grid on;

subplot(3,2,2)
plot((-0.5:1/400:0.5-1/400)*40,(abs(fftshift(fft(h2a,400)))),'linewidth',2)
grid on;
hold on
stem(6,max(abs(fftshift(fft(h1a,400)))),'r-o','linewidth',2)
stem(-6,max(abs(fftshift(fft(h1a,400)))),'r-o','linewidth',2)
hold off
title('H2A signal Frequency Spectrum')
xlabel('Frequency')
ylabel('Magnitude')

% h2b=sin(2*pi*(f1/fs)*(1:20));
subplot(3,2,3)
plot(h2b,'linewidth',1.5)
title('H2B signal')
xlabel('Time Index')
ylabel('Amplitude')
grid on

subplot(3,2,4)
plot((-0.5:1/400:0.5-1/400)*40,(abs(fftshift(fft(h2b,400)))),'linewidth',2)
hold on
stem(6,max(abs(fftshift(fft(h2a,400)))),'r-o','linewidth',2)
stem(-6,max(abs(fftshift(fft(h2a,400)))),'r-o','linewidth',2)
hold off
grid on
title('H2B signal Frequency Spectrum')
xlabel('Frequency')
ylabel('Magnitude')

%h2c=exp(-j*2*pi*(f1/fs)*(1:20));
subplot(3,2,5)
plot(1:20,h2c,'linewidth',1.5)
grid on
title('H2C signal')
xlabel('Time Index')
ylabel('Amplitude')
subplot(3,2,6)
plot((-0.5:1/400:0.5-1/400)*40,(abs(fftshift(fft(h2c,400)))),'linewidth',2)
grid on
hold on
stem(6,max(abs(fftshift(fft(h2c,400)))),'r-o','linewidth',2)
hold off
title('H2C signal Frequency Spectrum')
xlabel('Frequency')
ylabel('Magnitude')


%% PART D (POLE-ZERO plot of H1A H1B H1C H2A H2B H2C)
figure(3)

subplot(3,2,1)
plot(exp(j*2*pi*(0:0.01:1)),'linewidth',2)
hold on
plot(0,0,'x r','linewidth',2,'Markersize',12)
plot(roots(h1a),'o r','linewidth',2)
hold off
grid on
axis('square')
axis([-1.1 1.1 -1.1 1.1])
title('Pole-Zero Plot of H1A')
xlabel('Real Axis')
ylabel('Imaginary Axis')

subplot(3,2,3)
plot(exp(j*2*pi*(0:0.01:1)),'linewidth',2)
hold on
plot(0,0,'x r','linewidth',2,'Markersize',12)
plot(roots(h1b),'o r','linewidth',2)
hold off
grid on
axis('square')
axis([-1.1 1.1 -1.1 1.1])
title('Pole-Zero Plot of H1B')
xlabel('Real Axis')
ylabel('Imaginary Axis')

subplot(3,2,5)
plot(exp(j*2*pi*(0:0.01:1)),'linewidth',2)
hold on
plot(0,0,'x r','linewidth',2,'Markersize',12)
plot(roots(h1c),'o r','linewidth',2)
hold off
grid on
axis('square')
axis([-1.1 1.1 -1.1 1.1])
title('Pole-Zero Plot of H1C')
xlabel('Real Axis')
ylabel('Imaginary Axis')

subplot(3,2,2)
plot(exp(j*2*pi*(0:0.01:1)),'linewidth',2)
hold on
plot(0,0,'x r','linewidth',2,'Markersize',12)
plot(roots(h2a),'o r','linewidth',2)
hold off
grid on
axis('square')
axis([-1.1 1.1 -1.1 1.1])
title('Pole-Zero Plot of H2A')
xlabel('Real Axis')
ylabel('Imaginary Axis')

subplot(3,2,4)
plot(exp(j*2*pi*(0:0.01:1)),'linewidth',1.5)
hold on
plot(0,0,'x r','linewidth',2,'Markersize',12)
plot(roots(h2b),'o r','linewidth',2,'Markersize',12)
hold off
grid on
axis('square')
axis([-1.1 1.1 -1.1 1.1])
title('Pole-Zero Plot of H2B')
xlabel('Real Axis')
ylabel('Imaginary Axis')

subplot(3,2,6)
plot(exp(j*2*pi*(0:0.01:1)),'linewidth',2)
hold on
plot(0,0,'x r','linewidth',2,'Markersize',12)
plot(roots(h2c),'o r','linewidth',2)
hold off
grid on
axis('square')
axis([-1.1 1.1 -1.1 1.1])
title('Pole-Zero Plot of H2C')
xlabel('Real Axis')
ylabel('Imaginary Axis')

%% PART E (FSK Modulation)
dec=abs('yadavpratikvijay'); %Binary conversion of data Name:'yadav Pratik Vijay'
bin=dec2bin(dec);
data=reshape(bin',1,16*7);
FSK_MOD=[]                  % Frequency f0=2KHz for binary 1 and Frequency f1=6KHz for binary 0    
 for i= 1:1:105
     if (data(i)=='0')
         FSK_MOD=[FSK_MOD h1b];
     else FSK_MOD=[FSK_MOD h2b];
     end
 end
 
 %% PART F 
 figure(4)

subplot(2,1,1)

plot(FSK_MOD(1:201),'linewidth',2)
axis([ 0 210 -1.1 1.1])
grid on;
title('FSK Modulated Signal')
xlabel('Time')
ylabel('Amplitude')

subplot(2,1,2)
plot((-0.5:1/400:0.5-1/400)*40,(abs(fftshift(fft(FSK_MOD,400)))),'linewidth',2)
grid on

axis([-20 20 0 130])
title('Frequency response of FSK modulated signal')
xlabel('Frequency')
ylabel('Magnitude')

%% PART G (DETECTION USING BAND PASS FILTERS)
filter_two=0.1*filter(h1b,1,FSK_MOD);%filter for 2KHz
filter_six=0.1*filter(h2b,1,FSK_MOD);%filter for 6KHz

figure(5)

subplot(2,1,1)
plot(filter_two(1:400),'linewidth',2)
grid on
title('Signal Filtered through 2khz (sine wave)')
xlabel('Time')
ylabel('Amplitude')

subplot(2,1,2)
plot(filter_six(1:400),'linewidth',2)
grid on
title('Signal Filteres through 6khz (sine wave)')
xlabel('Time')
ylabel('Amplitude')

figure(6)

subplot(2,1,1)
avg_filter_two = conv(abs(filter_two)/12.63,ones(1,20));
plot(avg_filter_two,'linewidth',2)
grid on
title('Averaged output to unity peak for filter of 2KHz (sine wave)')
axis([0 500 -0.1 1.5])
xlabel('Real Axis')
ylabel('Imaginary Axis')

subplot(2,1,2)
avg_filter_six = conv(abs(filter_six)/12.63,ones(1,20));
plot(avg_filter_six,'linewidth',2)
grid on
title('Averaged output to unity peak for filter of 6KHz (sine wave)')
axis([0 500 -0.1 1.5])
xlabel('Real Axis')
ylabel('Imaginary Axis')

FSK_DEMOD=[];  %FSK demodulated data
for i=40:20:length(avg_filter_six)
    if (avg_filter_two(i)>avg_filter_six(i))
        FSK_DEMOD=[FSK_DEMOD 0]
    else FSK_DEMOD=[FSK_DEMOD 1]
    end
end

figure(7)
  
subplot(2,1,1)
plot(FSK_MOD(1:400),'linewidth',2)
title('FSK Modulated Signal Transmitted')
xlabel('Time Series')
ylabel('Amplitude')
grid on
  
subplot(2,1,2)
stairs(FSK_DEMOD(2:20),'linewidth',2)
title('Detected Signal from Band Pass Filters')
xlabel('Time Series')
ylabel('Amplitude')
axis([1 20 -0.2 1.2])
grid on


%% PART H(using complex sinusoids to avoid band pass filter for detection of FSK)
FSK_MOD1=[]  %FSK modulation using h1c and h2c
 for i= 1:1:105
     if (data(i)=='0') 
          FSK_MOD1=[FSK_MOD1 h1c];
     else FSK_MOD1=[FSK_MOD1 h2c];
     end
 end
 
filter_two1=1/20*filter(h1c,1,FSK_MOD1);
filter_six1=1/20*filter(h2c,1,FSK_MOD1);

figure(8)

subplot(2,1,1)
abs_filter_two1=abs(filter_two1);
plot(abs_filter_two1,'linewidth',2)
axis([0 500 -0.1 1.5])
grid on
title('FSK signal filtered with 2khz (complex sinusoid)')
xlabel('Time Index')
ylabel('Amplitude')

subplot(2,1,2)
abs_filter_six1=abs(filter_six1);
plot(abs_filter_six1,'linewidth',2)
axis([0 500 -0.1 1.5])
grid on
title('FSK signal filtered with 6khz (complex sinusoid)')
xlabel('Time Index')
ylabel('Amplitude')

FSK_DEMOD1=[];  %FSK demodulated data
for i=20:20:2100
    if (abs_filter_two1(i)>abs_filter_six1(i))
        FSK_DEMOD1=[FSK_DEMOD1 0];
    else FSK_DEMOD1=[FSK_DEMOD1 1];
    end
end


FSK_DEMOD2=[];  %FSK demodulated data
for i=40:20:length(avg_filter_six)
    if (avg_filter_two(i)>avg_filter_six(i))
        FSK_DEMOD2=[FSK_DEMOD 0]
    else FSK_DEMOD2=[FSK_DEMOD 1]
    end
end

figure(9)
  
subplot(2,1,1)
plot(FSK_MOD(1:400),'linewidth',2)
title('FSK Modulated Signal Transmitted')
xlabel('Time Index')
ylabel('Amplitude')
grid on
  
subplot(2,1,2)
stairs(FSK_DEMOD2(2:20),'linewidth',2)
title('Detected Signal from Band Pass Filters')
xlabel('Time Index')
ylabel('Amplitude')
axis([1 20 -0.2 1.2])
grid on

%%PART I (Using derivative function)
deriv = 1./ [1:10].*(-1).^(1:10);
ful_deriv=[-fliplr(deriv) 0 deriv].*hanning(21)';
imp_deriv=filter(ful_deriv,1,[1 zeros(1,100)]); %impulse response of derivative function

figure(10)

subplot(2,1,1)
plot(imp_deriv,'linewidth',2)
grid on;
axis([0 100 -1.1 1.1])
title('Impulse of derivative')
xlabel('Time Index')
ylabel('Amplitude')

subplot(2,1,2)
plot((-0.5:1/400:0.5-1/400),(abs(fftshift(fft(imp_deriv,400)))),'linewidth',2);% frequency response of derivative function
grid on;
title('Frequency Response of derivative function')
xlabel('Normalized Frequency')
ylabel('Amplitude')
filter_derv=filter(ful_deriv,1,FSK_MOD);

figure(11)

subplot(3,1,1)
plot(FSK_MOD(1:201),'linewidth',2);
axis([0 215 -1.2 1.2]);
grid on;
title('Input to Derivative function')
xlabel('Time Index')
ylabel('Amplitude')

subplot(3,1,2)
plot(filter_derv(1:201),'linewidth',2)
axis([0 215 -1.2 1.2]);
grid on;
title('Output of derivative filter')
xlabel('Time Index')
ylabel('Amplitude')

subplot(3,1,3)
FD=conv(abs(filter_derv)/2,ones(1,20));
plot(abs(filter_derv(1:201)),'linewidth',2)
hold on;
plot(FD(10:210)/6,'r','linewidth',2);
hold off;
axis([0 215 0 1.2]);
grid on;
title('Average of output by Box car')
xlabel('Time Index')
ylabel('Amplitude')

x=[]
for i=20:20:2100+19
    x=[x mean(FD(i:i+19))];
   end
FSK_DEMOD3=[];
  for k=1:105
      if ( x(k)>mean(FD))
       FSK_DEMOD3(k)=1;
          else FSK_DEMOD3(k)=0;
      end
  end
  
figure(12) %demod3 
  
subplot(2,1,1)
plot(FSK_MOD(1:400),'linewidth',2)
title('Transmitted FSK Modulated Signal ')
xlabel('Time Index')
ylabel('Amplitude')
grid on
  
subplot(2,1,2)
stairs(FSK_DEMOD3(2:20),'linewidth',2)
title('Detected Signal from Derivative Filter')
xlabel('Time Index')
ylabel('Amplitude')
axis([1 20 -0.2 1.2])
grid on


%% PART K,L (using Hilbert transform to detect FSK)
hilbert_func = (2/pi)*sin(pi*[-10:10]/2).^2./[-10:-1,1,1:10].*hanning(21)';
hilbert_filter = filter(hilbert_func,1,[1 zeros(1,100)]);

figure(13)%(Impulse and frequency of hilbert)

subplot(2,1,1)
plot (hilbert_filter,'linewidth',2);
axis([0 100 -0.8 0.8]);
grid on;
title('Impulse response of hilbert filter');
xlabel('Time Index');
ylabel('Amplitude')

subplot(2,1,2)
plot((-0.5:1/400:0.5-1/400),(abs(fftshift(fft(hilbert_filter,400)))),'linewidth',2);
grid on
title('Frequency Response of hilbert filter')
xlabel('Nomrmalized Frequency')
ylabel('Amplitude')

figure(14)%(FSK detection)

filter1=filter(hilbert_func,1,FSK_MOD);% pass FSK from hilbert function

subplot(3,1,1)
plot(FSK_MOD(1:201),'linewidth',2)
grid on;
axis([0 215 -1.2 1.2])
title('Input to Hilbert Filter')
xlabel('Time Index')
ylabel('Amplitude')

subplot(3,1,2)
plot(filter1(1:201),'linewidth',2)
grid on;
axis([0 215 -1.2 1.2])
title('Output of Hilbert Filter')
xlabel('Time Index')
ylabel('Amplitude')

complex_formed=[FSK_MOD+j*filter1]
ang=angle(complex_formed);
unw=unwrap(ang,pi);
filter2=filter(ful_deriv,1,unw);%unwrap derivative

subplot(3,1,3)
plot(abs(filter2(1:201)),'linewidth',2)
axis([0 215 0 2.5])
grid on;
title('Output of unwraped differentiator signal')
xlabel('Time Index')
ylabel('Amplitude')
x=[]
for i=20:20:2100-19
    x=[x mean(filter2(i:i+19))];
   end
FSK_DEMOD4=[];
  for k=1:104
      if ( x(k)> mean(filter2))
       FSK_DEMOD4(k)=0;
          else FSK_DEMOD4(k)=1;
      end
  end
  
figure(15) %demod4 
  
subplot(2,1,1)
plot(FSK_MOD(1:400),'linewidth',2)
title('Transmitted FSK Modulated Signal ')
xlabel('Time Index')
ylabel('Amplitude')
grid on
  
subplot(2,1,2)
stairs(FSK_DEMOD4(2:20),'linewidth',2)
title('Detected Signal from unwrapped Filter')
xlabel('Time Index')
ylabel('Amplitude')
axis([1 20 -0.2 1.2])
grid on


%% PART M (using arc tangent for FSK detection)
filter3=filter(ful_deriv,1,filter1);
arc_tan=[((FSK_MOD.*filter3)-(filter_derv.*filter1))./(FSK_MOD.^2 + filter1.^2)];

figure(16)

subplot(2,1,1)
plot(FSK_MOD(1:201),'linewidth',2)
grid on
axis([0 215 -1.1 1.1])
title('FSK Modulated Input')
xlabel('Time Index')
ylabel('Amplitude')

subplot(2,1,2)
plot(abs(arc_tan(1:201)),'linewidth',2)
grid on
axis([0 215 -0.1 1.5])
title('output of arc tangent')
xlabel('Time Index')
ylabel('Amplitude')


x=[]
for i=20:20:2100-19
    x=[x mean(arc_tan(i:i+19))];
   end
FSK_DEMOD5=[];
  for k=1:104
      if ( x(k)> 0.5)
       FSK_DEMOD5(k)=1;
          else FSK_DEMOD5(k)=0;
      end
  end
  
figure(17) %demod5 
  
subplot(2,1,1)
plot(FSK_MOD(1:400),'linewidth',2)
title('Transmitted FSK Modulated Signal ')
xlabel('Time Index')
ylabel('Amplitude')
grid on
  
subplot(2,1,2)
stairs(FSK_DEMOD5(1:20),'linewidth',2)
title('Detected Signal from Arc Tan Filter')
xlabel('Time Index')
ylabel('Amplitude')
axis([1 20 -0.2 1.2])
grid on
