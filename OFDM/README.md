This project is presented as a final exam for the course Modem Design (EE-655):

There are two Matlab software model implementation:

Part 1: OFDM modem
OFDM Modem is implemented as a 128 point transform with 54 occupied frequency bins. In particular, bins –27-to-1 and +1 to +27, (bin 0 empty). 
Each bin is modulated with 16-QAM. The time series formed at the output of the transform is separated by a guard band of 32 zero valued samples 
located at the beginning of the OFDM symbol. An OFDM packet composed of 50 OFDM symbol is formed and transmitted through a channel with impulse 
response and demodulated.

Part 2: SCOFDM modem
OFDM Modem is implemented as a 128 point transform with 65 occupied frequency bins. In particular, bins –27-to-1 and +1 to +27 (bin 0 empty). 
Each bin is modulated 16-QAM. The time series formed at the output of the transform is separated by a guard band of 32 zero valued samples located 
at the beginning of the OFDM symbol. An OFDM packet composed of 50 OFDM symbol is created. A single carrier OFDM signal is generated by 
cascade transforms with approximately constant spectral amplitude. This is achieved with a linear FM sweep known as CAZAC (constant amplitude, 
zero auto correlation). The final SC-OFDM signal is then passed through a channel to simulate transmission and demodulated.