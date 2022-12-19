% Example how to use multitaper methods for spectrum estimation.

clear all; close all;

[s,fs] = wavread('7_text8_1.wav');

framesize   = 240; 
frameshift  = 120;
frames      = enframe(s, framesize, frameshift);

numwindows = 8;
[thomson_tapers,thomson_weights]      = thomson(framesize, numwindows);
[multipeak_tapers, multipeak_weights] = multipeakwind(framesize, numwindows);
[SWCE_tapers, SWCE_weights]           = SWCE(framesize, numwindows);

%% Compute the spectra using different methods. For comparison, include
%% Hamming method as well.
NFFT = 512;
spec_hamming = multitaperspectra(frames, hamming(framesize), 1, NFFT);
spec_thomson = multitaperspectra(frames, thomson_tapers, thomson_weights, NFFT);
spec_multip  = multitaperspectra(frames, multipeak_tapers, multipeak_weights, NFFT);
spec_SWCE    = multitaperspectra(frames, SWCE_tapers, SWCE_weights, NFFT);
