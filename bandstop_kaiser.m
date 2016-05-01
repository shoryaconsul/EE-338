%Bandpass filter with Kaiser window
close all
clear all

%deriving filter specifications
number = 104;
if(number>75)
    number = number -75;
end

q = floor(number/10);
r = mod(number, 10);
BSL = 4 + 0.9*q + 2*r;
BSH = BSL+10;
BPL = BSL-2;
BPH = BSH + 2;

samp = 100; %sampling frequency in kHz

%normalized specifications
bpl = BPL/samp*2*pi;
bph = BPH/samp*2*pi;
bsl = BSL/samp*2*pi;
bsh = BSH/samp*2*pi;

%cutoff frequencies of low pass filters used
wch = (bph+bsh)/2;
wcl = (bpl+bsl)/2;

delta  = 0.15;
A = -20*log10(delta);
dwt = 2/100*2*pi;

N = ceil((A-8)/(2*2.285*dwt));
if(A<21) 
    alpha = 0;
elseif(A>=21 && A<=50)
    aplha = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else 
    alpha = 0.1102*(A-8.7);
end


N = N+6;            %increasing window length to reduce ripple
beta = alpha/N;

window = kaiser(2*N+1, beta);       %kaiser window
hid_up = ideal_lp(wch, 2*N+1);     %"Bigger" LPF
hid_low = ideal_lp(wcl, 2*N+1);    %"Smaller" LPF

hid_bp = hid_up - hid_low;         %getting ideal band pass impulse response (converse of band stop)

allpass = zeros(size(hid_bp));
allpass(N+1) = 1;                 %impulse response of all pass filter
hid = allpass - hid_bp;

filt = hid.*window';            %windowing impulse repsponse

tf(filt, [1], 0.1, 'variable','z^-1')
fvtool(filt,[1]);               %coefficients of z-transform given by terms of impulse response


