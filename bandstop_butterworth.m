%Designing a Buttweorth low pass filter
%Low pass specifications (analog domain): Stop band edge is at 1.4372
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

%corresponding analog filter specifications
wpl = tan(bpl/2);
wph = tan(bph/2);
wsl = tan(bsl/2);
wsh = tan(bsh/2);

%parameters for band stop to low pass transformation
omega_0 = sqrt(wpl*wph);
B = wph - wpl;

Omega_s1 = (B*wsl)/(omega_0^2 - wsl^2);
Omega_s2 = (B*wsh)/(omega_0^2 - wsh^2);

ws = min(abs(Omega_s1), abs(Omega_s2));

%Finding Butterworth filter parameters

D1 = 1/0.85^2 -1;
D2 = 1/0.15^2 -1;

N = ceil(log(sqrt(D2/D1))/log(ws));
wc = 0.5*(D1^(-1/(2*N))+ws*D2^(-1/(2*N)));

%N = 8;
%wc = 1.08;

poles = zeros(N,1);

for k=1: N
    poles(k) = 1.08i*exp(1i*(2*k-1)*pi/(2*N));
end


%displaying poles and zeros
zplane([],poles);
grid;
title('Pole - zero plot of Butterworth analog LPF');

dc = prod(poles);     %normalizing the lpf transfer function given by zp2tf

%analog LPF transfer function and its freq response
[b_lp,a_lp] = zp2tf([], poles, dc);
figure;
freqs(b_lp,a_lp);
axis([0.1 2 0.3 1.05]);
title('Frequency response of analog low pass filter');

%making coefficients real to avoid numerical inaccuracies
a_lp = real(a_lp);
b_lp = real(b_lp);

fprintf('The analog low pass transfer function is:');  
tf(b_lp,a_lp)

%{
%finding analog low pass transfer function
H_L = tf(b,a);
fprintf('The analog low pass transfer function is:');  
H_L
%}

syms s_l s;      %s_l is low pass variable, s is band pass variable
denom_poly = poly2sym(a_lp, s_l);
num_poly = poly2sym(b_lp, s_l);

lp_poly = num_poly/denom_poly;  %low pass transfer function in floating point numbers 

%fprintf('The analog low pass transfer function is:');  
%vpa(lp_poly, 5)

%getting bandstop transfer function
omega_0 = sqrt(wpl*wph);
B = wph - wpl;

syms lb;
lb = (B*s)/(s^2 + omega_0^2);           %low pass to band pass transformation
bp_poly = simplifyFraction(expand(subs(lp_poly, s_l, lb)));

%fprintf('The analog band stop transfer function is:');  
%vpa(bp_poly, 5)

%coverting band stop polynomial to vectors
[num_bp, denom_bp] = numden(bp_poly);
b_bp = sym2poly(num_bp);
a_bp = sym2poly(denom_bp);

fprintf('The analog band stop transfer function is:');  
b_bpf = b_bp./a_bp(1);
a_bpf = a_bp./a_bp(1);
tf(b_bpf, a_bpf)

%plotting frequency response of band pass transfer function
figure;
freqs(b_bp,a_bp);
axis([0.1 2 0.0005 1.05]);
title('Frequency response of analog band stop filter');
 
%bilinear transformation
syms sz zi;             %zi is 1/z
sz = (1-zi)/(1+zi);
bp_z = simplifyFraction(expand(subs(bp_poly, s, sz)));
%fprintf('The discrete time filter transfer function is:');  
%vpa(bp_z, 5)            %discrete time filter transfer function

[numz, denomz] = numden(bp_z);
bz = fliplr(sym2poly(numz));
az = fliplr(sym2poly(denomz));

%printing discrete filter transfer function
bzf = bz./az(1);
azf = az./az(1);
%fprintf('The discrete time filter transfer function is:');  
tf(bzf,azf,0.1,'variable','z^-1')

%freqz(bz, az, 100);
fvtool(bz,az);
