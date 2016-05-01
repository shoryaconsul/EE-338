function [ h ] = ideal_lp( wc, M )

a = (M-1)/2;
n = [0:1:M-1];
m = n - a;
if (mod((M-1), 2) == 0)
    m(a+1) = 1;
end

h = sin(wc*m)./(pi*m);

if (mod((M-1), 2) == 0)
    h(a+1) = wc/pi;
end

end

