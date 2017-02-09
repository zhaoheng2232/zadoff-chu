function [zc] = gen_zc_seq(r,L)

n = 0:L-1;
if mod(L,2) == 1
    zc = exp(1i*pi*r*n.*(n+1)/L);
else
    zc = exp(1i*pi*r*n.*n/L);
end
end