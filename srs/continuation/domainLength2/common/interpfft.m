function [u] = interpfft(u,n)
%INTERPFFT Summary of this function goes here
%   Detailed explanation goes here
u = u(:);
u = [u; flipud(u(1:end-1))];
n = 2*n-1;
fu = zeros(n, 1);
mid = floor(n/2);
offsetPlus = ceil(numel(u)/2);
offsetMinus = offsetPlus - numel(u) + 1;
fu(mid+offsetMinus:mid+offsetPlus) = n / numel(u) * fftshift(fft(u));
u = real(ifft(ifftshift(fu))); 
u = u(1:(n+1)/2);
end

