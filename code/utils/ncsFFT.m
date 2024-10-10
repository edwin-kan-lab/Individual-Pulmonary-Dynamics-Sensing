function [P1, f] = ncsFFT(ncs, fncs, show, range)

Y = fft(ncs); L = length(ncs); L = L-mod(L, 2);

P2 = abs(Y/L); 
P1 = P2(1:L/2+1); P1(2:end-1) = 2*P1(2:end-1);
f = fncs*(0:(L/2))/L;

if show
    figure; plot(f, P1); title('FFT of NCS'); xlim([range(1) range(2)]);
end

end