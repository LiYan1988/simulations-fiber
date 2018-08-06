function out=Dispersion_matrix(P)
% P is a structure containing the channel parameters
% out is a vector containing the diagonal enteries of the dispersion matrix
% the dispersion step would be ifft(out.*fft(signal))

x=(P.L/2-abs(P.L/2-[0:P.L-1]));
out=exp(P.betta2*1i/2*P.dz*(2*pi/P.T)^2*x.^2);

end
