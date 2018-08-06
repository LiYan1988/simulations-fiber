function h = gaussian_entropy_2D(y,sg2y)
%h = gaussian_entropy(y,sg2y) estimates the differential entropy of the
% variable y assuming a complex Gaussian distribution
%
% y: realizations of the variable
% sg2y: total variance of the distribution (per polarization)


N_campioni=length(y);
ususg2=1/sg2y;
h=2*log2(pi*sg2y)+log2(exp(1))*ususg2*(cumsum(abs(y(1,:)).^2+abs(y(2,:)).^2)./(1:N_campioni));

end

