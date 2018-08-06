%
% Genera N campioni di rumore gaussiano bianco complesso 
% mediante metodo di Box-Muller (polare)
%
%
% W(N) - double complex array
%   Output: gli N campioni di rumore generati 
%
% N - integer
%   Input: numero di campioni
%
% sgn - double precision
%    Input: deviazione standard di ciascuna componente I o Q di rumore
%
function W=gaussgen(N,sgn)

W = sgn*(randn(1,N)+1i*randn(1,N));