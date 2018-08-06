%% Funzione di trasferimento della fibra dovuta a dispersione e
%% dispersion-slope
%transfer function of chromatic dispersion in freq. 

function Fb = fiber(b2,b3,y)
      
Fb = exp(1i*wrapToPi(-b2*y.^2-b3*y.^3)).';