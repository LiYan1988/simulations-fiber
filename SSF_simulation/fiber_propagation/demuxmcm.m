
function Uf_demux=demuxmcm(Uif,Nsc,Nu,pol)

Uf_demux = zeros(Nu,Nsc,pol);

% Demux di tutte le sotto-portanti del canale osservato
for n=1:Nsc
   for n2=1:Nu
      ind = (n-1)*Nu+n2;
      Uf_demux(n2,n,:) = Uif(ind,:);
   end
end

end





