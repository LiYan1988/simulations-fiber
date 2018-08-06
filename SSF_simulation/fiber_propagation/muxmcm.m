
function Uf_mux=muxmcm(Uif,Nsc,Nu,Nch)

Uf_mux = zeros(Nu*Nsc,Nch,2);
        %uit(symbols , WDM channels , subcariere , Polarizaton)

for i=1:Nch
   for n1=1:Nsc
      for n2=1:Nu
         ind = (n1-1)*Nu+n2;
         Uf_mux(ind,i,:) = Uif(n2,i,n1,:);
      end
   end
end
%Uf_mux(sub_cr1,subcr2,..., WDM channel, Pol)
end





