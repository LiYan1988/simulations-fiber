function [fasi,pesi] = particle_resample(fasi_old,pesi_old)
%[fasi,pesi] = particle_resample(fasi_old,pesi_old) resample the particle set


N=length(fasi_old);
p=[0;cumsum(pesi_old)];
x=rand(N,1);
[~,ind]=histc(x,p);
fasi=fasi_old(ind,:);
uniweigth=1/N;
pesi=repmat(uniweigth,N,1);

end

