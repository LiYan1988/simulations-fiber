
function res=AIR_calc_aux6p(x_tmp,y_tmp,N_particle,eta,init_param)

adaptation = true;

res.parametri_iniziali = init_param;
res.parametri_stimati = init_param;

if adaptation
    [res.parametri,~] = channel_estimate_polrot_6p(x_tmp.',...
        y_tmp.',res.parametri_stimati,eta,N_particle);
    res.parametri_stimati = res.parametri(:,end);
end

%* Aggiornamento della struttura ch_aux
res.h0       = res.parametri_stimati(1);
res.sg2_awgn = res.parametri_stimati(2)^2;
res.sg2_R    = res.parametri_stimati(3:6).^2;
res.snr      = 10*log10(abs(res.h0)^2/res.sg2_awgn);

%* Stima della AIR
%* Estimate the conditional entropy h(Y|X)
[res.h_YdatoX,res.nr_YdatoX] = conditional_entropy_particle_polrot_6p(x_tmp.',...
    y_tmp.',N_particle,res.parametri_stimati);

%* Estimate the output entropy h(Y)
res.sg2y = abs(res.h0)^2+res.sg2_awgn;
[res.h_Y]  = gaussian_entropy_2D(y_tmp.',res.sg2y);

%* Stima dell'AIR
res.air = (res.h_Y-res.h_YdatoX).';     % stima basata sul metodo particle per l'entropia condizionata e su calcolo analitico per entropia di uscita

%**********************************************************************
%* Canale ausiliario AWGN
cxxs_tmp = sum(abs(x_tmp).^2);
cyys_tmp = sum(abs(y_tmp).^2);
cxys_tmp = sum(y_tmp.*conj(x_tmp));
xc_tmp = cxys_tmp./sqrt(cxxs_tmp.*cyys_tmp);
h0_tmp = cxys_tmp./cxxs_tmp;
z_tmp = y_tmp./h0_tmp;
xcam_tmp = abs(xc_tmp);
snr_tmp = xcam_tmp.^2./(1-xcam_tmp.^2);
res.AWGN.h0 = h0_tmp;
res.AWGN.snrdb = 10*log10(snr_tmp);

%* Stima dell'AIR
res.AWGN.air = sum(log2(1+snr_tmp)-snr_tmp/log10(2).*mean(abs(z_tmp-x_tmp).^2-abs(z_tmp).^2./(snr_tmp+1)));
