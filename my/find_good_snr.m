function [snr_good] = find_good_snr(snr_input, spacing_GHz, ...
    snr_QAM_th, snr_OOK_th)
% Find the powers of QAM and OOK that have good SNRs
% spacing_GHz: the channel spacing in GHz
% snr_QAM_th: SNR threshold of QAM
% snr_OOK_th: SNR threshold of OOK
snr_part = snr_input(snr_input{:, 'spacing_GHz'}==spacing_GHz, :);
idx = (snr_part{:, 'snr_OOK_mean'}>=snr_OOK_th) & (snr_part{:, 'snr_QAM'}>=snr_QAM_th);
snr_good = snr_part(idx, :);