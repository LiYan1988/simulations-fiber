function [signal, signal_complex_received_out, signal_complex_original_out] = ...
    align_signals(signal, sym_rm, param, cidx)
% There is a time shift between the received constellation signal and
% original signal. These signals are symbols/points on the constellation.
% This function finds the time offset between them.

% First truncate off the head and tail of original signal
u = param.data_mod_symbol_channel{cidx}(sym_rm-1:sym_rm+size(signal, 1)-2);
% Convert the received signal from an Nx2 real vector to an Nx1 complex
% vector
v = signal(:, 1)+signal(:, 2)*1i;
% cross correlation between v and u
[r, l] = xcorr(v, u);
% find the offset of the max correlation
[~, idx_xcorr] = max(r);
max_lag = l(idx_xcorr);

% if max_lag<0, v remove abs(max_lag) tail elements and u remove
% abs(max_lag) head elements
% if max_lag>0, v remove abs(max_lag) head elements and u remove
% abs(max_lag) tail elements
% if max_lag==0, keep them the same
if max_lag==0
    vv = v;
    uu = u;
elseif max_lag<0
    vv = v(1:end+max_lag);
    uu = u(1-max_lag:end);
elseif max_lag>0
    vv = v(1+max_lag:end);
    uu = u(1:end-max_lag);
end

signal_complex_received_out = vv;
signal_complex_original_out = uu;

signal = zeros(size(signal_complex_received_out, 1), 2);
signal(:, 1) = real(signal_complex_received_out);
signal(:, 2) = imag(signal_complex_received_out);

%% for testing
% figure;
% hold on;
% plot(l, abs(r))
% tmp = abs(sum(vv.*conj(uu)));
% 
% circx = zeros(size(vv));
% for n=1:size(vv, 1)
%     circx(n) = sum(circshift(vv, n).*conj(uu));
% end
% 
% figure;
% plot(abs(circx))

