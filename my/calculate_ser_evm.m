function param = calculate_ser_evm(param)
% Calculate SER and EVM for each channel

% SER
param.ser_channel = zeros(1, param.channel_number);
% EVM
param.evm_channel = zeros(1, param.channel_number);
% Q value
param.q_channel_1 = zeros(1, param.channel_number);
param.q_channel_2 = zeros(1, param.channel_number);

for cidx = 1:param.channel_number
    % rotate constellation back
    tx = param.signal_tx_complex{cidx};
    rx = param.signal_rx_complex{cidx};
    
    objfcn = @(v)v(1)+v(2)*rx - tx;
    
    opts = optimoptions(@lsqnonlin,'Display','off','UseParallel',false);
    x0 = (1+1i)*[1;1]; % arbitrary initial guess
    x_sol = lsqnonlin(objfcn,x0,[],[],opts);
    
    rx2 = x_sol(1) + x_sol(2)*rx;
    
    % Calculate SER
    tx_unique = unique(tx);
    a = abs(rx2 - tx_unique.');
    [u, pidx] = min(a, [], 2);
    v = abs(rx2 - tx);
    param.ser_channel(cidx) = sum(u<v)/length(u);
    
    % EVM
    a = abs(rx2-tx).^2;
    b = abs(tx).^2;
    param.evm_channel(cidx) = sqrt(mean(a)/mean(b));
    
    % calculate Q value for OOK channels
    if strcmp(param.channel_type{cidx}, 'ook')
        rx3 = abs(rx2).^2;
        s1 = std(rx3(tx==1));
        s0 = std(rx3(tx==0));
        i1 = mean(rx3(tx==1));
        i0 = mean(rx3(tx==0));
        param.q_channel_1(cidx) = (i1-i0)/(s1+s0);
        
        rx3 = abs(rx2);
        s1 = std(rx3(tx==1));
        s0 = std(rx3(tx==0));
        i1 = mean(rx3(tx==1));
        i0 = mean(rx3(tx==0));
        param.q_channel_2(cidx) = (i1-i0)/(s1+s0);
    end
end