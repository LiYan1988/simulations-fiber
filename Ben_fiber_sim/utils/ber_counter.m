function [ber, q_db, total_error,total_bit] = ber_counter(rx_bits,tx_bits)
% Calculates the BER and Q_BER (in dB)

% modified by Jianqiang Li on 4 August, 2011
% Modified by Benjamin Foo, 2018-01-30

    if (size(rx_bits, 1) ~=1 && size(rx_bits, 2) ~= 1) || (size(tx_bits, 1) ~=1 && size(tx_bits, 2) ~= 1)
        error('Expected vectors, not arrays')
    end

    data_len      = length(rx_bits);             % Total length of RX data set
    code_word_len = length(tx_bits);          % Length of TX data pattern
    N_code_words  = floor((data_len)/(code_word_len)); % Number of times TX pattern repeats in RX data set
    
    if N_code_words >= 1
        % Truncating and reshaping RX data set into matrix where each column is
        % one TX pattern
        rx_bit_array=reshape(rx_bits(1:code_word_len*N_code_words),[code_word_len,N_code_words]);
        test_pattern=rx_bit_array(:,1); % Using the first column as a test pattern for synchronization
    else
        % If only part of the pattern is received
        rx_bit_array = rx_bits; % If less than one pattern repitition, assign directly to output
        test_pattern = rx_bit_array;
    end
    
    % Synchronizing RX data to TX pattern
    [sync, ~] = xcov(tx_bits,test_pattern);
    sync_offset = find(sync==max(sync)) - code_word_len;
    tx_bits=circshift(tx_bits, -sync_offset); % Circular shift of TX pattern to match RX data    
%     [cov, lags] = xcov(tx_bits, test_pattern);
%     figure(11); plot(lags, cov); % Debugging
    
    if N_code_words % 1 or more patterns received
        
        tx_bit_array=zeros(length(tx_bits),N_code_words);
        for ind=1:N_code_words
            tx_bit_array(:,ind)= tx_bits;
        end
        total_bit = code_word_len*N_code_words;
        
    else
        % Less than 1 pattern received
        
%         sync_offset = find(sync==max(sync)) - code_word_len;
%         tx_bit_array = circshift(tx_bits, -sync_offset); % Synchronizing
        tx_bit_array = tx_bits(1:length(rx_bit_array)); % Truncating TX pattern to match RX length
        total_bit = data_len;
    end
    
    total_error=sum(sum(xor(rx_bit_array,tx_bit_array)));
    ber=total_error/total_bit;
    q_db=20*log10(sqrt(2)*erfcinv(2*ber));
    
end