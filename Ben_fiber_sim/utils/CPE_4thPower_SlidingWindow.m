function [CarrierRecovered_Signal, Phase_error_estimated] = CPE_4thPower_SlidingWindow(Complex_Samples,FreOffset,BlockSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency offset compensation and carrier phase estimation based on
% 4th-power method with sliding window for QPSK
%
% Input: 
%  Complex_Samples: complex samples at 1 sample/symbol
%  FreOffset: value of frequency offset in rad
%  BlockSize: block size for carrier phase estimation in sliding window mode
%
% Output:
%  CarrierRecovered_Signal: complex signal after carrier recovery
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Author: Jianqiang Li  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
symbol_num=length(Complex_Samples);
Complex_Samples=reshape(Complex_Samples,symbol_num,1); %To ensure column vector

%% Frequency offset compensation
if (FreOffset~=0)
    sym_offset=zeros(size(Complex_Samples));
    for Ind=2:1:symbol_num
        sym_offset(Ind)=mod((sym_offset(Ind-1)+FreOffset),2*pi);
    end
    Complex_Samples=Complex_Samples.*exp(1i.*(-sym_offset));
end

%% Carrier phase estimation
phase_error_4thPower=zeros(symbol_num,1);
for Ind=1:1:symbol_num
    current_block=get_block(Complex_Samples,Ind,BlockSize);
    block_4thPower=current_block.^4;
    phase_error_4thPower(Ind)=angle(sum(-block_4thPower));
end

Phase_error_estimated=unwrap(phase_error_4thPower)/4;
CarrierRecovered_Signal_temp=Complex_Samples.*exp(-1i.*Phase_error_estimated);
CarrierRecovered_Signal=CarrierRecovered_Signal_temp(floor(BlockSize/2)+1:end-floor(BlockSize/2)); % Removes symbols from front and back

end

%% Sub functions
function current_block=get_block(dataArray,dataIdx,blocksize)
	block_start=dataIdx-floor(blocksize/2);
	block_end=dataIdx+floor(blocksize/2);
	len=length(dataArray);
	if (block_start<1); block_start=1; end;
	if (block_end>len); block_end=len; end;
	current_block=dataArray(block_start:block_end);
end
