function [ydec, ymp] = demod(yr, imod)

if (imod == 5) %16-QAM demodulator
    ydec = qamdemod(yr,16,'gray','UnitAveragePower',true);
    ymp  = qammod(ydec,16,'gray','UnitAveragePower',true);
elseif (imod == 6) %64-QAM demodulator
    ydec = qamdemod(yr,64,'gray','UnitAveragePower',true);
    ymp  = qammod(ydec,64,'gray','UnitAveragePower',true);
elseif (imod == 4) %BPSK demodulator
    ydec = pskdemod(yr,2,0,'gray');
    ymp  = pskmod(ydec,2,0,'gray');
elseif (imod == 1) %QPSK demodulator
    ydec = pskdemod(yr,4,0,'gray');
    ymp  = pskmod(ydec,4,0,'gray');
end