function [BER11,BER14,BER22,BER41] = mimoOSTBCWithPCT(frLen, numPackets,EbNo)
% MIMOOSTBCWITHPCT Simulates BER calculations in parallel for different
% orthogonal space-time block coding systems of diversity order 4.
%
%   [BER11,BER14,BER22,BER41] = MIMOOSTBCWITHPCT(FRLEN, NUMPACKETS, EBNO)
%   calculates the BER of an orthogonal space-time block coding system over
%   a range of Eb/No values in parallel across the number of available
%   workers. Three systems with four transmit antennas, four receive
%   antennas and two transmit and two receive antennas respectively, are
%   considered since they all have a diversity order of 4. The performance
%   of these systems and a 1X1 system is compared with the BER values
%   obtained. FRLEN is the frame length, NUMPACKETS is the number of
%   packets simulated and EBNO is the Eb/No range.

% Copyright 2014 The MathWorks, Inc.

if isempty(gcp('nocreate'))
    parpool;
end
pool = gcp;
numWorkers = pool.NumWorkers;
% Simulation Parameters
lenEbNo = length(EbNo);     % length of the Eb/No vector
seed = 589*(1:numWorkers);  % seed for the random number generator

[errs11, errs14, errs22, errs41] = deal(zeros(numWorkers,lenEbNo));
[bits11, bits14, bits22, bits41] = deal(zeros(numWorkers,lenEbNo));

parfor n = 1:numWorkers
    for idx = 1:lenEbNo
        ber11 = mrc1m_pct(1,frLen,numPackets/numWorkers,EbNo,idx,seed(n));
        ber14 = mrc1m_pct(4,frLen,numPackets/numWorkers,EbNo,idx,seed(n));
        ber22 = ostbc2m_pct(2,frLen,numPackets/numWorkers,EbNo,idx,seed(n));
        ber41 = ostbc4m_pct(1,frLen,numPackets/numWorkers,EbNo,idx,seed(n));
        errs11(n,idx) = ber11(:,2);
        errs14(n,idx) = ber14(:,2);
        errs22(n,idx) = ber22(:,2);
        errs41(n,idx) = ber41(:,2);
        bits11(n,idx) = ber11(:,3);
        bits14(n,idx) = ber14(:,3);
        bits22(n,idx) = ber22(:,3);
        bits41(n,idx) = ber41(:,3);
    end
end

BER11 = sum(errs11,1)./sum(bits11,1);
BER14 = sum(errs14,1)./sum(bits14,1);
BER22 = sum(errs22,1)./sum(bits22,1);
BER41 = sum(errs41,1)./sum(bits41,1);

