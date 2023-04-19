function BER1m = mrc1m(M, frLen, numPackets, EbNo)
%MRC1M  Maximal-Ratio Combining for 1xM antenna configurations.
%
%   BER1M = MRC1M(M, FRLEN, NUMPACKETS, EBNOVEC) computes the bit-error rate 
%   estimates via simulation for a Maximal-Ratio Combined configuration using 
%	one transmit antenna and M receive antennas, where the frame length, number
%   of packets simulated and the Eb/No range of values are given by FRLEN, 
%	NUMPACKETS, EBNOVEC parameters respectively.
%
%   The simulation uses BPSK modulated symbols with appropriate receiver 
%   combining.
%
%   Suggested parameter values:
%       M = 1 to 4; FRLEN = 100; NUMPACKETS = 1000; EBNOVEC = 0:2:20;
%
%   Example:
%       ber12 = mrc1m(2, 100, 1000, 0:2:20);
%
%   See also OSTBC2M, OSTBC4M.

%   References:
%   [1] J. G. Proakis, "Digital Communications", McGraw Hill, New York, 
%		4th Ed., 2000. 
%
%   [2] D. G. Brennan, "Linear Diversity Combining Techniques", Proceedings of 
%		the IRE, vol. 47, June 1959, pp. 1075-1102.

%   Copyright 2006-2021 The MathWorks, Inc.

%% Simulation parameters
% Create comm.BPSKModulator and comm.BPSKDemodulator objects
P = 2; % modulation order
bpskMod = comm.BPSKModulator;
bpskDemod = comm.BPSKDemodulator('OutputDataType','double');

% Convert Eb/No values to SNR values. The output of the BPSK modulator
% generates unit power signals.
SNR = convertSNR(EbNo,"ebno","BitsPerSymbol",1);

% Create a comm.ErrorRate calculator System object to evaluate BER.
errorCalc = comm.ErrorRate;

%%  Pre-allocate variables for speed
ber_MaxRatio = zeros(3,length(EbNo));
h = waitbar(0, 'Percentage Completed');
set(h, 'name', 'Please wait...');
wb = 100/length(EbNo);

%% Loop over EbNo points
for idx = 1:length(EbNo) 
    reset(errorCalc);
    
    % Loop over the number of packets
    for packetIdx = 1:numPackets
        % generate data vector per user/channel
        data = randi([0 P-1], frLen, 1);     
        % modulate data
        tx = bpskMod(data);               
        % Repeat for all Rx antennas
        tx_M = tx(:, ones(1,M));                                                   
        % Create the Rayleigh channel response matrix
        H = (randn(frLen, M) + 1i*randn(frLen, M))/sqrt(2);
        % Received signal for each Rx antenna
        r = awgn(H.*tx_M,SNR(idx));
        % Combiner - assume channel response known at Rx
        z = sum(r .* conj(H), 2);
        % ML Detector (minimum Euclidean distance)
        demod1m = bpskDemod(z); % MR combined 
        % Determine and update BER for current EbNo value
        ber_MaxRatio(:,idx) = errorCalc(data, demod1m);
    end % end of FOR loop for numPackets

    str_bar = [num2str(wb) '% Completed'];
    waitbar(wb/100, h, str_bar);
    wb = wb + 100/length(EbNo);
end  % end of for loop for EbNo
BER1m = ber_MaxRatio(1,:);
close(h);

% [EOF]
