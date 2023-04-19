function BER1m = mrc1m_pct(M, frLen, numPackets, EbNo, idx, seed)
%MRC1M_PCT  Maximal-Ratio Combining for 1xM antenna configurations.
%
%   BER1M = MRC1M_PCT(M, FRLEN, NUMPACKETS, EBNO, IDX, SEED) computes the
%   bit-error rate estimates via simulation for a Maximal-Ratio Combined
%   configuration using one transmit antenna and M receive antennas for a
%   single Eb/No value, where the frame length, number of packets simulated
%   and the Eb/No range of values are given by FRLEN, NUMPACKETS, EBNO
%   parameters respectively. IDX is the index of the current Eb/No value
%   and SEED is the random number generator seed.
%
%   The simulation uses BPSK modulated symbols with appropriate receiver
%   combining.
%
%   Suggested parameter values:
%       M = 1 to 4; FRLEN = 100; NUMPACKETS = 1000; EBNO = 0:2:20;
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

% Set the global random stream for repeatability. Use the combined
% multiple recursive generator since it supports substreams.
s = RandStream.create('mrg32k3a', 'seed',seed);
prevStream = RandStream.setGlobalStream(s);

% Convert Eb/No values to SNR values. The output of the BPSK modulator
% generates unit power signals.
SNR = convertSNR(EbNo(idx),"ebno","BitsPerSymbol",1);

% Create a comm.ErrorRate calculator System object to evaluate BER.
errorCalc = comm.ErrorRate;

%  Pre-allocate variables for speed
ber_MaxRatio = zeros(length(EbNo),3);

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
    r = awgn(H.*tx_M,SNR);
    % Combiner - assume channel response known at Rx
    z = sum(r .* conj(H), 2);
    % ML Detector (minimum Euclidean distance)
    demod1m = bpskDemod(z); % MR combined
    % Determine and update BER for current EbNo value
    ber_MaxRatio(idx,:) = errorCalc(data, demod1m);
end % end of FOR loop for numPackets



BER1m = ber_MaxRatio(idx,:);
RandStream.setGlobalStream(prevStream);


% [EOF]
