function BER2m = ostbc2m_pct(M, frLen, numPackets, EbNo, idx, seed)
%OSTBC2M_PCT  Orthogonal space-time block coding for 2xM antenna
%configurations.
%
%   BER2M = OSTBC2M_PCT(M, FRLEN, NUMPACKETS, EBNO, IDX, SEED) computes the
%   bit-error rate estimates via simulation for an orthogonal space-time
%   block coded configuration using two transmit antennas and M receive
%   antennas for a single Eb/No value, where the frame length, number of
%   packets simulated and the Eb/No range of values are given by FRLEN,
%   NUMPACKETS, and EBNO parameters respectively. IDX is the index of the
%   current Eb/No value and SEED is the random number generator seed.
%
%   The simulation uses the full-rate Alamouti encoding scheme for BPSK
%   modulated symbols with appropriate receiver combining.
%
%   Suggested parameter values:
%       M = 1 or 2; FRLEN = 100; NUMPACKETS = 1000; EBNO = 0:2:20;
%
%   Example:
%       ber22 = ostbc2m(2, 100, 1000, 0:2:20);
%
%   See also MRC1M, OSTBC2M_E, OSTBC4M.

%   References:
%   [1] S. M. Alamouti, "A simple transmit diversity technique for wireless
%       communications", IEEE Journal on Selected Areas in Communications,
%       Vol. 16, No. 8, Oct. 1998, pp. 1451-1458.

%   Copyright 2006-2021 The MathWorks, Inc.

%% Simulation parameters
P = 2;      % Modulation order
N = 2;      % Number of transmit antennas
rate = 1;   % Space-time block code rate
blkLen = 2; % Space-time block code length

% Set the global random stream for repeatability. Use the combined
% multiple recursive generator since it supports substreams.
s = RandStream.create('mrg32k3a', 'seed',seed);
prevStream = RandStream.setGlobalStream(s);

% Create comm.BPSKModulator and comm.BPSKDemodulator System objects
bpskhMod = comm.BPSKModulator;
bpskDemod = comm.BPSKDemodulator('OutputDataType','double');

% Create comm.OSTBCEncoder and comm.OSTBCCombiner System objects
ostbcEnc = comm.OSTBCEncoder;
ostbcComb = comm.OSTBCCombiner('NumReceiveAntennas', M);

% Convert Eb/No values to SNR values. The output of the BPSK modulator
% generates unit power signals.
SNR = convertSNR(EbNo(idx),"ebno","BitsPerSymbol",1);

% Create a comm.ErrorRate calculator System object to evaluate BER.
errorCalc = comm.ErrorRate;

% Pre-allocate variables for speed
ber_Alamouti = zeros(length(EbNo),3);

reset(errorCalc);

% Loop over the number of packets
for packetIdx = 1:numPackets
    % Generate data vector per user/channel
    data = randi([0 P-1], frLen, 1);
    % Modulate data
    tx = bpskhMod(data);
    % Alamouti encoder
    txEnc = ostbcEnc(tx);
    % Create the Rayleigh channel response matrix
    H = (randn(frLen/rate/blkLen, N, M) + ...
        1i*randn(frLen/rate/blkLen, N, M))/sqrt(2);
    %   held constant for blkLen symbol periods
    H = H(kron((1:frLen/rate/blkLen), ones(blkLen,1)),:,:);
    % Received signal with power normalization
    chanOut = squeeze(sum(H .* repmat(txEnc,[1,1,M]), 2))/sqrt(N);
    % Add AWGN
    rxSig = awgn(chanOut,SNR);
    % Alamouti combiner
    rxDec = ostbcComb(rxSig, H);
    % ML Detector (minimum Euclidean distance)
    demod2m = bpskDemod(rxDec);
    % Determine and update BER for current EbNo value
    ber_Alamouti(idx,:) = errorCalc(data, demod2m);
end % End of FOR loop for numPackets


BER2m = ber_Alamouti(idx,:);
RandStream.setGlobalStream(prevStream);

% [EOF]
