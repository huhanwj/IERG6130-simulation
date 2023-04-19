function BER2m = ostbc2m(M, frLen, numPackets, EbNo)
%OSTBC2M  Orthogonal space-time block coding for 2xM antenna configurations.
%
%   BER2M = OSTBC2M(M, FRLEN, NUMPACKETS, EBNOVEC) computes the bit-error
%   rate estimates via simulation for an orthogonal space-time block coded
%	configuration using two transmit antennas and M receive antennas, where
%	the frame length, number of packets simulated and the Eb/No range of
%	values are given by FRLEN, NUMPACKETS, and EBNOVEC parameters
%	respectively.
%
%   The simulation uses the full-rate Alamouti encoding scheme for BPSK 
%   modulated symbols with appropriate receiver combining.
%
%   Suggested parameter values:
%       M = 1 or 2; FRLEN = 100; NUMPACKETS = 1000; EBNOVEC = 0:2:20;
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

% Create comm.BPSKModulator and comm.BPSKDemodulator System objects
bpskMod = comm.BPSKModulator;
bpskDemod = comm.BPSKDemodulator('OutputDataType','double');

% Create comm.OSTBCEncoder and comm.OSTBCCombiner System objects
ostbcEnc = comm.OSTBCEncoder;
ostbcComb = comm.OSTBCCombiner('NumReceiveAntennas', M);

% Convert Eb/No values to SNR values. The output of the BPSK modulator
% generates unit power signals.
SNR = convertSNR(EbNo,"ebno","BitsPerSymbol",1);

% Create a comm.ErrorRate calculator System object to evaluate BER.
errorCalc = comm.ErrorRate;

%% Pre-allocate variables for speed
ber_Alamouti = zeros(3,length(EbNo));
h = waitbar(0, 'Percentage Completed');
set(h, 'name', 'Please wait...');
wb = 100/length(EbNo);

%% Loop over EbNo points
for idx = 1:length(EbNo)   
    reset(errorCalc);
    
    % Loop over the number of packets
    for packetIdx = 1:numPackets
        % Generate data vector per user/channel
        data = randi([0 P-1], frLen, 1);
        % Modulate data
        tx = bpskMod(data); 
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
        rxSig = awgn(chanOut,SNR(idx));
        % Alamouti combiner
        rxDec = ostbcComb(rxSig, H);
        % ML Detector (minimum Euclidean distance)
        demod2m = bpskDemod(rxDec);
        % Determine and update BER for current EbNo value
        ber_Alamouti(:,idx) =  errorCalc(data, demod2m); 
    end % End of FOR loop for numPackets

    str_bar = [num2str(wb) '% Completed'];
    waitbar(wb/100, h, str_bar);
    wb = wb + 100/length(EbNo);
end  % End of for loop for EbNo

BER2m = ber_Alamouti(1,:);
close(h);
    
% [EOF]
