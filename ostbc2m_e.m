function BER2m_e = ostbc2m_e(M, frLen, numPackets, EbNo, pLen)
%OSTBC2M_E  Orthogonal space-time block coding with channel estimation for
%   2xM antenna configurations. 
%
%   BER2M_E = STBC2M_E(M, FRLEN, NUMPACKETS, EBNOVEC, PLEN) computes the
%   bit-error rate estimates via simulation for an orthogonal space-time
%   block coded configuration using two transmit antennas and M receive
%   antennas, where the frame length, number of packets simulated, Eb/No
%   range of values and the number of pilot symbols prepended per frame are
%   given by FRLEN, NUMPACKETS, EBNOVEC and PLEN parameters respectively.
%
%   The simulation uses the full-rate Alamouti encoding scheme for BPSK
%   modulated symbols with appropriate receiver combining. It uses the
%   pilot-aided Minimum-Mean-Square-Error (MMSE) method for estimating the
%   channel coefficients at the receiver. It is assumed the channel is
%   slowly fading (i.e. it remains constant for the whole frame of data and
%   changes independently from one frame to the other).
%
%   Suggested parameter values:
%       M = 1 or 2; FRLEN = 100; NUMPACKETS = 1000; EBNOVEC = 0:2:20, PLEN = 8;
%
%   Example:
%       ber22_e = ostbc2m_e(2, 100, 1000, 0:2:20, 8);
%
%   See also OSTBC2M, OSTBC4M, MRC1M.

%   References:
%   [1] A.F. Naguib, V. Tarokh, N. Seshadri, and A.R. Calderbank, "Space-time
%       codes for high data rate wireless communication: Mismatch analysis", 
%       Proceedings of IEEE International Conf. on Communications, 
%       pp. 309-313, June 1997.        
%
%   [2] S. M. Alamouti, "A simple transmit diversity technique for wireless 
%       communications", IEEE Journal on Selected Areas in Communications, 
%       Vol. 16, No. 8, Oct. 1998, pp. 1451-1458.

%   Copyright 2006-2021 The MathWorks, Inc.

%% Simulation parameters
P = 2;      % Modulation order
N = 2;      % Number of transmit antennas
rate = 1;   % Space-time block code rate

% Create comm.BPSKModulator and comm.BPSKDemodulator System objects
bpskMod = comm.BPSKModulator;
bpskDemod = comm.BPSKDemodulator( ...
    'OutputDataType','double');

% Create comm.OSTBCEncoder and comm.OSTBCCombiner System objects
ostbcEnc = comm.OSTBCEncoder;
ostbcComb = comm.OSTBCCombiner( ...
    'NumReceiveAntennas', M);

% Create a comm.MIMOChannel System object to simulate the 2xM spatially
% independent flat-fading Rayleigh channel
chanMIMO = comm.MIMOChannel( ...
    'MaximumDopplerShift', 0, ...
    'SpatialCorrelationSpecification', 'None', ...
    'NumReceiveAntennas', M);

% Convert Eb/No values to SNR values. The output of the BPSK modulator
% generates unit power signals.
SNR = convertSNR(EbNo,"ebno","BitsPerSymbol",1);

% Create a comm.ErrorRate calculator System object to evaluate BER.
errorCalc = comm.ErrorRate;

% Pilot sequences - orthogonal set over N
W = hadamard(pLen); % order gives the number of pilot symbols prepended/frame
pilots = W(:, 1:N); 

%%  Pre-allocate variables for speed
HEst = zeros(frLen/rate, N, M); 
ber_Estimate  = zeros(3,length(EbNo));
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
        % Prepend pilot symbols for each frame
        txSig = [pilots; txEnc];        
        % Pass through the 2xM channel
        reset(chanMIMO);
        chanOut = chanMIMO(txSig);
        % Add AWGN
        rxSig = awgn(chanOut,SNR(idx));
        % Channel Estimation
        %   For each link => N*M estimates
        HEst(1,:,:) = pilots(:,:).' * rxSig(1:pLen, :) / pLen;
        %   held constant for the whole frame
        HEst = HEst(ones(frLen/rate, 1), :, :);
        % Alamouti combiner
        rxDec = ostbcComb(rxSig(pLen+1:end,:), HEst);
        % ML Detector (minimum Euclidean distance)
        demod2m_e = bpskDemod(rxDec); 
        % Determine and update BER
        ber_Estimate (:,idx) = errorCalc(data, demod2m_e);
    end % End of FOR loop for numPackets

    str_bar = [num2str(wb) '% Completed'];
    waitbar(wb/100, h, str_bar);
    wb = wb + 100/length(EbNo);
end  % End of for loop for EbNo

BER2m_e = ber_Estimate(1,:);
close(h);
    
% [EOF]
