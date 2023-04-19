function BER4m = ostbc4m(M, frLen, numPackets, EbNo)
%OSTBC4M  Orthogonal space-time block coding for 4xM antenna configurations.
%
%   BER4M = OSTBC4M(M, FRLEN, NUMPACKETS, EBNOVEC) computes the bit-error
%   rate estimates via simulation for an orthogonal space-time block coded
%	configuration using four transmit antennas and M receive antennas,
%	where the frame length, number of packets simulated and the Eb/No range
%	of values are given by FRLEN, NUMPACKETS, and EBNOVEC parameters
%	respectively.
%
%   The simulation uses a half-rate orthogonal STBC encoding scheme with
%   QPSK modulated symbols to achieve a 1 bit/sec/Hz throughput for the
%   channel. Appropriate combining is performed at the receiver to account
%   for the multiple transmitter antennas.
%
%   Suggested parameter values:
%       M = 1 to 4; FRLEN = 100; NUMPACKETS = 1000; EBNOVEC = 0:2:20;
%
%   Example:
%       ber42 = ostbc4m(2, 100, 1000, 0:2:20);
%
%   See also MRC1M, OSTBC2M.

%   References:
%   [1] S. M. Alamouti, "A simple transmit diversity technique for wireless 
%       communications", IEEE Journal on Selected Areas in Communications, 
%       Vol. 16, No. 8, Oct. 1998, pp. 1451-1458.
%
%   [2] V. Tarokh, H. Jafarkhami, and A.R. Calderbank, "Space-time block 
%       codes from orthogonal designs", IEEE Transactions on Information 
%       Theory, Vol. 45, No. 5, Jul. 1999, pp. 1456-1467.
%
%   [3] V. Tarokh, H. Jafarkhami, and A.R. Calderbank, "Space-time block 
%       codes for wireless communications: Performance results", IEEE 
%       Journal on Selected Areas in Communications, Vol. 17,  No. 3, 
%       Mar. 1999, pp. 451-460.

%   Copyright 2006-2021 The MathWorks, Inc.

%% Simulation parameters
P = 4;      % Modulation order
N = 4;      % Number of transmit antennas
rate = 0.5; % Space-time block code rate
blkLen = 8; % Space-time block code length

% Create comm.QPSKModulator and comm.QPSKDemodulator System objects
qpskMod = comm.QPSKModulator('PhaseOffset', 0);
qpskDemod = comm.QPSKDemodulator('PhaseOffset', 0, 'OutputDataType', 'double');

% Convert Eb/No values to SNR values. The output of the BPSK modulator
% generates unit power signals.
SNR = convertSNR(EbNo,"ebno","BitsPerSymbol",1);

% Create a comm.ErrorRate calculator System object to compare decoded bits
% to the original transmitted bits.
errorCalc = comm.ErrorRate;

%% Pre-allocate variables for speed
txEnc = zeros(frLen/rate, N); 
H = zeros(frLen/rate, N, M);
z = zeros(frLen, M); z1 = zeros(frLen/N, M); z2 = z1; z3 = z1; z4 = z1;
ber_Estimate = zeros(3,length(EbNo));

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
        tx = qpskMod(data); 
        % Space-Time Block Encoder - G4, 1/2 rate
        %   G4Half = [s1 s2 s3 s4;-s2 s1 -s4 s3;-s3 s4 s1 -s2;-s4 -s3 s2 s1];
        %   G4 = [G4Half; conj(G4Half)];
        s1 = tx(1:N:end); s2 = tx(2:N:end); s3 = tx(3:N:end); s4 = tx(4:N:end);
        txEnc(1:blkLen:end, :) = [ s1  s2  s3  s4];
        txEnc(2:blkLen:end, :) = [-s2  s1 -s4  s3];
        txEnc(3:blkLen:end, :) = [-s3  s4  s1 -s2];
        txEnc(4:blkLen:end, :) = [-s4 -s3  s2  s1];
        for i = 1:4
            txEnc(i+4:blkLen:end, :) = conj(txEnc(i:blkLen:end, :));
        end          
        % Create the Rayleigh channel response matrix
        H(1:blkLen:end, :, :) = (randn(frLen/rate/blkLen, N, M) + ...
                              1i*randn(frLen/rate/blkLen, N, M))/sqrt(2);
        %   held constant for blkLen symbol periods
        for i = 2:blkLen
            H(i:blkLen:end, :, :) = H(1:blkLen:end, :, :); 
        end
        % Received signal with power normalization
        chanOut = squeeze(sum(H .* repmat(txEnc,[1,1,M]), 2))/sqrt(N);
        % Add AWGN
        r = awgn(chanOut,SNR(idx));
        % Combiner - assume channel response known at Rx
        hidx = 1:blkLen:length(H);
        for i = 1:M
            z1(:, i) = r(1:blkLen:end, i).* conj(H(hidx, 1, i)) + ...
                       r(2:blkLen:end, i).* conj(H(hidx, 2, i)) + ...
                       r(3:blkLen:end, i).* conj(H(hidx, 3, i)) + ...
                       r(4:blkLen:end, i).* conj(H(hidx, 4, i)) + ...
                       conj(r(5:blkLen:end, i)).* H(hidx, 1, i) + ...
                       conj(r(6:blkLen:end, i)).* H(hidx, 2, i) + ...
                       conj(r(7:blkLen:end, i)).* H(hidx, 3, i) + ...
                       conj(r(8:blkLen:end, i)).* H(hidx, 4, i);
    
            z2(:, i) = r(1:blkLen:end, i).* conj(H(hidx, 2, i)) - ...
                       r(2:blkLen:end, i).* conj(H(hidx, 1, i)) - ...
                       r(3:blkLen:end, i).* conj(H(hidx, 4, i)) + ...
                       r(4:blkLen:end, i).* conj(H(hidx, 3, i)) + ...
                       conj(r(5:blkLen:end, i)).* H(hidx, 2, i) - ...
                       conj(r(6:blkLen:end, i)).* H(hidx, 1, i) - ...
                       conj(r(7:blkLen:end, i)).* H(hidx, 4, i) + ...
                       conj(r(8:blkLen:end, i)).* H(hidx, 3, i);
    
            z3(:, i) = r(1:blkLen:end, i).* conj(H(hidx, 3, i)) + ...
                       r(2:blkLen:end, i).* conj(H(hidx, 4, i)) - ...
                       r(3:blkLen:end, i).* conj(H(hidx, 1, i)) - ...
                       r(4:blkLen:end, i).* conj(H(hidx, 2, i)) + ...
                       conj(r(5:blkLen:end, i)).* H(hidx, 3, i) + ...
                       conj(r(6:blkLen:end, i)).* H(hidx, 4, i) - ...
                       conj(r(7:blkLen:end, i)).* H(hidx, 1, i) - ...
                       conj(r(8:blkLen:end, i)).* H(hidx, 2, i);
    
            z4(:, i) = r(1:blkLen:end, i).* conj(H(hidx, 4, i)) - ...
                       r(2:blkLen:end, i).* conj(H(hidx, 3, i)) + ...
                       r(3:blkLen:end, i).* conj(H(hidx, 2, i)) - ...
                       r(4:blkLen:end, i).* conj(H(hidx, 1, i)) + ... 
                       conj(r(5:blkLen:end, i)).* H(hidx, 4, i) - ...
                       conj(r(6:blkLen:end, i)).* H(hidx, 3, i) + ...
                       conj(r(7:blkLen:end, i)).* H(hidx, 2, i) - ...
                       conj(r(8:blkLen:end, i)).* H(hidx, 1, i);
        end
        z(1:N:end, :) = z1; z(2:N:end, :) = z2;
        z(3:N:end, :) = z3; z(4:N:end, :) = z4;

        % ML Detector (minimum Euclidean distance)
        demod4m = qpskDemod(sum(z, 2));

        % Convert symbols to bits
        txBits = int2bit(data, log2(P));
        rxBits = int2bit(demod4m, log2(P));

        % Determine and update BER
        ber_Estimate(:,idx) = errorCalc(txBits, rxBits);
    end % End of FOR loop for numPackets

    str_bar = [num2str(wb) '% Completed'];
    waitbar(wb/100, h, str_bar);
    wb = wb + 100/length(EbNo);

end  % End of for loop for EbNo

BER4m = ber_Estimate(1,:);
close(h);

% [EOF]
