function [cfo,sco] = trackingPlotSCOCFOEstimates(cpe,peg,cfg)
%trackingPlotSCOCFOEstimates Featured example helper function
%
%   Plots PEG and CPE and estimates SCO and CFO using linear regression

%   Copyright 2016-2021 The MathWorks, Inc.

if isa(cfg,'wlanVHTConfig')
    ofdmInfo = wlanVHTOFDMInfo('VHT-Data',cfg);
elseif isa(cfg,'wlanHTConfig')
    ofdmInfo = wlanHTOFDMInfo('HT-Data',cfg);
else
    ofdmInfo = wlanNonHTOFDMInfo('NonHT-Data');
end
N = ofdmInfo.FFTLength;
Ns = N+ofdmInfo.CPLength; % Number of samples in an OFDM symbol
sr = wlanSampleRate(cfg);

figure;
subplot(211);
unwrappedPEG = (unwrap(peg*N)/N);
unwrappedPEG(peg==0) = 0; % Make sure 0 PEG is kept after unwrapping - signifies symbol not demodulated
plot(0:size(peg,1)-1,unwrappedPEG,'xb');
hold on;
x = (0:(size(peg,1)-1)).';
X = [ones(length(x),1) x];
b = X\unwrappedPEG; % Unwrap over skip/dup
plot(0:size(cpe,1)-1,X*b,'-');
sco = (b(2)/(2*pi*Ns/N))*1e6;
legend('Measured PEG','Linear fit','location','northeast');
title('Measured phase error gradient');
ylabel('Rad/subcarrier');
xlabel('OFDM symbol');
grid on;
dim = [0.15 0.625 0.22 0.07];
annotation('textbox',dim,'String',sprintf('%.1f PPM estimated SCO',sco),'FitBoxToText','on','BackgroundColor',[1 1 1]);
subplot(212);
plot(0:size(peg,1)-1,cpe,'or');
hold on;
x = (0:(size(cpe,1)-1)).';
X = [ones(length(x),1) x];
b = X\unwrap(cpe);
plot(0:size(peg,1)-1,X*b,'-r');
cfo = (b(2)/Ns)*sr*1/(2*pi);
legend('Measured CPE','Linear fit','location','northeast');
title('Measured common phase error');
ylabel('Rad');
xlabel('OFDM symbol');
grid on;
dim = [0.15 0.15 0.3 0.07];
annotation('textbox',dim,'String',sprintf('%.1f Hz estimated residual CFO',cfo),'FitBoxToText','on','BackgroundColor',[1 1 1]);

end
