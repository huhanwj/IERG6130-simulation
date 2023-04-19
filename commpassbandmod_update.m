function commpassbandmod_update
% This function will set the power on the nonlinearity to either 1 or
% 3. This effectively toggles the nonlinearity on or off.

% Copyright 2007-2016 The MathWorks, Inc.

nonlinOrderSource = [bdroot, '/Interferer/Constant'];
nonlinSubsys      = [bdroot, '/Interferer'];
BER_reset         = [bdroot, '/Compute BER/Constant'];

pow =  get_param(nonlinOrderSource, 'Value');
if pow == '1'
    pow = '3';
else
    pow = '1';
end

set_param(nonlinSubsys, 'nonlinOrder', pow);
set_param(nonlinOrderSource, 'Value', pow);

% Reset the BER
set_param(BER_reset,'Value','1');
pause(0.1)
set_param(BER_reset,'Value','0');
