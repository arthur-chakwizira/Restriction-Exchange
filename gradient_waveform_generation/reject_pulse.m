function reject = reject_pulse(pulse, waveform_params)
%Checks if waveform obeys constraints
reject = false;

slew_max = waveform_params.slew_max;
g_max = waveform_params.g_max;
dt = waveform_params.dt;

dgdt = diff(pulse)/dt;
if max(dgdt(:)) > slew_max; reject = true; end
% %Add more constraints...
if max(abs(pulse(:))) > g_max; reject = true; end
% 
if any(isnan((pulse(:)))); reject = true; end

end