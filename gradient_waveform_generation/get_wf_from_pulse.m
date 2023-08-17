function waveform_info = get_wf_from_pulse(waveform_params)
%   This function accepts the "waveform_params" structure which contains the 
%   gradient waveform parameters below.
%   Its output is a free waveform and the corresponding b-value and q-vector
if (nargin == 0)
waveform_params.T = 140e-3; %total encoding time in seconds
waveform_params.dt = 0.1e-3; %time-step in seconds
waveform_params.g_max = 80e-3; %maximum gradient amplitude in Tesla
waveform_params.slew_max = 70; %maximum slew rate in Tesla/seconds
waveform_params.t_start = 10e-3; %start-time of first gradient pulse in seconds
waveform_params.delta = 40e-3; %gradient pulse duration in seconds
waveform_params.interval = 20e-3; %time interval between end of pulse 1 and start of pulse 2 in seconds (interval + delta = DELTA)
waveform_params.time = linspace(0, waveform_params.T,round(waveform_params.T/waveform_params.dt)+1); %time vector = linspace(0, T, N_steps) in seconds
waveform_params.N_steps = length(waveform_params.time); %number of time-steps = round(T/dt)+1
waveform_params.gam_ma = 2*pi*42.6e6; %gamma constant for the 1H nucleus, 42.6e6 Hz/T
[~] = get_wf_from_pulse(waveform_params); 
end
%% Retrieve waveform parameters (all quantities are in SI units)
dt = waveform_params.dt; %time-step, eg. 0.1e-3
t_start = waveform_params.t_start; %start-time of first gradient pulse, eg. 1e-3
N_steps = waveform_params.N_steps; %number of time-steps = round(T/dt)+1
time = waveform_params.time; %time vector = linspace(0, T, N_steps)
gam_ma = waveform_params.gam_ma; %gamma constant for the 1H nucleus, 42.6e6*2pi
delta = waveform_params.delta; %gradient pulse duration, eg. 10e-3
interval = waveform_params.interval; %time interval between end of pulse 1 and start of pulse 2, eg. 6e-3
%% Form a pulse
waveform = zeros(N_steps, 1); %initialise waveform
N_pulse = round(waveform_params.delta/waveform_params.dt) +1; %pulse duration in terms of time points
pulse = waveform_params.pulse;

tol = 5; %will use this to compare floating points because "==" is problematic. We regard two points within tol*eps of each other as equal.
%sick and tired of pulse mispositioning errors
block1 = (time > t_start | ismembertol(time,t_start, tol*eps(t_start))) &...
 (time < t_start + delta | ismembertol(time,t_start + delta, tol*eps(t_start + delta)));
offset = N_pulse-sum(block1); last_true = find(block1, 1, 'last');
if offset > 0; block1(last_true+1:last_true+offset) = true; end
if offset < 0; block1(last_true:-1:last_true+offset+1) = false; end

block2 = (time > t_start+delta+interval | ismembertol(time,t_start+delta+interval, tol*eps(t_start+delta+interval))) &... %flip pulse 1 and place it in ...
 (time < t_start+2*delta+interval) | ismembertol(time,t_start+2*delta+interval, tol*eps(t_start+2*delta+interval));
offset = N_pulse-sum(block2); 
last_true = find(block2, 1, 'last');
if offset > 0; block2(last_true+1:last_true+offset) = true; end
if offset < 0; block2(last_true:-1:last_true+offset+1) = false; end

%% Place the pulse(s) in the correct positions in time
waveform(block1) = pulse; %place pulse 1 in the interval [t_start, t_start + delta]
waveform(block2) = -flip(pulse);  %[t_start+delta+DELTA, t_start+2*delta+DELTA]

%% Crop waveform to remove unnecessary zeros
% [~, waveform] = N_from_gwf(waveform);
%% Determine q-vector and b-value
q_vector = gam_ma*cumsum(waveform)*dt;
b_value = sum(q_vector.^2)*dt;
%% Pack output
waveform_info.waveform = waveform;
waveform_info.b_value = b_value;
waveform_info.q_vector = q_vector;

end
