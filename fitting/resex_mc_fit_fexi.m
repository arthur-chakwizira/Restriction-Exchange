function mfs = resex_mc_fit_fexi(s_simulated, xps)
%Fits restricted unified model representation

bf = xps.b_1;
bd = xps.b_2;
tm = xps.tm;
DELTA_1 = xps.DELTA_1;

% protocol.done = false;
%ADC, sig_ma, AXR
%S = S_0*exo(-(TE + TEf)/T2)*exp(-tm/T1)*exp(-bf*ADC)*exp(-bd*ADC*(1 - sig_ma*exp(-AXR*tm)))
x0 = [1,      0.3,  10];
lb = [eps,        eps,    0];
ub = [5,        1,    20];

fun = @(x)get_target(x); %objective
% lb = []; ub = [];

options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e4, ...
    'MaxIterations', 1e4, 'Display', 'off');
x = lsqnonlin(fun,x0,lb, ub, options);

%     A = []; B = []; Aeq = []; beq = []; nonlcon = [];
%     x = patternsearch(fun,x0,A,B,Aeq,beq,lb,ub,nonlcon, options);
%% Objective
  function target = get_target(x0)

%Model parameters
ADC = x0(1)*1e-9;
sig_ma = x0(2);
AXR = x0(3);

signal = exp(-bf*ADC).*exp(-bd*ADC.*(1 - sig_ma*exp(-AXR*tm)));

        target = (signal)-(s_simulated);
    end
%% Evaluate final results
% protocol.done = true;
mfs.signal = get_target(x)+s_simulated;
k = x(3)/(1 + x(3)*DELTA_1);
x(4) = k;
mfs.params = x;
mfs.par_names = {'ADC','Sigma', 'AXR', 'k'};

end
