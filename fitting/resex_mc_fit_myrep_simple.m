function mfs = resex_mc_fit_myrep_simple(s_simulated, xps)
%Fits unified full compartment model

unit_to_SI = [1e-15   1e-9    1    1]; % R, D2, f1, k

t_ub      = [inf    0.7        0.75        100];
t_lb      = [0      0.5        0.65        0  ];


% m_ub = (t_ub).*unit_to_SI;
% m_lb = (t_lb).*unit_to_SI;
% m_guess = msf_fit_random_guess(@resex_mc_signal_myrep, s_simulated, xps, m_lb, m_ub);
% t_0 = m_guess./unit_to_SI;

c = 7/1536;
R_0 = c*(5^4)/1; %d = 5, D0 = 1
t_0 = [R_0         2         0.5           5];


fun = @(x)get_target(x); %objective

% t_lb = []; t_ub = [];
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e6,  'MaxIterations', 1e6,...
    'Display', 'off');
t = lsqnonlin(fun,t_0,t_lb, t_ub, options);

%     A = []; B = []; Aeq = []; beq = []; nonlcon = [];
%     x = patternsearch(fun,x0,A,B,Aeq,beq,lb,ub,nonlcon, options);
%% Objective
    function target = get_target(t_0)
        m = t_0.*unit_to_SI;
        signal = resex_mc_signal_myrep_simple(m, xps);
               target = (signal)-(s_simulated);
    end
%% Evaluate final results
mfs.signal = get_target(t)+ s_simulated;
mfs.params = t;
mfs.par_names = {'R','D_2', 'f_1', 'k'};

end
