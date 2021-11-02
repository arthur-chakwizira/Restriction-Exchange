function mfs = resex_mc_fit_myrep(s_simulated, xps)
%Fits unified full compartment model
if (~iscolumn(s_simulated)); error("Fatal error: input signal has wrong dimensions"); end

unit_to_SI = [1e-9 1e-9 1 1 1e-6]; % D1, D2, f1, k, d

t_ub      = [3    3        1        100           100 ];
t_lb      = [0   0       0        0           0 ];

% t_lb = -t_ub;

% t_ub      = [0.6    1.6       0.8       100          6];
% t_lb      = [0.4   1.4    0.6      0         4];
% 
% m_ub = (t_ub).*unit_to_SI;
% m_lb = (t_lb).*unit_to_SI;
% m_guess = msf_fit_random_guess(@resex_mc_signal_myrep, s_simulated, xps, m_lb, m_ub);
% t_0 = m_guess./unit_to_SI;
% 
% 
t_0 = [0.5         1.5         0.5           5           5];


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
        signal = resex_mc_signal_myrep(m, xps);
               target = (signal)-(s_simulated);
    end
%% Evaluate final results
mfs.signal = get_target(t)+ s_simulated;
mfs.params = t;
mfs.par_names = {'D_1','D_2', 'f_1', 'k', 'd'};

end
