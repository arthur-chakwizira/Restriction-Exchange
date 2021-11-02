function mfs = resex_mc_fit_fm_any_geo(s_simulated, xps)
%Fits unified full compartment model

unit_to_SI = [1e-9 1e-9 1 1 1e-6, 1, 1]; % D1, D2, f1, k, d, M, N

t_ub      = [5     5       1         200           200    inf   inf];
t_lb      = [0       0        0         0             0     0     0];


% t_ub      = [0.1091     1.62        0.7046         200           200 ];
% t_lb      = [0.109       1.61        0.704         0             0];

mu = 1.8;
t_0 = [1         1         0.5           5           5    2*(mu^2)*(mu^2-1)   4*mu^2];


fun = @(x)get_target(x); %objective

% t_lb = []; t_ub = [];
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e4,  'MaxIterations', 1e4,...
    'Display', 'off');
t = lsqnonlin(fun,t_0,t_lb, t_ub, options);

%     A = []; B = []; Aeq = []; beq = []; nonlcon = [];
%     x = patternsearch(fun,x0,A,B,Aeq,beq,lb,ub,nonlcon, options);
%% Objective
    function target = get_target(t_0)
        m = t_0.*unit_to_SI;
        signal = resex_mc_signal_fm_any_geo(m, xps);
               target = (signal)-(s_simulated);
    end
%% Evaluate final results
mfs.signal = get_target(t)+ s_simulated;
mfs.params = t;
mfs.par_names = {'D_1','D_2', 'f_1', 'k', 'd', 'M', 'N'};

end
