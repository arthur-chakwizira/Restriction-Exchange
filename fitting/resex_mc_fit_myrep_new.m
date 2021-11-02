function mfs = resex_mc_fit_myrep_new(s_simulated, xps)
%Fits unified second-order representation
%Does not use Gamma and Vomega at all

unit_to_SI = [1e-9 1e-9 1 1 1e-6]; % D1, D2, f1, k, d

t_ub      = [3    3        1        100           100 ];
t_lb      = [0   0       0        0           0 ];

t_0 = [1         2         0.5           10           20];


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
        signal = resex_mc_signal_myrep_new(m, xps);
               target = (signal)-(s_simulated);
    end
%% Evaluate final results
mfs.signal = get_target(t)+ s_simulated;
mfs.params = t;
mfs.par_names = {'D_1','D_2', 'f_1', 'k', 'd'};

end
