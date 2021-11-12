function mfs = resex_mc_fit_ning(s_simulated, xps)
%Fits representatino from Ning et al. 2018

unit_to_SI = [1e-9 1e-9 1 1 1e-6]; % D1, D2, f1, k, d

t_ub      = [3.5     5        1         100         ];
t_lb      = t_ub*0;

t_0 = [1         2         0.5           10         ];


fun = @(x)get_target(x); %objective

% lb = []; ub = [];
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e4,  'MaxIterations', 1e4,...
    'Display', 'off');
t = lsqnonlin(fun,t_0,t_lb, t_ub, options);

%     A = []; B = []; Aeq = []; beq = []; nonlcon = [];
%     x = patternsearch(fun,x0,A,B,Aeq,beq,lb,ub,nonlcon, options);
%% Objective
    function target = get_target(t_0)
        m = t_0.*unit_to_SI;
        signal = resex_mc_signal_ning(m, xps);
               target = (signal)-(s_simulated);
    end
%% Evaluate final results
mfs.signal = get_target(t)+ s_simulated;
mfs.params = t;
mfs.par_names = {'D_1','D_2', 'f_1', 'k'};

end
