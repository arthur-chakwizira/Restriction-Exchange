function mfs = resex_mc_fit_km(s_simulated, xps)
%Fits modified Kärger model
if (~iscolumn(s_simulated)); error("Fatal error: input signal has wrong dimensions"); end

unit_to_SI = [1e-9 1e-9 1 1]; % D1, D2, f1, k

t_ub      = [1    2      1         20 ];
t_lb      = [0       0        0         0 ];



t_0 = [1         1         0.5           5];


fun = @(x)get_target(x); %objective

options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e4,  'MaxIterations', 1e4,...
    'Display', 'off');
t = lsqnonlin(fun,t_0,t_lb, t_ub, options);

%% Objective
    function target = get_target(t_0)
        m = t_0.*unit_to_SI;
        signal = resex_mc_signal_general_km(m, xps);
               target = (signal)-(s_simulated);
    end
%% Evaluate final results
mfs.signal = get_target(t)+ s_simulated;
mfs.params = t;
mfs.par_names = {'D_1','D_2', 'f_1', 'k',};

end
