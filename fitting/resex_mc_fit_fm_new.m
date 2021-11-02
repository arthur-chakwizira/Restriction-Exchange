function mfs = resex_mc_fit_fm_new(s_simulated, xps)
%Fits unified full compartment model

%we are about to do something bizarre
max_n = 0;
for l = 1:length(xps.gwf)
   max_n = max(max_n, numel(xps.gwf{l})); 
end

xps.max_n = max_n;

unit_to_SI = ones(max_n*2 + 3, 1)*1e-9;
t_ub = ones(1, max_n*2 + 3)*3.5;
t_lb = zeros(1, max_n*2 + 3);


unit_to_SI(end-2:end) = [1 1 1e-6]; % D1, D2, f1, k, d

t_ub(end-2:end)      = [1    100           100 ];

t_0 = [ones(1, max_n)*1   ones(1, max_n)*2    0.5           10           10];


fun = @(x)get_target(x); %objective

t_lb = []; t_ub = [];
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e4,  'MaxIterations', 1e4,...
    'Display', 'iter');
t = lsqnonlin(fun,t_0,t_lb, t_ub, options);

%     A = []; B = []; Aeq = []; beq = []; nonlcon = [];
%     x = patternsearch(fun,x0,A,B,Aeq,beq,lb,ub,nonlcon, options);
%% Objective
    function target = get_target(t_0)
        m = t_0.*unit_to_SI;
        signal = resex_mc_signal_fm_new(m, xps);
               target = (signal)-(s_simulated);
    end
%% Evaluate final results
mfs.signal = get_target(t)+ s_simulated;
mfs.params = t;
mfs.max_n = max_n;
mfs.par_names = {'D_1','D_2', 'f_1', 'k'};

end
