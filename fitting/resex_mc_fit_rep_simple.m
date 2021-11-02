function mfs = resex_mc_fit_rep_simple(s_simulated, xps)
%Fits unified model representation

%Bunching up V_R, C_DR and V_D into one parameter V
 unit_to_SI = [1e-9 1e-15 1e-18 1];  %E_D, E_R, V, k
             % E_D,     E_R,        V            k
t_ub      = [ 5-eps    10e4        1e9       100-eps];
t_lb =  zeros(size(t_ub));
%           []     [um^2/ms] [um^2*ms] [um^4/ms^2]

x0 =   [ 1,       40,        10,           5];


fun = @(x)get_target(x); %objective

% lb = []; ub = [];
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e4,  'MaxIterations', 1e4,...
    'Display', 'off');
x = lsqnonlin(fun,x0,t_lb, t_ub, options);

%     A = []; B = []; Aeq = []; beq = []; nonlcon = [];
%     x = patternsearch(fun,x0,A,B,Aeq,beq,lb,ub,nonlcon, options);
%% Objective
    function target = get_target(x0)                  
                m = x0.*unit_to_SI;
               signal = resex_mc_signal_rep_simple(m, xps);
               target = (signal)-(s_simulated);
   end
%% Evaluate final results
mfs.signal = get_target(x)+ s_simulated;
mfs.params = x;
mfs.par_names = {'E_D','E_R', 'V', 'k'};

end
