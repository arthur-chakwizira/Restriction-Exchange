function mfs = resex_mc_fit_monoexp(s_simulated, xps)
%Fits unified model representation

unit_to_SI = 1e-9;

ub      = 5;
lb      = 0;

x0 = 2;

fun = @(x)get_target(x); %objective

% lb = []; ub = [];
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e4,  'MaxIterations', 1e4,...
    'Display', 'off');
x = lsqnonlin(fun,x0,lb, ub, options);

%     A = []; B = []; Aeq = []; beq = []; nonlcon = [];
%     x = patternsearch(fun,x0,A,B,Aeq,beq,lb,ub,nonlcon, options);
%% Objective
    function target = get_target(x0)
                  
                m = x0.*unit_to_SI;
               signal = exp(-m.*xps.b);
               target = (signal)-(s_simulated);
   end
%% Evaluate final results
mfs.signal = get_target(x)+ s_simulated;
mfs.params = x;
mfs.par_names = {'D'};

end
