function mfs = resex_mc_fit_rep(s_simulated, xps)
%Fits unified model representation
if (~iscolumn(s_simulated)); error("Fatal error: input signal has wrong dimensions"); end

unit_to_SI = [1e-9 1e-15 1e-18 1e-24 1e-30  1];  %E_D, E_R, V_D, C_DR,  V_R, k

% protocol.done = false;
%             E_D,     E_R,        V_D,        C_DR,            V_R,            k
t_ub      = [5     1e5            1e9          1e9              1e9            100];
t_lb      = zeros(size(t_ub)) + eps;
%            [um^2/ms]  [um^2*ms]  [um^4/ms^2]      [um^4]      [um^4*ms^2]      [/s]


% 
% m_ub = (t_ub).*unit_to_SI;
% m_lb = (t_lb).*unit_to_SI;
% m_guess = msf_fit_random_guess(@resex_mc_signal_rep, s_simulated, xps, m_lb, m_ub);
% x0 = m_guess./unit_to_SI;


x0 = [1      20             9          500               1500            5];


fun = @(x)get_target(x); %objective

% lb = []; 
t_ub = [];
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e4,  'MaxIterations', 1e4,...
    'Display', 'off');
x = lsqnonlin(fun,x0,t_lb, t_ub, options);

%     A = []; B = []; Aeq = []; beq = []; nonlcon = [];
%     x = patternsearch(fun,x0,A,B,Aeq,beq,lb,ub,nonlcon, options);
%% Objective
    function target = get_target(x0)
                  
                m = x0.*unit_to_SI;
               signal = resex_mc_signal_rep(m, xps);
               target = (signal)-(s_simulated);
   end
%% Evaluate final results
mfs.signal = get_target(x)+ s_simulated;
mfs.params = x;
mfs.par_names = {'E_D','E_R', 'V_D', 'C_DR', 'V_R', 'k'};

end
