function dock_figures(swtch)
if nargin < 1%isstring(swtch)|| ischar(swtch)
    h = findobj('type', 'figure'); %all open figures
    for c_h = 1:numel(h)
        set(h(c_h), 'windowStyle', 'docked');
    end
else
    if swtch
        set(0,'DefaultFigureWindowStyle','docked')%dock all figures we open
    else
        set(0,'DefaultFigureWindowStyle','normal')%dock all figures we open
    end
end
end