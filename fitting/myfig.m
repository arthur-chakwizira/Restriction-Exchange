function ax = myfig(x,y)
%can generate single axes or subplots on x*y grid
%ax is a cell is x and y are supplied
if nargin < 1
    f = figure();
    set(f, 'Color', 'w')
    ax = axes(f, 'FontName', 'Helvetica');
    xlabel('', 'Interpreter', 'tex', 'fontsize', 40)
    ylabel('', 'Interpreter', 'tex', 'fontsize', 40)
    % legend([], 'Interpreter', 'latex', 'fontsize', 20)
    title('', 'Interpreter', 'tex', 'fontsize', 40)
    grid minor
    set(ax, 'Box', 'on', 'fontsize', 20, 'LineWidth', 1)
    set(ax,'tickdir', 'out', 'ticklength', [0.03 0.05], 'linewidth', 1.5, 'fontsize', 20, 'fontname', 'arial', 'box', 'off');
    axis(ax, 'square')
    hold(ax, 'on')
    
else
    f = figure();
    set(f, 'Color', 'w')
    ax = cell(x*y, 1);
    for c_ax = 1:x*y
        sub_ax = subplot(x,y,c_ax);
        xlabel('', 'Interpreter', 'tex', 'fontsize', 40)
        ylabel('', 'Interpreter', 'tex', 'fontsize', 40)
        % legend([], 'Interpreter', 'latex', 'fontsize', 20)
        title('', 'Interpreter', 'tex', 'fontsize', 40)
        grid minor
        set(sub_ax, 'Box', 'on', 'fontsize', 20, 'LineWidth', 1)
        set(sub_ax,'tickdir', 'out', 'ticklength', [0.03 0.05], 'linewidth', 1.5, 'fontsize', 16, 'fontname', 'arial', 'box', 'off');
%         axis(sub_ax, 'square')
        hold(sub_ax, 'on')
        
        ax{c_ax} = sub_ax;
    end
end


end