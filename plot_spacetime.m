function [fig_p, fig_v] = plot_spacetime(t_axis, x_axis, p, v, title_str)
%PLOT_SPACETIME Plots the spacetime solutions for pressure and velocity
%   [fig_p, fig_v] = PLOT_SPACETIME(t_axis, x_axis, p, v, title_str) takes
%   in the time axis t_axis, the space axis x_axis, the pressure solution p
%   and the velocity solution v, as well as an optional title string
%   title_str. It then plots the spacetime solutions for both pressure and
%   velocity, with time on the x-axis, space on the y-axis, and pressure or
%   velocity on the z-axis.
%
%   fig_p and fig_v are the figure handles for the pressure and velocity
%   plots, respectively.
%
%   If title_str is provided, it is concatenated before the words
%   'Pressure Solution' or 'Velocity Solution' in the title of the
%   respective plots.

    % Plot p
    fig_p = figure();
    surf(t_axis, x_axis, p,'EdgeColor','none','FaceColor','interp');
    xlabel('Time [s]');
    ylabel('Space [m]');
    zlabel('Pressure');
    if nargin == 5
        title(['Pressure Solution: ', title_str]);
    else
        title('Pressure Solution');
    end
    view(0, 90);  % set view to show from the top
    colorbar;
    
    % Plot v
    fig_v = figure();
    surf(t_axis, x_axis, v,'EdgeColor','none','FaceColor','interp');
    xlabel('Time [s]');
    ylabel('Space [m]');
    zlabel('Velocity');
    if nargin == 5
        title(['Velocity Solution: ', title_str]);
    else
        title('Velocity Solution');
    end
    view(0, 90);  % set view to show from the top
    colorbar;
end
