function plot_snapshot(x_axis,len_x,p,v,f,db_plot)
%PLOT_SNAPSHOT Summary of this function goes here
%   Detailed explanation goes here
    % Plot
    figure(f);

    if db_plot == false
        % Plot p
        subplot(2,1,1);
        plot(x_axis, p);
        title('Pressure');
        xlim([0,len_x]);
        ylim([-1,1]*2e-1);
        xlabel("x");
        ylabel("p");
    
        % Plot v
        subplot(2,1,2);
        plot(x_axis, v);
        title('Velocity');
        xlim([0,len_x]);
        ylim([-c0,c0]*5e-1);
        xlabel("x");
        ylabel("v");
    else
        % Plot p
        subplot(211);
        plot(x_axis, db(p));
        hold on;
        line([5 5], [-150 0], 'Color', 'red', 'LineStyle', '--');
        hold off;
        title('Pressure (dB)');
        grid on;
        xlim([0,len_x]);
        ylim([-150 0]);
        yticks(-150:10:0);
        xlabel("x");
        ylabel("p (dB)");
    
        % Plot v
        subplot(212);
        plot(x_axis, db(v));
        hold on;
        line([5 5], [-150 0], 'Color', 'red', 'LineStyle', '--');
        hold off;
        title('Velocity (dB)');
        grid on;
        xlim([0,len_x]);
        ylim([-150 0]);
        yticks(-150:10:0);
        xlabel("x");
        ylabel("v (dB)");
    end
end

