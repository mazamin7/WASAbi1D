clear all, close all, clc;

msg = "Choose the merge approach";
opts = ["Pre 7 points" "Post 7 points" "Borrel 3 points" "Borrel 7 points"];
choice = menu(msg, opts);

msg2 = "Choose the update for left";
opts2 = ["FDTD" "Fourier" "FEM"];
choice2 = menu(msg2, opts2);

msg3 = "Choose the update for right";
choice3 = menu(msg3, opts2);

msg4 = "Is left damped?";
opts4 = ["Yes" "No"];
choice4 = menu(msg4, opts4);

msg5 = "Is right damped?";
choice5 = menu(msg5, opts4);

if choice2 ~= 2
    msg6 = "Choose boundary condition for left";
    opts6 = ["N" "D"];
    choice6 = menu(msg6, opts6);

    boundCondLeft = opts6(choice6);
end

if choice3 ~= 2
    msg7 = "Choose boundary condition for right";
    choice7 = menu(msg7, opts6);

    boundCondRight = opts6(choice7);
end

dh = 1/128;
dt = 0.003;
c = 1;

assert(dt <= dh / sqrt(3) / c)

alpha_abs = 10;

len = 1;
dur = 5000;

dur_samples = floor(dur / dt);
N = len/dh;
x_axis = 1:N;

p_prev = zeros(N,1);
p_curr = p_prev * 0;
p_next = p_prev * 0;
force = p_prev * 0;

pulse_width = 1/16;
pulse_pos = 3/4;

pulse_width_x = pulse_width * N;
pulse_pos_x = pulse_pos * N;

pulse_axis = 1:pulse_width_x;
pulse = 1/2 - 1/2 * cos(2*pi*pulse_axis/pulse_width_x);

p_curr(pulse_pos_x-pulse_width_x/2+1:pulse_pos_x+pulse_width_x/2) = pulse;
p_prev(pulse_pos_x-pulse_width_x/2+1:pulse_pos_x+pulse_width_x/2) = pulse;

alpha = 1/90;
beta = -3/20;
gamma = 3/2;
delta = -49/18;

C = sparse(N,N);

C(N/2-2,N/2:N/2+1) = [-alpha, alpha];
C(N/2-1,N/2-1:N/2+2) = [-alpha, -beta, beta, alpha];
C(N/2,N/2-2:N/2+3) = [-alpha, -beta, -gamma, gamma, beta, alpha];
C(N/2+1,N/2-2:N/2+3) = -C(N/2,N/2-2:N/2+3);
C(N/2+2,N/2-1:N/2+2) = -C(N/2-1,N/2-1:N/2+2);
C(N/2+3,N/2:N/2+1) = -C(N/2-2,N/2:N/2+1);

for n = 1:dur_samples

    % Pre-merge
    if choice == 1
        force = (c / dh)^2 * C * p_curr;
    end
    
    % Update left
    if choice2 == 1
        p_next(1:N/2) = update_FDTD(p_curr(1:N/2), p_prev(1:N/2), c, dt, dh, force(1:N/2), choice4 == 1, alpha_abs, choice > 2, boundCondLeft, "N");
    elseif choice2 == 2
        p_next(1:N/2) = update_Fourier(p_curr(1:N/2), p_prev(1:N/2), c, dt, dh, force(1:N/2), choice4 == 1, alpha_abs);
    else
        p_next(1:N/2) = update_FEM(p_curr(1:N/2), p_prev(1:N/2), c, dt, dh, force(1:N/2), choice4 == 1, alpha_abs, boundCondLeft, "N");
    end
    
    % Update right
    if choice3 == 1
        p_next(N/2+1:N) = update_FDTD(p_curr(N/2+1:N), p_prev(N/2+1:N), c, dt, dh, force(N/2+1:N), choice5 == 1, alpha_abs, choice > 2, "N", boundCondRight);
    elseif choice3 == 2
        p_next(N/2+1:N) = update_Fourier(p_curr(N/2+1:N), p_prev(N/2+1:N), c, dt, dh, force(N/2+1:N), choice5 == 1, alpha_abs);
    else
        p_next(N/2+1:N) = update_FEM(p_curr(N/2+1:N), p_prev(N/2+1:N), c, dt, dh, force(N/2+1:N), choice4 == 1, alpha_abs, "N", boundCondRight);
    end
    
    % Post-merge
    if choice == 2
        p_next = p_next + (c * dt / dh)^2 * C * p_curr;
    elseif choice == 3
        % Compute residual parts
        res1 = c^2 * dt^2 / dh^2 * p_curr(N/2-1);
        res2 = c^2 * dt^2 / dh^2 * p_curr(N/2+2);

        % Remove the residual part
        p_next(N/2) = p_next(N/2) - res1;
        p_next(N/2+1) = p_next(N/2+1) - res2;

        % Transfer the removed residual part to the other domain
        p_next(N/2) = p_next(N/2) + res2;
        p_next(N/2+1) = p_next(N/2+1) + res1;
    elseif choice == 4
        % Compute residual parts
        res1 = c^2 * dt^2 / dh^2 * [alpha, beta, gamma] * p_curr(N/2-3:N/2-1);
        res2 = c^2 * dt^2 / dh^2 * [gamma, beta, alpha] * p_curr(N/2+2:N/2+4);

        % Remove the residual part
        p_next(N/2) = p_next(N/2) - res1;
        p_next(N/2+1) = p_next(N/2+1) - res2;

        % Transfer the removed residual part to the other domain
        p_next(N/2) = p_next(N/2) + res2;
        p_next(N/2+1) = p_next(N/2+1) + res1;
    end
    
    % Update
    p_prev = p_curr;
    p_curr = p_next;
    
    % Plot
    f = figure(1);
    f.Position = [100, 100, 1500, 400];
    plot(x_axis*dh, p_next);
    ylim([-1,1]);

    sgtitle(['instant [s]: ' num2str((n+1)*dt, '%4.3f') ' / ' num2str(dur, '%4.3f') ' ( ' num2str((n+1)/dur_samples*100, '%4.1f') '% )']);

end