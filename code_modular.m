clear all, close all, clc;

% User menu
msg = "Choose the merge approach";
opts = ["Pre 7 points" "Post 7 points"];
choice = menu(msg, opts);

msg2 = "Choose the update for left";
opts2 = ["FDTD 2ord" "FDTD 1ord" "Fourier 2ord" "Fourier 1ord" "PML"];
choice2 = menu(msg2, opts2);

msg3 = "Choose the update for right";
choice3 = menu(msg3, opts2);

% Simulation parameters
c0 = 1;
len_x = 10; % Domain length
T_sec = 2; % Simulation duration
alpha_abs_left = 0.05; % Absorption coefficient
alpha_abs_right = 0.05;
dh = 0.2;
dt = 0.02;
transmittivity = 1; % Transmittance of the middle boundary

assert(dt < dh / 2 / c0);

% Boundary conditions
bc_left = "D";
bc_right = "D";

% Is partition damped?
left_damped = false;
right_damped = false;
% left_damped = true;
% right_damped = true;

if left_damped == false
    alpha_abs_left = 0;
end

if right_damped == false
    alpha_abs_right = 0;
end

if (choice2 == 2) || (choice2 == 4)
    order_left = 1;
else
    order_left = 2;
end

if (choice3 == 2) || (choice3 == 4)
    order_right = 1;
else
    order_right = 2;
end

assert(~(choice2 == 3 && left_damped));
assert(~(choice3 == 3 && right_damped));

assert(order_left == order_right);
% otherwise domain decomposition doesn't work

% Source
freq_source = 1;
source_fun = @(n) sin(2*pi*freq_source*n*dt) * (n*dt <= 1/freq_source);
source_mu_ratio_x = 1/4;
source_sigma_ratio_x = 1/40;
mach_x = 0;

% Checking parameters validity
assert(dt <= dh / sqrt(3) / c0)

% Defining time and space axis
N_t = floor(T_sec / dt);
N_x = floor(len_x / dh);

% Building residue matrix
C = get_residue_matrix(N_x, 6);

% Initialize solution data
p = zeros(N_x,N_t);
v = zeros(N_x,N_t);

% Defining force spatial envelope
x_axis = linspace(0,len_x,N_x);
mu = len_x * source_mu_ratio_x;
sigma = len_x * source_sigma_ratio_x;       % standard deviation of force spatial envelope (Gaussian)
force_envelope = @(x,mu) 1/(sigma * sqrt(2 * pi)) * exp(-(x-mu).^2/(2*sigma^2)); % Gaussian function

% Initializing update methods
if choice2 <= 2 || choice2 == 5
    data_left = init_FDTD(len_x/2, c0, dt, dh, alpha_abs_left, bc_left, "N", choice2 == 5, order_left);
elseif choice2 >= 3
    data_left = init_Fourier(len_x/2, c0, dt, dh, order_left, alpha_abs_left);
end

if choice3 <= 2 || choice3 == 5
    data_right = init_FDTD(len_x/2, c0, dt, dh, alpha_abs_right, "N", bc_right, choice3 == 5, order_right);
elseif choice3 >= 3
    data_right = init_Fourier(len_x/2, c0, dt, dh, order_right, alpha_abs_right);
end

% Simulation loop
for n = 3:N_t

    % Update force position
    mu = mu + mach_x * c0 * dt;

    % Impose force
    force_envelope_temp = force_envelope(x_axis,mu);
    force = source_fun(n) * force_envelope_temp';

%     % Plot force envelope
%     figure(1);
%     plot(x_axis, force_envelope_temp)
%     xlabel('Position')
%     ylabel('Intensity')
%     title('Gaussian Force Envelope')

    % Impose non-homogeneous b.c.
    g1 = 0; % 1/2*0.3*sin(2*pi*4*n*dt) * (n <= 1/4 / dt);
    g2 = 0; % g1;
    % FOR NOW, ONLY FDTD SUPPORTS NON-HOMOGENEOUS DIRICHLET/NEUMANN

    % Residual calculation
    residual = (c0 / dh)^2 * C * p(:,n-1);

    % Pre-merge
    if choice == 1
        force = force + transmittivity^2 * residual;
    end
    
    % Update left
    if choice2 <= 2 || choice2 == 5
        [p(1:N_x/2,n),v(1:N_x/2,n)] = update_FDTD(data_left, p(1:N_x/2,n-1), p(1:N_x/2,n-2), force(1:N_x/2), v(1:N_x/2,n-1), v(1:N_x/2,n-2), g1, 0);
    elseif choice2 >= 3
        [p(1:N_x/2,n),v(1:N_x/2,n)] = update_Fourier(data_left, p(1:N_x/2,n-1), p(1:N_x/2,n-2), force(1:N_x/2), v(1:N_x/2,n-1), v(1:N_x/2,n-2));
    end
    
    % Update right
    if choice3 <= 2 || choice3 == 5
        [p(N_x/2+1:N_x,n),v(N_x/2+1:N_x,n)] = update_FDTD(data_right, p(N_x/2+1:N_x,n-1), p(N_x/2+1:N_x,n-2), force(N_x/2+1:N_x), v(N_x/2+1:N_x,n-1), v(N_x/2+1:N_x,n-2), 0, g2);
    elseif choice3 >= 3
        [p(N_x/2+1:N_x,n),v(N_x/2+1:N_x,n)] = update_Fourier(data_right, p(N_x/2+1:N_x,n-1), p(N_x/2+1:N_x,n-2), force(N_x/2+1:N_x), v(N_x/2+1:N_x,n-1), v(N_x/2+1:N_x,n-2));
    end

    % Post-merge
    if choice == 2
        if order_left == 1
            v = v + transmittivity^2 * 2 * dt * residual;
        elseif order_left == 2
            p = p + transmittivity^2 * dt*dt * residual;
        end
    end
    
    % Plot
    f = figure(2);
    f.Position = [100, 100, 1500, 900];
    sgtitle(['instant [s]: ' num2str(n*dt, '%4.3f') ' / ' ...
        num2str(T_sec, '%4.3f') ' ( ' num2str(n/N_t*100, '%4.1f') '% )']);

    % Plot p
    subplot(2,1,1);
    plot(x_axis, p(:,n));
    title('Pressure');
    xlim([0,len_x]);
    ylim([-1,1]*2e-1);

    % Plot v
    subplot(2,1,2);
    plot(x_axis, v(:,n));
    title('Velocity');
    xlim([0,len_x]);
    ylim([-c0,c0]*5e-1);

    % pause(0.1);

end

% Plot p
figure;
surf(p,'EdgeColor','none','FaceColor','interp');
xlabel('Time');
ylabel('Space');
zlabel('Pressure');
title('Pressure Solution');

% Plot v
figure;
surf(v,'EdgeColor','none','FaceColor','interp');
xlabel('Time');
ylabel('Space');
zlabel('Velocity');
title('Velocity Solution');


