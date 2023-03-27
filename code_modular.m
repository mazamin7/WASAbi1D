clear all, close all, clc;

% User menu
% msg = "Choose the merge approach";
% opts = ["Pre 7 points" "Post 7 points"];
% choice = menu(msg, opts);
choice = 1;

% msg2 = "Choose the update for left";
% opts2 = ["FDTD" "Fourier" "FEM" "PML"];
% choice2 = menu(msg2, opts2);
choice2 = 2;

% msg3 = "Choose the update for right";
% choice3 = menu(msg3, opts2);
choice3 = 2;

% if choice2 < 4
%     msg4 = "Is left damped?";
%     opts4 = ["Yes" "No"];
%     choice4 = menu(msg4, opts4);
% else
%     choice4 = 1;
% end
choice4 = 2;

% if choice3 < 4
%     msg5 = "Is right damped?";
%     choice5 = menu(msg5, opts4);
% else
%     choice5 = 1;
% end
choice5 = 2;

opts6 = ["N" "D"];

% if choice2 ~= 4 && choice2 ~= 2
%     msg6 = "Choose boundary condition for left";
%     choice6 = menu(msg6, opts6);
% 
%     boundCondLeft = opts6(choice6);
% else
    boundCondLeft = "N";
% end
choice6 = 1;

% if choice3 ~= 4 && choice3 ~= 2
%     msg7 = "Choose boundary condition for right";
%     choice7 = menu(msg7, opts6);
% 
%     boundCondRight = opts6(choice7);
% else
    boundCondRight = "N";
% end
choice7 = 1;

% Impose transmittance of the middle boundary
T = 1;

% Decide whether FDTD treats boundaries explicitly or not
explicitBoundariesFDTD = false; % WORKS BETTER IF SET TO FALSE
% false -> boundary at N_x/2
% true -> boundary at N_x

% Shift values for Borrel-merge (treat explicitly boundary)
shiftLeft = choice2 == 3 || (choice2 == 1 && explicitBoundariesFDTD == true);
shiftRight = choice3 == 3 || (choice3 == 1 && explicitBoundariesFDTD == true);

% Simulation parameters
c0 = 343;
len_x = 20; % Domain length
T_sec = 0.2; % Simulation duration
alpha_abs = 10; % Absorption coefficient
exact_damping = false; % exact damping only works without interfaces

% msg8 = "Choose the refinement";
% opts8 = ["Finest" "Finer" "Fine" "Coarse" "Custom"];
% choice8 = menu(msg3, opts3);
choice8 = 3;

% SIMULATION PARAMETERS
if choice8 == 1
    dh = 0.05;
elseif choice8 == 2
    dh = 0.1;
elseif choice8 == 3
    dh = 0.2;
elseif choice8 == 4
    dh = 0.5;
end

dh = 1/2^3;
dt = dh / 2 / c0;

disp(['Simulating with dh = ' num2str(dh) ', dt = ' num2str(dt)]);

% SOURCE
freq_source = 300;
source_pos_ratio_x = 5/10;
mach_x = 0;

% Defining time and space axis
dur_samples = floor(T_sec / dt);

% Checking parameters validity
assert(dt <= dh / sqrt(3) / c0)

% Initializing solution data
N_x = floor(len_x/dh);
N_x = 2 * floor(N_x/2);

p_prev = zeros(N_x,1);
p_curr = p_prev * 0;
p_next = p_prev * 0;

q_next_dct = zeros(N_x,1);

% Imposing initial conditions
pulse_width = 1/2^4;
pulse_pos = 1/4;

pulse_width_x = floor(pulse_width * N_x);
pulse_pos_x = floor(pulse_pos * N_x);

pulse_axis = 1:pulse_width_x;
pulse = 1/2 - 1/2 * cos(2*pi*pulse_axis/pulse_width_x);

% p_curr(pulse_pos_x-pulse_width_x/2+1:pulse_pos_x+pulse_width_x/2) = pulse;
% p_prev(pulse_pos_x-pulse_width_x/2+1:pulse_pos_x+pulse_width_x/2) = pulse;

% Defining force spatial envelope
x_axis = linspace(0,len_x,N_x);

mu = len_x * source_pos_ratio_x;           % mean of Gaussian
sigma = len_x/40;       % standard deviation of Gaussian

gauss = @(x) 1/(sigma * sqrt(2 * pi)) * exp(-(x-mu).^2/(2*sigma^2)); % Gaussian function

force_envelope = gauss(x_axis)';
force_envelope(force_envelope < 1e-3) = 0;

% Plot force envelope
figure(1);
plot(x_axis, force_envelope)
xlabel('Position')
ylabel('Intensity')
title('Gaussian Force Envelope')

% Building pre/post-merge matrix
alpha = 1/90;
beta = -3/20;
gamma = 3/2;
delta = -49/18;

C = sparse(N_x,N_x);

if shiftLeft == false
    C(N_x/2-2,N_x/2:N_x/2+1) = [-alpha, alpha];
    C(N_x/2-1,N_x/2-1:N_x/2+2) = [-alpha, -beta, beta, alpha];
    C(N_x/2,N_x/2-2:N_x/2+3) = [-alpha, -beta, -gamma, gamma, beta, alpha];
else
    C(N_x/2-2,N_x/2-1:N_x/2+2) = [-alpha, 0, 0, alpha];
    C(N_x/2-1,N_x/2-2:N_x/2+3) = [-alpha, -beta, 0, 0, beta, alpha];
    C(N_x/2,N_x/2-3:N_x/2+4) = [-alpha, -beta, -gamma, 0, 0, gamma, beta, alpha];
end

if shiftRight == false
    C(N_x/2+1,N_x/2-2:N_x/2+3) = -[-alpha, -beta, -gamma, gamma, beta, alpha];
    C(N_x/2+2,N_x/2-1:N_x/2+2) = -[-alpha, -beta, beta, alpha];
    C(N_x/2+3,N_x/2:N_x/2+1) = -[-alpha, alpha];
else
    C(N_x/2+1,N_x/2-3:N_x/2+4) = -[-alpha, -beta, -gamma, 0, 0, gamma, beta, alpha];
    C(N_x/2+2,N_x/2-2:N_x/2+3) = -[-alpha, -beta, 0, 0, beta, alpha];
    C(N_x/2+3,N_x/2-1:N_x/2+2) = -[-alpha, 0, 0, alpha];
end

% Initializing update methods
if choice2 == 1 || choice2 == 4
    FDTD_data_left = init_FDTD(len_x/2, c0, dt, dh, choice4 == 1, alpha_abs, explicitBoundariesFDTD == true, boundCondLeft, "N", choice2 == 4);
elseif choice2 == 2
    Fourier_data_left = init_Fourier(len_x/2, c0, dt, dh, choice4 == 1, alpha_abs);
else
    FEM_data_left = init_FEM(len_x/2, c0, dt, dh, choice4 == 1, alpha_abs, boundCondLeft, "N");
end

if choice3 == 1 || choice3 == 4
    FDTD_data_right = init_FDTD(len_x/2, c0, dt, dh, choice5 == 1, alpha_abs, explicitBoundariesFDTD == true, "N", boundCondRight, choice3 == 4);
elseif choice3 == 2
    Fourier_data_right = init_Fourier(len_x/2, c0, dt, dh, choice5 == 1, alpha_abs);
else
    FEM_data_right = init_FEM(len_x/2, c0, dt, dh, choice4 == 1, alpha_abs, "N", boundCondRight);
end

force = zeros(N_x,1);

% Simulation loop
for n = 1:dur_samples

    force = 8e5 * sin(2*pi*freq_source*n*dt) * force_envelope * (n*dt <= 1/freq_source);

    g1 = 0; % 1/2*0.3*sin(2*pi*4*n*dt) * (n <= 1/4 / dt);
    g2 = 0; % g1;

    % Pre-merge
    if choice == 1
        force = force + T^2 * (c0 / dh)^2 * C * p_curr;
    end
    
    % Update left
    if choice2 == 1 || choice2 == 4
        p_next(1:N_x/2) = update_FDTD(FDTD_data_left, p_curr(1:N_x/2), p_prev(1:N_x/2), force(1:N_x/2), g1, 0);
    elseif choice2 == 2
        [p_next(1:N_x/2),q_next_dct(1:N_x/2)] = update_Fourier(Fourier_data_left, p_curr(1:N_x/2), p_prev(1:N_x/2), force(1:N_x/2), q_next_dct(1:N_x/2), exact_damping);
    else
        p_next(1:N_x/2) = update_FEM(FEM_data_left, p_curr(1:N_x/2), p_prev(1:N_x/2), force(1:N_x/2));
    end
    
    % Update right
    if choice3 == 1 || choice3 == 4
        p_next(N_x/2+1:N_x) = update_FDTD(FDTD_data_right, p_curr(N_x/2+1:N_x), p_prev(N_x/2+1:N_x), force(N_x/2+1:N_x), 0, g2);
    elseif choice3 == 2
        [p_next(N_x/2+1:N_x),q_next_dct(N_x/2+1:N_x)] = update_Fourier(Fourier_data_right, p_curr(N_x/2+1:N_x), p_prev(N_x/2+1:N_x), force(N_x/2+1:N_x), q_next_dct(N_x/2+1:N_x), exact_damping);
    else
        p_next(N_x/2+1:N_x) = update_FEM(FEM_data_right, p_curr(N_x/2+1:N_x), p_prev(N_x/2+1:N_x), force(N_x/2+1:N_x));
    end

    % FOR NOW, ONLY FDTD SUPPORTS NON-HOMOGENEOUS DIRICHLET/NEUMANN
    
    % Post-merge
    if choice == 2
        p_next = p_next + T^2 * (c0 * dt / dh)^2 * C * p_curr;
    end

    % Fixing interfaces if explicit boundaries
    if shiftLeft == true && shiftRight == true
        p_next(N_x/2+1) = 0.5*(p_next(N_x/2) + p_next(N_x/2+1));
        p_next(N_x/2) = p_next(N_x/2+1);
    end

    % Update
    p_prev = p_curr;
    p_curr = p_next;
    
    % Plot
    f = figure(2);
    f.Position = [100, 100, 1500, 400];
    plot(x_axis, full(p_next));
    xlim([0,len_x])
    ylim([-1,1]);

    sgtitle(['instant [s]: ' num2str((n+1)*dt, '%4.3f') ' / ' num2str(T_sec, '%4.3f') ' ( ' num2str((n+1)/dur_samples*100, '%4.1f') '% )']);

    % pause(0.1);

end