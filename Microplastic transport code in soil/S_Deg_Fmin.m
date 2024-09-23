function calibrated_params = SL_Deg_Fmin()
    clear all; clc;
    start_time = tic;
    
    % Define initial guess for parameters
    x0 = [1, 0.2, 0.0005, 1, 0.43, 1500];

    % Set lower and upper bounds for parameters
    lb = [0, 0, 0, 0, 0.02, 500];  % Lower bounds for lamda, Ka, Kd, Kst, parametric_beta, Smax
    ub = [50, 0.01, 0.1, 2, 2.28, 10000];  % Upper bounds for lamda, Ka, Kd, Kst, parametric_beta, Smax

    % Define optimization options
    options = optimoptions('fmincon','Algorithm','interior-point', 'Display', 'iter', 'MaxIterations', 200, 'FunctionTolerance', 1e-6);
    % options = optimoptions('fmincon','Algorithm','sqp', 'Display', 'iter', 'MaxIterations', 200, 'FunctionTolerance', 1e-6);

    % Run optimization for parameter estimation
    calibrated_params = fmincon(@(params) objective_function(params, x0), x0, [], [], [], [], lb, ub, [], options);

    % Display the calibrated parameters with names
    disp('calibrated_params:');
    disp(calibrated_params);

    % Call the main function with calibrated parameters
    [t, C, x, S] = advection_dispersion_model(calibrated_params(1), calibrated_params(2), calibrated_params(3), calibrated_params(4), calibrated_params(5), calibrated_params(6));

    % Plot results
    plot_results(t, C, x, S);

    % Calculate errors
    [error, error_C, error_S] = objective_function(calibrated_params, x0);
    disp(['Error: ', num2str(error)]);
    disp(['Error_C (MP particles/ml of water): ', num2str(error_C)]);
    disp(['Error_S (MP particles/g of soil): ', num2str(error_S)]);

    % Measure the elapsed time
    elapsed_time = toc(start_time);
    disp(['Elapsed Time: ', num2str(elapsed_time), ' seconds']);
end

function [error, error_C, error_S] = objective_function(params, x0)
    % Call the main function with current parameter values
    [t, C, x, S] = advection_dispersion_model(params(1), params(2), params(3), params(4), params(5), params(6));

     % Provided Experimental Results for C and S
    t_exp = [0; 30; 48; 80; 112; 144; 184; 244];
    C_exp = [0; 13; 22; 28; 17; 11; 6; 4];

    x_exp = [0.16; 1.43; 3.81; 6.35; 8.89]; % Adjust as needed
    S_exp = [50; 13; 8; 4; 3]; % Define your experimental retention profile here

    % Interpolate model results to match experimental positions for S
    C_model = interp1(t, C(:,18), t_exp, 'linear', 'extrap');
    S_model = interp1(x, S(end,:), x_exp, 'linear', 'extrap');

    % Calculate the errors
    error_C = sqrt(mean((C_model - C_exp).^2));
    error_S = sqrt(mean((S_model - S_exp).^2));
    Weightage_C = 0.6;
    Weightage_S = 0.4;
    error = Weightage_C.*error_C + Weightage_S.*error_S;
end

function [t, C, x, S] = advection_dispersion_model(lamda, Ka, Kd, Kst, parametric_beta, S_max)
   
    %-------advection_dispersion_model function code here
    
%--- time and space discritization-------

dt = 0.1;                     % Time step
n = 18;                      % number of nodes 

%--- Soil column dimensions--------------

L = 10;                  % length of column (cm)
d = 5.08;                 % diameter of column (cm)
a = (22./7).*(d.^2)./4;  % area of the column
vol = a.*L;              % volume of column
dx = L./(n-1);             % space step

%--- soil, physcial and flow properties-----------------

d50 = 0.0025;     % mean diameter of porous media particles (cm)
theeta = 0.475;   % soil porosity (cm3/cm3)
rho = 1.39;     % soil medium bulk density (g/cc)

Q = 1;      % (ml/min)
V = Q./(theeta.*a); % pore velocity (cm/min)
Vv = vol.*theeta;        % Volume of voids (cm3)
Ci = 10667./(dx.*a.*theeta);    % concentration of colloids (particle/ml)

%--- matrix creation----------------

t = (0:dt:244);               % number of time steps
PV = Q.*t./Vv;
x = 0:dx:L;

CC = zeros(length(t),n);
C_star = zeros(length(t),n);
S = zeros(length(t),n);
CCC_star = zeros(length(t),n);
SSS_star = zeros(length(t),n);


D = lamda.*V; % Diffusion coefficient (cm2/min)

alpha = V.*dx./(2.*D);
beta = (dx.^2)./(D.*dt); % (1/beta) = Courant number (<1)

ls(1,1) = 0;                    % lower sub diagonal
us(1,1) = 1-alpha; us(n,1) = 0;       % uper sub diagonal
d(1,1) = ((1+alpha).*(1-(2.*alpha)))-(2+beta);      % main diagonal (Cauchy BC)
b(1,1) = Ci; b(n,1) = 0;     % boundary conditions
 
ls(2:n-1,1) =  (1+alpha);
us(2:n-1,1) =  (1-alpha);
d(2:n-1,1) = -(2+beta);
b(2:n-1,1) = 0;  % initial condition
% left hand side boundary condition
ls(n,1) = (2.*alpha) ;
d(n,1) = -((2.*alpha)+beta);


%   creating coefficient matrix
A = zeros(n);
for j = 1:n
    A(j,j) = d(j,1);
end
for j=1:n-1
    A(j,j+1) = us(j,1);
end
for j= 2:n
    A(j,j-1) = ls(j,1);
end

% BTC by implicit scheme for all time sereis 
C = zeros(length(t),n);
C(1,:) = b;
for i = 2:length(t)
         b = -beta.*C(i-1,:).';
       
        C(i,:) = A\b;
     C(i-1,:) = C(i,:);

    % %------Retardation (ODE) equation Solver by Runge-Kutta 4th order (RK4) scheme-----------

    for j = 1:n

    %--- defining the equation for slope-----------
      z(1) = 0;
    z(j+1) = z(j) + dx; % down gradient depth from inlet
    si_st = (d50./(d50+z(j))).^parametric_beta;  % depth dependent power law function for straining
    si_at = (1-(S(j)./S_max));

        dS_dt= @(t,C,S) (1./rho).*(theeta.*Ka.*(si_at).*C - Kd.*rho.*S + (theeta.*Kst.*(si_st).*C));
        dC_dt= @(t,C,S) (rho./theeta).*(1./rho).*((theeta.*Ka.*(si_at).*C - Kd.*rho.*S + (theeta.*Kst.*(si_st).*C)));
    %--- solving for average slope (using four slopes) to find concentration value at nest time step-----  

        KS1 = dS_dt(t(j),C(i-1,j),S(i,j));
        KC1 = dC_dt(t(j),C(i-1,j),S(i,j));
        SS_star = S(i-1,j) + KS1.*(dt/2);
        CC_star = C(i-1,j) - KC1.*(dt/2);

        KS2 = dS_dt(t(j)+(dt./2),CC_star,SS_star);
        KC2 = dC_dt(t(j)+(dt./2),CC_star,SS_star);
        SS_doublestar = S(i-1,j) + KS2.*(dt/2);
        CC_doublestar = C(i-1,j) - KC2.*(dt/2);

        KS3 = dS_dt(t(j)+(dt./2),CC_doublestar,SS_doublestar);
        KC3 = dC_dt(t(j)+(dt./2),CC_doublestar,SS_doublestar);
        SS_triplestar = S(i-1,j) + KS3.*(dt);
        CC_triplestar = C(i-1,j) - KC3.*(dt);

        KS4 = dS_dt(t(j)+(dt),CC_triplestar,SS_triplestar);
        KC4 = dC_dt(t(j)+(dt),CC_triplestar,SS_triplestar);

        KS_avg = (KS1+2.*KS2+2.*KS3+KS4)./6;
        KC_avg = (KC1+2.*KC2+2.*KC3+KC4)./6;
        SSS_star(i,j) = S(i-1,j) + KS_avg.*(dt);
        CCC_star(i,j) = C(i-1,j) - KC_avg.*(dt);
        S(i,j) = SSS_star(i,j);
        C(i,j) = CCC_star(i,j);

    end

end
end

function plot_results(t, C, x, S)
    % Plot Breakthrough Curve (t vs C)
    figure;
    subplot(2, 1, 1);  % Two rows, one column, first plot
    plot(t, C(:, 18), '-', 'DisplayName', 'Simulated');
    hold on;
    % Plot Experimental Results for C
    t_exp = [0; 30; 48; 80; 112; 144; 184; 244];
    C_exp = [0; 13; 22; 28; 17; 11; 6; 4];

    
    plot(t_exp, C_exp, 'o', 'DisplayName', 'Experimental');
    title('Breakthrough Curve');
    xlabel('Time (min)');
    ylabel('Concentration');
    legend('show');
    hold off;
    % Plot Retention Profile (x vs S)
    subplot(2, 1, 2);  % Two rows, one column, second plot
    plot(x, S(end, :), '-', 'DisplayName', 'Simulated');
    hold on;
    x_exp = [0.16; 1.43; 3.81; 6.35; 8.89]; % Adjust as needed
    S_exp = [50; 13; 8; 4; 3]; % Define your experimental retention profile here
    plot(x_exp, S_exp, 'o', 'DisplayName', 'Experimental');
    title('Retention Profile');
    xlabel('Distance (cm)');
    ylabel('Retention');
    legend('show');
    hold off;
end