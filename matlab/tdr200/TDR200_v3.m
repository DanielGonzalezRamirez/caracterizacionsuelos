%% TDR data processing
% This script reads the TDR output files (from gprMax and PCTDR) and
% obtains the complex relative permittivity through time-domain and
% frequency-domain (FFT) analysis of the signals.
% Author: 
%   Daniel Julián González Ramírez - dj.gonzalez1203@uniandes.edu.co
% june 2019
close all
clear, clc
%% References
% [1] Heimovaara, T. J., & Bouten, W. (1990). A Computer-Controlled
% 36-Channel Time Domain Reflectometry System for Monitoring Soil Water 
% Contents. Water Resources Research, 26(10), 2311–2316.
%% Anonymus functions
% Calculation of Havriliak Negami permittivity model:
Havriliak_Negami = @(eps_inf, eps_s, f, t_o, alpha, beta) eps_inf +...
    (eps_s - eps_inf)./ ((1 + (1i * 2 * pi * t_o .*f).^(1 - alpha)).^beta);
% Calculation of complex rho for the coaxial probe:
rho_calc = @(z,eps) (1 - z.*sqrt(eps))./(1 + z.*sqrt(eps));
% Calculation of complex gamma for the coaxial probe:
gamma_calc = @(f,eps,c) (1i*2*pi.*f.*sqrt(eps))./c;
% Calculation of the S11 parameter for the coaxial probe:
S11_calc = @(rho,gamma,l) (rho + exp(-2*l.*gamma))./(1 +...
    rho.*exp(-2*l.*gamma));
%% Constants
% The following are the values for the constant parameters for the coaxial
% TDR probe.
eps_o = 8.854187817e-12;    % Vacuum permittivity
mu_o = 1.256637061e-6;      % Vacuum permeability
c_o = 1/sqrt(eps_o*mu_o);   % Propagation velocity in vacuum [m/s]
Z_c = 50;                   % Impedance of the cable [Ohms]
low_freq = 0.01;            % Low frequency of TDR range [GHz]
high_freq = 1;              % High frequency of TDR range [GHz]
%% File load
% Prompt the user to select the source of the data to be used.
% Version 1 of TDR script will only accept gprMax files.
def_route = 'D:\OneDrive\OneDrive - Pontificia Universidad Javeriana\2019_10\TESIS\TDR200\*.dat';
[filename,pathname] = uigetfile(def_route, 'Select output file');
fullfilename = [pathname filename];
%% Retrival of the desired information from the output file
[header, Vinc, Vref, Vinc_prime, Vref_prime, S11, time, freq, thead]...
    = TDR200_read(fullfilename);
tsamp = time(2) - time(1); % Sampling time [ns]
fsamp = freq(2) - freq(1); % Sampling freq [GHz]
dt = tsamp * 10 ^-9;
%% Read header
Lmech = header.Lmech;   % Length of the probe [m]
a = header.a;           % Rod diameter [m]
b = header.b;           % Distance between outer rods [m]
Z_p = 60 * log(b/a);    % Impedance of the probe in vaccum [Ohms]
z = Z_c/Z_p;            % Impedance ratio
L_epoxy = header.Lhead; % Length of the epoxy head [m]
eps_epoxy = header.eps; % Relative permittivity of the epoxy
%% Load calibration
% The electrical length of the probe might differ from the actual length.
% In order to calibate the length to the electrical distance, an open-air
% signal is analyzed.
%load('calibration_2019-06-19-09-01-00.mat')
%load('calibration_2019-06-25-11-32-05.mat')
load('calibration_2019-07-16-12-07-08.mat') %Comercial
%% Time-Domain analysis
% Time-Domain analysis accounts for the determination of the electrical
% conductivity and the apparent relative permittivity of the medium.

% High-frequency components account for errors in the determination of the
% inflection points for the time-domain signal. Time-windowed smoothing is
% used to elimiate the undesirable frequency components.
cut_off_freq = 12e9; % Cut-off frequency [Hz]
cut_off_samples = floor(1 / (cut_off_freq * dt));
Vinc_smoothed = smooth(Vinc,cut_off_samples)';
Vref_smoothed = smooth(Vref,cut_off_samples)';

Vinc_prime_smooth = diff(Vinc_smoothed);
Vinc_prime_smooth = [Vinc_prime_smooth Vinc_prime_smooth(end)]./tsamp;
Vref_prime_smooth = diff(Vref_smoothed);
Vref_prime_smooth = [Vref_prime_smooth Vref_prime_smooth(end)]./tsamp;

% To find the reflection point in the TDR reflected signal one must
% locate the inflection point of the first reflection (S), through the
% first derivative of the reflected voltage.
S_loc = S_point_location(Vref_prime_smooth);

% After the location of the inflection point one must determine the time
% intervals for the reflection point analysis as in Heimovaara and Bouten
% [1]. This algorithm has a different approach than [1] for the calculation
% of the time intervals, although the current algotithm relies also on
% empirical parameters it has a systematic approach to the calculation
% based on percentages of the inflection point time.

dt_min = 0.1*time(S_loc);
time_min = time(S_loc) - dt_min;
min_loc = find(time >= time_min, 1);
time_min = time(min_loc:S_loc);

dt_base = 60*tsamp;
dt_infl = 20*tsamp;

% With the time intervals set now the linear approximation for the base and
% inflection intervals are needed. In this algorithm they are done using
% the first and last points of the interval to determine the slope and
% intercept of the lines. This approach showed better results than linear
% reggression and weigther linear regression.

line_base = line_calc(Vref,time,dt_base,time_min,min_loc,1);
line_infl = line_calc(Vref,time,dt_infl,time_min,S_loc,2);

% Apparent permittivity calculation. Using the interception point of both
% inflection line and base line as the reflection point (sto point), and 
% the first peak as the start point
start_time_loc = find(Vref ~= 0, 1);
line_difference = abs(line_infl.^2 - line_base.^2);
% stop_time_loc = find(line_infl >= line_base,1); % Location in time_min
[~,stop_time_loc] = min(line_difference);
stop_time_loc = find(time == time_min(stop_time_loc));

flight_time = (time(stop_time_loc) - time(start_time_loc))*1e-9; % [s]

apparent_permittivity = ((flight_time * c_o)/(L))^2;
display(apparent_permittivity)

% Electrical conductivity calculation. The electrical conductivity is 
% calculated using the reflected voltage at infinity-time and the geometry 
% of the probe.
rho_infty = Vref(end)/Vinc(end);
sigma_dc = ((eps_o*c_o)/(z*L))*((1-rho_infty)/(1+rho_infty));
display(sigma_dc)

% Response of an ideal medium. The area between the reflected waveforms of
% an ideal medium with eps = 1 -j*sigma/(eps_o*2*pi*f) and the actual
% reflected waveform gives the value of the static relattive permittivity
% eps_s.
%delay = 2*L/c_o;
delay = 0;
steps_delay = delay/dt;
%time_index = ceil(start_time_loc + steps_delay);
time_index = start_time_loc;

length_high = length(Vref) - time_index;
%first_nonzero = find(Vinc > 0,1);
first_nonzero = 1095;

Vref_ideal = [zeros(1,time_index) Vinc(first_nonzero:first_nonzero+length_high-1)];
%Vref_ideal = Vinc;
Vref_ideal = Vref_ideal*rho_infty;

wf_difference = Vref_ideal(1:12280) - Vref(1:12280);
area = trapz(wf_difference,time(1:12280).*10^-9);
eps_static = abs((area*c_o)/(g*(Vref(end))*z*L)) + 1;
%eps_static = 30.8998;

display(eps_static)
%% Graphics from time-domain analysis
% Incident and reflected voltage waveforms
figure('Name','Incident and reflected voltage')
subplot(2,1,1); plot(time, Vinc); grid minor; 
xlabel('time [ns]'); ylabel('Voltage [V]'); title('Incident voltage');
xlim([time(1) time(end)]); 
ylim([min([min(Vinc) min(Vref)]) max([max(Vinc) max(Vref)])]);
subplot(2,1,2); plot(time, Vref); grid minor; 
xlabel('time [ns]'); ylabel('Voltage [V]'); title('Reflected voltage');
xlim([time(1) time(end)]); 
ylim([min([min(Vinc) min(Vref)]) max([max(Vinc) max(Vref)])]);

% End point location and tangent lines
figure('Name','End reflection point')
hplot = plot(time, Vref); hold on; grid minor;
plot(time_min, line_base); hold on
plot(time_min, line_infl); hold on;
xlabel('time [ns]'); ylabel('Voltage [V]'); 
title('Reflected voltage end reflection');
xlim([time_min(1)-1 time_min(end)+1]);

dcm_obj = datacursormode();
dtip = createDatatip(dcm_obj, hplot);
dtip.Position = [time(stop_time_loc),Vref(stop_time_loc)];
text(time(stop_time_loc),Vref(stop_time_loc),'  End point')

clear hplot dtip dcm_obj

% Location of start and end point on the reflected voltage
figure('Name','Start and end point')
hplot = plot(time, Vref); hold on; grid minor;
xlabel('time [ns]'); ylabel('Voltage [V]'); 
title('Reflected voltage start and end points');
xlim([time(1) 1.5*time_min(end)]);

dcm_obj = datacursormode();
dtip1 = createDatatip(dcm_obj, hplot);
dtip1.Position = [time(start_time_loc),Vref(start_time_loc)];
text(time(stop_time_loc),Vref(stop_time_loc),'  End point')
dtip2 = createDatatip(dcm_obj, hplot);
dtip2.Position = [time(stop_time_loc),Vref(stop_time_loc)];
text(time(start_time_loc),Vref(start_time_loc),'  Start point')
text(time(stop_time_loc),Vref(stop_time_loc),'  End point')

dim = [.2 .5 .3 .3];
str = {['2-way flight time = ', num2str(flight_time), ' [s]'],...
    ['K_a = ',num2str(apparent_permittivity)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

clear hplot dtip1 dtip2 dcm_obj dim str

% Comparisson between incident amplitude and stable reflected amplitude
figure('Name','Incident and reflected voltage comparisson')
hplot1 = plot(time, Vinc); hold on; grid minor; 
hplot2 = plot(time, Vref);
xlabel('time [ns]'); ylabel('Voltage [V]'); 
title('Incident and reflected voltage');
legend('Incident', 'Reflected');
xlim([time(1) time(end)]); 

dcm_obj = datacursormode();
dtip1 = createDatatip(dcm_obj, hplot1);
dtip1.Position = [time(end),Vinc(end)];
dtip2 = createDatatip(dcm_obj, hplot2);
dtip2.Position = [time(end),Vref(end)];

dim = [.3 .5 .1 .1];
str = {['\rho_\infty = ', num2str(rho_infty)],...
    ['\sigma_{DC} = ',num2str(sigma_dc),' [S/m]']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

clear hplot1 hplot2 dtip1 dtip2 dcm_obj dim str

% Ideal reflected waveform and actual reflected waveform
figure('Name','Ideal and actual reflected waveforms')
plot(time, Vref_ideal); hold on; grid minor; 
plot(time, Vref);
xlabel('time [ns]'); ylabel('Voltage [V]'); 
title('Ideal and actual reflected waveforms');
legend('Ideal', 'Actual');
xlim([time(1) time(end)]);

dim = [.3 .5 .1 .1];
str = {['Area between waveforms = ', num2str(area), ' [Vs]'],...
    ['\epsilon_s = ',num2str(eps_static)],};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

clear dim str
%% Frequency-Domain analysis - Calibration
% Frequecy-Domain analysis accounts for the determination of the frequency 
% dependent complex relative permittivity of the medium.

S11_smoothed = fftshift(fft(Vref_prime_smooth))./fftshift(fft(Vinc_prime_smooth));

% To obtain the frequency dependent complex relative permittivity, the S11
% parameter of the measurement needs to be calibrated with the error
% functions for the system.
%S11_calibrated = (S11_smoothed - E_dir)./(E_smh.*(S11_smoothed - E_dir) + E_ftr);
%S11_calibrated = (S11_smoothed - E_dir)./(E_ftr);
S11_calibrated = smooth(S11_smoothed,27)';
figure('Name','Measured S11')
subplot(2,1,1);
plot(freq,real(S11_calibrated));
xlim([-high_freq high_freq]);
grid minor;
xlabel('Frequency [GHz]'); ylabel('re{S11}'); 
title('Measured S11');
subplot(2,1,2);
plot(freq,imag(S11_calibrated));
xlim([-high_freq high_freq]);
grid minor;
xlabel('Frequency [GHz]'); ylabel('im{S11}'); 

%% Recovery of the parameters
%  Then the S11 parameter will be used to solve the permittivity from the 
% coaxial line equation using particle swarm optimization (PSO) method.
% Each particle will move in three dimensions that correspond to the three
% variables not yet determined from a Cole-Cole dispersion model for
% permittivity: Relaxation frequency, and Eps_infty.

Nd = 2; % Dimensiones de las partículas
Qparticles = 100; % Cantidad de partículas a utilizar
rfreq_l = 10e6;
rfreq_h = 20e9;
epsinfty_l = 1;
epsinfty_h = 5;

low_index = find(freq >= -0.5,1);
high_index = find(freq >= 0.5,1);
freq_cut = freq(low_index:high_index);
S11_calibrated_cut = S11_calibrated(low_index:high_index);
S11_smoothed_cut = S11_smoothed(low_index:high_index);

QPSOruns = 10;
globalBests = zeros(QPSOruns,Nd);
globalFitness = zeros(QPSOruns,1);
%%
for global_runs = 1:QPSOruns
    particles = zeros(Qparticles,Nd); % Arreglo de particulas
    velocities = zeros(Qparticles,Nd); % Arreglo de velocidades de partículas

    % Particles are located randomly within constrains of their dimensions.
    for i = 1:Qparticles
        [particles(i,1), particles(i,2)] = ...
            random_particle_locator_v2(rfreq_l, rfreq_h, epsinfty_l, epsinfty_h);
    end

    % Stores the adjustment of particle produces on the current iteration:
    current_fitness = zeros(Qparticles, 1);
    % Stores the best adjustment each particle produces (particle best fitness):
    pBestFitness = zeros(Qparticles, 1);
    % Stores the location of each particle (particle best) where it best
    % adjusted the model:
    pBest = zeros(Qparticles, Nd);
    % Stores the best global adjustment the particles produced:
    gBestFitness = 0;
    % Stores the location of the particle that produced the best adjustment:
    gBest = zeros(1, Nd);

    % Velocity update constants:
    w = 0.5; c1 = 1.4; c2 = 1.4;

    % Exit condition is the number of iterations:
    Qiterations = 100;
    %% Iterative process
    save_1 = particles;
    % figure
    for c_iter = 1:Qiterations
    %     scatter(particles(:,1), particles(:,2))
    %     title(num2str(c_iter))
    %     xlim([rfreq_l rfreq_h]); ylim([0 10]);
    %     drawnow
        if c_iter == 33
            save_2 = particles;
        elseif c_iter == 66
            save_3 = particles;
        elseif c_iter == 100
            save_4 = particles;
        end
        for c_part = 1:Qparticles
            if particles(c_part,1) > rfreq_h || particles(c_part,1) < rfreq_l
                [particles(c_part,1), ~] = ...
            random_particle_locator_v2(rfreq_l, rfreq_h, epsinfty_l, epsinfty_h);
            end
            if particles(c_part,2) < epsinfty_l || particles(c_part,2) > eps_static
                [~, particles(c_part,2)] = ...
            random_particle_locator_v2(rfreq_l, rfreq_h, epsinfty_l, epsinfty_h);
            end

            frelax = particles(c_part,1);
            eps_inf = particles(c_part,2);
            alpha = 0;
            beta = 1;

            t_o = 1/(2*pi*frelax);
            eps_aux = Havriliak_Negami(eps_inf, eps_static, freq_cut.*10^9, t_o, alpha, beta);
            eps_recov = eps_aux - (1i*sigma_dc)./(2*pi*eps_o.*(freq_cut.*10^9));

            rho_recov = rho_calc(z,eps_recov);
            gamma_recov = gamma_calc(freq_cut.*10^9,eps_recov,c_o);
            S11_recov = S11_calc(rho_recov,gamma_recov,L);

            %real_fit = sum(real(S11_recov).^2 - real(S11_calibrated_cut).^2);
            %imag_fit = sum(imag(S11_recov).^2 - imag(S11_calibrated_cut).^2);
            %abs_fit = sum(abs(S11_recov).^2 - abs(S11_calibrated_cut).^2);

            peaks_recov = numel(findpeaks(real(S11_recov)));
            peaks_measu = numel(findpeaks(real(S11_calibrated_cut)));

            current_fitness(c_part) = abs(peaks_recov - peaks_measu);
            %current_fitness(c_part) = real_fit;

            if c_iter == 1
                pBestFitness(c_part) = current_fitness(c_part);
                pBest(c_part,:) = particles(c_part,:);
            else
                if current_fitness(c_part) < pBestFitness(c_part)
                    pBestFitness(c_part) = current_fitness(c_part);
                    pBest(c_part,:) = particles(c_part,:);
                end
            end
        end
        % The best fitness and particle location obtained among particles is
        % stored.
        [current_best, best_pos] = min(abs(pBestFitness));
        % The current best is compared against the former best. If the current
        % betters the former, it will be stored into gBestFitness and gBest.
        if c_iter == 1
            gBestFitness = current_best;
            gBest(:) = pBest(best_pos,:);
        else
            if current_best < gBestFitness
                gBestFitness = current_best;
                gBest(:) = pBest(best_pos,:);
            end
        end

        for c_part = 1:Qparticles
            velocities(c_part,:) = w.*velocities(c_part,:) +...
                c1 * rand .* (pBest(c_part,1:Nd) - particles(c_part,:)) +...
                c2 * rand .* (gBest(1:Nd) - particles(c_part,:));
            particles(c_part,:) = particles(c_part,:) + velocities(c_part,:);
        end
    end
    globalFitness(global_runs) = gBestFitness;
    globalBests(global_runs,:) = gBest;
end
[~,bestComp] = min(globalFitness);
% gBest(1) = mean(globalBests(:,1));
% gBest(2) = mean(globalBests(:,2));
gBest(1) = globalBests(bestComp,1);
gBest(2) = globalBests(bestComp,2);
%% Plot solution
% gBest = globalBests(5,:);
sol.frelax = gBest(1);
sol.eps_inf = gBest(2);
sol.eps_s = eps_static;
sol.sigma_dc = sigma_dc;
sol.t_o = 1/(2*pi*sol.frelax);
sol.alpha = 0;
sol.beta = 1;

eps_aux = Havriliak_Negami(sol.eps_inf, sol.eps_s, freq.*10^9, sol.t_o, sol.alpha, sol.beta);
complex_eps = eps_aux - (1i*sol.sigma_dc)./(2*pi*eps_o.*(freq.*10^9));

rho_sol = rho_calc(z,complex_eps);
gamma_sol = gamma_calc(freq.*10^9,complex_eps,c_o);
S11_sol = S11_calc(rho_sol,gamma_sol,L);


figure('Name','Recovered S11')
subplot(2,1,1);
plot(freq,real(S11_calibrated)); hold on;
plot(freq,real(S11_sol)); hold on
xlim([-high_freq high_freq]); legend('Measured', 'Recovered');
grid minor;
xlabel('Frequency [GHz]'); ylabel('re{S11}'); 
title('Recovered S11');
subplot(2,1,2);
plot(freq,imag(S11_calibrated)); hold on;
plot(freq,imag(S11_sol)); hold on
xlim([-high_freq high_freq]); legend('Measured', 'Recovered');
grid minor;
xlabel('Frequency [GHz]'); ylabel('im{S11}');

figure('Name', 'Recovered complex permittivity')
semilogx(freq, real(complex_eps), 'LineWidth', 2); hold on
semilogx(freq, abs(imag(complex_eps)), 'LineWidth', 2); grid minor
xlim([low_freq 12]);
xlabel('Frequency [GHz]'); ylabel('Relative permittivity');
legend('real', 'imag');
title('Complex permittivity measured with TDR');