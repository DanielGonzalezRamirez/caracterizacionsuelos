function [probe_info, Vinc, Vref, Vinc_prime, Vref_prime, S11, time, freq,...
    head_delay] = TDR200_read(fullfilename)
%% An�lisis de se�ales obtenidas con el TDR200
% El TDR200 de Campbell Scientific se conecta con el computador para hacer
% lecturas de la se�al reflejada por una sonda coaxial. A partir de estas
% lectura se puede determinar la permitividad compleja y la conductividad 
% del diel�ctrico donde est� insertada la sonda.
% En este script se estudia la onda reflejada para obtener la permitividad
% compleja y la conductividad de los suelos a partir de las mediciones.
%
% Daniel Juli�n Gonz�lez Ram�rez
% dj.gonzalez1203@uniandes.edu.co
% 2019-06-17
%% Constantes importantes para los c�lculos
eps_o = 8.854187817e-12;    % Permitividad del vac�o
mu_o = 1.256637061e-6;      % Permeabilidad del vac�o
c_o = 1/sqrt(eps_o*mu_o);   % Velocidad de propagaci�n en vac�o
%% Extracci�n del vector de tiempo y del coeficiente de reflexi�n
Qlines = 10;
[time, rho, header] = read_dat(fullfilename,Qlines);
tsamp = time(2) - time(1);
dt = tsamp * 10 ^-9;
%% Caracter�sticas de la sonda de medici�n
Qsamples = Qlines*header{3};
eps_head = header{8}; % Permitividad relativa de la cabeza de la sonda [m]
% La separaci�n entre las varillas y el di�metro de las varillas se conoce
% de las dimensiones f�sicas. Estas dimensiones son diferentes para la
% sonda Campbell Scientific CS610 y la fabricada.
if eps_head == 2.3
    probe_diameter = 0.00476; % Di�metro de las varillas [m]
    probe_spacing = 0.046; % Separaci�n entre las varillas externas [m]
    head_length = header{7}; % Longitud de la cabeza de la sonda [m]
    probe_length = header{6} - head_length; % Longitud enterrable de las varillas [m]
    cable_length = 2; % Longitud del cable coaxial que alimenta la sonda [m]
    cable_vp = 0.66; % Porcentaje de velocidad de propagaci�n del cable
    eps_head = 4.6;
elseif eps_head == 1.74
    probe_diameter = 0.0048; % Di�metro de las varillas [m]
    probe_spacing = 0.045; % Separaci�n entre las varillas externas [m]
    head_length = header{7}; % Longitud de la cabeza de la sonda [m]
    probe_length = header{6}; % Longitud enterrable de las varillas [m]
    cable_length = 4.8768 + .5; % Longitud del cable coaxial que alimenta la sonda [m]
    cable_vp = 0.87; % Porcentaje de velocidad de propagaci�n del cable
else
    probe_diameter = 0.0042; % Di�metro de las varillas [m]
    probe_spacing = 0.046; % Separaci�n entre las varillas externas [m]
end

probe_info.Lmech = probe_length;
probe_info.Lhead = head_length;
probe_info.a = probe_diameter;
probe_info.b = probe_spacing;
probe_info.eps = eps_head;

cable_delay = cable_length / (c_o * cable_vp) * 10^9; % Delay del cable [ns]
%% Creaci�n del perfil de subida de la se�al
% El pulso de TDR tiene un tiempo de subida de 85 ps de 0 a 250 mV, aunque
% el fabricante del equipo no especifica la forma del perfil de subida se
% puede asumir de forma exponencial. Asumir que el pulso tiene una forma
% exponencial de subida tiene un gran beneficio para el algoritmo de FDTD,
% ya que la se�al entra al dominio de manera suave evitando errores de
% dispersi�n num�rica.
delay = cable_delay; % retardo para el inicio del pulso [ns]
amplitude = 0.24; % Amplitud del pulso [V]
rise_time = 0.85; % Tiempo de subida del pulso [ns]
gaussian = amplitude .* exp(-((time-delay)./(rise_time/2)).^2);
%% Creaci�n del pulso en estado estable
% Cuando el pulso alcanza su amplitud m�xima debe conservar la misma
% amplitud por el resto de la simulaci�n. Para lograrlo se utiliza el
% perfil de subida concatenado con un vector contante en 250 mV.
iterations = length(time);
[~,peak_index] = max(gaussian);
Vinc = [gaussian(1:peak_index) amplitude.*ones(1,iterations-peak_index)];
%% Eliminaci�n de los efectos del cabezal de la sonda
head_delay = (head_length * sqrt(eps_head) / c_o) * 10^9;
total_delay = cable_delay + head_delay;
delay_index = find(time >= total_delay,1);
rho_new = [zeros(1,delay_index) rho(delay_index+1:end)];
%% Recuperaci�n del voltaje reflejado
Vref = rho_new .* Vinc;
cut_off_freq = 12e9; % Cut-off frequency [Hz]
cut_off_samples = floor(1 / (cut_off_freq * dt)); 
Vref_smoothed = smooth(Vref,cut_off_samples)';

Vref_prime_smooth = diff(Vref_smoothed);
Vref_prime = [Vref_prime_smooth Vref_prime_smooth(end)]./tsamp;
%% Frequency vector
fsamp = (1/dt)*1E-9;
freq = (-Qsamples/2:Qsamples/2-1)*fsamp/Qsamples;
%% S11 parameter
Vinc_smoothed = smooth(Vinc,cut_off_samples)';
Vinc_prime = diff(Vinc_smoothed);
Vinc_prime = [Vinc_prime Vinc_prime(end)]./tsamp;

S11 = fftshift(fft(Vref_prime))./fftshift(fft(Vinc_prime));
end