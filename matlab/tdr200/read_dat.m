% Lectura del archivo .dat de salida del TDR200
% El TDR200 de Campbell Scientific se conecta con el computador para hacer
% lecturas de la señal reflejada por una sonda coaxial. A partir de estas
% lectura se puede determinar la permitividad compleja y la conductividad 
% del dieléctrico donde está insertada la sonda.
% En esta función lee el archivo .dat de salida que contiene los datos
% de la lectura. 
% Argumentos:
%   fullfilename: Ruta del archivo .dat a leer
%   rows: Cantidad de filas del archivo a leer
% Retorno:
%   time: vector con la base de tiempo.
%   rho: vector con el valor del coeficiente de reflexión
%
% Daniel Julián González Ramírez
% dj.gonzalez1203@uniandes.edu.co
% 2019-01-30
function [time, rho, header] = read_dat(fullfilename, rows)
% Al leer manualmente el archivo .dat es posible identificar que el
% delimitador entre valore es una ',', que la cabecera de las mediciones es
% de 14 valores y que los datos de cada medición se guardan en una fila del
% archivo.
delimiter = ',';
header_lenght = 14;

% Primero se recupera la información de la cabecera, en especial para
% conocer cuantas muestras se tomaron en la lectura.
fileID = fopen(fullfilename,'r');
header = textscan(fileID, '%*q%*q%*q%*q%*q%f%f%f%f%f%f%f%f%*q', 1,...
    'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);

averages = header{1};
vp = header{2};
Qsamples = header{3};
cable_length = header{4};
window = header{5};
probe_length = header{6};
probe_offset = header{7};
eps_head = header{8};

eps_o = 8.854187817e-12;    % Permitividad del vacío
mu_o = 1.256637061e-6;      % Permeabilidad del vacío
c_o = 1/sqrt(eps_o*mu_o);   % Velocidad de propagación en vacío

aux_rho = zeros(1,Qsamples);
for i=1:Qsamples
    aux_rho(i) = cell2mat(textscan(fileID, '%f',1,'Delimiter', delimiter, ...
        'EmptyValue' ,NaN,'ReturnOnError', false));
end

time_cable = (cable_length)/c_o * 1e9; % [ns]
total_distance = cable_length + window; % [m]
total_time = (total_distance)/c_o * 1e9; % [ns]

aux_time = linspace(time_cable,total_time,Qsamples);

rho = zeros(1,rows*Qsamples);
time = zeros(1,rows*Qsamples);

rho(1:Qsamples) = aux_rho;
time(1:Qsamples) = aux_time;

for count_rows = 1:rows-1
    header = textscan(fileID, '%*q%*q%*q%*q%*q%f%f%f%f%f%f%f%f%*q', 1,...
    'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
    
    averages = header{1};
    vp = header{2};
    Qsamples = header{3};
    cable_length = header{4};
    window = header{5};
    probe_length = header{6};
    probe_offset = header{7};
    eps_head = header{8};
    
    for i=1:Qsamples
        aux_rho(i) = cell2mat(textscan(fileID, '%f',1,'Delimiter', delimiter, ...
            'EmptyValue' ,NaN,'ReturnOnError', false));
    end
    
    time_cable = (cable_length)/c_o * 1e9; % [ns]
    total_distance = cable_length + window; % [m]
    total_time = (total_distance)/c_o * 1e9; % [ns]

    aux_time = linspace(time_cable,total_time,Qsamples);
    
    curr_sample = Qsamples*count_rows;
    rho(curr_sample+1:curr_sample + Qsamples) = aux_rho;
    time(curr_sample+1:curr_sample + Qsamples) = aux_time;
end





