%% Algoritmo PSO-Lineal para aproximar modelo Cole-Cole a Debye
% Este script desarrolla el algoritmo propuesto en [1] para aproximar el
% comporatmiento de la permitividad compleja, descrita por el modelo
% Havriliak-Negami, a un modelo de Debye con suma de polos.
%
% El modelo Havriliak-Negami describe la permitividad de diferentes
% dieléctricos. Los modelos Cole-Cole y Debye son casos especiales del 
% modelo Havriliak-Negami.
%
% Modelo Havriliak-Negami:
% eps_r* = eps_inf + (eps_s - eps_inf)/(1+(j*2*pi*f*T_0)^(1-a))^b
% Cole-Cole: b = 1. Debye: a = 0 y b = 1.
%
% [1] Kelley, D. F., Destan, T. J., & Luebbers, R. J. (2007). 
% Debye function expansions of complex permittivity using a hybrid 
% particle swarm-least squares optimization approach. 
% IEEE Transactions on Antennas and Propagation, 55(7), 1999–2005. 
% https://doi.org/10.1109/TAP.2007.900230
% 
% Daniel Julian Gonzalez Ramirez
%
% 2019-02-13
clear, clc
close all
%% Funciones anónimas
% Función Havrialiak_Negami
Havriliak_Negami = @(eps_inf, eps_s, f, t_o, alpha, beta) eps_inf +...
    (eps_s - eps_inf)./ ((1 + (1i * 2 * pi * t_o .*f).^(1 - alpha)).^beta);
%% Definición del modelo Cole-Cole a aproximar e intervalo de frecuencias
% Parámetros Cole-Cole:
% t_o = 11e-9; % Tiempo de relajación [s]
% alpha = 0.25; % Constante del exponente de las frecuencias
% beta = 1; % Constante del exponente del denominador
% eps_s = 8.9; % Permitividad relativa en frecuencia 0
% eps_inf = 5.6; % Permitividad relativa en frecuencia infinita

t_o = 1/(2*pi*0.4e9); % Tiempo de relajación [s]
alpha = 0.11; % Constante del exponente de las frecuencias
beta = 1; % Constante del exponente del denominador
eps_s = 25; % Permitividad relativa en frecuencia 0
eps_inf = 4.7; % Permitividad relativa en frecuencia infinita

% Definición del intervalo de frecuencias: El cálculo del algoritmo
% requiere que el dominio de la frecuencia se encuentre distribuido
% logarítmicamente. Para distribuir el intervalo se entre que decadas de
% frecuencia se trabajará y se utiliza la función logspace para crear el
% arreglo de frecuencias.
% low_dec = 5; high_dec = 10; % Decadas de frecuencia a usar
low_dec = 7; high_dec = 12; % Decadas de frecuencia a usar
Qpoints = 4096; % Cantidad de puntos para el arreglo
freq = logspace(low_dec,high_dec, Qpoints); % Arreglo de frecuencia [Hz]

% A partir de la cantidad de décadas se calcula el número de polos a
% utilizar en la función Debye, el criterio es: #.dec <= Np <= 1.5*#.dec
%Np = floor(1.5*(high_dec - low_dec)); % Cantidad de polos a utilizar
Np = 3;

% Función Cole-Cole y gráfica:
eps_Cole_Cole = Havriliak_Negami(eps_inf, eps_s, freq, t_o, alpha, beta);
target_eps = abs(imag(eps_Cole_Cole));

figure('Name','Cole-Cole description of permittivity')
semilogx(freq, real(eps_Cole_Cole), 'LineWidth', 2); hold on
semilogx(freq, abs(imag(eps_Cole_Cole)), 'LineWidth', 2);
xlabel('Frequency [Hz]'); ylabel('Relative permitivitty'); grid minor
legend('Real', 'Imag'); title('Cole-Cole description of permittivity');

%% Condiciones inidicales para el algoritmo de optimización PSO-Lineal
% Se deben tener aproximadamente 4 muestras de frecuencia para cada polo 
% de la función. Las muestras de frecuencia preferiblemente deben ser
% igualmente espaciadas.
Nf = 4*Np;
sample_point = zeros(Nf, 1);
sample_freq = zeros(Nf, 1);

for i = 1:Nf
    sample_point(i) = round(i * (Qpoints / Nf));
    sample_freq(i) = freq(sample_point(i));
end

% Ubicación inicial de las partículas:
% Cada partícula se compone de una cantidad Np de frecuencias de
% relajación. Se puede entender cada polo como una dimensión para el
% dominio de la partícula, es decir: si se representa el modelo Debye como
% la combinación lineal de 3 polos, la partícula debe ubicarse en 3
% dimensiones (3 frecuencias de relajación); se se presenta el modelo a
% través de 7 polos, la partícula debe ubicarse en 7 dimensiones (7
% frecuencias de relajación).
% Se ubican las frecuencias de relajación a optimizar (partículas) de
% forma aleatoria en el rango de frecuencias definido anteriormente.
Qparticles = 30; % Cantidad de partículas a utilizar
particles = zeros(Qparticles, Np); % Arreglo de particulas
velocities = zeros(Qparticles,Np); % Arreglo de velocidades de partículas

for i = 1:Qparticles
    particles(i,:) = random_particle_locator(Np, low_dec, high_dec);
end

%% Optimización lineal para encontrar los pesos
% El algorítmo en cuestión encuentra los parámetros optimizados para la
% parte imaginaria de la permitividad. Esto es gracias a que si la parte
% imaginaria de la permitividad es descrita por la suma de funciones de
% Debye, entonces la parte real debe describirse por la parte real de la
% misma función. Esto se cumple, más una constante que se añade a la parte
% real [1]. 

% La ecuación que describe la parte imaginaria es:
% eps''(f) = (eps_s - eps_inf)*sum_{p=1}^N_p a_p (f/f_p) / (1 + (f/f_p)^2)
% Por facilidad, esto se puede escribir de forma matricial c = Da. Donde c
% es el vector de expresiones: eps''(f)/(eps_s - eps_inf). D es la matriz
% que contiene las expresiones de la forma: (f/f_p) / (1 + (f/f_p)^2). Y a
% es el vector que contiene los pesos de las partículas, a_p. Donde f es 
% una muestra en frecuencia y f_p es la frecuencia de relajación para una 
% partícula.

% Se crea el vector c, que es igual para todas las partículas:
c = c_vector(Nf, eps_Cole_Cole, sample_point, eps_s, eps_inf);

% Previamente se definen los vectores para almacenar el ajuste de cada
% partícula al modelo original:
% Almacena el ajuste actual de cada partícula:
current_fitness = zeros(Qparticles, 1);
% Almacena el mejor ajuste (particle best fitness) de cada partícula:
pBestFitness = zeros(Qparticles, 1);
% Almacena el valor de los parámetros, frecuencias y pesos, de mejor ajuste 
% para cada partícula (particle best):
pBest = zeros(Qparticles, 2*Np);
% Almacena el mejor ajuste general:
gBestFitness = 0;
% Almacena el valor de los parámetros, frecuencias y pesos, de mejor ajuste 
% general de todas las partículas:
gBest = zeros(1, 2*Np);

% Definición de las constantes de actualización de velocidad de las
% partículas:
w = 0.5; c1 = 1.4; c2 = 1.4;

% Desde este punto el procedimiento es iterativo para cada partícula. 
% La condición de salida que se utiliza es por cantidad de iteraciones:
Qiterations = 100;

for c_iter = 1:Qiterations
    for c_part = 1:Qparticles
        % Se crea la matriz D:
        D = D_matrix(Nf, Np, sample_freq, particles(c_part,:));

        % Los pesos se calculan a través de la expresión: 
        % a = (D'D + gammaI)^-1 D'c. Donde gamma es una constante arbitraria 
        % muy pequeña, que se incrementa iterativamente con el proposito que 
        % 'a' no tenga valores negativos.
        gamma = 0.01;
        % Se declaran dos condiciones de salida de la iteración que calcula los
        % pesos. La primera, es si todos los pesos hayados son positivos y por 
        % lo tanto se tiene ya una solución válida; está representada por la
        % bandera "flag_weights". La segunda, es si se han hecho más de cierto
        % número de intentos con diferentes valores iniciales de las
        % partículas; se representa a través de la variable "gamma_breached" la
        % cual no puede exceder un valor dado, en este caso 50.

        flag_weights = 0; gamma_breached = 0;
        while ~flag_weights && gamma_breached < 50
            weights = inv(D'*D + gamma.*eye(Np))*(D'*c);

            if isempty(find(weights < 0,1))
                flag_weights = 1;
            else
                gamma = gamma + 0.01;
                if gamma > 0.1
                    gamma_breached = gamma_breached + 1;

                    particles(c_part,:) =...
                        random_particle_locator(Np, low_dec, high_dec);
                    D = D_matrix(Nf, Np, sample_freq, particles(c_part,:));
                end
            end
        end

    %% Evaluación de "goodness of fit"
    % Para evaluar que tan buena es la solución encontrada se realiza
    % una prueba de ajuste ("goodness of fit").

        % Para hacer la prueba se debe construir la parte imaginaria
        % del modelo Debye con los parámetros hallados.

        debye_poles = zeros(Np, Qpoints);
        for i = 1:Np
            debye_poles(i,:) = weights(i).*((freq./particles(c_part,i))...
                ./(1 + (freq./particles(c_part,i)).^2));
        end
        debye_imag = (eps_s - eps_inf) .* sum(debye_poles);
        
        % La prueba de ajuste es la suma de las diferencias cuadradas entre
        % el modelo Cole-Cole y la función Debye obtenida.
        
        current_fitness(c_part) = sum(target_eps.^2 - debye_imag.^2);
        
        % Se comprueba si el ajuste actual es el mejor ajuste que ha tenido
        % la partícula. Si es el mejor ajuste, se almacena el valor del
        % ajuste, en el arreglo pBestFitness, y el valor de los parámetros 
        %(frecs. de relajación y pesos), en el arreglo pBest.
        if c_iter == 1
            pBestFitness(c_part) = current_fitness(c_part);
            pBest(c_part,:) = [particles(c_part,:) weights'];
        else
            if current_fitness(c_part) < pBestFitness(c_part)
                pBestFitness(c_part) = current_fitness(c_part);
                pBest(c_part,:) = [particles(c_part,:) weights'];
            end
        end
    end
    
    % Ahora se almacena el mejor ajuste de la partícula que tuvo el mejor
    % ajuste entre todas. Y se comprueba si el ajuste obtenido es mejor que
    % el anterior mejor ajuste general. Si el ajuste es mejor, se almacena
    % el valor del ajuste, en gBestFitness, y el valor de los parámetros en
    % el arreglo gBest.
    [current_best, best_pos] = min(pBestFitness);
    
    if c_iter == 1
        gBestFitness = current_best;
        gBest(:) = pBest(best_pos,:);
    else
        if current_best < gBestFitness
            gBestFitness = current_best;
            gBest(:) = pBest(best_pos,:);
        end
    end
    
    %% Cálculo de las nuevas velocidades y actualización de la posición
    % La velocidad que toman las partículas es función de los mejores
    % ajustes, de la partícula y en general. Como el método PSO se aplica
    % únicamente sobre las frecuencias de relajación y no sobre los pesos,
    % el único parámetro de optimización que se tiene en cuenta para la 
    % actualización de la velocidad es el valor delas frecuencias de 
    % relajación. También se utilizan las velocidades anteriores para la
    % velocidad actual, un generador de números aleatorios úniformemente
    % distribuidos entre [0,1] y los pesos w, c1, y c2 (definidos
    % previamente en el código).
    % Una vez se tienen las velocidades calculadas, se actualiza la
    % posición de las partículas. La nueva posición corresponde a la
    % posición actual más la velocidad.
    
    for c_part = 1:Qparticles
        velocities(c_part,:) = w.*velocities(c_part,:) +...
            c1 * rand .* (pBest(c_part,1:Np) - particles(c_part,:)) +...
            c2 * rand .* (gBest(1:Np) - particles(c_part,:));
        particles(c_part,:) = particles(c_part,:) + velocities(c_part,:);
    end 
end
%% Recuperación de los datos con mejor ajuste y gráfica de la función
best_freq = gBest(1:Np);
best_weights = gBest(Np + 1:end);

best_debye_poles = zeros(Np, Qpoints);

for i = 1:Np
    best_debye_poles(i,:) = best_weights(i).*((freq./best_freq(i))...
        ./(1 + (freq./best_freq(i)).^2));
end

best_debye_imag = (eps_s - eps_inf) .* sum(best_debye_poles);

figure('Name', 'Best Imaginary Debye fit')
semilogx(freq, abs(imag(eps_Cole_Cole)), 'LineWidth', 2); hold on
semilogx(freq, abs(best_debye_imag), 'LineWidth', 2); grid minor
xlabel('Frequency [Hz]'); ylabel('Relative permittivity');
legend('Original model', 'Debye approximation');
title('Imaginary permitivitty fit of the Debye approximation');

disp('Frecuencias de relajación')
disp(best_freq);
disp('Pesos')
disp(best_weights);

%% Cálculo de la función Debye compleja para la permitividad relativa
% Se calcula del factor de offset para la parte real de la función Debye:
delta_eps = calc_delta_eps(Nf, eps_Cole_Cole, sample_point, sample_freq, eps_s, eps_inf, best_weights, best_freq);

% Se calcula el valor de los polos Debye:
Debye_poles = zeros(Np, Qpoints);
for i = 1:Np
    Debye_poles(i,:) = best_weights(i)./(1 + 1i.*freq./best_freq(i));
end

% Se obtiene el modelo Debye completo:
eps_Debye = eps_inf + (eps_s - eps_inf) * sum(Debye_poles) + delta_eps;

figure('Name', 'Best Debye fit')
semilogx(freq, real(eps_Cole_Cole), 'LineWidth', 3); hold on
semilogx(freq, abs(imag(eps_Cole_Cole)), 'LineWidth', 3); hold on
semilogx(freq, real(eps_Debye), 'LineWidth', 2); hold on
semilogx(freq, abs(imag(eps_Debye)), 'LineWidth', 2); grid minor
xlabel('Frequency [Hz]'); ylabel('Relative permittivity');
legend('Original model real', 'Original model imag', 'Debye real', 'Debye imag');
title(['Permitivitty fit of the Debye approximation with ', num2str(Np), ' poles']);