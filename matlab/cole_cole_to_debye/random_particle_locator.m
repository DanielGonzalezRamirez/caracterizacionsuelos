function [relax_freqs] = random_particle_locator(Np, ld, hd)
%Ubica las frecuencias de relajación para cada partícula.
%Para ubicar las frecuencias garantiza que el número aleatorio seleccionado
%tenga distancia del número anterior.
%   Np: Número de frecuencias de relajación a ubicar (número de polos).
%   ld: Década de frecuencia inferior
%   hd: Década de frecuencia superior
    Qdec = hd - ld; % Cantidad de polos
    dec_per_pole = Qdec / Np; % Cantidad de décadas por polo
    
    relax_dec = zeros(1,Np);
    current_dec = ld;
    
    for i = 1:Np
        relax_dec(i) = current_dec + dec_per_pole .* rand;
        current_dec = current_dec + dec_per_pole;
    end

    relax_freqs = 10 .^ relax_dec;
end

