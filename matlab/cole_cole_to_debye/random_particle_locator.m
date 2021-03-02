function [relax_freqs] = random_particle_locator(Np, ld, hd)
%Ubica las frecuencias de relajaci�n para cada part�cula.
%Para ubicar las frecuencias garantiza que el n�mero aleatorio seleccionado
%tenga distancia del n�mero anterior.
%   Np: N�mero de frecuencias de relajaci�n a ubicar (n�mero de polos).
%   ld: D�cada de frecuencia inferior
%   hd: D�cada de frecuencia superior
    Qdec = hd - ld; % Cantidad de polos
    dec_per_pole = Qdec / Np; % Cantidad de d�cadas por polo
    
    relax_dec = zeros(1,Np);
    current_dec = ld;
    
    for i = 1:Np
        relax_dec(i) = current_dec + dec_per_pole .* rand;
        current_dec = current_dec + dec_per_pole;
    end

    relax_freqs = 10 .^ relax_dec;
end

