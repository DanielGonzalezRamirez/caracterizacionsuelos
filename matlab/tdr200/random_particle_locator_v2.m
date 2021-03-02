function [f, eps] = random_particle_locator_v2(freq_l, freq_h, eps_l, eps_h)
%Places the particles in a random location on each of its dimensions within
%the range given.
%   freq_l: lower limit of the relaxation frequency range.
%   freq_h: higher limit of the relaxation frequency range.
%   eps_l: lower limit of the eps_infty range.
%   eps_h: higher limit of the eps_infty range.
    
f = freq_l + (freq_l+freq_h)*rand;
eps = eps_l + (eps_l+eps_h)*rand;
end

