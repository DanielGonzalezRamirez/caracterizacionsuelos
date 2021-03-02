function [line] = line_calc(Vref,time,dt,time_min,loc,type)
%LINE_CALC Calculates the line that joins two points
%   This function calculates the slope and interception point for the lines
%   that describe the base and inflection intervals.
if type == 1
    t_low = time(loc) - dt;
    t_high = time(loc);
else
    t_low = time(loc) - dt/2;
    t_high = time(loc) + dt/2;
end
low_loc = find(time >= t_low, 1);
high_loc = find(time >= t_high, 1);

Vref_line = Vref(low_loc:high_loc);
time_line = time(low_loc:high_loc);

a = (Vref_line(end) - Vref_line(1)) / (time_line(end) - time_line(1));
b = Vref_line(end) - a*time_line(end);

line = a.*time_min + b;
end

