% Joel Lubinitsky
% Returns viscosity using Sutherland's rule (At SSL)
% 03/03/15

function visc = sutherlandVisc(temperature, units)

if     units == 'm'
    
    viscRef = 1.7894*10^-5; % kg/(m*s)
    tempRef = 288.16;       % K
    const   = 110;          % K
    
elseif units == 'e'
    
    viscRef = 3.7372*10^-7; % slug/(ft*s)
    tempRef = 518.688;      % R
    const   = 198;          % R
    
else
    
    error('Units must be Metric (m) or English (e)')
    
end

visc = viscRef * (temperature / tempRef) ^ 1.5 * ((tempRef + const) /...
                                              (temperature + const));