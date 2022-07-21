% **************************************************************************
%               T I L T  F U N C T I O N  f o r						
%        L U M P E D  B A R O R E F L E X  P O T S  M O D E L	
% **************************************************************************
%												
% License:										
% Copyright (C) 2022 J. R. Geddes, Johnny T. Ottesen, J. Mehlsen, and M. S. Olufsen
%
% Contact information:									
% Mette Olufsen (msolufse@ncsu.edu)
% North Carolina State University
% Raleigh, NC
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, and merge the Software subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% The authors shall be cited in any work resulting from the use of the 
% Software. The associated published article is
% https://doi.org/10.48550/arXiv.2109.14558 
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES 
% WHATSOEVER RESULTING FROM LOSS OF USE, OR DATA, WHETHER IN AN ACTION OF 
% CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN 
% CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE
%
%%********************************************************************************

%%********************************************************************************
% tiltftn.m (t,pars)
% t     - time
% pars  - parameter vector
%
% Tilt function for the closed loop 
% cardiovascular-baroreflex model outlined in the study 
% "Postural Orthostatic Tachycardia Syndrome (POTS) explained using 
% a baroreflex response model" by J. R. Geddes, Johnny T. Ottesen, J.
% Mehlsen, and M. S. Olufsen.
%
% Calculates the pressure contribution used in the model differential equations
% 
% Dependencies: 
%
% Output: Pressure due to gravity used in model.m
%*********************************************************************************

function [rhogh,arg] = tiltftn(t,pars)


% Define tilt function parameters
rho       = 1.06;   % Density of blood in g/cm^3
g         = 982;    % Gravitational acceleration in cm/s^2
                    % Smithsonian physical tables, Smithsonian Institution,
                    % Seventh Revised Edition, prepared by Frerick E Fowle,
                    % Vol 71, #4, 1921

ts =  pars(end-2);  % Time for start of tilt
td =  pars(end-1);  % Time for end of tilt
h  =  pars(end);    % height between compartments in question

theta     = 60;             % Maximal tilt-angle, in degrees.
a         = 60/14;          % Tilt-speed - 60 degrees over the course of 14 seconds
conv      = 1333.22;        % conversion cgs (cm,grams,sec) units to mmHg
  
if t < ts                    %If at rest argument is 0
   arg = 0;
elseif t < ts+14             %If tilting upward argument depends on angle at time t
   arg = a.*(t-ts);
elseif t < td                %If fully tilted argument depends on the max tilt angle
    arg = theta;
elseif t < td+14
    arg = theta - a*(t-td);  %Tilt down speed (not currently used)
else
    arg = 0;                 %If patient back at rest argument is 0
end


% Calculate the hydrostatic pressure in mmHg, 1Pa=0.0075mmHg
rhogh = rho*g*h*(sin(arg*pi/180)/conv); 


end

