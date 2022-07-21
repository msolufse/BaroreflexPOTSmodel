% **************************************************************************
%               E L A S T A N C E  F U N C T I O N  f o r						
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
% ElastanceBasic.m (t,Ts,Tr,Ed,Es)
% t  - time
% Ts - time for contraction of heart
% Tr - time for relaxation of heart
% Ed - End diastolic elastance of heart
% Es - End systolic elastance of heart
%
% Time varying elastance function for the closed loop 
% cardiovascular-baroreflex model outlined in the study 
% "Postural Orthostatic Tachycardia Syndrome (POTS) explained using 
% a baroreflex response model" by J. R. Geddes, Johnny T. Ottesen, J.
% Mehlsen, and M. S. Olufsen.
%
% Calculates the left heart elastance used in the model differential equations
% 
% Dependencies: 
%
% Output: Elastance of left heart used in model.m
%*********************************************************************************


function Elh = ElastanceBasic(t,Ts,Tr,Ed,Es)


for i = 1:length(t) %Solve for all values of t inputted
    if t(i)<=Ts     %if time is less than end contraction time follow upwards contraction equation
        Elh(i) = Ed + ((Es-Ed)/2).*(1-cos(pi.*t(i)/Ts));
    elseif (t(i)<=Ts+Tr) && (t(i)>Ts) %if time is less than end relaxation time but greater than end contraction time follow downwards equation
        Elh(i) = Ed + ((Es-Ed)./2)*(cos((pi./Tr).*(t(i)-Ts)) + 1);
    else %If past the end relaxation output a constant elastance equal to end diastolic elastance
        Elh(i)=Ed;
    end
end