% **************************************************************************
%                       M O D E L  f o r						
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
% model.m (t,y,pars,ts,T)
% t    - time
% y    - vector of states
% pars - parameter vector
% ts   - start of cardiac cycle
% T    - length of current cardiac cycle
%
% Model equations for the closed loop cardiovascular-baroreflex model outlined in
% the study "Postural Orthostatic Tachycardia Syndrome (POTS) explained using 
% a baroreflex response model" by J. R. Geddes, Johnny T. Ottesen, J.
% Mehlsen, and M. S. Olufsen.
%
% Calculates the right hand side of model differential equations
% 
% Dependencies: ElastanceBasic.m, load_global.m, tiltftn.m
%
% Output: Derivatives of model equations for ode15s to utilize
%*********************************************************************************

function xdot = model(t,y,pars,ts,T)

global Rav Rmv    % Declares global variables - atrial valve and mitrial valve resistance

% Define state variables
Vau    = y(1);   % Upper arterial volume
Vvu    = y(2);   % Upper venous volume
Val    = y(3);   % Lower arterial volume
Vvl    = y(4);   % Lower venous volume
Vlv    = y(5);   % Left ventricular volume (left heart)
pcm    = y(6);   % Mean carotid pressure
Raup   = y(7);   % Upper peripheral resistance
Ralp   = y(8);   % Lower peripheral resistance
Ed     = y(9);   % End diastolic elastance
Hc     = y(10);  % Heart rate

% Define resistances
Ral   = pars(1); % Resitance from upper to lower body arteries
Rvl   = pars(2); % Resistance from lower to upper body veins

% Define compliances
Cau = pars(3);   % Upper arterial compliance
Cal	= pars(4);   % Lower arterial compliance
Cvu = pars(5);   % Upper venous compliance

% Pressures
pau = Vau./Cau;  % Upper arterial pressure
pal = Val./Cal;  % Lower arterial pressure

VM_vl = pars(27);% Max lower body vein volume
m_vl = pars(28); % paramter linking volume and compliance
pvl = (1/m_vl) * log(VM_vl/(VM_vl - Vvl)); % Lower body vein pressure
pvu = Vvu./Cvu;  % Upper body vein pressure

% Heart parameters
Ts = 0.001.*(0.82/1.82)*(522-1.87*60/T); % Calculate Ts from formula
Tr = 0.001.*(1/1.82)*(522-1.87*60/T);    % Calculate Tr from formula 
Es     = pars(6); % Systolic elastance
Vd     = pars(7); % Unstressed left ventricular volume

% Tilt Function
rhogh  = tiltftn(t,pars); %Calculate pressure contribution from tilt

%Carotid pressure
tilde_pars = pars;      %Create new vector of paramters
tilde_pars(end) = 20;   %Distance from aorta and carotid
rhogh_tilde  = tiltftn(t,tilde_pars); %Pressure difference between aorta and carotid
pc = pau - rhogh_tilde; %Calculate carotid pressure

Elv  = ElastanceBasic(t-ts,Ts,Tr,Ed,Es); % Ventriular Elastance function 
plv  = Elv*Vlv; % Stressed ventricular volume & Vd(consant), ventricular volume at zero dyastolic pressure

%Calculate flows
if plv > pau  % "if" statement acts as a valve
   qav  = (plv - pau)/Rav; % Calculate flow using Ohm's law
else
   qav = 0;   % Don't allow backward flow
end

if pvu > plv  % "if" statement acts as a valve
  qmv  = (pvu - plv)/Rmv; % Calculate flow using Ohm's law
else
  qmv = 0;    % Don't allow backward flow
end

qal  = (pau - pal + rhogh)/Ral; % Calculate flow using Ohm's law
qaup = (pau - pvu)/Raup;        % Calculate flow using Ohm's law
qalp = (pal - pvl)/Ralp;        % Calculate flow using Ohm's law

if pvl > pvu   % "if" statement acts as a valve
    qvl  = (pvl - pvu - rhogh)/Rvl; % Calculate flow using Ohm's law
else
    qvl = 0;   % Don't allow backward flow
end

% Mean carotid pressure
tauPm = pars(29);       %Time constant for tracking mean carotid pressure
dpcm  = (pc-pcm)/tauPm; %ODE to track mean carotid pressure

% Resistance Control
kR   = pars(8); %Hill coefficient for resistance control
taur = pars(9); %Time scale for resistance control

RaupM= pars(10); %Max upper peripheral resistance control
Raupm= pars(11); %Min upper peripheral resistance control
p2Ru = pars(12); %Half saturation upper peripheral resistance control

RalpM= pars(13); %Max lower peripheral resistance control
Ralpm= pars(14); %Min lower peripheral resistance control
p2Ra = pars(15); %Half saturation lower peripheral resistance control

%Upper peripheral control Hill equation
Raupf = (RaupM - Raupm)*p2Ru^kR/(pcm^kR + p2Ru^kR) +Raupm;
dRaup = (-Raup+Raupf)/taur; % Upper peripheral control ODE

%Lower peripheral control Hill equation
Ralpf = (RalpM - Ralpm)*p2Ra^kR/(pcm^kR + p2Ra^kR) +Ralpm; 
dRalp = (-Ralp+Ralpf)/taur; % Lower peripheral control ODE



%Cardiac contractility (decreasing minimum elastance)
kE  = pars(16);% Hill coefficient for elastance control
tauE= pars(17);% Time constant for elastance control
EdM = pars(18);% Max response for elastance control
Edm = pars(19);% Min response for elastance control
p2E = pars(20);% Half saturation for elastance control

% End diastolic elastance control Hill equation
Edf = (EdM - Edm)*pcm.^kE./(pcm.^kE + p2E^kE) +Edm;
dEd = (-Ed + Edf)/tauE; %End diastolic elastance control ODE

%heart rate
kH  = pars(21);% Hill coefficient for heart rate control
tauH= pars(22);% Time constant for heart rate control
HM  = pars(23);% Max response for heart rate control
Hm  = pars(24);% Min response for heart rate control
p2H = pars(25);% Half saturation for heart rate control

% Heart rate control Hill function
Hf = (HM-(Hm))*p2H^kH./(pcm.^kH+p2H^kH) + Hm;
dHc = (-Hc + Hf)/tauH; %Heart rate control ODE


% Right hand side
xdot = [qav - qal - qaup; ...   %dVau 
        qvl + qaup - qmv; ...   %dVvu 
        qal - qalp; ...         %dVal
        qalp - qvl; ...         %dVvl
        qmv - qav ; ...         %dVlv
        dpcm; dRaup; dRalp; dEd; dHc];