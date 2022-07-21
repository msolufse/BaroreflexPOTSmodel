% **************************************************************************
%                       L O A D  G L O B A L  f o r						
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
% load_global.m (kH,kR,BV)
% kH - Hill coefficient for heart rate control
% kR - Hill coefficient for resistance control
% BV - Total blood volume
%
% Parameters for the closed loop cardiovascular-baroreflex model outlined in
% the study "Postural Orthostatic Tachycardia Syndrome (POTS) explained using 
% a baroreflex response model" by J. R. Geddes, Johnny T. Ottesen, J.
% Mehlsen, and M. S. Olufsen.
%
% Calculates the paramters used in the model differential equations
% 
% Dependencies: 
%
% Output: Paramaters used in model.m and initial conditions for
% differential equations
%*********************************************************************************


function [x0, Init] = load_global(kH,kR,BV)

global ODE_TOL    %Declare ODE_TOL as a global variable
global Rav Rmv    %Declare atrial and mitrial valve resistance as global variables

ODE_TOL  = 1e-8;  %Set ODE solver (ode15s) tolerance

TotalVol = BV;             % Let TotalVol be equal to the input BV
TotFlow = TotalVol/60;     % Assume blood circulates in one min Cardiac output [ml/sec]
COd     = TotFlow*60/1000; % Steady state only [ l/min ]

HI      = .96;             % initial heart rate [beats/sec]

% Flows (related to subject) 
qaup = TotFlow*0.80; % Upper body flow, arteries --> veins
qal  = TotFlow*0.20; % Upper body arteries --> lower body arteries
qvl  = qal;          % Lower body veins --> upper body veins
qalp = qal;          % Lower body flow, arteries --> veins

% Pressures (related to subject)     
pauD   = 80;         % Ideal diastolic pressure
pauS   = 120;        % Ideal systolic pressure
pau    = (2/3)*pauD+(1/3)*pauS;  % Mean systolic pressure data (Upper body arteries) - medical formula
pm     = pau;        % Initial mean pressure

pal    = pau*0.99;   % Pressure in lower body arteries
pvu    = 2.75;       % Pressure in upper body veins - reflects pressures in pulmonary veins
pvl    = 3.00;       % Pressure in lower body veins
plvD   = 2.5;        % Left ventricle diastolic pressure
Vd     = 10;         % Left ventricular unstressed (residual) volume
Vlvm    = 50-Vd;     % Minimum left ventricular volume 
VlvM   = 110-Vd;     % Maximum left ventricular volume

% Resistances (Ohm's law)
RaupI = (pau-pvu)/qaup;  % Initial upper body peripheral resistance
Ral   = (pau-pal)/qal;   % Upper body arteries --> lower body arteries
Rvl   = (pvl-pvu)/qvl;   % Lower body veins --> upper body veins
RalpI = (pal-pvl)/qalp;  % Initial lower body peripheral resistance
Rav   = 0.0001;          % Atrial valve resistance
Rmv   = 0.0001;          % Mitrial valve resistance

TotalVolSV = TotalVol*0.85; % 85 percent of blood volume in systemic circulation
Vart = TotalVolSV*.15;      % 15 percent arteries and arterioles
Vven = TotalVolSV*.85;      % 85 percent veins, capillaries, central veins

% Volume distribution (Beneken and deWit)
Vau = Vart*0.80*.30;     % 80 percent in upper body arteries of these 70 percent is unstressed
Val = Vart*0.20*.30;     % 20 percent in lower body arteries of these 70 percent is unstressed
Vvu = Vven*0.80*.075;    % 80 percent in upper body veins of these 92.5 percent is unstressed
Vvl = Vven*0.20*.075;    % 20 percent in lower body veins of these 92.5 percent is unstressed
                  
% Compliances, stressed volume percentages from Beneken are weighted
Cau = Vau/pauD;          % Upper body artery compliance 
Cvu = Vvu/pvu;           % Upper body venous compliance
Cal  = Val/pal;          % Lower body artery compliance 

% Lover body volumes
VMvl = 4*Vvl;                       %Choosen based on BP drop without controls 
mvl  = log(VMvl/(VMvl - Vvl))/pvl;  %Solving for desired initial value

% Ventricle parameters 
Tsfrac = 0.12;     % Time fraction for systolic phase, increasing elasticity 
Trfrac = 0.14;     % Time fraction for Systolic phase, decreasing elasticity
% % % % % % % % % EmI    = 0.0273;   % Minimum elastance of the heart [ Em = pveins / (Max Vlv - Vd) = 3 / (120-10) ]
% % % % % % % % % EM     = 3.6667;   % Maximum elastance of the heart [ EM = part / (min Vlv - Vd) = 120/(40-10) ]  % 120 is max arterial pressure
EdI     = plvD/VlvM;  % Initial minimum elastance of the heart [ Ed = pveins / (Max Vlv - Vd) 
Es     = pauS/Vlvm;   % Maximum elastance of the heart [ EdM = part / (min Vlv - Vd)

% Set initial conditions from calculated values above
VauI = Vau; 
VvuI = Vvu;
ValI = Val;
VvlI = Vvl;
VlvI = VlvM;

%Resistance Upper Body Pars
taur   = 12.5;     %Resistance time constant 
RaupM = 3*RaupI;   %Upper body resistance max value
Raupm = 0.2*RaupI; %Upper body resistance min value
alpha = (RaupI-Raupm)/(RaupM-Raupm); % alpha used to calculate p2Ru
%Calculate upper half saturation to ensure proper initial placement on response curve
p2Ru  = (pm^kR * alpha/(1 - alpha))^(1/kR); 

%Resistance Lower Body Pars
RalpM  = 3*RalpI;   %Lower body resistance max value
Ralpm  = 0.2*RalpI; %Lower body resistance min value
alpha  = (RalpI-Ralpm)/(RalpM-Ralpm); % alpha used to calculate p2Ra
%Calculate lower half saturation to ensure proper initial placement on response curve
p2Ra   = (pm^kR * alpha/(1 - alpha))^(1/kR); 

%Contractility (Em)
EdM  = 1.25*EdI; % End diastolic elastance max value
Edm  = 0.01*EdI; % End diastolic elastance min value
tauE = 12.5;     % End diastolic elastance time constant
kE = 7; %        % End diastolic elastance Hill coefficient 
alpha= (EdI-Edm)/(EdM-Edm);% alpha used to calculate p2E
%Calculate half saturation to ensure proper initial placement on response curve
p2E  = (pm^kE*(1 - alpha)/alpha)^(1/kE); 


% Heart rate control parameters 
HM   = 200/60; %Max HR 
Hm   = .3;     %Min HR
tauH = 6.25;   %Time constant for HR
alpha= (HI-Hm)/(HM-Hm);% alpha used to calculate p2H
%Calculate half saturation to ensure proper initial placement on response curve
p2H  =(pm^kH * alpha/(1 - alpha))^(1/kH);

tauP = 2.5; %Time constant to track mean carotid pressure

Init = [VauI VvuI ValI VvlI VlvI pm RaupI RalpI EdI HI]; %Initial conditions



% Parameter vector
x0 = [Ral Rvl ...                %1-2
      Cau Cal Cvu...       %3-5
      Es Vd ...    %6-7
      kR taur RaupM Raupm p2Ru RalpM Ralpm p2Ra ... %8-15 
      kE tauE EdM Edm p2E...     %16-20  
      kH tauH HM Hm p2H,TotalVol,... %21-26
      VMvl, mvl,tauP]; %27-29
  
  
