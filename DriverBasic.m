
% **************************************************************************
%                       D R I V E R  f o r						
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
% DriverBasic.m ()
%
% Simulates the closed loop cardiovascular-baroreflex model outlined in the
% study "Postural Orthostatic Tachycardia Syndrome (POTS) explained using 
% a baroreflex response model" by J. R. Geddes, Johnny T. Ottesen, J.
% Mehlsen, and M. S. Olufsen.
%
% Solves the proposed model, beat by beat, inside a while loop. 
% Model is written in the file "model.m"
% Paramters are calculated in called from file "load_global.m"
% Time varying Elastance function is calculated in the file
%   "ElastanceBasic.m"
% Pressure induced by head-up Tilt is calculated in the file "tiltftn.m"
%
% Dependencies: ElastanceBasic.m, load_global.m, model.m, tiltftn.m
%
% Output: Graph of simulation
%*********************************************************************************

clear
global ODE_TOL % Declare global paramter "ODE_TOL"
dt   = .01;    % Time step to evalulate solution at
eps  = 1e-6;   % Small amount added to last index to insure solving for full time
height = 25;   % Height between upper and lower body compartments (cm)

kR = 25;       % Input value of kR (Hill coefficient of Resistance Hill eq)
kH = 25;       % Input value of kH (Hill coefficient of heart rate Hill eq)
BV = 4500;     % Input value of blood volume (ml)
tup  = 200;    % Input length of rest to solve for (when to tilt upwards)
tend = tup + 100; % Input length of time for entire simulation

[x0, Init] = load_global(kH,kR,BV); % Load in parameters from file "load_global.m"
pars = [x0 tup tend height]; % Augment parameter vector with inputed quantities


% Initialize empty vectors 
timeS = [];
SolS  = [];

H     = Init(end);        %Initial Heart rate 
T     = round(1/H/dt)*dt; %Calculate length of first cardiac cycle

k1    = 1;              % Initialize index
k2    = round(T/dt)+k1; % Index of end of cardiac cycle
UP = 0;
while k2 < tend/dt+1+eps % Continue to solve model until time step reaches "tend"

    tdc = [k1-1:k2-1]*dt; % Desired time span to solve for cardiac cycle

    options = odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL); % Set options for ode15s
    sol     = ode15s(@model,[tdc(1) tdc(end)],Init,options,pars,tdc(1),T); % Solve the ODEs
    sols    = deval(sol,tdc); % Evaluate ODEs at desired times
    
    % Add solutions for current cycle to the previous solutions matricies
    timeS   = [timeS tdc(1:end-1)];
    SolS    = [SolS  sols(:,1:end-1)];

    Init    = sols(:,end); % Use the solutions at the last time step as the
    % initial conditions for the next cardiac cycle
    
    %Add noise to cardiac cycle length
    T       = round(1./sols(10,end)*(1+unifrnd(-1,1)*.02)/dt)*dt; 
    
    k1      = k2;        % Sets last index of this loop as first index for next loop
    k2 = k2+round(T/dt); % Redefining the last index
end

  timeS     = [timeS tdc(end)];    % Add last point to the vector
  SolS      = [SolS  sols(:,end)]; % Add last point to the vector
  
  % Calculate desired quantities from solution matrix
  Cau = pars(3);     % Upper arterial compliance (scalar paramter)
  VauS = SolS(1,:);  % Upper arterial volume (vector)
  pauS = VauS./Cau;  % Calculate upper arterial pressure (vector)
  HcS = SolS(10,:);  % Extract HR time series (vector)
 

%% Plotting

fontsize = 15; % Define desired fontsize
co = 'k';      % Choose line color

ylh = [.7,2.1];   % y-axis limits for heart rate plot
ylp = [50,140];   % y-axis limits for blood pressure plot

tbefore = 75;          % How much of simulation before tilt to plot
tafter = 75;           % How much of simulation after tilt to plot
tstartS = tup-tbefore; % Time value to start the plot
tendS = tup + tafter;  % Time value to end the plot

tilt_indm = find(abs(timeS-tup) == min(abs(timeS-tup)));    %index for tilt
ind1 = find(abs(timeS-tstartS) == min(abs(timeS-tstartS))); %index for start
ind2 = find(abs(timeS-tendS) == min(abs(timeS-tendS)));     %index for end
sm1 = ind1:ind2; %index span for simulation plot


figure(1)
clf
subplot(2,1,1)
hold on
plot(timeS(sm1),HcS(sm1),co,'linewidth',2) %Plot model heart rate
plot(ones(2,1).*timeS(tilt_indm),ylh,'--k','linewidth',2)  %Dashed line to denote tilt
xticks([])                   % Remove xticks on axis
ylabel('H (bps)')            % y-axis label
set(gca,'fontsize',fontsize) % Set font size
xlim([tstartS,tendS])        % limits on x-axis
ylim(ylh)                    % limits on y-axis

subplot(2,1,2)
hold on
plot(ones(2,1).*timeS(tilt_indm),ylp,'--k','linewidth',2) %Dashed line to denote tilt
plot(timeS(sm1),pauS(sm1),'r')                %Plot upper arterial pressure
plot(timeS(sm1),SolS(6,sm1),'k','linewidth',2)  %Plot mean carotid blood pressure
ylabel('P_{au} (mmHg)')       % y-axis label
xlabel('Time (s)')            % x-axis label
set(gca,'fontsize',fontsize)  % Set font size
xlim([tstartS,tendS])         % limits on x-axis
ylim(ylp)                     % limits on y-axis








