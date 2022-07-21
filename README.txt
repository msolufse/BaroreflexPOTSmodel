README File for UPEMD analysis of Baroreflex

Code is based on the study:
"Postural Orthostatic Tachycardia Syndrome (POTS) explained using a baroreflex response model" by J. R. Geddes, Johnny T. Ottesen, J. Mehlsen, and M. S. Olufsen. DOI: https://doi.org/10.48550/arXiv.2109.14558 

This software simulates the closed-loop model outlined in the above publication.

Usage:
The script "DriverBasic.m" is the main file that solves the model. The file loads in parameters from the file "load_global.m" and then solves the model beat by beat in a while loop by calling the function "model.m" inside the ordinary differential equation solver "ode15s" (MATLAB function). The equations in "model.m" depend on the functions "ElastanceBasic.m", which calculates the elastance of the left heart at a given time point in the cardiac cycle, and "tiltftn.m" which calculates the pressure contribution from gravity for a given tilt angle. The second part of "DriverBasic.m" plots the heart rate, upper arterial pressure, and mean carotid pressure of the simulation. Running the entire file "DriverBasic.m" will solve the model during rest and head-up tilt and plot the heart rate, upper arterial pressure, and mean carotid pressure.


-DriverBasic-
Input:
Output: Graph of heart rate, upper arterial pressure, and mean carotid pressure of the simulation.

-model-
Input: time (t), vector of states (y), vector of parameters (pars), the start of the current cardiac cycle (ts), and length of the current cardiac cycle (T).
Output: Right-hand side of model differential equations

-load_global-
Input: Hill coefficient for heart rate control (kH), Hill coefficient for resistance control (kR), total blood volume (BV)
Output: Parameters used in the model.m and initial conditions for differential equations

-ElastanceBasic-
Input: time (t), time for contraction of the heart (Ts), time for relaxation of the heart (Tr), end-diastolic elastance (Ed), end-systolic elastance (Es)
Output: Value of elastance of left ventricle (left heart) at a given time

-tiltftn-
Input: time (t), parameter vector (pars)
Output: Gravity pressure contribution (rhogh) and tilting argument (arg)




License:										
Copyright (C) 2022 J. R. Geddes, Johnny T. Ottesen, J. Mehlsen, and M. S. Olufsen

Contact information:									
Mette Olufsen (msolufse@ncsu.edu)
North Carolina State University
Raleigh, NC
 
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, and merge the Software subject to the following conditions: 

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

The authors shall be cited in any work resulting from the use of the Software. The associated published article is https://doi.org/10.48550/arXiv.2109.14558 

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, OR DATA, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE

