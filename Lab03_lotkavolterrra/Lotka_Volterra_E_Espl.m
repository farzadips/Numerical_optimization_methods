% function [xn,yn,x_eq,y_eq,tn] = Lotka_Volterra_E_Espl(a,b,c,d,...
%    x0,y0,t0,T,N)
%
% Function that return the time serieses of the Prey-Predator model of 
% Lotka-Volterra through the explicit Euler Method.
%
%   INPUTS:
%   a = preys' reproduction coefficient;
%   b = preys' death coefficient due to predators' hunting;
%   c = predators' death coefficient;
%   d = predators' survival % reproduction coefficient hunting preys;
%   x0 = starting number of prey population;
%   y0 = starting number of prey population;
%   t0 = starting time of integration;
%   T = final time of integration;
%   N = number of intervals for splitting [t0,T]
%
%   OUTPUTS:
%   xn = time series of prey population;
%   yn = time series of predator population;
%   x_eq = x-axis of the stable equilibrium point (if a>0, otherwise 0);
%   y_eq = y-axis of the stable equilibrium point (if a>0, otherwise 0);
%   tn = total number of instants characterizing the time serieses.
%
%

function [xn,yn,x_eq,y_eq,tn] = Lotka_Volterra_E_Espl(a,b,c,d,...
    x0,y0,t0,T,N)

if a>0
    x_eq=c/d; 
    y_eq=a/b;
else
    x_eq=0;
    y_eq=0;
end

h=(T-t0)/N;
tn=linspace(t0,T,N+1);

xn=zeros (1,N+1);
yn=zeros (1,N+1);
xn(1)=x0; yn(1)=y0;

for i =1: N
xn(i+1)= xn(i)+(a-b*yn(i))*xn(i)*h;
yn(i+1)= yn(i)+(-c+d*xn(i))*yn(i)*h;
end

end