%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                       Ecology 16-Feb-20                               %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exercise 1: Logistic growth (Verhulst model)
% Find a locla maximum: dN/dt =0 + d2N/dndt <0: 
% here a polynomial function order 2 with -aN2 then there is one maximum we
% only need to solve dN/dt = 0

% solving using symbolic language:  
syms r N K
dNdt = r*N*(1-N/K);
Sol = solve(dNdt == 0, N)

%  going a bit further: ---------------------------------------------------
% plot the diff eq: 
N = linspace(0,100,100); %ind
param.r = 1.5; % Ind per day
param.K = 50; % ind
dNdt = param.r*N.*(1-N/param.K);

subplot(1,2,1)
plot(N, dNdt)
hold on 
plot([0, 100], [0, 0], 'k--', [50, 50], [-160, 0])
ylabel('dN/dt')
xlabel('N')

% solve the diff equation: 
tspan = [0, 10];
param.r = [0.5, 1 1.5]'; % 3 different growth rate
N0 = [4, 4, 4]; 
[t,y] = ode45(@(t, N) solve_pop(t, N, param), tspan, N0); % See definition 

subplot(1,2,2)
plot(t, y)
leg = legend(num2str(param.r))
title(leg, 'r')
ylabel('N')
xlabel('t (days)')

%% Exercice2: Life Table
delta = 0.1; % d-1
kappa = 5; %d
m = 25; % eggs
x = linspace(0, 100, 1000); % age axis 
dx = x(2:end) - x(1:end-1); % dx (variation of age)
mx = x; mx(x<kappa) = 0; mx(x>=kappa) = m;

% Survival rate:
subplot(2,2,1)
lx = exp(-delta*x); 
plot(x, lx)
xlabel('Age (day)')
ylabel('Survival prob')

% egg production: 
subplot(2,2,2) 
plot(x, mx)
xlabel('Age (day)')
ylabel('Egg production')

% potential egg production: 
subplot(2,2,3)
plot(x, lx.*mx)
xlabel('Age (day)')
ylabel('Potential egg production')

% solving R0
R0 = lx(1:end-1).*mx(1:end-1) * dx' % # of offspring (Integration as a product of matrix )
% solving T: 
T = (x(1:end-1).*lx(1:end-1).*mx(1:end-1) * dx')/R0 % day
r = log(R0)/T % log(egg)/day

%% Exercice 3: Lotka?Volterra
% 2. Solve equations
syms r N c P a m
dNdt = r*N - c*N*P
dPdt = a*c*N*P - m*P

Neq = solve(dPdt == 0, N)
Peq = solve(dNdt == 0, P)

% diametre: 1 micron = 10-6 m = 10-4cm
% radius = 0.5*10-4cm
volume = (4/3)*pi*(5*10^(-4))^3; % cm^3 <=> mL
clear = 24 * volume * 10^5 % cm^3/hour <=> mL/hour
PredgrowthEff = 2/3000;

sym2poly(subs(Neq, [a, c, m], [PredgrowthEff, clear, 1]))
sym2poly(subs(Peq, [r, c], [1, clear]))

[t,y] = ode45(@(t,y) solve_ex3(t, y), [1, 100], [4,4])

subplot(1,2,1)
semilogy(t, y)
legend('Prey', 'Predator')
ylabel('Concentration (#/V)')
xlabel('time (day)')
title('Concentration over time')

subplot(1,2,2)
plot( y(:,2), y(:,1))
ylabel('Prey (#/V)')
xlabel('Predator (#/V)')
title('Prey as a function of Pred')

%% function definitions
function dNdt = solve_pop(t, N0, param)
% Initial condition: 
N = N0;

% ODE: 
dNdt = param.r.*N.*(1-N/param.K);
end

function y = solve_ex3(t, y0)
% parameters: 
c = 0.0013; 
r = 1; 
m = 0.3; 
a = 2/3000; 

% Initial condition 
N = y0(1);
P = y0(2); 

% ODE: 
dNdt = r*N - c*N*P;
dPdt = a*c*N*P - m*P;

y = [dNdt, dPdt]'
end





