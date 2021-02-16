%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                       Ecology 16-Feb-20                               %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exercise 1: Logistic growth (Verhulst model)
% Find a locla maximum: d(dN/dt)/dN = 0 
% here dN/dt is a polynomial function of order 2 with -aN2 then there is one maximum
% exemple with a k = 50
clf
figure
% solving using symbolic language:  
syms r N K
G = r*N*(1-N/K); % logistic growth 
dGdN = r - 2*r*N/K % derivative of G
Sol = solve(dGdN == 0, N) % max of G

%  going a bit further: ---------------------------------------------------
% solve the diff equation: 
tspan = [0, 10];
param.r = [0.5, 1 1.5]'; % 3 different growth rate
N0 = [4, 4, 4]; 
[t,y] = ode45(@(t, N) solve_pop(t, N, param), tspan, N0); % See definition 

subplot(2,2,1)
plot(t, y)
leg = legend(num2str(param.r))
title(leg, 'r')
ylabel('N population size')
xlabel('t (days)')
title('Population dynamic')

% plot the diff eq: 
N = linspace(0,100,100); % Ind
param.r = 1.5; % Ind per day
param.K = 50; % ind
dNdt = param.r*N.*(1-N/param.K);

subplot(2,2,2)
plot(N, dNdt)
hold on 
plot([0, 100], [0, 0], 'k--', [50, 50], [-160, 0])
ylabel('G = dN/dt')
xlabel('N population size')
title('Fisheries yield G')

% plot the derivative of G:
dGdN = param.r - 2*param.r*N/param.K;
subplot(2,2,3)
semilogx(N, dGdN)
hold on 
plot([1, 100], [0, 0], 'k--')
hold on 
plot([25, 25], [-10, 0])
ylabel('dG/dN = d(dN/dt)/dN')
xlabel('N population size')
title('Derivative of Fisheries yield dG/dN')

%% Exercice2: Life Table
figure
delta = 0.1; % d-1
kappa = 5; %d
m = 25; % eggs
x = linspace(0, 100, 1000); % age axis 
dx = x(2:end) - x(1:end-1); % dx (variation of age)
mx = x; mx(x<kappa) = 0; mx(x>=kappa) = m;

% Survival rate:
subplot(2,2,1)
lx = exp(-delta*x); 
plot(x, lx, 'Linewidth', 1.5)
xlabel('Age (day)')
ylabel('Survival prob')
title('l(x) = exp(- \delta x) ')

% egg production: 
subplot(2,2,2) 
plot(x, mx, 'Linewidth', 1.5)
xlabel('Age (day)')
ylabel('Egg production')
title('Fecondity')

% potential egg production: 
subplot(2,2,3)
plot(x, lx.*mx, 'Linewidth', 1.5)
xlabel('Age (day)')
ylabel('R_0')
title('Reproductive output, R_0 = \int{ l_x m_x dx}')

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

[t,y] = ode45(@(t,y) solve_ex3(t, y), [1, 100], [4,4]);

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





