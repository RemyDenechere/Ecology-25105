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







