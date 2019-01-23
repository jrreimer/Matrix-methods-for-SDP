%%% matrix version of Chan and Godfray's (1993) resource pool model
close all
clear

tic

%%%%% set SDP parameters %%%%%
T = 1000; % terminal time
N = 75;  % number of states considered (i.e., x_max = N)
eta = 0.2; % probability of encountering a host
mu = 0.0125; % instantaneous mortality rate
alpha = 30; % energy gained from host feeding
beta = 4; % energy required to mature an egg 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MATRIX METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pimax = ones(N,1); % initial policy guess; all 1's (i.e., all parasitize)
Pmax = zeros(N,N); % empty P matrix
Gmax = zeros(N,1); % G column vector for initial policy guess

for xval = 1:N % Set up initial P and G (for policy of all 1s)
    if xval ~= 1
        Pmax(xval,xval-1) = (1-eta)*exp(-mu);
    end
    if xval-1-beta > 0
        Gmax(xval,1) = eta;
        Pmax(xval,xval-1-beta) = eta*exp(-mu);
    end    
end

Fstar = sum((inv(eye(N) - Pmax))*Gmax); % calculate F* = (I-P)^(-1) * G and then sum over components

for xval = 1:N % change pi through each xval successively, and see if it improves F*
    
    Pi = pimax;
    P = Pmax;
    G = Gmax;
    Pi(xval,1) = 2;
    
    % 1 means parsitize; 2 means host feed
    P(xval,:) = zeros(1,N);
    G(xval,1) = 0;

    if xval ~= 1 % don't find a host, so x -> x-1
        P(xval,xval-1) = (1-eta)*exp(-mu);
    end 
    if xval-1+alpha <= N % find a host and host feed; x -> x-1+alpha
        P(xval, xval-1+alpha) = eta*exp(-mu);
    else % if x-1+alpha > N, it x maxes out at N
        P(xval, N) = eta*exp(-mu);
    end
        
    FstarNew = sum((inv(eye(N) - P))*G); % calculate new F* = (I-P)^(-1) * G
    if FstarNew > Fstar % if the new F* is bigger than the previous max, make this the new baseline policy
        Fstar = FstarNew;
        pimax = Pi;
        Pmax = P;
        Gmax = G;
    end
end
toc

disp('predicted policy is:')
disp(pimax')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% BACKWARDS INDUCTION (with matrices) %%%%%%%%%%%%%%%%%%%%%
tic
f = zeros(N,T); % fitness values for all time steps and states
d = zeros(N,T-1); % optimal decisions for all time steps and states
% terimnal condition is zero for all states, so don't need to define it.

for t = T-1:-1:1 % work backwards in time
    pimax2 = ones(N,1); % initial policy guess; all ones (all parasitize)
    Pmax = zeros(N,N); % initialize optimal P matrix
    Gmax = zeros(N,1); % initialize optimal G vector
    for xval = 1:N % Set up initial P 
        if xval ~= 1 % if don't find host, x -> x-1
            Pmax(xval,xval-1) = (1-eta)*exp(-mu);
        end
        if xval-1-beta > 0 % if find host and parasitize, x -> x-1-beta
            Gmax(xval,1) = eta;
            Pmax(xval,xval-1-beta) = eta*exp(-mu);
        end
    end
    Ftemp = Gmax + Pmax*f(:,t+1); % calculate f(:,t) for policy of all 1s

    for xval = 1:N
        Pi = pimax2;
        Pi(xval,1) = 2; % 1 means parsitize; 2 means host feed
        P = Pmax;
        G = Gmax;
        P(xval,:) = zeros(1,N);
        G(xval,1) = 0; 
        if xval ~= 1 % if don't find host, x -> x-1
            P(xval,xval-1) = (1-eta)*exp(-mu);
        end
        if xval-1+alpha <= N % if find host and host feed, x -> x-1+alpha
            P(xval, xval-1+alpha) = eta*exp(-mu);
        else % if x-1+alpha > N, then x -> N
            P(xval, N) = eta*exp(-mu);
        end
        
        % check if new F* is bigger than current Ftemp, if so, keep this
        % policy
        if sum(G + P*f(:,t+1)) > sum(Ftemp)
            Ftemp = G + P*f(:,t+1);
            pimax2 = Pi; 
            Gmax = G;
            Pmax = P;
        end       
        f(:,t) = Ftemp;
        d(:,t) = pimax2;   
    end
end % for t = T-1:-1:1

% note: pimax2 should match pimax, calculated above using the analytic
% method; if it doesn't, display the following message:
if (sum(pimax - pimax2 ~= 0))
    disp('not as stationary decisions yet')
end
toc

% make figure of optimal decisions at each time step, using backwards
% induction (for beta = 4; below, do the same thing for beta = 16)
figure
colormap([0 0.408 0.545; 0.8 0.8 0.8]); 
set(gcf, 'Position', [10, 50, 550, 650])
subplot(2,1,1)
imagesc(d)
ylabel('state')
set(gca,'fontsize',14) 


%%%%%%%%%%%%%%%% now for beta = 16 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the exact same as the above backwwards induction code for beta =4
tic 
beta = 16;
f = zeros(N,T);
d = zeros(N,T-1);
% terimnal condition is zeros, so don't need to define it.

for t = T-1:-1:1
    pimax2 = ones(N,1); % initial policy guess; all ones (all parasitize)
    Pmax = zeros(N,N);
    Gmax = zeros(N,1);
    for xval = 1:N % Set up initial P 
        if xval ~= 1
            Pmax(xval,xval-1) = (1-eta)*exp(-mu);
        end
        if xval-1-beta > 0
            Gmax(xval,1) = eta;
            Pmax(xval,xval-1-beta) = eta*exp(-mu);
        end
    end
    Ftemp = Gmax + Pmax*f(:,t+1);

    for xval = 1:N
        Pi = pimax2;
        Pi(xval,1) = 2; % 1 means parsitize; 2 means host feed
        P = Pmax;
        G = Gmax;
        P(xval,:) = zeros(1,N);
        G(xval,1) = 0; 
        if xval ~= 1
            P(xval,xval-1) = (1-eta)*exp(-mu);
        end
        if xval-1+alpha <= N
            P(xval, xval-1+alpha) = eta*exp(-mu);
        else
            P(xval, N) = eta*exp(-mu);
        end
                 
        if sum(G + P*f(:,t+1)) > sum(Ftemp)
            Ftemp = G + P*f(:,t+1);
            pimax2 = Pi; 
            Gmax = G;
            Pmax = P;
        end       
        f(:,t) = Ftemp;
        d(:,t) = pimax2;   
    end
end % for t = T-1:-1:1
% note: pimax2 should match pimax, calculated above using the analytic
% method
if (sum(pimax - pimax2 ~= 0))
    disp('not as stationary decisions yet')
end
toc
colormap([0 0.408 0.545; 0.8 0.8 0.8]); 
set(gcf, 'Position', [10, 50, 550, 650])
subplot(2,1,2)
imagesc(d)
ylabel('state')
xlabel('time')
set(gca,'fontsize',14) 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Markov chains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tmax = 15; % how long to run the Markov chains for
X = zeros(N,Tmax); % matrix to save probability of being in each state at time t
x0 = poisspdf(1:N,floor(N/2)); % initial condition: poisson distribution with mean N/2
X(:,1) = x0; % start with initial condition in first column of X
P = zeros(N,N,Tmax); % matrices corresponding to optimal decisions at each time

for c = 1:Tmax
    Pi = d(:,c); % get optimal policy at each time step c
    
    % define P(pi) and G(pi)
    for xval = 1:N % for each value of x
        if xval ~= 1 % if don't find host, x -> x-1
            P(xval,xval-1,c) = (1-eta)*exp(-mu);
        end 
        if Pi(xval) == 1 % if parasitize
            if xval-1-beta > 0
                G(xval,1) = eta; 
                P(xval,xval-1-beta,c) = eta*exp(-mu);
            end
        else % if host feed
            if xval-1+alpha <= N
                P(xval, xval-1+alpha,c) = eta*exp(-mu);
            else
                P(xval, N,c) = eta*exp(-mu);
            end
        end
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    X(:,c+1) = transpose(P(:,:,c))*X(:,c); % from Markov chain equation in paper
end
% hack: if we want to include the zero state, we can take "1 - column sum" of all non-zero
% states
zro = 1-sum(X);
Xfull = [zro; X];

figure
subplot(3,1,2)
colormap(flipud(gray))
imagesc(Xfull)
set(gca,'YDir','normal')
ylabel('state')
set(gca,'fontsize',14) 
colorbar('east','Ticks',0.03:0.05:0.17);
ylim([0,75])
xlim([1,Tmax])
yticks([0 20 40 60])

% if want to only show probabilities, conditional on being alive,
% normalize the column vectors of x
Xnorm = X./sum(X);
subplot(3,1,3)
colormap(flipud(gray))
imagesc(Xnorm)
set(gca,'YDir','normal')
xlabel('time')
ylabel('state')
set(gca,'fontsize',14) 
colorbar('east','Ticks',0.01:0.02:0.07);

ylim([1,75])
xlim([1,Tmax])
yticks([0 20 40 60])


%%%%%%%%%% Monte Carlo simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsim = 20; % number of Monte Carlo simulations to perform
X_MC = zeros(nsim,Tmax); % store x values; each row is one simulation
for c = 1:nsim % do the following for each simulation
    X_MC(c,1) = poissrnd(floor(N/2)); % get random initial condition drawn from poisson distribution
    for t = 2:Tmax % for each time step until Tmax
        dx = d(X_MC(c,t-1),t-1); % get optimal decision for realized state
        if rand < eta % encounter host with probability eta
            if rand < exp(-mu) % survive with probability exp(-mu)
                if dx == 1 % parasitize, if optimal
                    X_MC(c,t) = max(X_MC(c,t-1)-1-beta,0);
                else % host feed, if optimal
                    X_MC(c,t) = min(X_MC(c,t-1)-1+alpha,N);
                end
            else % don't survive
                break
            end
        else % don't encounter host 
            X_MC(c,t) = max(X_MC(c,t-1) - 1,0);
        end
        if X_MC(c,t) == 0 % if dead, stop simulation
            break
        end
    end % end t = 1:Tmax
end % end simulation loop

subplot(3,1,1)
plot(X_MC.','k')
set(gca,'fontsize',14) 
ylabel('state')

ylim([0,75])
xlim([1,Tmax])









