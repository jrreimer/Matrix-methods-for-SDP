%%% matrix version of Chan and Godfray's (1993) resource pool model
%%% vary beta to look at sensitivity of threshold to this value
%%% in their paper, they looked at beta = 4, and beta = 16; we can look 
%%% at a much wider range for this threshold. This code is nearly identical
%%% to that found in Parasitic_wasp_example.m; here we are just varying
%%% parameter values and plotting the stationary decision. 

close all
clear

% SDP parameters
N = 75;  % number of states considered 
eta = 0.2;
mu = 0.0125;
alpha = 30;
betavals = 1:20; % range of beta values considered
Pi_thresholds = zeros(N,length(betavals));

for beta = betavals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SENSITIVITY TO BETA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define all corresponding P and G matrices
Pi = ones(N,1); % initial policy guess; all ones (all parasitize)
P = zeros(N,N);
G = zeros(N,1);
for xval = 1:N % Set up P
    if xval ~= 1
        P(xval,xval-1) = (1-eta)*exp(-mu);
    end
    G(xval,1) = eta; 
    if xval-1-beta > 0
        P(xval,xval-1-beta) = eta*exp(-mu);
    end     
end

valmax = sum((inv(eye(N) - P))*G);
pimax = Pi;

for pi = 1:N % change pi through each xval
    
    Pi = pimax;
    Pi(pi,1) = 2;
    
    % 1 means parsitize; 2 means host feed
    P = zeros(N,N);
    G = zeros(N,1);

    for xval = 1:N % for each value of x
        if xval ~= 1
            P(xval,xval-1) = (1-eta)*exp(-mu);
        end
        if Pi(xval) == 1 % parasitize
            G(xval,1) = eta; % should this be here? Or inside the if statement?
            if xval-1-beta > 0
                P(xval,xval-1-beta) = eta*exp(-mu);
            end
        else % host feed
            if xval-1+alpha <= N
                P(xval, xval-1+alpha) = eta*exp(-mu);
            else
                P(xval, N) = eta*exp(-mu);
            end
        end
    end
    vals = sum((inv(eye(N) - P))*G);
    if vals > valmax
        valmax = vals;
        pimax = Pi;
    end
end

Pi_thresholds(:,beta) = pimax;
end

figure
subplot(1,2,2)
colormap([0 0.408 0.545; 0.8 0.8 0.8]); 
set(gcf, 'Position', [100, 100, 550, 350])
imagesc(Pi_thresholds)
xticklabels = 0:5:betavals(length(betavals));
xticks = 0:5:length(betavals);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xlabel('\beta')
set(gca,'fontsize',14) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SENSITIVITY TO ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SDP parameters
N = 75;  % number of states considered 
eta = 0.2;
mu = 0.0125;
alphavals = 20:40;
beta = 10; % range of beta values considered
Pi_thresholds = zeros(N,length(alphavals));

for alpha = alphavals
% define all corresponding P and G matrices
Pi = ones(N,1); % initial policy guess; all ones (all parasitize)
P = zeros(N,N);
G = zeros(N,1);
for xval = 1:N % Set up P
    if xval ~= 1
        P(xval,xval-1) = (1-eta)*exp(-mu);
    end
    G(xval,1) = eta; 
    if xval-1-beta > 0
        P(xval,xval-1-beta) = eta*exp(-mu);
    end
     
end

valmax = sum((inv(eye(N) - P))*G);
pimax = Pi;

for pi = 1:N % change pi through each xval
    
    Pi = pimax;
    Pi(pi,1) = 2;
    
    % 1 means parsitize; 2 means host feed
    P = zeros(N,N);
    G = zeros(N,1);

    for xval = 1:N % for each value of x
        if xval ~= 1
            P(xval,xval-1) = (1-eta)*exp(-mu);
        end
        if Pi(xval) == 1 % parasitize
            G(xval,1) = eta; % should this be here? Or inside the if statement?
            if xval-1-beta > 0
                P(xval,xval-1-beta) = eta*exp(-mu);
            end
        else % host feed
            if xval-1+alpha <= N
                P(xval, xval-1+alpha) = eta*exp(-mu);
            else
                P(xval, N) = eta*exp(-mu);
            end
        end
    end
    vals = sum((inv(eye(N) - P))*G);
    if vals > valmax
        valmax = vals;
        pimax = Pi;
    end
end

Pi_thresholds(:,alpha-alphavals(1)+1) = pimax;
end

subplot(1,2,1)
colormap([0 0.408 0.545; 0.8 0.8 0.8]); 
imagesc(Pi_thresholds)
xticklabels = alphavals(1):5:alphavals(length(alphavals));
xticks = 1:5:length(alphavals);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
ylabel('state')
xlabel('\alpha')
set(gca,'fontsize',14) 