%%% illustrative minimum example; this is a 2 patch model with 5 possible
%%% states; following the structure of the first canonical equation for 
%%% patch choice

close all
clear

% define SDP parameters
T = 20; % terminal time
f = zeros(5,T); % matrix of fitness values
d = zeros(5,T-1); % matrix of optimal decisions at each time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RowOpts = zeros(2, 5, 5); % storage of all possible rows; need one 2x5 
% matrix for each of the 5 states

% for oscillating decisions, use these parameters
p1 = 0.3;
q1 = 0.99; 

% for no oscillations, use these parameters
% p1 = 0.4;
% q1 = 0.99;

p_vals = [p1 p1*2]; % probability of finding prey in each patch
q_vals = [q1 q1*0.9]; % probability of surviving predation in each patch
% note: you are twice as likely to find prey in patch 2, but 10% less likely
% to survive predation (1st patch is low risk, low reward; 2nd patch is
% higher risk, higher reward)

% First, define matrix where it is only decision 1 or decision 2,
% regardless of state
for c1 = 1:2 % 2 possible actions for each state; patch 1 or 2
    p = p_vals(c1); % get probabiliy of finding prey for a given patch
    q = q_vals(c1); % get probability of surviving in a given patch
    
    % define matrix P 
    P = [0 0 q*p 0 0;
    q*(1-p) 0 0 q*p 0;
    0 q*(1-p) 0 0 q*p;
    0 0 q*(1-p) 0 q*p;
    0 0 0 q*(1-p) q*p];
    
    for c2 = 1:5 % for each of 5 possible states
        RowOpts(c1,:,c2) = P(c2,:); % save rows into RowOpts to be called later
    end
end

% Now, construct all possible P matrices from the component rows in RowOpts
P = zeros(5,5,2^5); % create zero matrix
patches = zeros(5,2^5); % create zero patches placeholder for each policy

% go through all permutation of 1s and 2s for each possible policy
% save policy in matrix "patches" and corresponding matrix P
ct = 1;
for c1 = 1:2
    for c2 = 1:2
        for c3 = 1:2
            for c4 = 1:2
                for c5 = 1:2
                    patches(:,ct) = [c1 c2 c3 c4 c5]';
                    P(1,:,ct) = RowOpts(c1,:,1);
                    P(2,:,ct) = RowOpts(c2,:,2);
                    P(3,:,ct) = RowOpts(c3,:,3);
                    P(4,:,ct) = RowOpts(c4,:,4);
                    P(5,:,ct) = RowOpts(c5,:,5);
                    ct = ct+1; 
                end
            end
        end
    end    
end

% calculate dominant eigenvalue of P for each policy
eigs = zeros(1,2^5); 
for ct = 1:2^5
    eigs(ct) = max(abs(eig(P(:,:,ct))));
end

stationary = find(eigs == max(eigs)); % find largest dominant eigenvalue 
% then P(:,:,stationary) is the matrix with biggest e-value

% now get corresponding right and left eigenvectors
[V,D,W] = eig(P(:,:,stationary)); % V is right eigenvectors, W is left
D = diag(D);
[out,ord] = sort(abs(D),'descend');
D = D(ord); lambda1 = D(1); lambda2 = abs(D(2)); % eigenvalues
rho = lambda1/lambda2; % damping ratio
V = V(:,ord); V = V./sum(V,1); V1 = V(:,1); % right eigenvectors; stable structure is V1
W = W(:,ord); W = W./W(1,:); W = W(:,1); % left eigenvectors




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now solve using backwards induction (on matrices) %%%%%%%%%%%%%%%%%%%%%%%

% set terminal condition
f(:,T) = ones(5,1);
 
% set up SDP model
for t = T-1:-1:1
    V = []; % value functions at each step    
    for j = 1:2^5
        V(j) = sum(P(:,:,j)*f(:,t+1)); % summing matrix multiplication of P*f(t+1) 
    end   
    allSame = find(V == max(V)); % find maximum value(s); call it allSame
    d(t) = allSame(randi(length(allSame))); % if more than one maximum 
    % value, choose randomly (this will usually not be an issue)
    f(:,t) = P(:,:,d(t))*f(:,t+1); % update fitness function according to optimal decision
    decision(:,t) = patches(:,d(t)); % save policy which was optimal at time t
end

%disp(decision) % display all optimal decisions
%disp(eig(P(:,:,stationary))) % display eigenvalues of stationary decision matrix
%disp(stationary) % matrix number of stationary matrix 
%disp(eigs) % display dominant eigenvalue of each matrix (for interest)
disp(patches(:,stationary)) % display stationary policy
dampingRaiot = D(1)/abs(D(2)); % calculate the damping ratio
disp(dampingRaiot) % display the damping ratio


%%%%%%%%%%%%%%%% Plot fitness and decisions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf, 'Position', [10, 50, 550, 750])
subplot(3,1,1)
plot(f','k','LineWidth',2) % following optimal decisions
xlabel('time')
ylabel('F(t)')
yticks(0.2:0.4:1.0)
ylim([0.1,1])
yticklabels({'0.2'; '0.6'; '1.0'})
set(gca,'fontsize',14) 
xticks(5:5:T-5)
xlim([1,T])

%figure
subplot(3,1,2)
asympt = ones(5,30);
for ct = 1:5
    asympt(ct,:) = V1(ct)*asympt(ct,:);
end
plot(asympt','--','Color',[0.7 0.7 0.7],'LineWidth',2) % following rule of thumb (i.e., the asymptotic) decisions
hold on
plot((f./sum(f,1))','k','LineWidth',2) % following optimal decisions
xlabel('time')
ylabel('normalized F(t)')
yticks(0.10:.1:0.3)
set(gca,'fontsize',14) 
xticks(5:5:T-5)
xlim([1,T])
ylim([0.1,0.3])

subplot(3,1,3)
colormap([0 0.408 0.545; 0.8 0.8 0.8]); 
imagesc(decision)
xlabel('time')
ylabel('state')
stateNames = {'x_1'; 'x_2'; 'x_3'; 'x_4'; 'x_5'};
set(gca, 'YTick', 1:1:5, 'YTickLabel', stateNames) 
set(gca,'fontsize',14) 
xticks(0:5:T-1)

