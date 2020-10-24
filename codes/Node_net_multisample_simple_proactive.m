function [S_total,E_total,I_total,U_total,Q_total,R_total] ...
                            = Node_net_multisample_simple_proactive(S,E,I,U,Q,R,b,T,n,gamma,Dq)
% This function handle solves the modified ODE system (Proactive Measure) 
% i.e., unreported ratio = 0 and a quarantined compartment is included

% S_total: The time series of susceptible compartment
% E_total: The time series of exposed compartment
% I_total: The time series of reported compartment
% U_total: The time series of unreported compartment
% Q_total: The time series of quarantined compartment
% R_total: The time series of removed compartment
% S: The initial susceptible compartment
% E: The initial exposed compartment
% I: The initial reported compartment
% U: The initial unreported compartment
% Q: The initial quarantined compartment
% R: The initial removed compartment
% b: The transmission rate
% r: The reported ratio
% Dq: The quarantine rate
% T: Terminal time (Day)
% n: transportation matrix
% gamma: relative infectivity -- presymptomatic/symptomatic

n_state = size(n,1); Nsample = size(S,2);

Dc = 2.3; % Dc: the average duration of infection for critical cases. 
Dl = 6; % Dl: the average duration of infection for mild cases.
cI = 0.1; % cI: proportion of critical cases among all reported cases. 
cU = 0.2; % cU: proportion of critical cases among all unreported cases.

dt = 1/24; N = T/dt; % Time step & number of time steps

n_out = sum(n)'; n_out_state = n_out*ones(1,Nsample);

%%%%%%%%%%%%%%%%%
% initialization
S_total = zeros(n_state,Nsample,N+1); E_total = zeros(n_state,Nsample,N+1);
I_total = zeros(n_state,Nsample,N+1); U_total = zeros(n_state,Nsample,N+1);
Q_total = zeros(n_state,Nsample,N+1); R_total = zeros(n_state,Nsample,N+1);

S_total(:,:,1) = S; E_total(:,:,1) = E;
Q_total(:,:,1) = Q; U_total(:,:,1) = U;
I_total(:,:,1) = I; R_total(:,:,1) = R;

%%%%%%%%%%%%%%
% evolution
for step = 1:N
   
    P = S + E + U; P(P==0) = eps;
    
    dSdt = (b.*(gamma*E+U).*S)./P;      
    dIdt = E/Dq;
    
    
    S_n = S - dt*dSdt - dt*n_out_state.*(S./P) + dt*n*(S./P);
    
    E_n = E + dt*dSdt - dt*dIdt - dt*n_out_state.*(E./P) + dt*n*(E./P);
    
    I_n = I - cI*dt*I/Dc - (1-cI)*dt*I/Dl;
    
    A_n = U - cU*dt*U/Dc - (1-cU)*dt*U/Dl ...
        - dt*n_out_state.*(U./P) + dt*n*(U./P);

    R_n = R + cI*dt*I/Dc + (1-cI)*dt*I/Dl + cU*dt*U/Dc + (1-cU)*dt*U/Dl;
    Q_n = Q + dt*dIdt;

    % avoid non-physical quanities
    l = find(S_n<0,n_state*Nsample);
    if ~isempty(l)
        S_n(l) = 0;
    end
    k = find(E_n<0,n_state*Nsample);
    if ~isempty(k)
        E_n(k) = 0;
    end
    m = find(A_n<0,n_state*Nsample);
    if ~isempty(m)
        A_n(m) = 0;
    end
    
    %%%%%%%%%%%%%%
    % update the quanities
    S = S_n; U = A_n; E = E_n; Q = Q_n; I = I_n; R = R_n;
    
    %%%%%%%%%%%%%%
    % record the history
    S_total(:,:,step+1) = S; E_total(:,:,step+1) = E; 
    Q_total(:,:,step+1) = Q; U_total(:,:,step+1) = U;
    I_total(:,:,step+1) = I; R_total(:,:,step+1) = R;
   
    
end


end
