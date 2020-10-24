% The main program for data analysis

load('data_traffic/top_state51.mat','N','n_tr','state_name'); % load total population N; traffic matrix n_tr; names of states
alpha_tr = 0.5;   % this factor is used to modify the transport matrix and accouts for the travel decline during the beginning period of epidemic compared to the data set
n_tr = alpha_tr*n_tr;
n_state = size(N,1);

load('data_medical/EKI_T_Itrue_51state_Mar1_Mar20.mat','I_true','T');

beta = 0;     % No quarantined compartment
gamma = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial distribution %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_range = [0.1,0.3];
b_range = [1,1.5];
alpha_range = [0,500];
U_range = [0,200];

Nsample = 2000; sigma = 10;

rng('shuffle');

r = rand(n_state,Nsample); 
r = (r_range(2)-r_range(1))*r + r_range(1);
r = cat(3,r,zeros(n_state,Nsample,T));

b = rand(n_state,Nsample); 
b = (b_range(2)-b_range(1))*b + b_range(1);
b = cat(3,b,zeros(n_state,Nsample,T));

alpha = rand(n_state,Nsample);
alpha = (alpha_range(2)-alpha_range(1))*alpha + alpha_range(1);

U0 = rand(n_state,Nsample);
U0 = (U_range(2)-U_range(1))*U0 + U_range(1);



%%%  Initialize S-E-I-U-R
S_total_sample = zeros(n_state,Nsample,T+1);
E_total_sample = zeros(n_state,Nsample,T+1);
I_total_sample = zeros(n_state,Nsample,T+1);
U_total_sample = zeros(n_state,Nsample,T+1);
R_total_sample = zeros(n_state,Nsample,T+1);


I = repmat(I_true(:,1),1,Nsample);
E = alpha;
Q = zeros(size(E)); 
R = zeros(size(E));
U = U0;
S = N*ones(1,Nsample) - E - I - U;

S_total_sample(:,:,1) = S;
E_total_sample(:,:,1) = E;
I_total_sample(:,:,1) = I;
U_total_sample(:,:,1) = U;
R_total_sample(:,:,1) = R;

%%% begin analyzing for each day sucessively
for t = 1:T
        
        % Run the ODE for 1 day 
        [S_total,E_total,I_total,U_total,~,R_total] ...
                            = Node_net_multisample_simple(S,E,I,U,Q,R,b(:,:,t),r(:,:,t),beta,1,n_tr,gamma);
        Se = S_total(:,:,end); Se = reshape(Se,n_state,Nsample);
        Ee = E_total(:,:,end); Ee = reshape(Ee,n_state,Nsample);
        Ie = I_total(:,:,end); Ie = reshape(Ie,n_state,Nsample);
        Ue = U_total(:,:,end); Ue = reshape(Ue,n_state,Nsample);
        Re = R_total(:,:,end); Re = reshape(Re,n_state,Nsample);
        
        
        % Ensemble Kalman Filter
        uE = [Se;Ee;Ie;Ue;Re;b(:,:,t);r(:,:,t)];
        uE_mean = mean(uE,2); uE_mean2 = uE_mean*ones(1,Nsample);
        
        MuE = Ie; 
        MuE_mean = mean(MuE,2); MuE_mean2 = MuE_mean*ones(1,Nsample);
        
        Cup = (uE-uE_mean2)*(MuE-MuE_mean2)'/Nsample;
        Cpp = (MuE-MuE_mean2)*(MuE-MuE_mean2)'/Nsample;
        
        d = I_true(:,t+1);
        
        uA = uE + Cup*( (Cpp+sigma^2*eye(n_state))\(d*ones(1,Nsample) + sigma*rand(n_state,Nsample) - MuE) );

        
        % Allocate S/E/I/U/R/b/r and avoid the non-physical quanities
        S = uA(1:n_state,:); 
        S(S<0) = 0; 
        S = min(S,S_total_sample(:,:,t)); % May replace it with S = min(S,N*ones(1,Nsample)); 
        
        E = uA(n_state+1:2*n_state,:); 
        E(E<0) = 0;
        
        I = uA(2*n_state+1:3*n_state,:); 
        I(I<0) = 0; 
        I = max(I,I_total_sample(:,:,t)); % May delete this line
        
        U = uA(3*n_state+1:4*n_state,:); 
        U(U<0) = 0;
        
        R = uA(4*n_state+1:5*n_state,:); 
        R(R<0) = 0; 
        R = max(R,R_total_sample(:,:,t)); % May delete this line
        
        b(:,:,t+1) = uA(5*n_state+1:6*n_state,:); 
        b(b<0) = 0;
        
        r(:,:,t+1) = uA(6*n_state+1:end,:); 
        r(r>1) = 1; 
        r(r<0.01) = 0.01; % May replace it with r(r<0) = eps;
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Record the history
        
        S_total_sample(:,:,t+1) = S;
        E_total_sample(:,:,t+1) = E;
        I_total_sample(:,:,t+1) = I;
        U_total_sample(:,:,t+1) = U;
        R_total_sample(:,:,t+1) = R;
        
end

mkdir(['gamma_',num2str(gamma),'/DA']);
save(['gamma_',num2str(gamma),'/DA/DA_sigma_',num2str(sigma),'.mat']);
