% This code plots the predicted Re - time curve (without proactive measure)

load('gamma_0.5/DA/DA_sigma_10.mat','S','E','I','U','Q','R','b','r','n_tr',...
                                      'n_state','Nsample','gamma','beta','sigma')
                                  
Date_i = datetime(2020,3,20); % initial date for prediction
T_predict = 15; % days to predict
time_day = Date_i + caldays(0:T_predict); % array for dates

dt_inv = 24; % 1/dt in the ODE solver (used to extract data for each day)

N_plot = 5; % Number of states to plot (the top N_plot states with largest I_true intially)




% Plot for each triple (alpha_r,alpha_b,alpha_tr)
ratio_ntr = 1;

for ratio_r = [0,0.3,0.6,1]
for ratio_b = [0.1,0.3,0.6,1]


 
n_tr_scaled = ratio_ntr*n_tr;
r_scale = ones(n_state,Nsample) - (ones(n_state,Nsample)-r)*ratio_r;
b_scale = ratio_b*b;

%%%%%%%%%%%%%%%%%%%%%
%%%   Prediction

% use the mean in data analysis step as intial data
Sp = mean(S,2);
Ep = mean(E,2);
Ip = mean(I,2);
Up = mean(U,2);
Qp = mean(Q,2); 
Rp = mean(R,2); 
bp = b_scale(:,:,T+1); bp = mean(bp,2); 
rp = r_scale(:,:,T+1); rp = mean(rp,2); 

[S_predict_sample,E_predict_sample,I_predict_sample,U_predict_sample,...
                          ~,R_predict_sample]...
    = Node_net_multisample_simple(Sp,Ep,Ip,Up,Qp,Rp,bp,rp,beta,T_predict,n_tr_scaled,gamma);

% extra data for each day
S_predict_sample = S_predict_sample(:,:,1:dt_inv:end);
E_predict_sample = E_predict_sample(:,:,1:dt_inv:end);
I_predict_sample = I_predict_sample(:,:,1:dt_inv:end);
U_predict_sample = U_predict_sample(:,:,1:dt_inv:end);
R_predict_sample = R_predict_sample(:,:,1:dt_inv:end);

%%% Plot Re
I_predict = reshape(I_predict_sample,size(I_predict_sample,1),size(I_predict_sample,3));
U_predict = reshape(U_predict_sample,size(U_predict_sample,1),size(U_predict_sample,3));
E_predict = reshape(E_predict_sample,size(E_predict_sample,1),size(E_predict_sample,3));
S_predict = reshape(S_predict_sample,size(S_predict_sample,1),size(S_predict_sample,3));
R_predict = reshape(R_predict_sample,size(R_predict_sample,1),size(R_predict_sample,3));

cA = 0.1; % cI = 0.19;
Dc = 2.3; Dl = 6; De = 5.3;
coef_E = gamma*De; coef_A = 1/( cA/Dc + (1-cA)/Dl);

[~,n_plot] = maxk(I_true(:,end),N_plot);

for i_plot = n_plot 
    
    EA = E_predict(i_plot,:) + U_predict(i_plot,:); EA(EA==0) = eps;
    
    Re = bp(i_plot)./EA .*( coef_E*E_predict(i_plot,:) + coef_A*U_predict(i_plot,:) );
    
    
    figure(2061)
    plot(time_day,Re,'-x','LineWidth',1); hold on;
    ylabel('R_e');
    
    pause(0.01);
    
end
legend(state_name(n_plot));
hold off;

print(['gamma_',num2str(gamma),'/predict51_b_',num2str(ratio_b),'_r_',num2str(ratio_r),'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/Re_top',int2str(N_plot),'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');



end
end
