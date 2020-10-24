% This code plots the predicted (S,E,I,U,R) curve (without proactive
% measure)


load('gamma_0.5/DA/DA_sigma_10.mat','S','E','I','U','Q','R','b','r','n_tr',...
                                      'n_state','Nsample','gamma','beta','sigma');

Date_i = datetime(2020,3,20); % initial date for prediction
T_predict = 15; % days to predict
time_day = Date_i + caldays(0:T_predict); % array for dates

dt_inv = 24; % 1/dt in the ODE solver (used to extract data for each day)


% Plot for each triple (alpha_r,alpha_b,alpha_tr)

ratio_ntr = 1;

for ratio_r = [0,0.3,0.6,1]
for ratio_b = [0.1,0.3,0.6,1]



n_tr = ratio_ntr*n_tr;
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
    = Node_net_multisample_simple(Sp,Ep,Ip,Up,Qp,Rp,bp,rp,beta,T_predict,n_tr,gamma);


% extra data for each day
S_predict_sample = S_predict_sample(:,:,1:dt_inv:end);
E_predict_sample = E_predict_sample(:,:,1:dt_inv:end);
I_predict_sample = I_predict_sample(:,:,1:dt_inv:end);
U_predict_sample = U_predict_sample(:,:,1:dt_inv:end);
R_predict_sample = R_predict_sample(:,:,1:dt_inv:end);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Plot (S,E,I,U,R) 
I_predict = reshape(I_predict_sample,size(I_predict_sample,1),size(I_predict_sample,3));
U_predict = reshape(U_predict_sample,size(U_predict_sample,1),size(U_predict_sample,3));
E_predict = reshape(E_predict_sample,size(E_predict_sample,1),size(E_predict_sample,3));
S_predict = reshape(S_predict_sample,size(S_predict_sample,1),size(S_predict_sample,3));
R_predict = reshape(R_predict_sample,size(R_predict_sample,1),size(R_predict_sample,3));


mkdir(['gamma_',num2str(gamma),'/predict51_b_',num2str(ratio_b),'_r_',num2str(ratio_r),'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/log_small']);
mkdir(['gamma_',num2str(gamma),'/predict51_b_',num2str(ratio_b),'_r_',num2str(ratio_r),'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/log_median']);
mkdir(['gamma_',num2str(gamma),'/predict51_b_',num2str(ratio_b),'_r_',num2str(ratio_r),'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/log_large']);

for i_plot = 1:n_state   
    figure(102)
    subplot(3,2,1); 
    plot(time_day,max(log10(S_predict(i_plot,:)),0),'-x'); ylabel('log_{10}S'); title('Susceptible(S)');
    subplot(3,2,2); 
    plot(time_day,max(log10(E_predict(i_plot,:)),0),'-x'); ylabel('log_{10}E'); title('Latent(E)');
    subplot(3,2,3); 
    plot(time_day,max(log10(U_predict(i_plot,:)),0),'-x'); ylabel('log_{10}U'); title('Unreported(U)');
    subplot(3,2,4); 
    plot(time_day,max(log10(I_predict(i_plot,:)),0),'-x'); ylabel('log_{10}I'); title('Reported(I)');
    subplot(3,2,5); 
    plot(time_day,max(log10(R_predict(i_plot,:)),0),'-x'); ylabel('log_{10}R'); title('Resolved(R)');
    
    sgtitle(state_name{i_plot});

    if I_predict(i_plot,1)<50
        print(['gamma_',num2str(gamma),'/predict51_b_',num2str(ratio_b),'_r_',num2str(ratio_r),'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/log_small/',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    elseif I_predict(i_plot,1)<200
        print(['gamma_',num2str(gamma),'/predict51_b_',num2str(ratio_b),'_r_',num2str(ratio_r),'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/log_median/',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    else
        print(['gamma_',num2str(gamma),'/predict51_b_',num2str(ratio_b),'_r_',num2str(ratio_r),'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/log_large/',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');    
    end

    hold off;
    pause(0.01);
    
end

save(['gamma_',num2str(gamma),'/predict51_b_',num2str(ratio_b),'_r_',num2str(ratio_r),'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/predict51_b_',num2str(ratio_b),'_r_',num2str(ratio_r),'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'.mat']);


end
end
