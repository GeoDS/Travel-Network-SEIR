clear;

% This code creates the box plot of threshold Dq, using maximal
% slope, for all the states

load('gamma_0.5/DA/DA_sigma_10.mat'); 
load('gamma_0.5/predict51_fig2/Index_min_S_end_I_S0_51.mat','Index_min_S_end_I','S0');

Date_i = datetime(2020,3,20);
% time_day = Date_i + caldays(0:T);
T_predict = 40;
time_day = Date_i + caldays(0:T_predict);
% load('data_medical/EKI_T_Itrue_51state_Mar20.mat','I_true','T');


ratio_ntr = 1; n_tr = ratio_ntr*n_tr;
ratio_b = 1; b = ratio_b*b;
 

%%%%%%%%%%%%%%%%%%%%%
%%%   Prediction

Sp = mean(S,2);
Ep = mean(E,2);
Ip = mean(I,2);
Ap = mean(A,2);
Qp = mean(Q,2); 
Rp = mean(R,2);
bp = b(:,:,T+1); bp = mean(bp,2); %bp = reshape(bp,[],1,1);

Dq_range = 0.5:0.05:5.5;

S_end = zeros(n_county,length(Dq_range));
E_end = zeros(n_county,length(Dq_range));
A_end = zeros(n_county,length(Dq_range));
I_end = zeros(n_county,length(Dq_range));
Q_end = zeros(n_county,length(Dq_range));
R_end = zeros(n_county,length(Dq_range));

i = 1;
for Dq = Dq_range
    
    [S_predict,E_predict,I_predict,A_predict,Q_predict,R_predict]...
        = Node_net_multisample_simple_predict(Sp,Ep,Ip,Ap,Qp,Rp,bp,T_predict,n_tr,gamma,Dq);
    
    S_predict = reshape(S_predict,size(S_predict,1),size(S_predict,3));
    E_predict = reshape(E_predict,size(E_predict,1),size(E_predict,3));
    I_predict = reshape(I_predict,size(I_predict,1),size(I_predict,3));
    A_predict = reshape(A_predict,size(A_predict,1),size(A_predict,3));
    Q_predict = reshape(Q_predict,size(Q_predict,1),size(Q_predict,3));
    R_predict = reshape(R_predict,size(R_predict,1),size(R_predict,3));
    
    S_end(:,i) = S_predict(:,end);
    E_end(:,i) = E_predict(:,end);
    I_end(:,i) = I_predict(:,end);
    A_end(:,i) = A_predict(:,end);
    Q_end(:,i) = Q_predict(:,end);
    R_end(:,i) = R_predict(:,end);

    i = i+1;
end



log_S = max(log10(S_end./S0),0);
dlog_S = log_S(:,2:end) - log_S(:,1:end-1);
[~,idx] = min(dlog_S,[],2);

Dq_max_slp = Dq_range(idx);


figure(108)
boxplot(Dq_max_slp); % ylim([0,Inf]);
ylabel('$D_q$','Interpreter','latex');

% print(['gamma_',num2str(gamma),'/predict51_fig2/S','_top',num2str(length(Index_min_S_end_I)),'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');



% save(['gamma_',num2str(gamma),'/predict51_fig2/predict51_new_b_1_Dq_',num2str(length(Index_min_S_end_I)),'.mat']);
