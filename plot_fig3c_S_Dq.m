% This code plots the fig.3c (predicted S as a function of Dq in 
% the proactive measure) in the main text

gamma = 0.5; n_top = 15;
T_predict = 15;

load(['gamma_',num2str(gamma),'/DA/DA_sigma_10.mat'],'S','E','I','U','Q','R','b','r',...
                                                   'N','T','n_tr','state_name','n_state','sigma','beta'); 
load(['gamma_',num2str(gamma),'/predict51_fig3_T_',int2str(T_predict),'/Index_top_',num2str(n_top),'.mat'],'Index_S_unimproved','S0');

ratio_ntr = 1; n_tr = ratio_ntr*n_tr;
ratio_b = 1; b = ratio_b*b;

% specify range of Dq
Dq_range = 0.5:0.1:5.5;



%%%%%%%%%%%%%%%%%%%%%
%%%   Prediction

% use the mean in data analysis step as intial data

Sp = mean(S,2);
Ep = mean(E,2);
Ip = mean(I,2);
Up = mean(U,2);
Qp = mean(Q,2); 
Rp = mean(R,2);
bp = b(:,:,T+1); bp = mean(bp,2);

S_end = zeros(n_state,length(Dq_range));
E_end = zeros(n_state,length(Dq_range));
U_end = zeros(n_state,length(Dq_range));
I_end = zeros(n_state,length(Dq_range));
Q_end = zeros(n_state,length(Dq_range));
R_end = zeros(n_state,length(Dq_range));

i = 1;
for Dq = Dq_range
    
    [S_predict,E_predict,I_predict,U_predict,Q_predict,R_predict]...
        = Node_net_multisample_simple_proactive(Sp,Ep,Ip,Up,Qp,Rp,bp,T_predict,n_tr,gamma,Dq);
    
    S_predict = reshape(S_predict,size(S_predict,1),size(S_predict,3));
    E_predict = reshape(E_predict,size(E_predict,1),size(E_predict,3));
    I_predict = reshape(I_predict,size(I_predict,1),size(I_predict,3));
    U_predict = reshape(U_predict,size(U_predict,1),size(U_predict,3));
    Q_predict = reshape(Q_predict,size(Q_predict,1),size(Q_predict,3));
    R_predict = reshape(R_predict,size(R_predict,1),size(R_predict,3));
   
    S_end(:,i) = S_predict(:,end);
    E_end(:,i) = E_predict(:,end);
    I_end(:,i) = I_predict(:,end);
    U_end(:,i) = U_predict(:,end);
    Q_end(:,i) = Q_predict(:,end);
    R_end(:,i) = R_predict(:,end);

    i = i+1;
end




%%%%%%%%%%%%%%%%%%%
%%% Plot Fig 3c
i = 1; 
for i_plot = Index_S_unimproved
    
    if mod(i,3) == 1
        marker = '-v'; 
    elseif mod(i,3) == 2
        marker = '-o';
    else
        marker = '-s';
    end
    
    figure(105)
    gca = plot(Dq_range,log10(S_end(i_plot,:)/S0(i_plot)),marker,'LineWidth',1); 
    hold on;
    ylabel('$\log_{10}\frac{S}{S_0(\alpha_r = 1)}$','Interpreter','latex'); 
    xlim([Dq_range(1),Dq_range(end)]);
    xlabel('$D_q$','Interpreter','latex'); 
    color = get(gca,'Color');
    set(gca,'MarkerFaceColor',color,'MarkerSize',5);
    
    
    figure(106)
    gca = plot(Dq_range,log10(S_end(i_plot,:)/N(i_plot)),marker,'LineWidth',1); hold on;
    ylabel('$\log_{10}\frac{S}{N}$','Interpreter','latex'); 
    xlim([Dq_range(1),Dq_range(end)]);
    xlabel('$D_q$','Interpreter','latex'); 
    set(gca,'MarkerFaceColor',color,'MarkerSize',5);
    
    i = i+1;
    pause(0.01);
    
end

figure(105)
hold off; legend(state_name(Index_S_unimproved),'Location','northeast');
figure(106)
hold off; legend(state_name(Index_S_unimproved),'Location','southwest');



% Save the graphs
mkdir(['gamma_',num2str(gamma),'/predict51_fig3_T_',int2str(T_predict)]);
figure(105)
print(['gamma_',num2str(gamma),'/predict51_fig3_T_',int2str(T_predict),'/fig3c_S_S0_top',int2str(n_top),'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
figure(106)
print(['gamma_',num2str(gamma),'/predict51_fig3_T_',int2str(T_predict),'/fig3c_S_N_top',int2str(n_top),'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
