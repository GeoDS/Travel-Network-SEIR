% This code creates the box plot for each state in the data
% analysis step (Box plot for (S,E,U,R) and I and b and r)

load('gamma_0.5/DA/DA_sigma_10.mat');
load('data_medical/EKI_T_Itrue_51state_Mar1_Mar20.mat','I_true','T');

% the 'small/median/large' folders are for states with small/median/large 
% infected population
mkdir(['gamma_',num2str(gamma),'/DA/Box/small']);
mkdir(['gamma_',num2str(gamma),'/DA/Box/median']);
mkdir(['gamma_',num2str(gamma),'/DA/Box/large']);



% The following four parts are independent, and can be commented out
% separately.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Plot (S,E,I,U) 

for i_plot = 1:n_state
    
    
    S_sample_plot = reshape(S_total_sample(i_plot,:,:),Nsample,[]);
    E_sample_plot = reshape(E_total_sample(i_plot,:,:),Nsample,[]);
    U_sample_plot = reshape(U_total_sample(i_plot,:,:),Nsample,[]);
    I_sample_plot = reshape(I_total_sample(i_plot,:,:),Nsample,[]);
    R_sample_plot = reshape(R_total_sample(i_plot,:,:),Nsample,[]);
    
    figure(1)
    subplot(2,2,1); 
    boxplot(S_sample_plot,'symbol',''); 
    xticks(1:2:20); xticklabels(1:2:20);
    xlabel('March');
    title('Susceptible(S)');
    
    subplot(2,2,2); 
    boxplot(E_sample_plot,'symbol',''); ylim([0,Inf]);
    xticks(1:2:20); xticklabels(1:2:20);
    xlabel('March');
    title('Latent(E)');
    
    subplot(2,2,3); 
    boxplot(U_sample_plot,'symbol',''); ylim([0,Inf]);
    xticks(1:2:20); xticklabels(1:2:20);
    xlabel('March');
    title('Unreported(U)');
    
    subplot(2,2,4); 
    boxplot(R_sample_plot,'symbol',''); ylim([0,Inf]);
    xticks(1:2:20); xticklabels(1:2:20);
    xlabel('March');
    title('Resolved(R)');
    
    sgtitle(state_name{i_plot});
    pause(0.01);
    
    m = median(I_sample_plot(:,end));
    if m < 50
        print(['gamma_',num2str(gamma),'/DA/Box/small/SEUR',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    elseif m < 200
        print(['gamma_',num2str(gamma),'/DA/Box/median/SEUR',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    else
        print(['gamma_',num2str(gamma),'/DA/Box/large/SEUR_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    end

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Plot reported cases
for i_plot = 1:n_state
    
    I_sample_plot = reshape(I_total_sample(i_plot,:,:),Nsample,[]);
    
    figure(100)
    boxplot(max(log10(I_sample_plot),0),'symbol',''); 
    hold on;
    plot(1:20,max(log10(I_true(i_plot,:)),0),'x','LineWidth',1); 
    hold off;
    xticks(1:2:20); xticklabels(1:2:20);
    xlabel('March');
    ylim([0,Inf]);
    
    
    title(state_name{i_plot});
    pause(0.01);
    ylabel('No. of reported cases (log_{10}I)');
    
    m = median(I_sample_plot(:,end));
    if m < 50
        print(['gamma_',num2str(gamma),'/DA/Box/small/reported_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    elseif m < 200
        print(['gamma_',num2str(gamma),'/DA/Box/median/reported_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    else
        print(['gamma_',num2str(gamma),'/DA/Box/large/reported_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    end


end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot transmission rate b
for i_plot = 1:n_state
    
    b_plot = reshape(b(i_plot,:,:),Nsample,[]);

    figure(2000)
    boxplot(b_plot,'symbol',''); ylim([0,2]);
    xticks(1:2:20); xticklabels(1:2:20); xlabel('March');
    ylabel('b');
    title(state_name{i_plot}); pause(0.01);
    
    I_sample_plot = reshape(I_total_sample(i_plot,:,:),Nsample,[]);
    m = median(I_sample_plot(:,end));
    if m < 50
        print(['gamma_',num2str(gamma),'/DA/Box/small/b_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    elseif m < 200
        print(['gamma_',num2str(gamma),'/DA/Box/median/b_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    else
        print(['gamma_',num2str(gamma),'/DA/Box/large/b_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot report ratio r
for i_plot = 1:n_state
    
    r_plot = reshape(r(i_plot,:,:),Nsample,[]);

    figure(5000)
    boxplot(r_plot,'symbol',''); ylim([0,0.5]);
    xticks(1:2:20); xticklabels(1:2:20); xlabel('March');
    ylabel('r');
    title(state_name{i_plot}); pause(0.01);
    
    I_sample_plot = reshape(I_total_sample(i_plot,:,:),Nsample,[]);
    m = median(I_sample_plot(:,end));
    if m < 50
        print(['gamma_',num2str(gamma),'/DA/Box/small/r_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    elseif m < 200
        print(['gamma_',num2str(gamma),'/DA/Box/median/r_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    else
        print(['gamma_',num2str(gamma),'/DA/Box/large/r_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    end

end


