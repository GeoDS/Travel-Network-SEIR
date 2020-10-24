% This code creates the box plots of Re in the data analysis step

load('gamma_1.5/DA/DA_sigma_10.mat');

mkdir(['gamma_',num2str(gamma),'/DA/Box/small']);
mkdir(['gamma_',num2str(gamma),'/DA/Box/median']);
mkdir(['gamma_',num2str(gamma),'/DA/Box/large']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Plot Re

cU = 0.1; % cI = 0.19;
Dc = 2.3; Dl = 6; De = 5.3;
coef_E = gamma*De; coef_A = 1/( cU/Dc + (1-cU)/Dl);

for i_plot = 1:n_state
    
    
    E_sample_plot = reshape(E_total_sample(i_plot,:,:),Nsample,[]);
    U_sample_plot = reshape(U_total_sample(i_plot,:,:),Nsample,[]);
    b_plot = reshape(b(i_plot,:,:),Nsample,[]);
    
    EA = E_sample_plot + U_sample_plot; EA(EA==0) = eps;
    
    Re = b_plot./EA.* ( coef_E*E_sample_plot + coef_A*U_sample_plot );
    
    
    figure(2060)
    boxplot(Re,'symbol',''); ylim([0,Inf]);
    xticks(1:2:20); xticklabels(1:2:20); 
    xlabel('March'); ylabel('R_e');
    title(state_name{i_plot}); pause(0.01);
    
    I_sample_plot = reshape(I_total_sample(i_plot,:,:),Nsample,[]);
    m = median(I_sample_plot(:,end));
    if m < 50
        print(['gamma_',num2str(gamma),'/DA/Box/small/Re_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    elseif m < 200
        print(['gamma_',num2str(gamma),'/DA/Box/median/Re_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    else
        print(['gamma_',num2str(gamma),'/DA/Box/large/Re_',state_name{i_plot},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');
    end

end




