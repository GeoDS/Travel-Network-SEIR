% This code plots the heat map of the predicted Re at day 
% T_predict as a function of transmission rate b and reported rate r. 
% No proactive measure is taken.

load('gamma_0.5/DA/DA_sigma_10.mat','S','E','I','U','Q','R','b','r','n_tr',...
                                      'n_state','Nsample','gamma','beta','sigma')
Date_i = datetime(2020,3,20); % initial date for prediction
T_predict = 15; % days to predict
time_day = Date_i + caldays(0:T_predict); % array for dates

dt_inv = 24; % 1/dt in the ODE solver (used to extract data for each day)

% parameters used to compute the Re
cA = 0.2; % cI = 0.1;
Dc = 2.3; Dl = 6; De = 5.3;
coef_E = gamma*De; coef_A = 1/( cA/Dc + (1-cA)/Dl);


N_plot = 5;    % Number of states to be plotted
[~,n_plot] = maxk(I_true(:,end),N_plot); % Index for the states to be plotted


%%%%%%%%%%%%%%%%%%
%%% Compute Re
ratio_ntr = 1;
ratio_r_range = 0:0.05:1;
ratio_b_range = 0.1:0.05:1;

Re = zeros(length(ratio_r_range),length(ratio_b_range),N_plot);

for i = 1:length(ratio_b_range)
    for j = 1:length(ratio_r_range)
        
        ratio_b = ratio_b_range(i);
        ratio_r = ratio_r_range(j);
        
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
        
        
        [~,E_predict_sample,~,U_predict_sample,~,~]...
            = Node_net_multisample_simple(Sp,Ep,Ip,Up,Qp,Rp,bp,rp,beta,T_predict,n_tr_scaled,gamma);
        
        % extra data for each day
        E_predict_sample = E_predict_sample(:,:,1:dt_inv:end);
        U_predict_sample = U_predict_sample(:,:,1:dt_inv:end);
        
        
        U_predict = reshape(U_predict_sample,size(U_predict_sample,1),size(U_predict_sample,3));
        E_predict = reshape(E_predict_sample,size(E_predict_sample,1),size(E_predict_sample,3));
        
        
        %%%%%%%%%%%%
        %%% compute Re
        EA = E_predict(n_plot,end) + U_predict(n_plot,end); EA(EA==0) = eps;
        
        Re_temp = bp(n_plot)./EA .*( coef_E*E_predict(n_plot,end) + coef_A*U_predict(n_plot,end) );
        
        Re_temp = reshape(Re_temp,1,1,[]);
        
        Re(j,i,:) = Re_temp;
        
        
    end

    
end



%%%%%%%%%%%%%%%%%%%%%55
%%% Plot the heat map

[B,R] = meshgrid(ratio_b_range,ratio_r_range);

for i_plot = 1:N_plot   
    
    figure(2065+i_plot)
    hold on; 
    
    pcolor(B,R,Re(:,:,i_plot)); shading flat; % plot the heat map of Re
    contour(B,R,Re(:,:,i_plot),1,'LineWidth',2,'LineColor','r'); % plot the Re = 1 level set
    
    title(state_name{n_plot(i_plot)}); 
    xlabel('\alpha_b'); ylabel('\alpha_r'); colorbar;
    hold off;
    
    print(['gamma_',num2str(gamma),'/Re_end_',state_name{n_plot(i_plot)},'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'_T_',int2str(T_predict),'.pdf'],'-dpdf');
    
    pause(0.01);
    
end


