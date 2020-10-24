% This code plots the fig.3a in the main text

load('gamma_0.5/DA/DA_sigma_10.mat','S','E','I','U','Q','R','b','r',...
                                   'N','gamma','beta','n_state','state_name',...
                                   'T','I_true','n_tr');
T_predict = 15; % days to predict


ratio_ntr = 1; n_tr = ratio_ntr*n_tr;
ratio_b = 1; b = ratio_b*b;
 

% Number of states to be plotted
n_top = 15;

% specify range of alpha_r
unreport_range = 0:0.1:1;


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
rp = r(:,:,T+1); rp = mean(rp,2);


i = 1;S_end = zeros(n_state,length(unreport_range));
for ratio_r = unreport_range
    
    r_scaled = ones(n_state,1) - (ones(n_state,1)-rp)*ratio_r;
    
    [S_predict,E_predict,~,~,~,~] ...
        = Node_net_multisample_simple(Sp,Ep,Ip,Up,Qp,Rp,bp,r_scaled,beta,T_predict,n_tr,gamma);
    
    S_predict = reshape(S_predict,size(S_predict,1),size(S_predict,3));
    
    S_end(:,i) = S_predict(:,end);
    
    i = i+1;
end



%%%%%%%%%%%%%%%%%%%%
%%% Plot Fig 3a

S0 = S_end(:,end);

% Plot the states with max I_true
[~,Index_top] = maxk(I_true(:,end),n_top);

% compute S/S0
S_S0 = log10(S_end./S0); 
S_S0 = S_S0(Index_top,:);
S_S0 = S_S0(:);

% compute S/N
S_N = log10(S_end./N); 
S_N = S_N(Index_top,:);
S_N = S_N(:);

% find the states that are S0./N < threshold_saturate (used to plot fig3c)
threshold_saturate = 1;
S_N_test = S_end(:,1)./N; S_N_test = S_N_test(Index_top);
Index_S_unimproved = find(S_N_test<threshold_saturate,n_state);
Index_S_unimproved = Index_top(Index_S_unimproved);



% Plot fig3a
r_ratio_plot = repmat(unreport_range,n_top,1); r_ratio_plot = r_ratio_plot(:);

grp = state_name(Index_top); grp = repmat(grp,1,length(unreport_range)); grp = grp(:);

sym = ''; maker_size = [];
for i = 1:n_top
    if mod(i,3) == 1
        sym = strcat(sym,'v'); maker_size = [maker_size,5];
    elseif mod(i,3) == 2
        sym = strcat(sym,'o'); maker_size = [maker_size,5];
    else
        sym = strcat(sym,'s'); maker_size = [maker_size,6];
    end
end
cmap = hsv(n_top);

figure(103)
gca = gscatter(r_ratio_plot,S_S0,grp,cmap,sym,maker_size);
xlabel('$\alpha_r$','Interpreter','latex'); 
ylabel('$\log_{10}\frac{S}{S(\alpha_r = 1)}$','Interpreter','latex');
% ylim([0,1.5]);

for i = 1:length(gca)
  set(gca(i), 'MarkerFaceColor', cmap(i,:));
end




figure(104)
gca = gscatter(r_ratio_plot,S_N,grp,cmap,sym,maker_size);
xlabel('$\alpha_r$','Interpreter','latex'); 
ylabel('$\log_{10}\frac{S}{N}$','Interpreter','latex');
legend('Location','eastoutside');
% ylim([-3,0]);

for i = 1:length(gca)
  set(gca(i), 'MarkerFaceColor', cmap(i,:));
end




% Save the graphs
mkdir(['gamma_',num2str(gamma),'/predict51_fig3_T_',int2str(T_predict)]);

figure(103)
print(['gamma_',num2str(gamma),'/predict51_fig3_T_',int2str(T_predict),'/fig3a_S_S0_top',int2str(n_top),'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');

figure(104)
print(['gamma_',num2str(gamma),'/predict51_fig3_T_',int2str(T_predict),'/fig3a_S_N_top',int2str(n_top),'_sigma',num2str(sigma),'_beta',num2str(beta),'_gamma',num2str(gamma),'.pdf'],'-dpdf');

% Save the index to plot fig3c
Index_S_unimproved = Index_S_unimproved';
save(['gamma_',num2str(gamma),'/predict51_fig3_T_',int2str(T_predict),'/Index_top_',num2str(n_top),'.mat'],'Index_S_unimproved','S0');

