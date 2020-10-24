clear;

% This code generates *.csv files given the prediction data

T_predict = 15;

gamma = 0.5; 
ratio_ntr = 1;

for b_plot = [0.1,0.3,0.6,1]
    for alpha_r_plot = [0,0.3,0.6,1]

        b_str = num2str(b_plot);
        r_str = num2str(alpha_r_plot);
        
        load(['gamma_',num2str(gamma),'/predict51_b_',b_str,'_r_',r_str,...
            '_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),...
            '/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),...
            '_T_',int2str(T_predict),'.mat'],'S_predict','U_predict','I_predict','E_predict','R_predict',...
            'state_name','time_day');
        
        time_day = datestr(time_day); time_day = cellstr(time_day);
        
        UIER_predict = U_predict + I_predict + E_predict + R_predict;
        
        
        m_U = array2table(U_predict);
        m_S = array2table(S_predict);
        m_E = array2table(E_predict);
        m_I = array2table(I_predict);
        m_R = array2table(R_predict);
        m_UIER = array2table(UIER_predict);
        
        
        
        m_U.Properties.RowNames = state_name;
        m_U.Properties.VariableNames = time_day;
        
        writetable(m_U,['gamma_',num2str(gamma),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'_A.csv'],'WriteRowNames',true);
        
        
        m_E.Properties.RowNames = state_name;
        m_E.Properties.VariableNames = time_day;
        
        writetable(m_E,['gamma_',num2str(gamma),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'_E.csv'],'WriteRowNames',true);
        
        
        m_S.Properties.RowNames = state_name;
        m_S.Properties.VariableNames = time_day;
        
        writetable(m_S,['gamma_',num2str(gamma),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'_S.csv'],'WriteRowNames',true);
        
        
        m_I.Properties.RowNames = state_name;
        m_I.Properties.VariableNames = time_day;
        
        writetable(m_I,['gamma_',num2str(gamma),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'_I.csv'],'WriteRowNames',true);
        
        
        m_R.Properties.RowNames = state_name;
        m_R.Properties.VariableNames = time_day;
        
        writetable(m_R,['gamma_',num2str(gamma),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_R.csv'],'WriteRowNames',true);
        
        
        m_UIER.Properties.RowNames = state_name;
        m_UIER.Properties.VariableNames = time_day;
        
        writetable(m_UIER,['gamma_',num2str(gamma),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'/predict51_b_',b_str,'_r_',r_str,'_travel_',num2str(ratio_ntr),'_T_',int2str(T_predict),'_AIER.csv'],'WriteRowNames',true);
        
        

    end
end
