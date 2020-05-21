function [Pi_i_taus] = prof_taus(q_dstars)
    global eta gamma a b c_upn c_o2 f c_o   
    P_n_tau = b - eta*(q_dstars(1,1) + q_dstars(2,1)) - gamma*(q_dstars(3,1) + q_dstars(4,1))
    P_o_tau = a - eta*(q_dstars(3,1) + q_dstars(4,1)) - gamma*(q_dstars(1,1) + q_dstars(2,1))

    digitsold=digits(10);   
    
    Pi_i_tau = (P_n_tau - c_upn)*q_dstars(1,1) + (P_o_tau - c_o2)*q_dstars(3,1) 
    
    P_n_tau1 = b - eta*(q_dstars(5,1) + q_dstars(6,1)) - gamma*(q_dstars(7,1) + q_dstars(8,1))
    P_o_tau1 = a - eta*(q_dstars(7,1) + q_dstars(8,1)) - gamma*(q_dstars(5,1) + q_dstars(6,1))   
    
    Pi_i_tau1 = (P_n_tau1 - c_upn*f(q_dstars(1,1)))*q_dstars(5,1) + (P_o_tau1 - c_o)*q_dstars(7,1) 
    
    Pi_i_taus = Pi_i_tau + Pi_i_tau1
end