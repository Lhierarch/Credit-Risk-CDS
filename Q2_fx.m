function s=Q2_fx(St)
    global C;
    global t;
    global S;
    % global P;
    
       S(t+1)=St;
    rf=0.3;
    %R=0.4;
    R=0.4;
    Premiumleg=0;
    Protectionleg=0;
   
    for i=1:t
        % Premiumleg=Premiumleg+C(t)*P(t)*S(i)+(C(t)/2)*P(t-0.5)*(S(i-1)-S(i));
        % Premiumleg=Premiumleg+C(t)*exp(-rf*i)*S(i)+...
        % (C(t)/2)*exp(-rf*(i-0.5))*(S(i-1)-S(i));
        Premiumleg=Premiumleg+C(t)*exp(-rf*i)*S(i+1)+...
            (C(t)/2)*exp(-rf*(i-0.5))*(S(i)-S(i+1));
    
    end  
    
    for i=1:t
        Protectionleg=Protectionleg+(1-R)*exp(-rf*(i-0.5))*(S(i)-S(i+1));
    end

    s=(Protectionleg-Premiumleg)^2;

end



