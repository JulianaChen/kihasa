%% simulate husband wages

% Drawing Types

for n=1:1:G.n_pop
        seed(n)=rand;
        if seed(n)<0.16
            type(n,1)=1;
            abi(n,1)=1;
            edu(n,1)=1;
        elseif seed(n)<0.189 && seed(n)>=0.16
            type(n,1)=2;
            abi(n,1)=1;
            edu(n,1)=2;
        elseif seed(n)<0.215 && seed(n)>=0.189
            type(n,1)=3;
            abi(n,1)=1;
            edu(n,1)=3;
        elseif seed(n)<0.579 && seed(n)>=0.215
            type(n,1)=4;
            abi(n,1)=2;
            edu(n,1)=1;
        elseif seed(n)<0.73 && seed(n)>=0.579
            type(n,1)=5;
            abi(n,1)=2;
            edu(n,1)=2;
        elseif seed(n)<=1 && seed(n)>=0.73
            type(n,1)=6;
            abi(n,1)=2;
            edu(n,1)=3;
        end
end

for t=1:G.n_period-1
    
    age = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + t;

    for n=1:G.n_pop
    
    % Calculate husband wages

    wh_mean(n,t) = eta01 + eta02*(abi(n)==2) + eta11*(edu(n)==2) + eta12*(edu(n)==3) + eta2*age(t); %function of women's age, education and ability types
    wh_sd(n,t) = eta03 + eta04*(abi(n)==2) + eta21*(edu(n)==2) + eta22*(edu(n)==3) + eta3*age(t); %same as above
    wh(n,t) = normrnd(wh_mean(n,t), wh_sd(n,t));

    end 
    
end 

for t=1:G.n_period-1
    
    age = 18*(edu==1) + 20*(edu==2) + 22*(edu==3) + t;

    for n=1:G.n_pop
    
    % Calculate child investment level 

    Inv_mean(n,t) = iota01 + iota02*(abi(n)==2) + iota11*(edu(n)==2) + iota12*(edu(n)==3) + iota2*age(t);
    Inv_sd(n,t) = iota03 + iota04*(abi(n)==2) + iota21*(edu(n)==2) + iota22*(edu(n)==3) + iota3*age(t);
    Inv(n,t) = normrnd(Inv_mean(n,t),Inv_sd(n,t)); 

    end 
    
end 