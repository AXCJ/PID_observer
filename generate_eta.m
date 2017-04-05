function eta = generate_eta(A,B,C,D,n,p)

flag=0; k=0;
    while(~flag) 
        eta=(1*randn(p,p));
    %   [eta] = boundchang(eta, 1e0);

        C_tilde = eta*C;
        D_tilde = eta*D;
        Zeros_tilde=sort(TZOCS(A, B, C_tilde, D_tilde));
        rank_obs_=rank(obsv(A,C_tilde));
        realzers=real(Zeros_tilde );
        if(~isempty(realzers))
        if (sum(realzers>=0)==0)&&(rank_obs_==n)
            flag=1;
            save ('eta value', 'eta')

        end
        else
          k=k+1;
        end
    end



end

