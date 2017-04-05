function De = generate_De(A,C_tilde,n,p)
flag=0; k=0;
    while(~flag)
        De = randn(n,p);
        zero_m_e = tzero(A, De, C_tilde, zeros(p));
        realzers = real(zero_m_e);
        imagzers = imag(zero_m_e);
        if((imagzers==0))
            if(~isempty(realzers))
                if(sum(realzers>=0)==0)
                    flag = 1;
                    save('De_conference','De')
                    break
                end
            else
                  k=k+1;
            end
        end
    end
end

