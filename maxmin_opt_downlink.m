function[R_cf_opt_user] = maxmin_opt_downlink(M,L1,K,BETAA,tau,Pp,Pd,Phii_cf,etak_opt,I_temp)
%% Calculate UL Rate base on optimal power control of pilot sequences
etak = etak_opt;
%-------------------------------------------------------------------------
%% Create new Beta matrix
Beta_matrix = zeros(M,K);
for k=1:K
    Beta_matrix(I_temp(:,k),k) = BETAA(I_temp(:,k),k); % All inactive channels = 0
end
    
 %% Create Gamma matrix (eq.8)
Gammaa = zeros(M,K); % No changing Index
Gammaa_nom = zeros(M,K); % Gamma/Beta
mau=zeros(M,K);
for m=1:M
    for k=1:K
        mau(m,k)=norm(((Beta_matrix(m,:).*etak').^(1/2)).*(Phii_cf(:,k)'*Phii_cf))^2;
    end
end

for m=1:M
    for k=1:K
        Gammaa(m,k)=tau*Pp*etak(k,1)*Beta_matrix(m,k)^2/(tau*Pp*mau(m,k) + 1);
        Gammaa_nom(m,k)=tau*Pp*etak(k,1)*Beta_matrix(m,k)/(tau*Pp*mau(m,k) + 1);
    end
end

%% Compute etaa(m): (each AP transmits equal power to K terminals)
etaa=zeros(M,1);
for m=1:M
    if sum(Gammaa(m,:))==0
        etaa(m)=0;
    else
        etaa(m)=1/(L1*sum(Gammaa(m,:)));
    end
end
%% Compute Rate
SINR=zeros(1,K);
Ratestep=zeros(1,K);

%Pilot contamination
PC = zeros(K,K);
for ii=1:K
    for k=1:K
        PC(ii,k)=sum((etaa.^(1/2)).*((Gammaa_nom(:,ii)).*Beta_matrix(:,k)))*Phii_cf(:,ii)'*Phii_cf(:,k);
    end
end
PC1=(abs(PC)).^2;

for k=1:K
    num=0;
    for m=1:M
        num=num + (etaa(m)^(1/2))*Gammaa(m,k);
    end
    SINR(k) = L1^2*Pd*num^2/(1 + Pd*sum(Beta_matrix(:,k)) + L1^2*Pd*sum(PC1(:,k)) - L1^2*Pd*PC1(k,k));
    %Rate:
    Ratestep(k) = log2(1+ SINR(k));
end      
R_cf_min=min(Ratestep);

%% 2) Max-Min power allocations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmin=2^R_cf_min-1;
tmax=2^(2*R_cf_min+2.5)-1; %case when L1=5
epsi=max(tmin/5,0.01);

PhiPhi=zeros(K,K);
for ii=1:K
    for k=1:K
        PhiPhi(ii,k)=norm(Phii_cf(:,ii)'*Phii_cf(:,k));
    end
end
BETAAn=Beta_matrix*Pd;
Gammaan=Gammaa*Pd;
Xopt=ones(M,K);
yopt=sqrt(Pd)*ones(M,1);
Zopt=ones(K,K);
%cvx_solver sedumi
count = 0;
cvx_quiet true
            while( tmax - tmin > epsi)
            count = count + 1;
            tnext = (tmax+tmin)/2; 
           cvx_begin sdp
              variables X(M,K) y(M,1) Z(K,K)
              minimize(0)
              subject to
                for k=1:K
                     norm([sqrt(tnext)*[Z((1:(k-1)),k);Z(((k+1):K),k)].*[PhiPhi((1:(k-1)),k);PhiPhi(((k+1):K),k)];...  
                         sqrt(tnext/L1)*(BETAAn(:,k).^(1/2)).*y ; sqrt(tnext*Pd)/L1]) <= (Gammaan(:,k))'*X(:,k) ;
                end
                for m=1:M
                    norm(((Gammaan(m,:)).^(1/2)).*X(m,:)) <= y(m); 
                    y(m)<=(sqrt(Pd)/L1);
                end
                
                for k=1:K
                    for ii=1:K
                        sum( ((Gammaa_nom(:,ii)).*Beta_matrix(:,k)).*X(:,ii)  ) <= Z(ii,k);
                    end
                end
                
                for m=1:M
                 for k=1:K
                    X(m,k)>=0;
                 end
                end
               
            cvx_end


            % bisection
            if strfind(cvx_status,'Solved') % feasible
            fprintf(1,'Problem is feasible ',tnext);
            tmin = tnext;
            Xopt=X;
            yopt=y;
            Zopt=Z;
            else % not feasible
            fprintf(1,'Problem not feasible ',tnext);
            tmax = tnext;   
            end

            end
            R_cf_opt_user = zeros(1,K);
            for k=1:K
                muk=((Gammaan(:,k))'*Xopt(:,k))^2/( norm([[Zopt((1:(k-1)),k);Zopt(((k+1):K),k)].*[PhiPhi((1:(k-1)),k);PhiPhi(((k+1):K),k)];...
                    1/sqrt(L1)*(BETAAn(:,k).^(1/2)).*yopt ; sqrt(Pd)/L1])^2 );
                R_cf_opt_user(k) = log2(1+muk);
            
            end   