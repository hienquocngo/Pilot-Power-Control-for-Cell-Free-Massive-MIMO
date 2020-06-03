function[R_cf_opt_user] = maxmin_opt_uplink(L1,K,tau,Pp,Pu,Phii_cf,etak_opt,PhiPhi,Num,Betak)
%% Calculate UL Rate base on optimal power control of pilot sequences
etak = etak_opt;
%-------------------------------------------------------------------------
 %% Create Gamma matrix (eq.8)
Gammaa = zeros(Num,K);
mau=zeros(Num,K);
for m=1:Num
    for k=1:K
        mau(m,k)=norm(((Betak(m,:,k).*etak').^(1/2)).*(Phii_cf(:,k)'*Phii_cf))^2;
    end
end

for m=1:Num
    for k=1:K
        Gammaa(m,k)=tau*Pp*etak(k,1)*Betak(m,k,k)^2/(tau*Pp*mau(m,k) + 1);
    end
end


%% %% Compute UL Rate (consider data PL power control = 1) (eq.28_trang's paper) 

SINR=zeros(1,K);


%Pilot contamination
PC = zeros(K,K);
for ii=1:K
    for k=1:K
        PC(ii,k) = sum((Gammaa(:,k)./Betak(:,k,k)).*Betak(:,ii,k))*((etak(ii)/etak(k))^(1/2))*Phii_cf(:,k)'*Phii_cf(:,ii); %the first part of denominator
    end
end
PC1=(abs(PC)).^2;

Rate = zeros(1,K);

for k=1:K
    deno1=0;
    for m=1:Num
        deno1=deno1 + Gammaa(m,k)*sum(Betak(m,:,k));
    end
    SINR(1,k) = L1^2*Pu*(sum(Gammaa(:,k)))^2/(L1*sum(Gammaa(:,k)) + L1*Pu*deno1 + L1^2*Pu*sum(PC1(:,k)) - L1^2*Pu*PC1(k,k));
    %Rate: 
    Rate(1,k) = log2(1+ SINR(1,k));
end
R_cf_min = min(Rate);

%% 2) Max-Min power allocations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmin=2^R_cf_min-1;
tmax=2^(2*R_cf_min+2)-1; %case when L1=5
epsi=max(tmin/5,0.01);

% BETAAn=BETAA*Pu;
Te1 =zeros(K,K);
Te2 =zeros(K,K);

x_cf_opt=ones(K,1);


for ii=1:K
    for k=1:K
        Te1(ii,k)=L1*Pu*sum(Betak(:,ii,k).*Gammaa(:,k)); 
        Te2(ii,k) = L1^2*Pu*(sum((Gammaa(:,k)./Betak(:,k,k)).*Betak(:,ii,k))*((etak(ii)/etak(k))^(1/2)))^2*PhiPhi(k,ii)^2; 
    end
end

%cvx_solver sedumi
cvx_quiet true
    while( tmax - tmin > epsi)

    tnext = (tmax+tmin)/2; 
   cvx_begin %sdp
      variables x(K,1); 
      minimize(0)
      subject to
        for k=1:K
            Te1(:,k)'*x(:,1) + [Te2(1:(k-1),k); Te2((k+1):K,k) ]'*[x(1:(k-1),1); x((k+1):K,1)] +...
                L1*sum(Gammaa(:,k)) <= (1/tnext)*L1^2*Pu*(sum(Gammaa(:,k)))^2*x(k,1) ;
        end              
        for k=1:K
            x(k,1)>=0;
            x(k,1)<=1;
        end

    cvx_end


            % bisection
            if strfind(cvx_status,'Solved') % feasible
            fprintf(1,'Problem is feasible ',tnext);
            tmin = tnext; % tmin =: tnext
            x_cf_opt=x;  % x_cf_opt =: x
            else % not feasible
            fprintf(1,'Problem not feasible ',tnext);
            tmax = tnext;   % tmax =: tnext
            end

     end
  R_cf_opt_user = zeros(K,1);
  mu1k = zeros(K,1);
            for k=1:K
                mu1k(k,1) = L1^2*Pu*(sum(Gammaa(:,k)))^2*x_cf_opt(k,1)/( Te1(:,k)'*x_cf_opt(:,1) +...
                    [Te2(1:(k-1),k); Te2((k+1):K,k) ]'*[x_cf_opt(1:(k-1),1); x_cf_opt((k+1):K,1)] + L1*sum(Gammaa(:,k)) );
                R_cf_opt_user(k,1) = log2(1+mu1k(k,1));
            end             