function [Num,I_temp,iter,Phii_cf,etak_opt,PhiPhi,t,Betak]= PilotPC(Phii,K,BETAA,tau,Pp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%       Pilot power control Optimization              %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Phii_cf = Phii; % pilot set of cell-free systems

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose Num best Beta

%%%%%%%%%%%%%%%%%%%%%%%
 [BETAA_temp1,I_temp1] = sort(BETAA,'descend');

%%%%%%%%%%%%%%%%%%%   Option 1: Choose active number of APs   %%%%%%%%%
 Num = 3;

%%%%  Option 2: Choose Num based on the number time decrease of B_mk %%%%%%%

% p_dec = 0.1; %Power decrease
% Ma_num = zeros(K,1);
% for k=1:K
%     iter =1; epsilon = 1;
%     while (epsilon >= p_dec)
%         iter=iter+1;
%         epsilon = BETAA_temp1(iter,k)/BETAA_temp1(1,k);
%     end
%     Ma_num(k,1)=iter;
% end
% Num = max(Ma_num);  

%%%%%% Option 3: Choose Num based on percentage of power %%%%%%%%%%%%%%%%%%%
% per_p = 0.3; %Get Percentage of power (1-per_p)
% Ma_num = zeros(K,1);
% for k=1:K
%     iter =1; epsilon = 1;
%     while (epsilon >= per_p)
%         iter=iter+1;
%         sum_p = 0;
%         for i=1:iter
%             sum_p = sum_p + BETAA_temp1(i,k);
%         end
%         epsilon = (sum(BETAA_temp1(:,k))-sum_p)/sum(BETAA_temp1(:,k));       
%     end
%     Ma_num(k,1)=iter;
% end
% Num = max(Ma_num);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_temp = zeros(Num,K);
for n = 1:Num
    for k = 1:K
        I_temp(n,k) = I_temp1(n,k);
    end
end
%%  Create PhiPhi(KxK) norm of pilot vector.
PhiPhi = zeros(K,K);
for k=1:K
    for i=1:K
        PhiPhi(k,i)=norm(Phii_cf(:,k)'*Phii_cf(:,i));
    end
end
%-------------------------------------------------------------------------
Betak = zeros(Num,K,K);
for k=1:K
    Betak(:,:,k)= BETAA(I_temp(:,k),:); %create beta matrix correspond with Num largest large-scale coefficient; each user k has one Beta matrix
end
Betap = tau*Pp*Betak;
T=zeros(Num,K,K);
for k=1:K
    for m =1:Num
        for i = 1:K
            T(m,i,k) = Betap(m,i,k)*PhiPhi(k,i)^2; %MS of k-th user
        end
    end
end

%% Iteration
% Initial setup
etak_ini = 0.5*ones(K,1);
iter = 0; epsilon = 1; 
ref_etak = zeros(K,1);
while (epsilon >= 5e-3)
    iter = iter + 1;
  
    if (iter == 1)   
       etak_n = etak_ini;
    else    
      etak_n = ref_etak;
    end
 %-----------------------------------------------------------------
 %Initial parameter
a_n = zeros(Num,K);
b_n = zeros(Num,K);
c_n = zeros(Num,K);
for m = 1:Num
    for k = 1:K
        a_n(m,k) = 3*(Betap(m,k,k)*etak_n(k,1))/(T(m,:,k)*etak_n(:,1)+1);
        b_n(m,k) = (Betap(m,k,k)*etak_n(k,1))^2/(T(m,:,k)*etak_n(:,1)+1);
        c_n(m,k) = (Betap(m,k,k)*etak_n(k,1))/(T(m,:,k)*etak_n(:,1)+1)^2;

    end
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 x=zeros(Num,K);
 y=zeros(Num,K);
% ref_etak = zeros(K,1);
 cvx_quiet true
          cvx_begin %sdp'
              variable etak(K,1);
              variable t;
              expression x(Num,K);
              expression y(Num,K);
             for n = 1:Num
                  for k=1:K
                      x(n,k) = inv_pos(Betap(n,k,k)*etak(k,1));
                     y(n,k) = T(n,:,k)*etak(:,1)+1;
                  end
             end
              minimize(t);
              subject to
              
             for k = 1:K
                 t >= Num - sum(a_n(:,k)) + sum(b_n(:,k).*x(:,k)) + sum(c_n(:,k).*y(:,k));  
             end            
              for k = 1:K
                0.01 <= etak(k,1);
                 etak(k,1) <= 1;
              end

           cvx_end
    
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ref_etak = etak;
 ref_t(iter) = t;
    if (iter >= 2)
        epsilon = abs(ref_t(iter) - ref_t(iter-1))/abs(ref_t(iter-1));
    end    

    if (iter >= 15)
        epsilon = 0;
    end
end
etak_opt=etak;



