clc
clear all
close all

%Consider a rectangular area with DxD m^2
%M distributed APs serves K terminals, they all randomly located in the
%area
K=40; %number of terminals
M=100; %number of access points 
L1=3; %Number of antennas per AP
B=20; %Mhz
D=1; %in kilometer
tau=5; % uplink training duration (in samples)
[U,S,V]=svd(randn(tau,tau));%U(tau x tau) is unitary matrix includes tau orthogonal sequences(both columns and rows) 

Hb = 15; % Base station height in m
Hm = 1.65; % Mobile height in m
f = 1900; % Frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8); 
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;      

power_f=0.1; 
noise_p = 10^((-203.975+10*log10(20*10^6)+9)/10); %noise power
Pu = power_f/noise_p;%nomalized receive SNR
Pp = Pu;%pilot power
Pd=2*Pu; %downlink power: 200 mW

d0=0.01;%km
d1=0.05;%km

N=1; %Number of iterations

Rate_greedy = zeros(N,K);
R_1 = zeros(N,K);
% number_iter = zeros(N,1);
 for i=1:N
%Randomly locations of M APs:
AP=unifrnd(-D/2,D/2,M,2); % Mx2 numbers in [-D/2, D/2]

%Wrapped around (8 neighbor cells)
D1=zeros(M,2);
D1(:,1)=D1(:,1)+ D*ones(M,1);
AP1=AP+D1;

D2=zeros(M,2);
D2(:,2)=D2(:,2)+ D*ones(M,1);
AP2=AP+D2;

D3=zeros(M,2);
D3(:,1)=D3(:,1)- D*ones(M,1);
AP3=AP+D3;

D4=zeros(M,2);
D4(:,2)=D4(:,2)- D*ones(M,1);
AP4=AP+D4;

D5=zeros(M,2);
D5(:,1)=D5(:,1)+ D*ones(M,1);
D5(:,2)=D5(:,2)- D*ones(M,1);
AP5=AP+D5;

D6=zeros(M,2);
D6(:,1)=D6(:,1)- D*ones(M,1);
D6(:,2)=D6(:,2)+ D*ones(M,1);
AP6=AP+D6;

D7=zeros(M,2);
D7=D7+ D*ones(M,2);
AP7=AP+D7;

D8=zeros(M,2);
D8=D8- D*ones(M,2);
AP8=AP+D8;

%Randomly locations of K terminals:
Ter=unifrnd(-D/2,D/2,K,2);
sigma_shd=8; %in dB

%Create an MxK large-scale coefficients beta_mk
BETAA = zeros(M,K);
dist=zeros(M,K);
betadB = zeros(M,K);
for m=1:M
    for k=1:K
    dist(m,k) = min([norm(AP(m,:)-Ter(k,:)), norm(AP1(m,:)-Ter(k,:)),norm(AP2(m,:)-Ter(k,:)),norm(AP3(m,:)-Ter(k,:)),norm(AP4(m,:)-Ter(k,:)),norm(AP5(m,:)-Ter(k,:)),norm(AP6(m,:)-Ter(k,:)),norm(AP7(m,:)-Ter(k,:)),norm(AP8(m,:)-Ter(k,:)) ]); %distance between Terminal k and AP m
    % Terminal chooses the nearest AP.
    if dist(m,k)<d0
         betadB=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
         betadB= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
    else
    betadB = -L - 35*log10(dist(m,k)) + sigma_shd*randn(1,1); %large-scale in dB 
    end
    
    BETAA(m,k)=10^(betadB/10); % eq.52 change from dB to non_dB
    end
end

%% Pilot Asignment: (random assign tau pilots for K users) 
% Create Phii(tau x K) pilot matrix. 
Phii=zeros(tau,K);
for k=1:K 
   Point=randi([1,tau]); 
   Phii(:,k)=U(:,Point);
end 

%% Calculate eta based on proposed algorithm
[Num,I_temp,iter,Phii_cf,etak_opt,PhiPhi,t,Betak]= PilotPC(Phii,K,BETAA,tau,Pp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Downlink
%  Rate_0(i,:) = opt_downlink(M,L1,K,BETAA,tau,Pp,Pd,Phii_cf,etak_opt,I_temp); % Data power control
%  R_1(i,:) = maxmin_opt_downlink(M,L1,K,BETAA,tau,Pp,Pd,Phii_cf,etak_opt,I_temp);  % Maxmin + Data power control
%% Uplink
 Rate_0(i,:) =  opt_uplink(L1,K,tau,Pp,Pu,Phii_cf,etak_opt,PhiPhi,Num,Betak);  % Data power control
 R_1(i,:) =  maxmin_opt_uplink(L1,K,tau,Pp,Pu,Phii_cf,etak_opt,PhiPhi,Num,Betak); % Maxmin + Data power control 
 end
 % Plot
Y=linspace(0,1,N*K);
Rate_greedy_new = reshape(Rate_0,N*K,1);
Rate_greedy_maxmin=reshape(R_1,N*K,1);
hold on
plot(9*sort(Rate_greedy_new),Y(:),'--r')
plot(9*sort(Rate_greedy_maxmin),Y(:),'b')
legend('Data PC','Pilot + Data PC')
