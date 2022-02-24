clear; close all; clc;

Nt = 8;                 % Number of TX antennas
Nr = 4;                 % Number of RX antennas
Nris = 15^2;            % Number of RIS elements
       
K = 1;                  % Rician factor    
D = 500;                % TX-RX distance
dist_ris = 40;          % RIS distance from TX
f = 2e9;                % Frequency

lt = 20;                % TX position 
lr = 100;               % RX position 
Pt = 1;                 % Transmit power in Watts
N0 = -120;              % Noise power in dB
SNR = db2pow(-N0);      % SNR            
no_mat = 10;            % Number of channel realizations 

no_iter = 500;          % Number of iterations 
alpha_dir = 3;          % FSPL exponent of the direct link

% Generate channel matrices (WORNING: It is only works for Nris that is a square number)
[Hdirt,H1t,H2t] = chan_mat_RIS_surf_univ_new(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,alpha_dir);

Cpgm = zeros(1,no_iter+1);

for i = 1:no_mat 
    Hdir = Hdirt{i}; H1 = H1t{i}; H2 = H2t{i};                  
    
    % Scaling factor
    c = sqrt(norm(Hdir)/norm(H2*H1))*max(sqrt(Pt),1)/sqrt(Pt)*10;     
    
    % Initial Q matrix and RIS phase shifts    
    Qinit = eye(Nt)*(Pt/Nt);
    omega_init = ones(1,Nris);
    
    % PGM iterative optimization 
    [dCpgm,~] = PGM_opt(Pt,Hdir*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
    Cpgm = Cpgm+dCpgm;
    
end

semilogx(1:length(Cpgm),Cpgm/no_mat,'r','DisplayName','PGM');
xlabel('Iteration number'); ylabel('Achievable rate [bit/s/Hz]');
xlim([0 no_iter]); 
legend('show','Location','SouthEast');
print('../results/Achievable_Rate', '-dpdf')

