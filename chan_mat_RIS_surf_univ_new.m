function [Hdir,H1,H2] = chan_mat_RIS_surf_univ_new(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,varargin)

lambda = 3e8/f;     % Wavelength
dt = lambda/2;      % TX antenna space
dr = lambda/2;      % RX antenna space
dris = lambda/2;    % RIS element space
k = 2*pi/lambda;    % Wavenumber

% Geometrical placement
% x, y and z axis 
% TX antenna array
tx_arr(1,:) = zeros(1,Nt); 
tx_arr(2,:) = (sort(0:Nt-1,'descend')-(Nt-1)/2)*dt+lt; 
tx_arr(3,:) = zeros(1,Nt); 
% RX antenna array
rx_arr(1,:) = D*ones(1,Nr);
rx_arr(2,:) = (sort(0:Nr-1,'descend')-(Nr-1)/2)*dr+lr; 
rx_arr(3,:) = zeros(1,Nr);
% RIS 
center = [dist_ris 0]; % RIS center position 
N1 = sqrt(Nris);
N2 = N1;                                    % Number of RIS elements in two dimensions N1 and N2
ris_pos = RISPosition(N1,N2,dris,center);   % RIS elements' coordinates 
a = repmat(ris_pos{1},N1,1);                % Placing RIS elements in proper coordinates
ris_arr(1,:) = a(:)';        
ris_arr(2,:) = zeros(1,Nris);
ris_arr(3,:) = repmat(ris_pos{2},1,N2); 

if isempty(varargin)                        % Load the FSPL of the direct link
    alpha = 2;
else
    alpha = varargin{1};
end

% direct TX-RX paths/channel matrix
for i1 = 1:Nr                                        % Distance between the TX and RX antennas                                           
    for j1 = 1:Nt
        d(i1,j1) = norm(rx_arr(:,i1)-tx_arr(:,j1));
    end 
end
Hdir_los = exp(-1i*k*d);                             % Direct link, LOS matrix exponents 
tx_rx_dist = sqrt(D^2+(lt-lr)^2);                    % TX-RX distance   
FSPL_dir = (lambda/(4*pi))^2/tx_rx_dist^alpha(1);    % Inversion of the FSPL of the direct link 
Hdir = Rician(Hdir_los,sqrt(FSPL_dir),no_mat,K);     % Direct link channel matrix           


% indirect paths (TX-RIS-RX)
for l1 = 1:Nris                                      % Distance between the RIS elements and the RX antennas                                                            
    for r1 = 1:Nr  
        d2(r1,l1) = norm(rx_arr(:,r1)-ris_arr(:,l1)); 
    end
    for t1 = 1:Nt                                    % Distance between the RIS elements and the TX antennas                  
        d1(l1,t1) = norm(tx_arr(:,t1)-ris_arr(:,l1));   
    end
end

tx_ris_dist = sqrt(dist_ris^2+lt^2);                 % TX-RIS distance
ris_rx_dist = sqrt((D-dist_ris)^2+lr^2);             % RIS-RX distance   

FSPLindir = lambda^4/(256*pi^2)*...                  % Inversion of the FSPL of the indirect link 
           ((lt/tx_ris_dist+lr/ris_rx_dist)^2)*...
           1/(tx_ris_dist*ris_rx_dist)^2;

% TX-RIS channel matrix
H1_los = exp(-1i*k*d1);                             % TX-RIS link, LOS matrix exponents  
FSPL_1 = sqrt(FSPLindir);                           % FSPL of the indirect link is embedded in the TX-RIS channel matrix 
H1 = Rician(H1_los,FSPL_1,no_mat,K);

% RIS-RX channel matrix
H2_los = exp(-1i*k*d2);                             % RIS-RX link, LOS matrix exponents 
FSPL_2 = 1;
H2 = Rician(H2_los,FSPL_2,no_mat,K);
end

function pos = RISPosition(N1,N2,dist,center)                 % Determine positions of RIS elements
d1 = (0:N1-1)-(N1-1)/2;
d2 = (0:N2-1)-(N2-1)/2;
pos{1} = center(1)+d1*dist;
pos{2} = center(2)+d2*dist;
end 

function Hout = Rician(Hlos,FSPL,no_mat,K)                     % Create the Rician channel matices
Hlos = repmat(Hlos,no_mat,1);
Hnlos = sqrt(1/2)*(randn(size(Hlos))+1i*randn(size(Hlos)));
Htot = FSPL/sqrt(K+1)*(Hlos*sqrt(K)+Hnlos);
dim = size(Hlos,1)/no_mat;
for ind = 1:no_mat
   Hout{ind} = Htot((ind-1)*dim+1:ind*dim,:); 
end
end
