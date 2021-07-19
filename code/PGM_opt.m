% Function that implements the PGM iterative optimization

function [Cout,iter_time] = PGM_opt(Pt,Hdir,H1,H2,maxIter,Qinit,myomegaunitball,c)

% Load the initial covariance matrix and RIS phase shifts  
myomegaunitcirc = myomegaunitball.';
Q = Qinit; 

% Line-search parameters 
delta = 1e-5;       
rho = 0.5;

iIter = 0;         
stepsize = 10000;   % Inital step size 


Cout = [RIScap(Hdir,H2,H1,myomegaunitcirc,Q)]; % Initial achievable rate 
Cprev = Cout;

iter_time = 0;
tic
while iIter<maxIter
    iIter = iIter+1;
    q_Q = gracov(Hdir,H2,H1,myomegaunitcirc,Q);     % gradient w.r.t. Q
    g_RIS = gradRIS(Hdir,H2,H1,myomegaunitcirc,Q);  % gradient w.r.t. RIS
   
    for iLineSearch = 0:30

        % Updated covariance (Q) matrix, line 3 of Algorithm 1 
        Qnew = Q + stepsize*q_Q;
        Qnew = cov_mat_proj_modified(Qnew,Pt*c^2);
        
        % Updated RIS phase shifts, line 4 of Algorithm 1
        yunitcirc = myomegaunitcirc + stepsize*g_RIS;
        myomegaunitcircnext = projectontounitcircle(yunitcirc)/c;
        
        % New achievable rate 
        Cnew = RIScap(Hdir,H2,H1,myomegaunitcircnext,Qnew);
        
        % Line-search procedure
        if (((Cnew-Cprev) >= delta*(norm(Qnew-Q)^2+ ...
                norm(myomegaunitcircnext-myomegaunitcirc)^2)) ...
            || (stepsize<1e-4) )
            % If the step size satisfies (35c) OR it is too small:    
            % Obtained Q matrix and RIS phase shifts
            myomegaunitcirc = myomegaunitcircnext;
            Q = Qnew;
            Cprev = Cnew;
            break
        else 
            % Reduce the step size 
            stepsize=stepsize*rho;
        end
    end
    % New achievable rate
    Cout = [Cout Cnew];
    % Exection time of an iteration
    iter_time = [iter_time toc];
end
end

% Calculation of gradient w.r.t. Q
function y = gracov(Hdir,H2,H1,myomega,Q)
Z = Hdir+H2*diag(myomega)*H1;
Nr = size(H2,1);
y = Z'*inv(eye(Nr)+Z*Q*Z')*Z;
end

% Calculation of gradient w.r.t. RIS
function y = gradRIS(Hdir,H2,H1,myomega,Q)
Z = Hdir+H2*diag(myomega)*H1;
Nr = size(H2,1);
y = diag(H2'*inv(eye(Nr)+Z*Q*Z')*Z*Q*H1');
end

% Achievable rate calculation
function y = RIScap(Hdir,H2,H1,myomega,Q)
Z = Hdir+H2*diag(myomega)*H1;
Nr = size(H2,1);
y = real(log(det(eye(Nr)+Z*Q*Z')))/log(2);
end

% RIS projection
function y = projectontounitcircle(x)
y = x./abs(x);
end

% Covarinace matrix projection
function Qnew = cov_mat_proj_modified(Qold,Pt)
    [U,D] = eig(Qold);
    Dnew = water_fill(Pt,real(diag(D))).';
    Qnew = U*diag(Dnew)*U';
end 

% Water-filling algorithm
function vect_out = water_fill(Pt,vect_in)
    vect_in = vect_in.';
    [sort_val,sort_idx] = sort(vect_in,'descend');
    for n = length(vect_in):-1:1
        water_level = (sum(sort_val(1:n))-Pt)/n;
        di = sort_val(1:n)-water_level;
        if di(:)>=0
            break
        end   
    end
    vect_out = zeros(1,length(vect_in));
    vect_out(sort_idx(1:n)) = di;
end 




