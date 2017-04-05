
% ======================================================================= %
%                               2015/01/29 
%     Compute the Transmission Zeros of Continuous-Time MIMO Systems
% ======================================================================= %
% 
%  Limitation : m > p
%
% dx =  Ax + Bu         Open-loop system (Nonsquare system)
%  y =  Cx + Du
%
%  u = -Kx + Er         Control law
%
% dx = (A-BK)x + BEr    Closed-loop system (Square system) 
%  y = (C-DK)x + DEr
%
% Transmission Zeros of Continuous-Time Systems ( TZOCS )
%
% Transmission Zeros = TZOCS( Continuous-Time Systems )
%
% Syntax : Z = TZOCS(A,B,C,D)

function Z = TZOCS(A,B,C,D)

m = size(B,2); % No. of inputs
p = size(C,1); % No. of outputs

%% ******************************************************************** %%
%             Linear Quadratic State-Feedback Tracker Design                          
% *********************************************************************** %

% ~~~ Weighting matrices {Qc,Rc} , high ratio of Qc to Rc

Qc = eye(p);  
Rc = 1e-5*eye(m);
R_bar = Rc + D'*Qc*D;
N     = C'*Qc*D;

% ~~~ State-Feedback and Forward Gains {Kc,Ec} 

[Kc,Pc] =  lqr(A, B, C'*Qc*C, R_bar, N);
 Ec     = -inv(R_bar)*( ( C - D*Kc )*inv( A - B*Kc )*B - D )'*Qc;


% ~~~ State-Feedback and Forward Gains {Kc,Ec} 
%% ******************************************************************** %%
%                      Compute the Transmission Zeros  
%    ( The invariamt zeros are not changed by constant state feedback )
% *********************************************************************** %

Z = tzero(A-B*Kc, B*Ec, C-D*Kc, D*Ec);
z_p=eig(A-B*Kc);
%% ******************************************************************** %%
%                      Compute the Transmission Zeros  
%                               ( Check rank )
% *********************************************************************** %

% pole_cl = eig(A-B*Kc);
% 
% Sc = [ A-B*Kc , B*Ec ;
%        C-D*Kc , D*Ec ];
% 
% Sz = [ pole_cl(6,:)*eye(n) , zeros(n,p)  ;
%            zeros(p,n)       , zeros(p,p) ];  
%      
% [U,S,V] = svd(Sc - Sz);
% 
% min( diag(S) )

end
