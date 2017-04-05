
clc; clear all; close all;
format short g
% ******************************************************************** %%
%                       Simulation Parameters Setting                           
% *********************************************************************** %
dt      = 0.01;                 % Simulation interval
T_stop  = 3;                    % Final simulation time
Sim_t = 0:dt:T_stop;
T_length = length(Sim_t);


% ******************************************************************** %%
%            Original Plant( non-square , non-minimal phase )
% *********************************************************************** %
% load nonminimal_system    % { A , B , C , D}
% load De_term

% A = [   -2.1     0.3     0.5     0.9;
%          0.2     0.3       0     0.7;
%          0.5    -1.1    -2.3    -0.2;
%          0.8     1.7    -1.2    -0.3];
% B = [   1   -0.4;
%      -1.7   -1.7;
%      -0.2    1.2; 
%       0.7    0.1];
% C = [2  -2  -1   4;
%      2   3   4   4];

A = [-31.31 0 -2.833e4;0 -10.25 8001;1 -1 0];
B = [28.06;0;0];
Bd = [0;7.210;0];
C = [1 0 0];
D = zeros(size(C,1));

De = [0.997218557152221;0.869829858547889;7.02606509638271e-05];

n = size(A,1);               % No. of states
m = size(B,2);               % No. of inputs
p = size(C,1);               % No. of outputs
D = zeros(p,m);
Zeros_ori = TZOCS(A, B, C, D)
Qc = eye(p);
Rc = 1e-6*eye(m);
R_bar = Rc + D'*Qc*D;
N     = C'*Qc*D;
[Kc,~] =  lqr(A, B, C'*Qc*C, R_bar, N);
 Ec    = -inv(R_bar)*( ( C - D*Kc )*inv( A - B*Kc )*B - D )'*Qc;
Zeros_origin = tzero(A-B*Kc, B*Ec, C-D*Kc, D*Ec);
z_p = eig(A-B*Kc);
% ~~~ Check whether the original plant is controllable and observable

Mc_ori = ctrb(A, B);           % Controllability matrix 
Mo_ori = obsv(A, C);           % Observability matrix

Rc_ori = rank(Mc_ori);         % Check whether the plant is controllable (the system is controllable if Mc has full rank n )
Ro_ori = rank(Mo_ori);         % Check whether the plant is observable   (the system is observable   if Mo has full rank n )

% generate eta
% eta = generate_eta(A,B,C,D,p);

% load 'eta_value'                % eta ( trial and error  , randomly generated )
% eta=10*eta;
eta = eye(1);
C_tilde = eta * C;
D_tilde = eta * D;
Zeros_tilde = TZOCS(A, B, C_tilde, D_tilde);

Qc = eye(p);
Rc = 1e-6*eye(m);
R_bar = Rc + D_tilde'*Qc*D_tilde;
N = C_tilde'*Qc*D_tilde;
[Kc,~] =  lqr(A, B, C_tilde'*Qc*C_tilde, R_bar, N);
Ec = -inv(R_bar)*( ( C_tilde - D_tilde*Kc )*inv( A - B*Kc )*B - D_tilde )'*Qc;
Zeros_tilde = tzero(A-B*Kc, B*Ec, C_tilde-D_tilde*Kc, D_tilde*Ec);
z_p = eig(A-B*Kc);

% ~~~ Check whether the plant is observable
Mo_tilde = obsv(A, C_tilde);  % Observability matrix
Ro_tilde = rank(Mo_tilde);    % Check whether the plant is observable ( the system is observable if Mo has full rank n )
if Ro_tilde ~= size(A,1) || Ro_tilde ~= size(A,1)
    disp(['[' 8 '[err] Rank deficient]' 8])
    break
end

% generate De
% De = generate_De(A,C_tilde,n,p);

%%
% ~~~ Check whether the plant is observable
Rc_obs = rank(ctrb(A, De));         % Check whether the plant is controllable (the system is controllable if Mc has full rank n )
Ro_obs = rank(obsv(A, C_tilde));         % Check whether the plant is observable
if Rc_obs ~= size(A,1) || Ro_obs ~= size(A,1)
    disp(['[' 8 '[err] Rank deficient]' 8])
    break;
end
%% ******************************************************************** %%
%                        Create a Reference Input 
%                            ( discontinuous )
% *********************************************************************** % 

i = 0;
for ts = 0 : dt : T_stop
    
    i = i + 1;
    ts_c(:,i) = ts;            % time span (ts)   
    ref_c(:,i) = ref_func(ts);  % reference 
    
end
ref_c(2,:) = [];
ref_eta = eta*ref_c;
% x_aug(:,1) = zeros(n+p,1);             % Initial condition
% xc(:,1) = 0.001*[ 10 -2 -1 -5 1 2];
xc(:,1) = 0.001*[ 10 -2 -1 ];
x_eta(:,1) =(C_tilde*xc(:,1)-ref_eta(:,1))*dt; %%%%% PID controller----------------
x_aug(:,1) = [xc(:,1) ; x_eta(:,1) ];
xe(:,1) = pinv(C)*C*xc(:,1);                  % estimator state Initial condition
xe_aug(:,1) = [ xe(:,1) ; x_eta(:,1)];        % Augumented system PID Initial condition

%% ******************************************************************** %%
%                        Model-Following Approach 
%          ( square , minimal phase , controllable , observable )  
% ********************************************************************** %

A_r = [ -1  -1   0   0  ;
         1  -1   0   0  ;
         0   0  -2   6  ;
         0   0  -6  -2 ];

B_r = [ 1   ;
        0   ;
        0   ;
        0  ];

C_r = [ 1  0  0  0];

nr   = size(A_r,1);     % No. of states 
mr   = size(B_r,2);     % No. of inputs
pr   = size(C_r,1);     % No. of outputs

D_r  = zeros(pr,mr);    % Without a direct feed-through term

x_r(:,1) = zeros(nr,1); % Initial condition

% ~~~ Check whether the model is controllable and observable

Mc_r = ctrb(A_r, B_r);  % Controllability matrix 
Mo_r = obsv(A_r, C_r);  % Observability matrix

Rc_r = rank(Mc_r);      % Check whether the plant is controllable (the system is controllable if Mc has full rank nr )
Ro_r = rank(Mo_r);      % Check whether the plant is observable   (the system is observable   if Mo has full rank nr )

% ~~~ Linear Quadratic State-Feedback Tracker Design

Qr = eye(pr);
Rr = 1e-6*eye(mr);

[Kr,Pr] = lqr(A_r, B_r, C_r'*Qr*C_r, Rr);
Er = -inv(Rr)*B_r'*inv( A_r - B_r*Kr )'*C_r'*Qr;

% ~~~ Perform Simulation

opts_sim = simset('Solver', 'ode45', 'FixedStep', dt, 'InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-7);

[T_r, X_r, U_r, Y_r] = sim('Model_Following_Approach', ts_c, opts_sim, [ ts_c' , ref_eta' ]);

x_r = X_r';             % state
u_r = U_r';             % control input
y_r = Y_r';             % step response
e_r = ref_eta - y_r;      % tracking error
% plot(C_r*( A_r*x_r + B_r*u_r ))

%% ******************************************************************** %%
%                 Step 2 : PID Filter ( assigning zeros )                     
% *********************************************************************** %
% %  The assigned zeros 
PID_roots=[-1e2 -1e2];

i = 1;
for ii = 1:p
    S(ii,ii)=PID_roots(ii,i)+PID_roots(ii,i+1);
    P(ii,ii)=PID_roots(ii,i)*PID_roots(ii,i+1);
end

cof = 1e0;

Kd = cof*eye(p) % Derivative   gain (small value)
Ki = P*Kd
Kp = -S*Kd

Kp_1=  Kp(1,1);
% Kp_2=  Kp(2,2);
Ki_1 = Ki(1,1);               % Integral gain   
% Ki_2 = Ki(2,2);
Kd_1 = Kd(1,1);               % Integral gain   
% Kd_2 = Kd(2,2);


%% ******************************************************************** %%
%                         Step 3 : Augmented Plant   
%              ( minimal phase , controllable , observable )  
% *********************************************************************** %
% ----------------------------PID Controller-----------------------------
A_aug = [    A    , zeros(n,p) ;
          C_tilde , zeros(p,p) ];

B_aug = [     B       ;
          zeros(p,m) ];

C_aug = [ Kp*C_tilde + Kd*C_tilde*A , Ki];
 
D_aug = Kd*C_tilde*B;

% ~~~ Check whether the augmented plant is controllable and observable
Mc_aug = ctrb(A_aug, B_aug);           % Controllability matrix 
Mo_aug = obsv(A_aug, C_aug);           % Observability matrix

Rc_aug = rank(Mc_aug);                 % Check whether the plant is controllable ( the system is controllable if Mc has full rank n+p )
Ro_aug = rank(Mo_aug);                 % Check whether the plant is observable   ( the system is observable   if Mo has full rank n+p )
% ******************************************************************** %%
%             Linear Quadratic State-Feedback Tracker Design                          
% *********************************************************************** %

% ~~~ Weighting matrices {Qc,Rc} , high ratio of Qc to Rc

Qc = eye(p);
Rc = 1e-7*eye(m);

R_bar = Rc + D_aug'*Qc*D_aug;
N     = C_aug'*Qc*D_aug;

% ~~~ Stata-Feedback and Forward Gains {Kc,Ec} 

[Kc,Pc] =  lqr(A_aug, B_aug, C_aug'*Qc*C_aug, R_bar, N);
 Ec     = -inv(R_bar)*( ( C_aug - D_aug*Kc )*inv( A_aug - B_aug*Kc )*B_aug - D_aug )'*Qc;

Kc1 = Kc(:,1:n);
Kc2 = Kc(:,n+1:end);

% ******************************************************************** %%
%                  Compute the d(t) , s(t) , and Cc(t)                                    
% *********************************************************************** %

i = 0;
for ts = 0 : dt : T_stop
    
    i = i + 1;
     
    d_au(:,i) = [  zeros(n,1) ; -ref_eta(:,i)];
    s_au(:,i) = -Kp*ref_eta(:,i) - Kd*C_r*( A_r*x_r(:,i) + B_r*u_r(:,i) ); 
    Cc(:,i) = -Ec*s_au(:,i) + inv(R_bar)*B_aug'*inv( ( A_aug - B_aug*Kc )' )*Pc*d_au(:,i);
     
end

Ip = [ zeros(p,n)  eye(p) ];
Ip3= [ eye(n,n)  zeros(n,p) ];
Ip2 = [ zeros(n,p) ; eye(p)];
C_bar = [ C_tilde  zeros(p,p)];  

% ~~~ Check whether the plant is observable

Rc_obs = rank(ctrb(A, De));         % Check whether the plant is controllable (the system is controllable if Mc has full rank n )
Ro_obs = rank(obsv(A, C_tilde));         % Check whether the plant is observable
 

%% ---------------PID-filter Observer-------------------

% %  The assigned zeros 
PID_roots=[-1e2 -1e2];

i=1;
for ii=1:p
    S(ii,ii)=PID_roots(ii,i)+PID_roots(ii,i+1);
    P(ii,ii)=PID_roots(ii,i)*PID_roots(ii,i+1);
end

cof=1;

Ked=cof*eye(p); % Derivative   gain (small value)
Kei=P*Ked;
Kep=-S*Ked;

Kep_1=  Kep(1,1);
Kei_1 = Kei(1,1);               % Integral gain
Ked_1 = Ked(1,1);               % Integral gain

%%
%-------------------PID filter augmented plant_Observer----------------------
Ae_aug = [ A-De*inv(eye(p)+Ked*C_tilde*De)*Ked*C_tilde*A     De*inv(eye(p)+Ked*C_tilde*De)*Kei  ;
                        -C_tilde                                  zeros(p,p) ];

Be_aug = [ De*inv(eye(p)+Ked*C_tilde*De)*Kep;
                         zeros(p,p)         ];

Ce_aug = [C_tilde  zeros(p,p)];
De_aug = zeros(p,p);
 
% ~~~ Check whether the augmented plant is controllable and observable
Rec_aug = rank(ctrb(Ae_aug, Be_aug));                 % Check whether the plant is controllable ( the system is controllable if Mc has full rank n+m )
Reo_aug = rank(obsv(Ae_aug, Ce_aug));                 % Check whether the plant is observable   ( the system is observable   if Mo has full rank n+m )
%%
Qe = eye(p);
Re = 1e-9*eye(p);
% ~~~ State-Feedback and Forward Gains {Ke,Ee_aug}

[Ke,Pe_aug] =  lqr(Ae_aug - 0*eye(size(Ae_aug,1)), Be_aug, Ce_aug'*Qe*Ce_aug, Re);
Ee      = -inv(Re)*( ( Ce_aug )*inv( Ae_aug - Be_aug*Ke )*Be_aug )'*Qe;
Ke1 = Ke(:,1:n);
Ke2 = Ke(:,n+1:end);
i = 0;
%
tzero(Ae_aug, Be_aug, Ce_aug, De_aug)
tzero(Ae_aug-Be_aug*Ke, Be_aug*Ee, Ce_aug-De_aug*Ke, De_aug*Ee)
% eigAe_aug_Be_aug_Ke = eig(Ae_aug-Be_aug*Ke)

eigA_DeKe = eig(A-De*Ke1)

J = De*inv(eye(p)+Ked*C_tilde*De)*Ked;

de_aug_a = J*C_r*( A_r*x_r + B_r*u_r);
de_aug_b = J*C_tilde*B; % de_aug_b * u = de_aug
Ced=inv(Re)*Be_aug'*inv((Ae_aug-Be_aug*Ke)')*Pe_aug; % Ced* d = Ce

%% UIE
h = 0; % shift the y-axis
Qod = 1e8*eye(n);
Rod = 1e0*eye(p);
% H: shift the y-axis
Lt = lqr((A-De*Ke1)'+h*eye(n), C_tilde', Qod, Rod)'; % regulator
Kd = pinv(B)*Lt;
eigAt_LtC = eig((A-De*Ke1)-B*Kd*C_tilde)
eigAt_LtC = eig((A-De*Ke1)-Lt*C_tilde);
eigAt = eig(A-De*Ke1)
% rank(ctrb(ss(A-De*Ke1, B, C_tilde,0)));
Ccd = inv(R_bar)*B_aug'*inv( ( A_aug - B_aug*Kc )' )*Pc*B_aug;
opts_sim = simset('FixedStep', dt, 'InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-7);

% subsystem of disturbance
Bfi = 8e0;
Cfi = 1.25e2;

Af = -Bfi*Cfi*eye(size(B,2));
Bf = Bfi*eye(size(B,2));
Cf = Cfi*eye(size(B,2));

de = 10*disturb2(Sim_t);
% set(0,'ShowHiddenHandles','on'); set(gcf,'menubar','figure')



% Kc1 = lqr(A, B, C'*Qc*C, Rc); % regulator
L = lqr(A' + 30*eye(n), C_tilde', Qod, Rod)'; % regulator
eigAt_LtC = eig(A-L*C_tilde)

[Tc, xx, uc, xc, yc, y_xi, ye_xi, yf, de_hat] = ...
    sim('plant3_conference', ts_c, opts_sim);
% [Tc, xx, uc, xc, yc, y_xi, ye_xi, yf] = sim('plant3_Conference', ts_c, opts_sim);


figure % r vs y
% subplot(211), 
plot(ts_c, ref_c(1,:), 'b-', ts_c, yc(:,1)', 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} r_{c1} vs. y_{c1}')
legend('\fontsize{10} r_{c1}','\fontsize{10} y_{c1}','Location','NorthWest','Orientation','Horizontal')
%  axis([0 T_stop -100 100])
% subplot(212), plot(ts_c, ref_c(2,:), 'b-', ts_c, yc(:,2)', 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} r_{c2} vs. y_{c2}')
% legend('\fontsize{10} r_{c2}','\fontsize{10} y_{c2}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -120 170]) 

figure % de
% subplot(211),
plot(ts_c, de(:,1), 'b-', ts_c, de_hat(:,1)', 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} d_{e} vs. d_{e}hat')
legend('\fontsize{10} d_{e}','\fontsize{10} d_{e}hat','Location','NorthWest','Orientation','Horizontal')
% subplot(212), plot(ts_c, de(:,2), 'b-', ts_c, de_hat(:,2)', 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} r_{c2} vs. y_{c2}')
% legend('\fontsize{10} r_{c2}','\fontsize{10} y_{c2}','Location','NorthWest','Orientation','Horizontal')



