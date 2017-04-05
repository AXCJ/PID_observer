
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
A = [   -2.1     0.3     0.5     0.9;
         0.2     0.3       0     0.7;
         0.5    -1.1    -2.3    -0.2;
         0.8     1.7    -1.2    -0.3];
B = [   1   -0.4;
     -1.7   -1.7;
     -0.2    1.2; 
      0.7    0.1];
C = [2  -2  -1   4;
     2   3   4   4];
load nonminimal_system    % { A , B , C , D}
load De_term

n = size(A,1);               % No. of states
m = size(B,2);               % No. of inputs
p = size(C,1);               % No. of outputs
D = zeros(p,m);
% Zeros_ori = TZOCS(A, B, C, D)
Qc = eye(p);
Rc = 1e-6*eye(m);
R_bar = Rc + D'*Qc*D;
N     = C'*Qc*D;
[Kc,~] =  lqr(A, B, C'*Qc*C, R_bar, N);
 Ec    = -inv(R_bar)*( ( C - D*Kc )*inv( A - B*Kc )*B - D )'*Qc;
Zeros_origin = tzero(A-B*Kc, B*Ec, C-D*Kc, D*Ec)
z_p = eig(A-B*Kc);
% ~~~ Check whether the original plant is controllable and observable

Mc_ori = ctrb(A, B);           % Controllability matrix 
Mo_ori = obsv(A, C);           % Observability matrix

Rc_ori = rank(Mc_ori);         % Check whether the plant is controllable (the system is controllable if Mc has full rank n )
Ro_ori = rank(Mo_ori);         % Check whether the plant is observable   (the system is observable   if Mo has full rank n )

load 'eta_value'                % eta ( trial and error  , randomly generated )
% eta=10*eta;
% generate eta
% eta = generate_eta(A,B,C,D);

% eta = eye(2);
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
    break;
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

ref_eta = eta*ref_c;
% x_aug(:,1) = zeros(n+p,1);             % Initial condition
xc(:,1) = 0.001*[ 10 -2 -1 -5 1 2];
% xc(:,1) = 0.01*[ 10 -2 -1 -5];
x_eta(:,1) =(C_tilde*xc(:,1)-ref_eta(:,1))*dt; %%%%% PID controller----------------
x_aug(:,1) = [xc(:,1) ; x_eta(:,1) ];
xe(:,1) = pinv(C)*C*xc(:,1);                   % estimator state Initial condition
xe_aug(:,1) = [ xe(:,1) ; x_eta(:,1)];        % Augumented system PID Initial condition

%% ******************************************************************** %%
%                        Model-Following Approach 
%          ( square , minimal phase , controllable , observable )  
% ********************************************************************** %

A_r = [ -1  -1   0   0  ;
         1  -1   0   0  ;
         0   0  -2   6  ;
         0   0  -6  -2 ];

B_r = [ 1  0  ;
        0  0  ;
        0  1  ;
        0  0 ];

C_r = [ 1  0  0  0  ;
        0  0  1  0 ];

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

%% ******************************************************************** %%
%                 Step 2 : PID Filter ( assigning zeros )                     
% *********************************************************************** %
% %  The assigned zeros 
PID_roots=[-5e3 -5e0; -5e3 -5e0];

 i=1;
for ii=1:p
 S(ii,ii)=PID_roots(ii,i)+PID_roots(ii,i+1);
 P(ii,ii)=PID_roots(ii,i)*PID_roots(ii,i+1);
end

cof=1e-2;

Kd=cof*eye(p); % Derivative   gain (small value)
Ki=P*Kd;
Kp=-S*Kd;

Kp_1=  Kp(1,1);
Kp_2=  Kp(2,2);
Ki_1 = Ki(1,1);               % Integral gain   
Ki_2 = Ki(2,2);
Kd_1 = Kd(1,1);               % Integral gain   
Kd_2 = Kd(2,2);


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
Rc = 1e-4*eye(m);

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
simin = [ts_c' ref_eta'];
% ******************************************************************** %%
%                             Perform Simulation                           
% *********************************************************************** %
% ref_f=zeros(p,T_length);
% [Tc, Xc, Uc, Yf] = sim('Augmented_Plant', ts_c, opts_sim, [ ts_c' , ref_eta' , Cc' , d_au' , s_au' ]);
[Tc, xx, Xc, Uc, Yf] = sim('plant1_noObserver', ts_c, opts_sim);
 
u_c = Uc';            % control input
x_c = Xc(:,1:n)';     % states
y_c = C*x_c;          % response

x_eta = Xc(:,n+1:end)';
y_eta = C_tilde*x_c;
y_f = Yf';            % response

e_tra = ref_eta - y_f;  % tracking error

% ******************************************************************** %%
%                             Simulation Results                           
% *********************************************************************** %

figure
starbit = 10;
subplot(211), plot(ts_c(starbit:end), ref_eta(1,starbit:end), 'b-', ts_c(starbit:end), y_f(1,starbit:end), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} r_{\xi,1} vs. y_{f,1}')
legend('\fontsize{10} r_{\xi,1}','\fontsize{10} y_{f,1}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -100 100])
subplot(212), plot(ts_c(starbit:end), ref_eta(2,starbit:end), 'b-', ts_c(starbit:end), y_f(2,starbit:end), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} r_{\xi,2} vs. y_{f,2}')
legend('\fontsize{10} r_{\xi,2}','\fontsize{10} y_{f,2}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -120 170])

% figure
% subplot(211), plot(ts_c, ref_eta(1,:), 'b-', ts_c, y_e2ta(1,:), 'r:', 'LineWidth', 1.5)
% ylabel('\fontsize{10} r_{\xi,1} vs. y_{\xi,1}')
% legend('\fontsize{10} r_{\xi,1}','\fontsize{10} y_{\xi,1}','Location','NorthWest','Orientation','Horizontal')
% %  axis([0 T_stop -100 100])
% subplot(212), plot(ts_c, ref_eta(2,:), 'b-', ts_c, y_eta(2,:), 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} r_{\xi,2} vs. y_{\xi,2}')
% legend('\fontsize{10} r_{\xi,2}','\fontsize{10} y_{\xi,2}','Location','NorthWest','Orientation','Horizontal')
% % axis([0 T_stop -120 170])

figure
subplot(211), plot(ts_c, ref_c(1,:), 'b-', ts_c, y_c(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} r_{c,1} vs. y_{c,1}')
legend('\fontsize{10} r_{c,1}','\fontsize{10} y_{c,1}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -22 20])
subplot(212), plot(ts_c, ref_c(2,:), 'b-', ts_c, y_c(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} r_{c,2} vs. y_{c,2}')
legend('\fontsize{10} r_{c,2}','\fontsize{10} y_{c,2}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -12 12])
% 
% figure
% subplot(311), plot(ts_c, u_c(1,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,1}')
% % legend('\fontsize{10} u_{c,1}','Location','Northeast','Orientation','Horizontal')
% %  axis([0 T_stop -1.5e4 1.5e4])
% subplot(312), plot(ts_c, u_c(2,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,2}')
% % legend('\fontsize{10} u_{c,2}','Location','Southeast','Orientation','Horizontal')
% % axis([0 T_stop -1.5e4 1.7e4])
% subplot(313), plot(ts_c, u_c(3,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,3}')
% % legend('\fontsize{10} u_{c,3}','Location','Southeast','Orientation','Horizontal')
% % axis([0 T_stop -1.7e4 1.7e4])
% xlabel('\fontsize{10} Time ( sec. )')
% 
% figure
% subplot(211), plot(ts_c, x_c(1,:), 'b-', 'LineWidth', 1.5)
% legend('\fontsize{10} x_{c,1}')
% % axis([0 T_stop -0.1 0.1])
% subplot(212), plot(ts_c, x_c(2,:), 'b-', 'LineWidth', 1.5)
% legend('\fontsize{10} x_{c,2}')
% xlabel('\fontsize{10} Time ( sec. )')
% % axis([0 T_stop -0.1 0.1])

%% *****************ESTIMATOR********************************************
% *************************************************************************
%---------------------------PID Observer---------------------------------------------

%%
Ip = [ zeros(p,n)  eye(p) ];
Ip3= [ eye(n,n)  zeros(n,p) ];
Ip2 = [ zeros(n,p) ; eye(p)];
C_bar = [ C_tilde  zeros(p,p)];  

% ~~~ Check whether the plant is observable

Rc_obs = rank(ctrb(A, De));         % Check whether the plant is controllable (the system is controllable if Mc has full rank n )
Ro_obs = rank(obsv(A, C_tilde));         % Check whether the plant is observable
 

%% ---------------PID-filter Observer-------------------

% %  The assigned zeros 
PID_roots=[-5e3 -2; -5e3 -2];

 i=1;
for ii=1:p
 S(ii,ii)=PID_roots(ii,i)+PID_roots(ii,i+1);
 P(ii,ii)=PID_roots(ii,i)*PID_roots(ii,i+1);
end

cof=5e-3;

Ked=cof*eye(p); % Derivative   gain (small value)
Kei=P*Ked;
Kep=-S*Ked;

Kep_1=  Kep(1,1);
Kep_2=  Kep(2,2);
Kei_1 = Kei(1,1);               % Integral gain   
Kei_2 = Kei(2,2);
Ked_1 = Ked(1,1);               % Integral gain   
Ked_2 = Ked(2,2);
% % 
% % 

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
Re = 1e-6*eye(p);
% ~~~ State-Feedback and Forward Gains {Ke,Ee_aug}

[Ke,Pe_aug] =  lqr(Ae_aug - 0*eye(size(Ae_aug,1)), Be_aug, Ce_aug'*Qe*Ce_aug, Re);
Ee      = -inv(Re)*( ( Ce_aug )*inv( Ae_aug - Be_aug*Ke )*Be_aug )'*Qe;
Ke1 = Ke(:,1:n);
Ke2 = Ke(:,n+1:end);
i = 0;
%
% eigAe_aug_Be_aug_Ke = eig(Ae_aug-Be_aug*Ke)

eigA_DeKe = eig(A-De*Ke1)

J = De*inv(eye(p)+Ked*C_tilde*De)*Ked;

de_aug_a = J*C_r*( A_r*x_r + B_r*u_r);
de_aug_b = J*C_tilde*B; % de_aug_b * u = de_aug
Ced=inv(Re)*Be_aug'*inv((Ae_aug-Be_aug*Ke)')*Pe_aug; % Ced* d = Ce

%% UIE
h = 40; % shift the y-axis
Qod = 1e7*eye(n);
Rod = 1e0*eye(p);
% H: shift the y-axis
Lt = lqr((A-De*Ke1)'+h*eye(n), C_tilde', Qod, Rod)'; % regulator
Kd = pinv(B)*Lt;
eigAt_LtC = eig((A-De*Ke1)-B*Kd*C_tilde)    
eigAt_LtC = eig((A-De*Ke1)-Lt*C_tilde)
% rank(ctrb(ss(A-De*Ke1, B, C_tilde,0)));
Ccd = inv(R_bar)*B_aug'*inv( ( A_aug - B_aug*Kc )' )*Pc*B_aug;
opts_sim = simset('FixedStep', dt, 'InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-7);

% subsystem of disturbance
Bfi = 8e0;
Cfi = 1.25e2;

Af = -Bfi*Cfi*eye(size(B,2));
Bf = Bfi*eye(size(B,2));
Cf = Cfi*eye(size(B,2));

de = 50*disturb2(Sim_t);
% set(0,'ShowHiddenHandles','on'); set(gcf,'menubar','figure')

[Tc, xx, uc, xc, yc, y_xi, ye_xi, yf, de_hat] = sim('plant4_IEEE', ts_c, opts_sim);
% [Tc, xx, uc, xc, yc, y_xi, ye_xi, yf] = sim('plant3_Conference', ts_c, opts_sim);

% figure
% subplot(211), plot(ts_c, y_xi(:,1), 'b-', ts_c, ye_xi(:,1), 'r:', 'LineWidth', 1.5)
% ylabel('\fontsize{10} y_{\xi,1} vs. ye_{\xi,1}')
% legend('\fontsize{10} y_{\xi,1}','\fontsize{10} ye_{\xi,1}','Location','NorthWest','Orientation','Horizontal')
% %  axis([0 T_stop -100 100])
% subplot(212), plot(ts_c, y_xi(:,2), 'b-', ts_c, ye_xi(:,2), 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} y_{\xi,2} vs. ye_{\xi,2}')
% legend('\fontsize{10} y_{\xi,2}','\fontsize{10} ye_{\xi,2}','Location','NorthWest','Orientation','Horizontal')
% % axis([0 T_stop -120 170])
% 
% figure
% subplot(211), plot(ts_c, ref_eta(1,:), 'b-', ts_c, y_xi(:,1)', 'r:', 'LineWidth', 1.5)
% ylabel('\fontsize{10} r_{\xi1} vs. y_{\xi1}')
% legend('\fontsize{10} r_{\xi1}','\fontsize{10} y_{\xi1}','Location','NorthWest','Orientation','Horizontal')
% %  axis([0 T_stop -100 100])
% subplot(212), plot(ts_c, ref_eta(2,:), 'b-', ts_c, y_xi(:,2)', 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} r_{\xi2} vs. y_{\xi2}')
% legend('\fontsize{10} r_{\xi2}','\fontsize{10} y_{\xi2}','Location','NorthWest','Orientation','Horizontal')
% % axis([0 T_stop -120 170])

figure
subplot(211), plot(ts_c, ref_c(1,:), 'b-', ts_c, yc(:,1)', 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} r_{c1} vs. y_{c1}')
legend('\fontsize{10} r_{c1}','\fontsize{10} y_{c1}','Location','NorthWest','Orientation','Horizontal')
%  axis([0 T_stop -100 100])
subplot(212), plot(ts_c, ref_c(2,:), 'b-', ts_c, yc(:,2)', 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} r_{c2} vs. y_{c2}')
legend('\fontsize{10} r_{c2}','\fontsize{10} y_{c2}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -120 170])


figure
for i = 1:min(size(ref_c))
subplot(min(size(ref_c))*100+10+i), 
plot(ts_c, ref_c(i,:), 'b-', ts_c, yc(:,i)', 'r:', 'LineWidth', 1.5)
ylabel(['\fontsize{10} r_{c' num2str(i) '} vs. y_{c' num2str(i) '}'])
legend(['\fontsize{10} r_{c' num2str(i) '}'],['\fontsize{10} y_{c' num2str(i) '}'],'Location','NorthWest','Orientation','Horizontal')
% xlabel('\fontsize{10} Time ( sec. )')
end



figure % de
% subplot(211),
plot(ts_c, de(:,1), 'b-', ts_c, de_hat(:,1)', 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} d_{e} vs. d_{e}hat')
legend('\fontsize{10} d_{e}','\fontsize{10} d_{e}hat','Location','NorthWest','Orientation','Horizontal')
% subplot(212), plot(ts_c, de(:,2), 'b-', ts_c, de_hat(:,2)', 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} r_{c2} vs. y_{c2}')
% legend('\fontsize{10} r_{c2}','\fontsize{10} y_{c2}','Location','NorthWest','Orientation','Horizontal')
% figure % de
% plot(ts_c, de(:,2), 'b-', ts_c, de_hat(:,2)', 'r:', 'LineWidth', 1.5)
% 

