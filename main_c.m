
% =========================================================== %%
%                               2016/01/06
%                           Frequency Shaping 
%        ( continuous-time system , PID filter , m greater than p )
% ======================================================================= %

clc; clear all; close all;
format short g
% ******************************************************************** %%
%                       Simulation Parameters Setting                           
% *********************************************************************** %
dt      = 0.01;                 % Simulation interval
T_stop  = 5;                    % Final simulation time
Sim_t = 0:dt:T_stop;
T_length = length(Sim_t);

w_min  = -3;                    % Minimum frequency 10^w_min 
w_max  =  6;                    % Maximum frequency 10^w_max 
pts    = 1e4;                   % points
w = logspace(w_min,w_max,pts);  % Specifies the frequency range

% ******************************************************************** %%
%            Original Plant( non-square , non-minimal phase )
% *********************************************************************** %
load nonminimal_system    % { A , B , C , D}
load De_term

n   = size(A,1);               % No. of states 
m   = size(B,2);               % No. of inputs
p   = size(C,1);               % No. of outputs 
D=zeros(p,m);
% Zeros_ori = TZOCS(A, B, C, D)
Qc = eye(p);
Rc = 1e-6*eye(m);
R_bar = Rc + D'*Qc*D;
N     = C'*Qc*D;
[Kc,~] =  lqr(A, B, C'*Qc*C, R_bar, N);
 Ec     = -inv(R_bar)*( ( C - D*Kc )*inv( A - B*Kc )*B - D )'*Qc;
Zeros_ori = tzero(A-B*Kc, B*Ec, C-D*Kc, D*Ec)
z_p=eig(A-B*Kc)
% ~~~ Check whether the original plant is controllable and observable

Mc_ori = ctrb(A, B);           % Controllability matrix 
Mo_ori = obsv(A, C);           % Observability matrix

Rc_ori = rank(Mc_ori);         % Check whether the plant is controllable (the system is controllable if Mc has full rank n )
Ro_ori = rank(Mo_ori);         % Check whether the plant is observable   (the system is observable   if Mo has full rank n )
%% -----------Step_Response     % Double check whether the plant has unstable zeros (on RHP)
% ref=ones(p,T_length);
% Q = eye(p);  
% R = 1e-6*eye(m);
% 
% [Kcs,P1] = lqr(A, B, C'*Q*C, R);
%  Ecs     = -inv(R)*B'*inv( A - B*Kcs )'*C'*Q;
%  % ~~~ Perform Simulation
% x(:,1)=zeros(n,1);
% opts_sim = simset('Solver', 'ode45', 'FixedStep', dt, 'InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-7);
% [T_s, X_s, Us, Ys] = sim('Original_Plant', Sim_t, opts_sim, [ Sim_t' , ref' ]);
% 
% xs = X_s';             % state
% us = Us';             % control input
% ys = Ys';             % step response
% figure
% plot(Sim_t, ref(1, :),'b-',Sim_t, ys(1, :),'r:.', 'LineWidth', 1.5);
% legend('\fontsize{10} r_{1}','\fontsize{10} y_{c,1}','Location','SouthWest','Orientation','Horizontal')
% xlabel('Time (sec)');ylabel('r_{1} vs. y_{c,1}');
% % axis([0 T_stop  -2.2 1.7])
% figure
% plot(Sim_t, ref(2, :),'b-',Sim_t, ys(2, :),'r:.', 'LineWidth', 1.5);
% legend('\fontsize{10} r_{2}','\fontsize{10} y_{c,2}','Location','SouthWest','Orientation','Horizontal')
% xlabel('Time (sec)');ylabel('r_{2} vs. y_{c,2}');
% % axis([0 T_stop  -1.2 1.2])

%% ******************************************************************** %%
%% Step 1 :transform the non-square  non-minimal phase system to a minimal-phase and observable system 
% *********************************************************************** %
load 'eta_value'                % eta ( trial and error  , randomly generated )
% eta=[0.8937   -1.7223;
%      2.5072    2.5898];

%  flag=0; k=0;
% while(~flag) 
%     eta=(1*randn(p,p));
% %   [eta] = boundchang(eta, 1e0);
% 
%     C_tilde = eta*C;
%     D_tilde = eta*D;
%     Zeros_tilde=sort(TZOCS(A, B, C_tilde, D_tilde));
%     rank_obs_=rank(obsv(A,C_tilde));
%     realzers=real(Zeros_tilde );
%     if(~isempty(realzers))
%     if (sum(realzers>=0)==0)&&(rank_obs_==n)
%         flag=1;
% %         save ('eta value', 'eta')
% 
%     end
%     else
%       k=k+1;
%     end
% end
eta=10*eta;
C_tilde =eta *C;
D_tilde = eta*D;
Zeros_tilde=TZOCS(A, B, C_tilde, D_tilde)

Qc = eye(p);  
Rc = 1e-6*eye(m);
R_bar = Rc + D_tilde'*Qc*D_tilde;
N = C_tilde'*Qc*D_tilde;
[Kc,~] =  lqr(A, B, C_tilde'*Qc*C_tilde, R_bar, N);
Ec = -inv(R_bar)*( ( C_tilde - D_tilde*Kc )*inv( A - B*Kc )*B - D_tilde )'*Qc;
Zeros_tilde = tzero(A-B*Kc, B*Ec, C_tilde-D_tilde*Kc, D_tilde*Ec);
z_p=eig(A-B*Kc);
%%
% ~~~ Check whether the plant is observable

Mo_tilde = obsv(A, C_tilde);  % Observability matrix
Ro_tilde = rank(Mo_tilde);    % Check whether the plant is observable ( the system is observable if Mo has full rank n )

% -----------Step_Response_tilde     % Double check whether the plant has unstable zeros (on RHP)
% ref=ones(p,T_length);
% Q = eye(p);  
% R = 1e-6*eye(m);
% 
% [Kcs,P1] = lqr(A, B, C_tilde'*Q*C_tilde, R);
%  Ecs     = -inv(R)*B'*inv( A - B*Kcs )'*C_tilde'*Q;
%  % ~~~ Perform Simulation
% x(:,1)=zeros(n,1);
% opts_sim = simset('Solver', 'ode45', 'FixedStep', dt, 'InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-7);
% [T_s, X_s, Us, Ys] = sim('Original_Plant_tilde', Sim_t, opts_sim, [ Sim_t' , ref' ]);
% 
% xs = X_s';             % state
% us = Us';             % control input
% ys = Ys';             % step response
% figure
% plot(Sim_t, ref(1, :),'b-',Sim_t, ys(1, :),'r:.', 'LineWidth', 1.5);
% legend('\fontsize{10} r_{1}','\fontsize{10} y_{c,1}','Location','SouthWest','Orientation','Horizontal')
% xlabel('Time (sec)');ylabel('r_{1} vs. y_{c,1}');
% % axis([0 T_stop  -2.2 1.7])
% figure
% plot(Sim_t, ref(2, :),'b-',Sim_t, ys(2, :),'r:.', 'LineWidth', 1.5);
% legend('\fontsize{10} r_{2}','\fontsize{10} y_{c,2}','Location','SouthWest','Orientation','Horizontal')
% xlabel('Time (sec)');ylabel('r_{2} vs. y_{c,2}');
% % axis([0 T_stop  -1.2 1.2])

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
% flag=0;
% while(~flag)
% A_r=randn(n,n);
% B_r=randn(n,p);
% C_r=randn(p,n);
% D_r=zeros(p);         % Without a direct feed-through term
% [A_r] = boundchang(A_r, 1e0);
% [B_r] = boundchang(B_r, 1e0);
% [C_r] = boundchang(C_r, 1e0);

% [D_r] = boundchang(D_r, 1e0);
% 
% rank_obs_r=rank(obsv(Ar,Cr));

% rank_ctr_r=rank(ctrb(Ar,Br));
% pole_r=eig(Ar);
% zero_r = tzero(Ar,Br,Cr,Dr);
% realzer_r=real(zero_r );
% if (sum(realzer_r>=0)==0)&&(rank_obs_r==n)&&(rank_ctr_r==n)
%    flag=1;
% %    save('r(t) model', 'A_r','B_r','C_r','D_r')
% end
% end

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

% ~~~ Simulation Results
% 
% figure
% subplot(211), plot(ts_c, ref_eta(1,:), 'b-', ts_c, y_r(1,:), 'r:', 'LineWidth', 1.5)
% legend('\fontsize{10} r_{\xi,1}','\fontsize{10} y_{r,1}','Location','SouthWest','Orientation','Horizontal')
%  axis([0 T_stop -100 100])
% subplot(212), plot(ts_c, ref_eta(2,:), 'b-', ts_c, y_r(2,:), 'r:', 'LineWidth', 1.5)
% legend('\fontsize{10} r_{\xi,2}','\fontsize{10} y_{r,2}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -120 170])
% xlabel('\fontsize{10} Time ( sec. )')
% 
% figure
% subplot(211), plot(ts_c, u_r(1,:), 'b-', 'LineWidth', 1.5)
% legend('\fontsize{10} u_{r,1}')
% axis([0 T_stop -0.7e4 0.5e4])
% subplot(212), plot(ts_c, u_r(2,:), 'b-', 'LineWidth', 1.5)
% legend('\fontsize{10} u_{r,2}')
% xlabel('\fontsize{10} Time ( sec. )')
% axis([0 T_stop -2.5e4 0.5e4])

%% Traditional LQAT
Q = eye(p);  
R = 1e-6*eye(m);

[K,P1] = lqr(A, B, C'*Q*C, R);
 E     = -inv(R)*B'*inv(( A - B*K )')*C'*Q;
% %  % ~~~ Perform Simulation
% x(:,1)=xc(:,1);
% opts_sim = simset('Solver', 'ode45', 'FixedStep', dt, 'InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-7);
% [T_tr, X_tr, U_tr, Y_tr] = sim('model_plant', ts_c, opts_sim, [ ts_c' , ref_c' ]);
% % 
% x_tr = X_tr';      
% % state
% u_tr = U_tr';             % control input
% y_tr = Y_tr';             % step response
% figure
% subplot 211
% plot(ts_c, ref_c(1, :),'b-',ts_c, y_tr(1, :),'r:', 'LineWidth', 1.5);
% legend('\fontsize{10} r_{c,1}','\fontsize{12} y_{c,1}','Location','SouthWest','Orientation','Horizontal')
% xlabel('Time (sec)');ylabel('r_{c,1} vs. y_{c,1}');
% axis([0 T_stop -22 20])
% 
% subplot 212
% plot(ts_c, ref_c(2, :),'b-',ts_c, y_tr(2, :),'r:', 'LineWidth', 1.5);
% legend('\fontsize{10} r_{c,2}','\fontsize{12} y_{c,2}','Location','SouthWest','Orientation','Horizontal')
% xlabel('Time (sec)');ylabel('r_{c,2} vs. y_{c,2}');
% axis([0 T_stop  -12 12])
% % 
% figure
% subplot(311), plot(ts_c, u_tr(1,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,1}')
% legend('\fontsize{10} u_{c,1}','Location','SouthWest')
%  axis([0 T_stop -500 500])
% subplot(312), plot(ts_c, u_tr(2,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,2}')
% legend('\fontsize{10} u_{c,2}','Location','NorthWest')
%  axis([0 T_stop -500 500])
% subplot(313), plot(ts_c, u_tr(3,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,3}')
% legend('\fontsize{10} u_{c,3}','Location','SouthWest')
%  axis([0 T_stop -500 500])
% xlabel('\fontsize{10} Time ( sec. )')
% 
% figure
% subplot 311
%  plot(ts_c, x_tr (1,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,1} vs. x_{e,1}')
% legend('\fontsize{10} x_{c,1}','\fontsize{10} x_{e,1}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -22 20])
% subplot 312
% plot(ts_c, x_tr (3,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,3} vs. x_{e,3}')
% legend('\fontsize{10} x_{c,3}','\fontsize{10} x_{e,3}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -dt dt])
% subplot 313
%  plot(ts_c, x_tr (5,:), 'b-','LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,5} vs. x_{e,5}')
% legend('\fontsize{10} x_{c,5}','\fontsize{10} x_{e,5}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -20 20])
% figure
% subplot 311
%  plot(ts_c, x_tr (2,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,1} vs. x_{e,1}')
% legend('\fontsize{10} x_{c,1}','\fontsize{10} x_{e,1}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -22 20])
% subplot 312
% plot(ts_c, x_tr (4,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,3} vs. x_{e,3}')
% legend('\fontsize{10} x_{c,3}','\fontsize{10} x_{e,3}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -dt dt])
% subplot 313
%  plot(ts_c, x_tr (6,:), 'b-','LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,5} vs. x_{e,5}')
% legend('\fontsize{10} x_{c,5}','\fontsize{10} x_{e,5}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -20 20])
% 
eig(A-B*K)
%% Traditional Observer
% 
% Ro = 1e0*eye(p);
% Qo = 1e6*eye(n);
% [Ko,Po,Eo] = lqr(A',C',Qo,Ro);
% Lo = Po*C'*inv(Ro);
% 
% [To,Xo,U_otr,Y_otr,X_otr,Y_O,X_O] = sim('Traditional_Obsver', ts_c, opts_sim, [ ts_c' , ref_c']);   
% u_otr=U_otr';
% y_otr=Y_otr';
% x_otr=X_otr';
% y_etr=Y_O';
% x_etr=X_O';
% figure
% plot(ts_c, ref_c(1,:), 'b-', ts_c, y_etr(1,:), 'r:', 'LineWidth', 1.5)
% ylabel('\fontsize{10} r_{c,1} vs. y_{c,1}')
% legend('\fontsize{10} r_{c,1}','\fontsize{10} y_{e,1}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -22 20])
% figure
% plot(ts_c, ref_c(2,:), 'b-', ts_c, y_etr(2,:), 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} r_{c,2} vs. y_{c,2}')
% legend('\fontsize{10} r_{c,2}','\fontsize{10} y_{e,2}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -12 12])
% 
% figure
% subplot 211
%  plot(ts_c, ref_c(1,:), 'b-', ts_c, y_otr(1,:), 'r:', 'LineWidth', 1.5)
% ylabel('\fontsize{10} r_{c,1} vs. y_{c,1}')
% legend('\fontsize{10} r_{c,1}','\fontsize{10} y_{c,1}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -25 17])
% subplot 212; plot(ts_c, ref_c(2,:), 'b-', ts_c, y_otr(2,:), 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} r_{c,2} vs. y_{c,2}')
% legend('\fontsize{10} r_{c,2}','\fontsize{10} y_{c,2}','Location','SouthWest','Orientation','Horizontal')
%  axis([0 T_stop -20 12])
% 
% figure
% subplot(211), plot(ts_c, y_otr(1,:), 'b-', ts_c, y_etr(1,:), 'r:', 'LineWidth', 1.5)
% ylabel('\fontsize{10} y_{c,1} vs. y_{e,1}')
% legend('\fontsize{10} y_{c,1}','\fontsize{10} y_{e,1}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop ])
% subplot(212), plot(ts_c, y_otr(2,:), 'b-', ts_c, y_etr(2,:), 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} y_{c,2} vs. y_{e,2}')
% legend('\fontsize{10} y_{c,2}','\fontsize{10} y_{e,2}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop ])
% 
% figure
% subplot(311), plot(ts_c, u_otr(1,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,1}')
% legend('\fontsize{10} u_{c,1}','Location','SouthWest')
% %  axis([0 T_stop -5000 5000])
% subplot(312), plot(ts_c, u_otr(2,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,2}')
% legend('\fontsize{10} u_{c,2}','Location','NorthWest')
% %  axis([0 T_stop -5000 5000])
% subplot(313), plot(ts_c, u_otr(3,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,3}')
% legend('\fontsize{10} u_{c,3}','Location','SouthWest')
% %  axis([0 T_stop -5000 5000])
% xlabel('\fontsize{10} Time ( sec. )')
% 
% figure
% subplot 311
%  plot(ts_c, x_otr(1,:), 'b-', ts_c, x_etr(1,:), 'r:', 'LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,1} vs. x_{e,1}')
% legend('\fontsize{10} x_{c,1}','\fontsize{10} x_{e,1}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -22 20])
% subplot 312
% plot(ts_c, x_otr(3,:), 'b-', ts_c, x_etr(3,:), 'r:', 'LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,3} vs. x_{e,3}')
% legend('\fontsize{10} x_{c,3}','\fontsize{10} x_{e,3}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -dt dt])
% subplot 313
%  plot(ts_c, x_otr(5,:), 'b-', ts_c, x_etr(5,:), 'r:', 'LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,5} vs. x_{e,5}')
% legend('\fontsize{10} x_{c,5}','\fontsize{10} x_{e,5}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -20 20])
% 
% figure, 
% subplot 311
% plot(ts_c, x_otr(2,:), 'b-', ts_c, x_etr(2,:), 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} x_{c,2} vs. x_{e,2}')
% legend('\fontsize{10} x_{c,2}','\fontsize{10} x_{e,2}','Location','Southeast','Orientation','Horizontal')
% subplot 312
% plot(ts_c, x_otr(4,:), 'b-', ts_c, x_etr(4,:), 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} x_{c,4} vs. x_{e,4}')
% legend('\fontsize{10} x_{c,4}','\fontsize{10} x_{e,4}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -35 40])
% subplot 313
%  plot(ts_c, x_otr(6,:), 'b-', ts_c, x_etr(6,:), 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} x_{c,6} vs. x_{e,6}')
% legend('\fontsize{10} x_{c,6}','\fontsize{10} x_{e,6}','Location','SouthWest','Orientation','Horizontal')
% 
% figure
% subplot(311), plot(ts_c, u_otr(1,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,1}')
% legend('\fontsize{10} u_{c,1}','Location','SouthWest')
% %  axis([0 T_stop -2 2])
% subplot(312), plot(ts_c, u_otr(2,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,2}')
% legend('\fontsize{10} u_{c,2}','Location','NorthWest')
% % axis([0 T_stop -2 2])
% subplot(313), plot(ts_c, u_otr(3,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,3}')
% legend('\fontsize{10} u_{c,3}','Location','SouthWest')
% %  axis([0 T_stop -2 2])
% xlabel('\fontsize{10} Time ( sec. )')
% 

%% Traditional LQAT_C_tilde
% Q = eye(p);  
% R = 1e-8*eye(m);
% 
% [K,P1] = lqr(A, B, C_tilde'*Q*C_tilde, R);
%  E     = -inv(R)*B'*inv(( A - B*K )')*C_tilde'*Q;
%  % ~~~ Perform Simulation
% x(:,1)=xc(:,1);
% opts_sim = simset('Solver', 'ode45', 'FixedStep', dt, 'InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-7);
% [T_tr, X_tr, U_tr, Y_tr] = sim('model_plant_tilde', ts_c, opts_sim, [ ts_c' , ref_eta' ]);
% 
% x_tr = X_tr';      
% % state
% u_tr = U_tr';             % control input
% y_tr = Y_tr';             % step response
% figure
% subplot 211
% plot(ts_c, ref_eta(1, :),'b-',ts_c, y_tr(1, :),'r:', 'LineWidth', 1.5);
% legend('\fontsize{10} r_{c,1}','\fontsize{12} y_{c,1}','Location','SouthWest','Orientation','Horizontal')
% xlabel('Time (sec)');ylabel('r_{c,1} vs. y_{c,1}');
% % axis([0 T_stop -22 20])
% 
% subplot 212
% plot(ts_c, ref_eta(2, :),'b-',ts_c, y_tr(2, :),'r:', 'LineWidth', 1.5);
% legend('\fontsize{10} r_{c,2}','\fontsize{12} y_{c,2}','Location','SouthWest','Orientation','Horizontal')
% xlabel('Time (sec)');ylabel('r_{c,2} vs. y_{c,2}');
% % axis([0 T_stop  -12 12])
% % 
% figure
% subplot(311), plot(ts_c, u_tr(1,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,1}')
% legend('\fontsize{10} u_{c,1}','Location','SouthWest')
% %  axis([0 T_stop -2 2])
% subplot(312), plot(ts_c, u_tr(2,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,2}')
% legend('\fontsize{10} u_{c,2}','Location','NorthWest')
% % axis([0 T_stop -2 2])
% subplot(313), plot(ts_c, u_tr(3,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} u_{c,3}')
% legend('\fontsize{10} u_{c,3}','Location','SouthWest')
% %  axis([0 T_stop -2 2])
% xlabel('\fontsize{10} Time ( sec. )')
% 
% figure
% subplot 311
%  plot(ts_c, x_tr (1,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,1} vs. x_{e,1}')
% legend('\fontsize{10} x_{c,1}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -22 20])
% subplot 312
% plot(ts_c, x_tr (3,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,3} vs. x_{e,3}')
% legend('\fontsize{10} x_{c,3}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -dt dt])
% subplot 313
%  plot(ts_c, x_tr (5,:), 'b-','LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,5} vs. x_{e,5}')
% legend('\fontsize{10} x_{c,5}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -20 20])
% figure
% subplot 311
%  plot(ts_c, x_tr (2,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,2} vs. x_{e,1}')
% legend('\fontsize{10} x_{c,2}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -22 20])
% subplot 312
% plot(ts_c, x_tr (4,:), 'b-', 'LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,4} vs. x_{e,3}')
% legend('\fontsize{10} x_{c,4}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -dt dt])
% subplot 313
%  plot(ts_c, x_tr (6,:), 'b-','LineWidth', 1.5)
% ylabel('\fontsize{10} x_{c,6} vs. x_{e,5}')
% legend('\fontsize{10} x_{c,6}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -20 20])
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

% Kp_1 = 50;                % Proportional gain         
% Kp_2 = 50;     
% 
% Ki_1 = 500;                 % Integral     gain      
% Ki_2 = 500;
% 
% Kd_1 = 0.0001;            % Derivative   gain      
% Kd_2 = 0.0001;
% 
% Kp = diag([ Kp_1 Kp_2 ]); % Proportional gain matrix
% Ki = diag([ Ki_1 Ki_2 ]); % Integral     gain matrix
% Kd = diag([ Kd_1 Kd_2 ]); % Derivative   gain matrix
% 

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

% ******************************************************************** %%
%                             Perform Simulation                           
% *********************************************************************** %
% ref_f=zeros(p,T_length);
[Tc, Xc, Uc, Yf] = sim('Augmented_Plant', ts_c, opts_sim, [ ts_c' , ref_eta' , Cc' , d_au' , s_au' ]);
   
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
subplot(211), plot(ts_c, ref_eta(1,:), 'b-', ts_c, y_f(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} r_{\xi,1} vs. y_{f,1}')
legend('\fontsize{10} r_{\xi,1}','\fontsize{10} y_{f,1}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -100 100])
subplot(212), plot(ts_c, ref_eta(2,:), 'b-', ts_c, y_f(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} r_{\xi,2} vs. y_{f,2}')
legend('\fontsize{10} r_{\xi,2}','\fontsize{10} y_{f,2}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -120 170])

figure
subplot(211), plot(ts_c, ref_eta(1,:), 'b-', ts_c, y_eta(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} r_{\xi,1} vs. y_{\xi,1}')
legend('\fontsize{10} r_{\xi,1}','\fontsize{10} y_{\xi,1}','Location','NorthWest','Orientation','Horizontal')
%  axis([0 T_stop -100 100])
subplot(212), plot(ts_c, ref_eta(2,:), 'b-', ts_c, y_eta(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} r_{\xi,2} vs. y_{\xi,2}')
legend('\fontsize{10} r_{\xi,2}','\fontsize{10} y_{\xi,2}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -120 170])

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
figure
subplot(311), plot(ts_c, u_c(1,:), 'b-', 'LineWidth', 1.5)
ylabel('\fontsize{10} u_{c,1}')
% legend('\fontsize{10} u_{c,1}','Location','Northeast','Orientation','Horizontal')
%  axis([0 T_stop -1.5e4 1.5e4])
subplot(312), plot(ts_c, u_c(2,:), 'b-', 'LineWidth', 1.5)
ylabel('\fontsize{10} u_{c,2}')
% legend('\fontsize{10} u_{c,2}','Location','Southeast','Orientation','Horizontal')
% axis([0 T_stop -1.5e4 1.7e4])
subplot(313), plot(ts_c, u_c(3,:), 'b-', 'LineWidth', 1.5)
ylabel('\fontsize{10} u_{c,3}')
% legend('\fontsize{10} u_{c,3}','Location','Southeast','Orientation','Horizontal')
% axis([0 T_stop -1.7e4 1.7e4])
xlabel('\fontsize{10} Time ( sec. )')

% figure
% % subplot(211), plot(ts_c, e_tra(1,:), 'b-', 'LineWidth', 1.5)
% % legend('\fontsize{10} r_{1,1} -  y_{2,1}','Location','SouthEast')
% % axis([0 T_stop -0.01 0.01])
% % subplot(212), plot(ts_c, e_tra(2,:), 'b-', 'LineWidth', 1.5)
% % legend('\fontsize{10} r_{1,2} -  y_{2,2}','Location','SouthEast')
% % axis([0 T_stop -0.01 0.01])
% % xlabel('\fontsize{10} Time ( sec. )')
% 
figure
subplot(211), plot(ts_c, x_eta(1,:), 'b-', 'LineWidth', 1.5)
legend('\fontsize{10} x_{\xi,1}')
% axis([0 T_stop -0.1 0.1])
subplot(212), plot(ts_c, x_eta(2,:), 'b-', 'LineWidth', 1.5)
legend('\fontsize{10} x_{\xi,2}')
xlabel('\fontsize{10} Time ( sec. )')
% axis([0 T_stop -0.1 0.1])

figure
subplot(211), plot(ts_c, x_c(1,:), 'b-', 'LineWidth', 1.5)
legend('\fontsize{10} x_{c,1}')
% axis([0 T_stop -0.1 0.1])
subplot(212), plot(ts_c, x_c(2,:), 'b-', 'LineWidth', 1.5)
legend('\fontsize{10} x_{c,2}')
xlabel('\fontsize{10} Time ( sec. )')
% axis([0 T_stop -0.1 0.1])



% ******************************************************************** %%
%                             Pole-Zero Map                           
% *********************************************************************** %

% ~~~ Original plant ( non-square , non-minimal phase )

Poles_ori = eig(A);            % the poles of the original plant
% tzero(A, B, C, D)            % the transmission zeros of the original plant
Zeros_ori = TZOCS(A, B, C, D);


Sys_ori  = ss(A, B, C, D);

Re_p_ori = real(Poles_ori);
Im_p_ori = imag(Poles_ori);

Re_z_ori = real(Zeros_ori);
Im_z_ori = imag(Zeros_ori);

% figure
% % pzmap(Sys_ori)
% plot(Re_p_ori, Im_p_ori, 'LineWidth', 2, 'MarkerSize', 9, 'LineStyle', 'none', 'Marker', 'x', 'Color' , [0 0 1]), hold on
% plot(Re_z_ori, Im_z_ori, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle', 'none', 'Marker', 'o', 'Color' , [0 0 1]), hold off
% axis([ -6 11 -3 3 ])
% xlabel('\fontsize{10} Real')
% ylabel('\fontsize{10} Imag.')
% title('')
% 
% ~~~ Squaring the Plant ( square , minimal phase , observable)  

% Zeros_tilde = TZOCS(A, B, C_tilde, D_tilde);

Sys_tilde  = ss(A, B, C_tilde, D_tilde);

Re_p_tilde = real(Poles_ori);
Im_p_tilde = imag(Poles_ori);

Re_z_tilde = real(Zeros_tilde);
Im_z_tilde = imag(Zeros_tilde);

% figure
% % pzmap(Sys_tilde)
% plot(Re_p_tilde, Im_p_tilde, 'LineWidth', 2, 'MarkerSize', 9, 'LineStyle', 'none', 'Marker', 'x', 'Color' , [0 0 1]), hold on
% plot(Re_z_tilde, Im_z_tilde, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle', 'none', 'Marker', 'o', 'Color' , [0 0 1]), hold off
% axis([ -6 11 -3 3 ])
% xlabel('\fontsize{10} Real')
% ylabel('\fontsize{10} Imag.')
% title('')

% ~~~ Augmented plant ( square , minimal phase , controllable , observable)  
Zeros_aug = [ roots([ Kd_1 , Kp_1 , Ki_1 ]) ; roots([ Kd_2 , Kp_2 , Ki_2 ]) ; Zeros_tilde ];
Poles_aug = eig(A_aug);                % the poles of the augmented plant
Sys_aug  = ss(A_aug, B_aug, C_aug, D_aug);

Re_p_aug = real(Poles_aug);
Im_p_aug = imag(Poles_aug);

Re_z_aug = real(Zeros_aug);
Im_z_aug = imag(Zeros_aug);

% figure
% % pzmap(Sys_aug)
% plot(Re_p_aug, Im_p_aug, 'LineWidth', 2, 'MarkerSize', 9, 'LineStyle', 'none', 'Marker', 'x', 'Color' , [0 0 1]), hold on
% plot(Re_z_aug, Im_z_aug, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle', 'none', 'Marker', 'o', 'Color' , [0 0 1]), hold off
% % axis([ -1.1e5 2e4 -3 3 ])
% axis([ -6 6 -3 3 ])
% xlabel('\fontsize{10} Real')
% ylabel('\fontsize{10} Imag.')
% title('')

% ~~~ Closed-loop plant ( the closed-loop poles approach to open-loop zeros )

A_cl = A_aug - B_aug*Kc;
B_cl = B_aug*Ec;
C_cl = C_aug - D_aug*Kc;
D_cl = D_aug*Ec;

Sys_cl = ss(A_cl, B_cl, C_cl, D_cl);
eig(A_aug - B_aug*Kc)
% TZOCS(A_aug,B_aug,C_aug,D_aug)

Qc = eye(p);  
Rc = 1e-5*eye(m);
R_bar = Rc + D_aug'*Qc*D_aug;
N     = C_aug'*Qc*D_aug;
[K1,P1] =  lqr(A_aug, B_aug, C_aug'*Qc*C_aug, R_bar, N);
 E1     = -inv(R_bar)*( ( C_aug - D_aug*K1 )*inv( A_aug - B_aug*K1 )*B_aug - D_aug )'*Qc;
Z = tzero(A_aug-B_aug*K1, B_aug*E1, C_aug-D_aug*K1, D_aug*E1);
z_p=eig(A_aug-B_aug*K1)

% figure
% pzmap(Sys_cl)
% axis([ -6 1 -1 1 ])
% xlabel('\fontsize{10} Real')
% ylabel('\fontsize{10} Imag.')
% title('')

s   = tf('s');

LX = ss(A, B, Kc1, zeros(m,m)) + ss(A, B, Kc2*C_tilde, zeros(m,m))/s;
LY = ss(A - B*Kc1, B*Kc2, C_tilde, zeros(p,p))/s;
% % 
% % LYe = ss(A, De, Ke1, zeros(p,p)) + ss(A, De, Ke2*C_tilde, zeros(p,p))/s;
% % LXe = ss(A - De*Ke1, De*Ke2, C_tilde, zeros(p,p))/s;
% % 
% %  
% % 
% % 
% % figure
% % sigma(LXe,'b-')
% % xlabel('\fontsize{10} Freq.')
% % ylabel('\fontsize{10} Singular Values')
% % legend('\fontsize{10} Uncertainty loop')
% % axis([ 1e-6 1e6 -400 100])
% % title('')
% % 
% % figure
% % sigma(LYe,'b-')
% % xlabel('\fontsize{10} Freq.')
% % ylabel('\fontsize{10} Singular Values')
% % legend('\fontsize{10} Sensor estimator loop')
% % axis([ 1e-6 1e6 -50 160])
% % title('')
% % 
figure
sigma(LX,'b-')
xlabel('\fontsize{10} Freq.')
ylabel('\fontsize{10} Singular Values')
legend('\fontsize{10} Control loop')
axis([ 1e-6 1e6 -200 300])
title('')
plotboundary(0, 'y', '-.', 'k')

figure
sigma(LY,'b-')
xlabel('\fontsize{10} Freq.')
ylabel('\fontsize{10} Singular Values')
legend('\fontsize{10} Reference loop')
axis([ 1e-6 1e6 -200 300])
title('')
plotboundary(0, 'y', '-.', 'k')
% % %% ******************************************************************** %%
% % %                           Singular Value Plots                           
% % % *********************************************************************** %
figure
sigma(Sys_ori, 'r:', Sys_aug, 'b-')
legend('\fontsize{10} G_{ol}\rm (s)','\fontsize{10} G_{aug}\rm (s)')
legend('\fontsize{10} Original plant','\fontsize{10} Augmented plant')
xlabel('\fontsize{10} Freq.')
ylabel('\fontsize{10} Maximum and Minimum Singular Values')
title('')
axis([ 1e-6 1e6 -200 300])
plotboundary(0, 'y', '-.', 'k')
%% *****************ESTIMATOR********************************************
% *************************************************************************
%---------------------------PID Observer---------------------------------------------

%%
Ip = [ zeros(p,n)  eye(p) ];
Ip3= [ eye(n,n)  zeros(n,p) ];
Ip2 = [ zeros(n,p) ; eye(p)];
C_bar = [ C_tilde  zeros(p,p)];  

%%
% flag=0; k=0;
% while(~flag)
% De = randn(n,p);
%  zero_m_e = tzero(A, De, C_tilde, zeros(p))
%     realzers=real(zero_m_e );
%     imagzers = imag(zero_m_e);
%     if((imagzers==0))
%     if(~isempty(realzers))
%         if (sum(realzers>=0)==0)
%         flag=1;
% %         save('De_','De')
%         break
%         end
%     else
%           k=k+1;
%     end
%     end
% end
% ~~~ Check whether the plant is observable

Rc_obs = rank(ctrb(A, De));         % Check whether the plant is controllable (the system is controllable if Mc has full rank n )
Ro_obs = rank(obsv(A, C_tilde));         % Check whether the plant is observable
 

%% ---------------PID-filter Observer-------------------

% Kep_1 = 25;                % Proportional gain       
% Kep_2 = 25;     
% 
% Kei_1 = 50;                 % Integral     gain   
% Kei_2 = 50;
% 
% Ked_1 = 0.005;            % Derivative   gain       
% Ked_2 = 0.005;
% Kep = diag([ Kep_1 Kep_2 ]); % Proportional gain matrix
% Kei = diag([ Kei_1 Kei_2 ]); % Integral     gain matrix
% Ked = diag([ Ked_1 Ked_2 ]); % Derivative   gain matrix


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

[Ke,Pe_aug] =  lqr(Ae_aug, Be_aug, Ce_aug'*Qe*Ce_aug, Re);
Ee      = -inv(Re)*( ( Ce_aug )*inv( Ae_aug - Be_aug*Ke )*Be_aug )'*Qe;
Ke1 = Ke(:,1:n);
Ke2 = Ke(:,n+1:end);
i = 0;
%
eig(A-De*Ke1)

for ts = 0 : dt : T_stop
    i = i + 1;
    
    de_aug_e(:,i) = [  De*inv(eye(p)+Ked*eta*C_tilde*De)*Ked*C_r*( A_r*x_r(:,i) + B_r*u_r(:,i)) ; zeros(p,1)]; %errrrrrrorrrrr 多一個Xi
end
C_a=[De*inv(eye(p)+Ked*eta*C_tilde*De)*Ked*C_tilde*B; zeros(p,m)]; %errrrrrrorrrrr 要加負號

Ce=inv(Re)*Be_aug'*inv((Ae_aug-Be_aug*Ke)')*Pe_aug;
%% ******************************************************************** %%
%             test of Estimator Design                          
% *********************************************************************** %

% [Te,X,Xe_aug,Ye_aug,V2,X_e, Y_e] = sim('Augmented_Observer_',ts_c,opts_sim,[ts_c' , y_eta' , de_aug_e']);   
% 
% xe_aug = Xe_aug';
% ye_aug = Ye_aug';
% v2 = V2';
% x_e=X_e';
% y_e=Y_e';
% %
% figure
% subplot(211), plot(ts_c, ref_eta(1,:), 'b-', ts_c, ye_aug(1,:), 'r:', 'LineWidth', 1.5)
% ylabel('\fontsize{10} r_{\xi,1} vs. y_{\xie,1}')
% legend('\fontsize{10} r_{\xi,1}','\fontsize{10} y_{\xie,1}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -2.2 1.7])
% subplot(212), plot(ts_c, ref_eta(2,:), 'b-', ts_c, ye_aug(2,:), 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} r_{\xi,2} vs. y_{\xie,2}')
% legend('\fontsize{10} r_{\xi,2}','\fontsize{10} y_{\xie,2}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -1.2 1.2])
% 
% figure
% subplot(211), plot(ts_c, y_eta(1,:), 'b-', ts_c, ye_aug(1,:), 'r:', 'LineWidth', 1.5)
% ylabel('\fontsize{10} y_{\xi,1} vs. y_{\xie,1}')
% legend('\fontsize{10} y_{\xi,1}','\fontsize{10} y_{\xie,1}','Location','SouthWest','Orientation','Horizontal')
% % axis([0 T_stop -2.2 1.7])
% subplot(212), plot(ts_c, y_eta(2,:), 'b-', ts_c, ye_aug(2,:), 'r:', 'LineWidth', 1.5)
% xlabel('\fontsize{10} Time ( sec. )')
% ylabel('\fontsize{10} y_{\xi,2} vs. y_{\xie,2}')
% legend('\fontsize{10} y_{\xi,2}','\fontsize{10} y_{\xie,2}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -1.2 1.2])

%% By ---Elca--- UIE
h = 0; % shift the y-axis
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
Bfi = 8e2;
Cfi = 1.25e2;

Af = -Bfi*Cfi*eye(size(B,2));
Bf = Bfi*eye(size(B,2));
Cf = Cfi*eye(size(B,2));

de = disturb2(Sim_t);
% set(0,'ShowHiddenHandles','on'); set(gcf,'menubar','figure')
 %% PI State Estimate Tarcker Design 

[Te,X, Ue,X_aug ,Yef,Y_eta,Xe_aug,Ye_aug,Xe,Ye_eta] = sim('State_Estimate_Tarcker_IEEE_disturb', ts_c, opts_sim, [ ts_c' , ref_eta' , Cc' , d_au' , s_au' , de_aug_e'  ]);   

ue_c = Ue';            % control input
xe_c = X_aug(:,1:n)';     % states (xc)
yc= C*xe_c;          % response
% 
% xe_1 = Xe1(:,n+1:end)'; % states (x_eta)
% y_1 = C_tilde*x_c;    % response (y_eta)
y_eta = Y_eta';             % y_eta
y_eta_e =Ye_eta';
ye_f = Yef';             % yf
xe_aug = Xe_aug';
ye_aug = Ye_aug';
xe = Xe';
ye=C*xe;
                %y_eta_e
% y_error = Y_error';
% e_tra2 = ref_eta - y_2;  % tracking error

%                            Simulation Results                           
% *********************************************************************** %
close all;
figure
eyy = yc-ye;
subplot(211), plot(ts_c, eyy(1,:), 'b-', 'LineWidth', 1.5)
% % axis([0 T_stop -100 100])
subplot(212), plot(ts_c, eyy(2,:), 'b-', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')



figure
subplot(211), plot(ts_c, y_eta(1,:), 'b-', ts_c, y_eta_e(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} y_{\xi,1} vs. y_{\xi,1}')
legend('\fontsize{10} y_{\xi,1}','\fontsize{10} y_{\xie,1}','Location','NorthWest','Orientation','Horizontal')
% % axis([0 T_stop -100 100])
subplot(212), plot(ts_c, y_eta(2,:), 'b-', ts_c, y_eta_e(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} y_{\xi,2} vs. y_{\xi,2}')
legend('\fontsize{10} y_{\xi,2}','\fontsize{10} y_{\xie,2}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -120 170])

figure
subplot(211), plot(ts_c, ref_eta(1,:), 'b-', ts_c, y_eta_e(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} r_{\xi,1} vs. y_{\xi,1}')
legend('\fontsize{10} r_{\xi,1}','\fontsize{10} y_{\xi,1}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -100 100])
subplot(212), plot(ts_c, ref_eta(2,:), 'b-', ts_c, y_eta_e(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} r_{\xi,2} vs. y_{\xi,2}')
legend('\fontsize{10} r_{\xi,2}','\fontsize{10} y_{\xi,2}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -120 170])

figure
subplot(211), plot(ts_c, yc(1,:), 'b-', ts_c, ye(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} y_{c,1} vs. y_{e,1}')
legend('\fontsize{10} y_{c,1}','\fontsize{10} y_{e,1}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -22 20])
subplot(212), plot(ts_c, yc(2,:), 'b-', ts_c, ye(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} y_{c,2} vs. y_{e,2}')
legend('\fontsize{10} y_{c,2}','\fontsize{10} y_{e,2}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -12 12])


figure
subplot(211), plot(ts_c, ref_c(1,:), 'b-', ts_c, yc(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} r_{c,1} vs. y_{e,1}')
legend('\fontsize{10} r_{c,1}','\fontsize{10} y_{c,1}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -22 20])
subplot(212), plot(ts_c, ref_c(2,:), 'b-', ts_c, yc(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} r_{c,2} vs. y_{e,2}')
legend('\fontsize{10} r_{c,2}','\fontsize{10} y_{c,2}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -12 12])

figure
subplot 311
 plot(ts_c, xe_c(1,:), 'b-', ts_c, xe(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} x_{c,1} vs. x_{e,1}')
legend('\fontsize{10} x_{c,1}','\fontsize{10} x_{e,1}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -1.5e3 2e3])
subplot 312
plot(ts_c, xe_c(3,:), 'b-', ts_c, xe(3,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} x_{c,3} vs. x_{e,3}')
legend('\fontsize{10} x_{c,3}','\fontsize{10} x_{e,3}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -6e3 6e3])
subplot 313
 plot(ts_c, xe_c(5,:), 'b-', ts_c, xe(5,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} x_{c,5} vs. x_{e,5}')
legend('\fontsize{10} x_{c,5}','\fontsize{10} x_{e,5}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -1.2e3 1.5e3])

figure, 
subplot 311
plot(ts_c, xe_c(2,:), 'b-', ts_c, xe(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} x_{c,2} vs. x_{e,2}')
legend('\fontsize{10} x_{c,2}','\fontsize{10} x_{e,2}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -8e2 8e2])
subplot 312
plot(ts_c, xe_c(4,:), 'b-', ts_c, xe(4,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} x_{c,4} vs. x_{e,4}')
legend('\fontsize{10} x_{c,4}','\fontsize{10} x_{e,4}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -3e3 3e3])
subplot 313
 plot(ts_c, xe_c(6,:), 'b-', ts_c, xe(6,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} x_{c,6} vs. x_{e,6}')
legend('\fontsize{10} x_{c,6}','\fontsize{10} x_{e,6}','Location','NorthWest','Orientation','Horizontal')
% axis([0 T_stop -2e3 3e3])


figure
subplot(311), plot(ts_c, ue_c(1,:), 'b-', 'LineWidth', 1.5)
ylabel('\fontsize{10} u_{c,1}')
% legend('\fontsize{10} u_{c,1}','Location','SouthWest')
% axis([0 T_stop -1.22e4 1.2e4])
subplot(312), plot(ts_c, ue_c(2,:), 'b-', 'LineWidth', 1.5)
ylabel('\fontsize{10} u_{c,2}')
% legend('\fontsize{10} u_{c,2}','Location','NorthWest')
% axis([0 T_stop -2e4 2e4])
subplot(313), plot(ts_c, ue_c(3,:), 'b-', 'LineWidth', 1.5)
ylabel('\fontsize{10} u_{c,3}')
% legend('\fontsize{10} u_{c,3}','Location','SouthWest')
% axis([0 T_stop -2e4 2e4])
xlabel('\fontsize{10} Time ( sec. )')
% 
 %% ************************************Robustness**************************
% % % ************************************************************************
per1 = [mean(ref_c(1,:)) 0 ; 0 mean(ref_c(2,:))]*ones(p,T_length)+[std(ref_c(1,:)) 0;0 std(ref_c(2,:))]*randn(p,T_length);
per=eta*per1;
mean(per1(1,:))
mean(per1(2,:))
std(per1(1,:))
std(per1(2,:))
% per = [mean(ref_eta(1,:)) 0 ; 0 mean(ref_eta(2,:))]*ones(p,T_length)+[std(ref_eta(1,:)) 0;0 std(ref_eta(2,:))]*randn(p,T_length);

[Te,X_p, Ue_p,X_aug_p, Yef_p,Y_eta_p,Xe_aug_p,Ye_aug_p,Xe_p,Ye_eta_p] = sim('State_Estimate_Tarcker_p', ts_c, opts_sim, [ ts_c' , ref_eta' , Cc' , d_au' , s_au' , de_aug_e'  , per' ]);   
ue_c_p = Ue_p';            % control input

xe_c_p = X_aug_p(:,1:n)';     % states (xc)
ye_c_p = C*xe_c_p;          % response

ye_eta_p = Y_eta_p';             % y_eta
ye_f_p = Yef_p';             % yf     
xe_aug_p = Xe_aug_p';
ye_aug_p = Ye_aug_p';
xe_p = Xe_p';
ye_p = Ye_eta_p';                %y_eta_e

figure
subplot(211),plot(Sim_t,per1(1,:))
% legend('\fontsize{10} perturbation_{1}','Location','SouthWest','Orientation','Horizontal')
ylabel('\fontsize{10} n_{1} ')
% axis([0 T_stop -25 30])
% plotboundary(min(ref_c(1,:)), 'y', '-.', 'k')
% plotboundary(max(ref_c(1,:)), 'y', '-.', 'k')
subplot(212),plot(Sim_t,per1(2,:))
% legend('\fontsize{10} perturbation_{2}','Location','SouthWest','Orientation','Horizontal')
ylabel('\fontsize{10} n_{2} ')
% plotboundary(min(ref_c(2,:)), 'y', '-.', 'k')
% plotboundary(max(ref_c(2,:)), 'y', '-.', 'k')

figure
subplot(211),plot(Sim_t,per(1,:))
% legend('\fontsize{10} perturbation_{1}','Location','SouthWest','Orientation','Horizontal')
ylabel('\fontsize{10} n_{1} ')
% axis([0 T_stop -25 30])
plotboundary(min(ref_eta(1,:)), 'y', '-.', 'k')
plotboundary(max(ref_eta(1,:)), 'y', '-.', 'k')
subplot(212),plot(Sim_t,per(2,:))
% legend('\fontsize{10} perturbation_{2}','Location','SouthWest','Orientation','Horizontal')
ylabel('\fontsize{10} n_{2} ')
plotboundary(min(ref_eta(2,:)), 'y', '-.', 'k')
plotboundary(max(ref_eta(2,:)), 'y', '-.', 'k')

figure
subplot(211), plot(ts_c, ref_eta(1,:), 'b-', ts_c, ye_f_p(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} r_{f,1} vs. y_{f,1}')
legend('\fontsize{10} r_{f,1}','\fontsize{10} y_{f,1}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -2.2 1.7])
subplot(212), plot(ts_c, ref_eta(2,:), 'b-', ts_c, ye_f_p(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} r_{f,2} vs. y_{f,2}')
legend('\fontsize{10} r_{f,2}','\fontsize{10} y_{f,2}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -1.2 1.2])

figure
subplot(211), plot(ts_c, ref_eta(1,:), 'b-', ts_c, ye_eta_p(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} r_{\xi,1} vs. y_{\xi,1}')
legend('\fontsize{10} r_{\xi,1}','\fontsize{10} y_{\xi,1}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -2.2 1.7])
subplot(212), plot(ts_c, ref_eta(2,:), 'b-', ts_c, ye_eta_p(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} r_{\xi,2} vs. y_{\xi,2}')
legend('\fontsize{10} r_{\xi,2}','\fontsize{10} y_{\xi,2}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -1.2 1.2])

figure
subplot(211), plot(ts_c, ref_c(1,:), 'b-', ts_c, ye_c_p(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} r_{c,1} vs. y_{c,1}')
legend('\fontsize{10} r_{c,1}','\fontsize{10} y_{c,1}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -25 17])
subplot(212), plot(ts_c, ref_c(2,:), 'b-', ts_c, ye_c_p(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} r_{c,2} vs. y_{c,2}')
legend('\fontsize{10} r_{c,2}','\fontsize{10} y_{c,2}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -20 12])

figure
subplot(211), plot(ts_c, ye_eta_p(1,:), 'b-', ts_c, ye_p(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} y_{\xi,1} vs. y_{\xi,e1}')
legend('\fontsize{10} y_{\xi,1}','\fontsize{10} y_{\xi,e1}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop ])
subplot(212), plot(ts_c, ye_eta_p(2,:), 'b-', ts_c, ye_p(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} y_{\xi,2} vs. y_{\xi,e2}')
legend('\fontsize{10} y_{\xi,2}','\fontsize{10} y_{\xi,e2}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop ])

figure
subplot 311
 plot(ts_c, xe_c_p(1,:), 'b-', ts_c, xe_p(1,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} x_{c,1} vs. x_{e,1}')
legend('\fontsize{10} x_{c,1}','\fontsize{10} x_{e,1}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -22 20])
subplot 312
plot(ts_c, xe_c_p(3,:), 'b-', ts_c, xe_p(3,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} x_{c,3} vs. x_{e,3}')
legend('\fontsize{10} x_{c,3}','\fontsize{10} x_{e,3}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -dt dt])
subplot 313
 plot(ts_c, xe_c_p(5,:), 'b-', ts_c, xe_p(5,:), 'r:', 'LineWidth', 1.5)
ylabel('\fontsize{10} x_{c,5} vs. x_{e,5}')
legend('\fontsize{10} x_{c,5}','\fontsize{10} x_{e,5}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -20 20])

figure, 
subplot 311
plot(ts_c, xe_c_p(2,:), 'b-', ts_c, xe_p(2,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} x_{c,2} vs. x_{e,2}')
legend('\fontsize{10} x_{c,2}','\fontsize{10} x_{e,2}','Location','Southeast','Orientation','Horizontal')
subplot 312
plot(ts_c, xe_c_p(4,:), 'b-', ts_c, xe_p(4,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} x_{c,4} vs. x_{e,4}')
legend('\fontsize{10} x_{c,4}','\fontsize{10} x_{e,4}','Location','SouthWest','Orientation','Horizontal')
% axis([0 T_stop -35 40])
subplot 313
 plot(ts_c, xe_c_p(6,:), 'b-', ts_c, xe_p(6,:), 'r:', 'LineWidth', 1.5)
xlabel('\fontsize{10} Time ( sec. )')
ylabel('\fontsize{10} x_{c,6} vs. x_{e,6}')
legend('\fontsize{10} x_{c,6}','\fontsize{10} x_{e,6}','Location','SouthWest','Orientation','Horizontal')

figure
subplot(311), plot(ts_c, ue_c_p(1,:), 'b-', 'LineWidth', 1.5)
ylabel('\fontsize{10} u_{c,1}')
% legend('\fontsize{10} u_{c,1}','Location','SouthWest')
% axis([0 T_stop -2e5 2e5])
subplot(312), plot(ts_c, ue_c_p(2,:), 'b-', 'LineWidth', 1.5)
ylabel('\fontsize{10} u_{c,2}')
% legend('\fontsize{10} u_{c,2}','Location','NorthWest')
% axis([0 T_stop -2e5 2e5])
subplot(313), plot(ts_c, ue_c_p(3,:), 'b-', 'LineWidth', 1.5)
ylabel('\fontsize{10} u_{c,3}')
% legend('\fontsize{10} u_{c,3}','Location','SouthWest')
% axis([0 T_stop -2e5 2e5])
xlabel('\fontsize{10} Time ( sec. )')


%% ******************************************************************** %%
%                        Control and Reference Loops                           
% *********************************************************************** %


LYe = ss(A, De, Ke1, zeros(p,p)) + ss(A, De, Ke2*C_tilde, zeros(p,p)) / s;
LXe = ss(A - De*Ke1, De*Ke2, C_tilde, zeros(p,p)) / s;


figure
sigma(LXe,'b-')
xlabel('\fontsize{10} Freq.')
ylabel('\fontsize{10} Singular Values')
legend('\fontsize{10} Uncertainty loop')
axis([ 1e-6 1e6 -400 100])
title('')

figure
sigma(LYe,'b-')
xlabel('\fontsize{10} Freq.')
ylabel('\fontsize{10} Singular Values')
legend('\fontsize{10} Sensor estimator loop')
axis([ 1e-6 1e6 -50 160])
title('')



PID_roots_e=[roots([Ked_1 Kep_2 Kei_2]);roots([Ked_1 Kep_2 Kei_2])];
TZOCS(A, De, C_tilde, zeros(p))
Syse_aug  = ss(Ae_aug, Be_aug, Ce_aug, De_aug);
Poles_eaug = eig(Ae_aug);
tzero(Ae_aug, Be_aug, Ce_aug, De_aug)
Zeros_eaug = TZOCS(Ae_aug, Be_aug, Ce_aug, De_aug);
eig(A-De*Ke1)

