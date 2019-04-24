%========================================================================
%  Function: Dispatch of Zero-Carbon-Emisssion Micro Energy Internet with
%            seperate operation mode
%  Author: Li Rui
%  Version: 1.5
%  Data: 2016/04/19 PDN
%        2016/05/26 PDN + OLTC
%        2016/06/07 heat network + power network
%        2016/06/12 add CAES with power
%        2016/06/15 add CAES with heat
%        2016/06/15 add Power from Power grid
%        2016/06/16 add peak off price
%        2016/06/18 change the capacity of NSF-CAES system
%        2016/06/29 fix the active power equation and objective function
%========================================================================
clc
close all
clear all

format short

NT = 24;                % total dispatch period
Nc = 2;                 % total stages of compressor
Ne = 2;                 % total stages of turbine

Sb = 10;                % Base power MW
Vb= 12.66;              % Base Voltage kV
Zb = Vb^2/Sb;           % Base impedence
Ib = Sb/(sqrt(3)*Vb);   % kA

%% Power Bus
    %  No.(1)|Type(2)|Pd(3)|Qd(4)  
powerbus = [
    1	3	0	    0;	   % MW  MVar 
	2	1	0.100   0.060;	
	3	1	0.090	0.040;	
	4	1	0.120	0.080;
	5	1	0.060	0.030;
	6	1	0.060	0.020;
	7	1	0.200   0.100;
	8	1   0.200   0.100;
	9	1	0.060	0.020;
	10	1	0.060	0.020;
	11	1	0.045	0.030;
	12	1	0.060	0.035;
	13	1	0.060	0.035;
	14	1	0.120	0.080;
	15	1	0.060	0.010;
	16	1	0.060	0.020;
	17	1	0.060	0.020;
	18	1	0.090	0.040;
	19	1	0.090	0.040;
	20	1	0.090	0.040;
	21	1	0.090	0.040;
	22	1	0.090	0.040;
	23	1	0.090	0.050;
	24	1	0.420	0.200;
	25	1	0.420	0.200;
	26	1	0.060	0.025;
	27	1	0.060	0.025;
	28	1	0.060	0.020;
	29	1	0.120	0.070;
	30	1	0.200	0.100;
    31	1	0.150	0.070;
    32	1	0.210	0.100;
    33	1	0.060	0.040;
];
N_bus1 = size(powerbus,1);
Pd_ratio = powerbus(:,3)/sum(powerbus(:,3));    % Active load ratio
Qd_ratio = powerbus(:,4)/sum(powerbus(:,4));    % Reactive load ratio
Pd0 = [63 62 60 58 59 65 72 85 95 99 100 99 93 92 90 88 90 92 96 98 96 90 80 70]/15-1.5;% MW
Qd0 = [18 16 15 14 15.5 15 16 17 18 19 20 20.5 21 20.5 21 19.5 20 20 19.5 19.5 18.5 18.5 18 18]/10;  % system load MVar
Pd = Pd_ratio * Pd0; 
Qd = Qd_ratio * Qd0;
Pd = Pd/Sb; 
Qd = Qd/Sb; % p.u 

U2_min = 0.95^2; 
U2_max = 1.05^2;

%% Compensator
% Location(1)|Max(2)|Min(3)|Step(4)
ComCap = [                     
    5   0.2  0  0.05;
    10  0.2  0  0.05;
    13  0.2  0  0.05;
    17  0.2  0  0.05;
    20  0.2  0  0.05;
    23  0.2  0  0.05;
    30  0.2  0  0.05;
];
v = 2;                          % Step number for linearization
N_ComCap = size(ComCap,1);      % Number of compensator
Ind_ComCap  = ComCap(:,1);
S = ComCap(:,4); 

%% SVG 
SVG = [ % Mvar   % Location(1)|Max(2)|Min(3)
    4  0.1 0;
    9  0.1 0;
    14 0.1 0;
    ];
Ind_SVG = SVG(:,1);
SVG(:,2:3) = SVG(:,2:3)/Sb; % SVG p.u

%% Heat Node
%  No(1)|Hd(2)|Pr_SR_min(3)|tao_S_max(4)|tao_S_min(5)|tao_R_max(6)|tao_R_min(7)|mass flow(8)
heatnode = [ %kW      %par   %℃
    1	0      50000   120   90   80  60  0;
	2	0      50000   120   90   80  60  0;
	3	0      50000   120   90   80  60  0;
	4	0      50000   120   90   80  60  0;
	5	250    50000   120   90   80  60  2;
	6	250    50000   120   90   80  60  2;
	7	250    50000   120   90   80  60  2;
	8	500    50000   120   90   80  60  4;
];
N_bus2 = size(heatnode,1);
H_ratio = heatnode(:,2)/sum(heatnode(:,2)); 

H_hd0 = [1250*ones(1,4), 1150*ones(1,4), 1000*ones(1,4), 800*ones(1,4), 1150*ones(1,4), 1250*ones(1,4)]; % kW
H_Hd = H_ratio * H_hd0;

Nd_Hd = find(heatnode(:,2)>0); 
m_Hd = heatnode(Nd_Hd,8);   % Mass flow ratio of heat load

tao_NS_max = heatnode(:,4);
tao_NS_min = heatnode(:,5);
tao_NR_max = heatnode(:,6);
tao_NR_min = heatnode(:,7);

%% Power Gen
Ind_gen = [2 7 19 26];
% Ind_gen = [2];
% Wg1 = [11.7  11.3  11.3  12.3 13.5 14.9 16.4 17.2 17.7 18 17.9 17.4 ...   
%      16.3 16.1 16.2 16.6 16.8 16.9 16.8 16.6 16.4 16.5 16.6 16.8]/10;% Wind Gen #2(MW)  1.8 MW
% Wg1 = [58.27 82.12 89.22 84.73 77.25 65.13 75.91 71.55 73.40 49.11 30.71 18.09 ...
%      10.80 12.50 15.00 21.62 15.00 10.88 14.50 12.54 16.00 28.41 30.34 37.10]/50; % Wind Gen #1(MW)
% Wg1 = [58.27 82.12 89.22 84.73 77.25 75.13 75.91 71.55 73.40 49.11 30.71 28.09 ...
%      20.80 22.50 25.00 21.62 25.00 20.88 24.50 22.54 26.00 28.41 40.34 47.10]/30; % Wind Gen #1(MW)
Wg1 = [6.88  7.08  7.20  7.16  6.96  6.52  6.44  5.98  5.72  5.54  5.36  5.12 ...
          4.64  4.56  4.60  4.64  4.52  4.52  4.92  5.40  5.96  6.56  6.68  6.72]-4.2;% Wind Gen #1(MW) MW  3MW
Wg2 = [16.6 16.4 16.5 16.6 16.8 11.7 11.3  11.3  12.3 13.5 14.9 16.4 17.2 17.7 18 17.9 17.4 ...   
     16.3 16.1 16.2 16.6 16.8 16.9 16.8]/40;% Wind Gen #2(MW)  0.6 MW
Wg3 = Wg2;
Wg4 = Wg2;
Wg1 = Wg1/Sb; % p.u
Wg2 = Wg2/Sb;
Wg3 = Wg3/Sb;
Wg4 = Wg4/Sb;
Wg = [Wg1;Wg2;Wg3;Wg4];

%% Heat Gen
% Location(1)|Hg_max(2)|Hg_min(3)|C_A(4)|C_B(5)|C_(6)|Mass flow(7) %
heatgen = [% kW                  % kg/s
    2     2500  0   0.05  20 0 10
];
N_gen = size(heatgen,1);
Nd_HS = heatgen(:,1);  % Node of heat station equipped with heat pump
m_HS = heatgen(:,7);
Hg_min = heatgen(:,3);
Hg_max = heatgen(:,2);

%% CAES Hub
yita_comp = [0.80, 0.75];
yita_turb = [0.86, 0.86];
beta_comp = [11.6,8.15];
pi_turb = [8.9,8.9];

Pcomp_min = zeros(Nc,1); 
Pcomp_max = 500*ones(Nc,1);    %kW
Pturb_min = zeros(Ne,1); 
Pturb_max = 1000*ones(Ne,1);    %kW

Vst = 2000;    % m^3

% Pcomp_min = zeros(Nc,1); 
% Pcomp_max = 1500*ones(Nc,1);    %kW
% Pturb_min = zeros(Ne,1); 
% Pturb_max = 3000*ones(Ne,1);    %kW
% 
% Vst = 10000;    % m^3

k = 1.4;        % adiabatic exponent

Rg = 0.297;     % KJ(kg.K)
cp_a = 1.007;   % KJ(kg.K) 25℃
cp_w = 4.2;     % KJ(kg.K) 25℃
cp_s = 2.5;     % KJ(kg.K) 25℃。

tao_am = 15;    % ambient temperature
tao_K = 273.15; 
tao_am = tao_am + tao_K; 
tao_str = 40;
tao_str = tao_str + tao_K;
tao_salt_min = 60 + tao_K;
tao_salt_max = 320 + tao_K;

pr_am = 0.101*1e3; % Kpa  ambient pressure
pr_st_min = 8.4*1e3;   % Kpa
pr_st_max = 9.0*1e3;   % Kpa

qm_comp_min = 0;
qm_comp_max = 2.306/3.6; % kg/s
qm_turb_min = 0;
qm_turb_max = 8.869/3.6; % kg/s
H_str_min = 0.2*1e3; %kW
H_str_max = 3.0*1e3; %kW

% H_str0 = 1200; % Initial heat storage

% qm_comp_min = 0;
% qm_comp_max = 9.3/3.6; % kg/s
% qm_turb_min = 0;
% qm_turb_max = 43.7/3.6; % kg/s
% H_str_min = 1e3; %kW.h
% H_str_max = 20*1e3;%kW.h

pr_comp_in1 = pr_am*ones(1,NT); % Fix pressure of each stage 
pr_comp_out1 = beta_comp(1)*pr_comp_in1;
pr_comp_in2 = pr_comp_out1 ;
pr_comp_out2 = beta_comp(2)*pr_comp_in2;
y_comp1 = (beta_comp(1))^((k-1)/k);
y_comp2 = (beta_comp(2))^((k-1)/k);
 
pr_turb_in1 = pr_st_min*ones(1,NT);
pr_turb_out1 = pr_turb_in1/pi_turb(1); 
pr_turb_in2 = pr_turb_out1;
pr_turb_out2 = pr_turb_in2/pi_turb(2);
y_turb1 = (pi_turb(1))^(-(k-1)/k);
y_turb2 = (pi_turb(2))^(-(k-1)/k);

tao_comp_in1 = tao_am*ones(1,NT); % Fix pressure of each stage
tao_comp_in2 = (40 + tao_K)*ones(1,NT);
tao_comp_out1 = tao_comp_in1/yita_comp(1).*(y_comp1-1+yita_comp(1));
tao_comp_out2 = tao_comp_in2/yita_comp(2).*(y_comp2-1+yita_comp(2));
tao_cold_s_in1 = tao_comp_out1;
tao_cold_s_out1 = 90 + tao_K ; % Fix salt heat exchanger output temperature 90℃
tao_cold_w_in1 = tao_cold_s_out1;
tao_cold_w_out1 = 40 + tao_K ; % water heat exchanger output temperature 40℃
tao_cold_s_in2 = tao_comp_out2;
tao_cold_s_out2 = 90 + tao_K; 
tao_cold_w_in2 = tao_cold_s_out2;
tao_cold_w_out2 = 40 + tao_K ; 

tao_turb_in1  = (280 + tao_K)*ones(1,NT);
tao_turb_in2  = (280 + tao_K)*ones(1,NT);
tao_turb_out1 = tao_turb_in1*yita_turb(1).*(y_turb1-1+1/yita_turb(1));
tao_turb_out2 = tao_turb_in2*yita_turb(2).*(y_turb2-1+1/yita_turb(2));
tao_heat_in1 = tao_str;
tao_heat_in2 = tao_turb_out1;
tao_heat_out1 = tao_turb_in1;
tao_heat_out2 = tao_turb_in2;

CAES_ind = 2;
%% Power Line
% No.(1)|From bus(2)|To bus(3)|r(4)|x(5)|P_line_max(6)|P_line_min(7) %
branch = [
    1   1   2	0.0922	0.0470	9.9	 0; 
	2   2	3	0.4930	0.2512	9.9	 0;  
	3   3	4	0.3661	0.1864	9.9	 0;
    4   4	5	0.3811	0.1941	9.9  0;  
	5   5	6	0.8190	0.7070	9.9	 0; 
	6   6	7	0.1872	0.6188	9.9	 0; 
	7   7	8	0.7115	0.2351	9.9  0; 
	8   8	9	1.0299	0.7400	9.9  0;  
	9   9	10	1.0440	0.7400	9.9	 0; 
	10  10	11	0.1967	0.0651	9.9	 0; 
	11  11	12	0.3744	0.1298	9.9  0; 
	12  12	13	1.4680	1.1549	9.9	 0; 
	13  13	14	0.5416	0.7129	9.9	 0;
	14  14	15	0.5909	0.5260	9.9	 0;
	15  15	16	0.7462	0.5449	9.9	 0;
	16  16	17	1.2889	1.7210	9.9	 0;
	17  17	18	0.7320	0.5739	9.9	 0;
	18  2	19	0.1640	0.1565	9.9	 0;
	19  19	20	1.5042	1.3555	9.9  0;
	20  20	21	0.4095	0.4784	9.9	 0;
	21  21	22	0.7089	0.9373	9.9	 0;
	22  3	23	0.4512	0.3084	9.9	 0;
	23  23	24	0.8980	0.7091	9.9	 0;
	24  24	25	0.8959	0.7071	9.9	 0;
	25  6	26	0.2031	0.1034	9.9	 0;
	26  26	27	0.2842	0.1447	9.9	 0;
	27  27	28	1.0589	0.9338	9.9	 0;
	28  28	29	0.8043	0.7006	9.9	 0;
	29  29  30	0.5074	0.2585	9.9	 0;
	30  30	31	0.9745	0.9629	9.9	 0;
	31  31	32	0.3105	0.3619	9.9	 0;
	32  32	33	0.3411	0.5302	9.9	 0;
];
N_line = size(branch,1);
line_i = branch(:,2);
line_j = branch(:,3);
r = branch(:,4)/Zb;
x = branch(:,5)/Zb;
Pmax = branch(:,6)/Sb; 
Pmin = branch(:,7)/Sb;

%% OLTC
% Line No.(1)|K_max(2)|K_min(3)|K_Step(4)|%
OLTC = [
    1   1.05  0.95  0.01;
    18  1.05  0.95  0.01;
    22  1.05  0.95  0.01;
    25  1.05  0.95  0.01;
]; 
t_OLTC = 0.95:0.01:1.05; % available tap value 
T_OLTC = repmat(t_OLTC',1,NT);
n_OLTC= length(t_OLTC); % num of total tap value
N_OLTC = size(OLTC,1);% num of OLTC 
Ind_OLTC = OLTC(:,1); 

Ind_subline = zeros(N_line,2); % Index of Children line of power netwrok 
for i = 1: N_line
    temp = find(line_i == line_j(i));
    if  ~isempty(temp)
        Ind_subline(i,1:length(temp)) = temp;
    end
end

Pg_min = zeros(N_bus1,NT); 
Pg_max = zeros(N_bus1,NT);

%% Pipe
% No.(1)|From node(2)|To node(3)|L(4)|u_p(5)|u_T(6)|ms_max(7)|ms_min(8)|mr_max(9)|mr_min(10)
pipe = [
    1   1   2	1000  0.04  0.5*1e-3 1e6 0  10 10;
	2   2	3	1000  0.04  0.5*1e-3 1e6 0  8  8;
	3   3	4	1000  0.04  0.5*1e-3 1e6 0  6  6;
    4   2	5	1000  0.04  0.5*1e-3 1e6 0  2  2;
	5   3	6	1000  0.04  0.5*1e-3 1e6 0  2  2;
	6   4	7	1000  0.04  0.5*1e-3 1e6 0  2  2;
	7   4	8	1000  0.04  0.5*1e-3 1e6 0  4  4;
];
N_pipe = size(pipe,1);
pipe_i = pipe(:,2);
pipe_j = pipe(:,3);
L_pipe = pipe(:,4);
% miu_pipe = pipe(:,5);
lamada_pipe = pipe(:,6);
m_pipe_max = pipe(:,7);
m_pipe_min = pipe(:,8);
ms_pipe = pipe(:,9);
mr_pipe = pipe(:,10);

% Index
S_pipe_F = zeros(N_bus2,2);
S_pipe_T = zeros(N_bus2,2);

%% Variables
% % PDN
P = sdpvar(N_line,NT,'full'); % active power on each line
Q = sdpvar(N_line,NT,'full'); % reactive power on each line
Psub = sdpvar(N_line,NT,'full'); %
Qsub = sdpvar(N_line,NT,'full'); % 
I2 = sdpvar(N_line,NT,'full'); %  square of line current amplititude
U2 = sdpvar(N_bus1,NT,'full'); %  square of bus voltage amplititude
Pg = sdpvar(N_bus1,NT,'full'); % Injected generator active power to each bus 
% Qg = sdpvar(N_bus1,NT,'full'); % injected generator reactice power to each bus
Qc = sdpvar(N_bus1,NT,'full'); % SVG
Pgrid = sdpvar(1,NT); % add 2016/06/15 power bought from tranmission network
Qgrid = sdpvar(1,NT);
% %  DHN
tao_PS_F = sdpvar(N_pipe,NT,'full'); %'From' side temperature of supply pipe 
tao_PS_T = sdpvar(N_pipe,NT,'full'); %'To' side temperature of supply pipe 
tao_PR_F = sdpvar(N_pipe,NT,'full'); %'From' side temperature of return pipe
tao_PR_T = sdpvar(N_pipe,NT,'full'); %'To' side temperature of return pipe
tao_NS = sdpvar(N_bus2,NT,'full');   % Node temperature of supply network； 
tao_NR = sdpvar(N_bus2,NT,'full');   % Node temperature of return network； 
Hg_HP = sdpvar(N_bus2,NT,'full');    % Heat power of heat pump equipped with CAES

% %　CAES
on_comp = binvar(1,NT);  % on/off of comp.
on_turb = binvar(1,NT);  % on/off of turb.
Pcomp1 = sdpvar(1,NT);   % power cons. of comp1.
Pcomp2 = sdpvar(1,NT); 
Pturb1= sdpvar(1,NT);
Pturb2= sdpvar(1,NT);
Pcaes_d = sdpvar(1,NT);
Pcaes_g = sdpvar(1,NT);
pr_st = sdpvar(1,NT);  % 储气罐压强
pr_st0 = sdpvar(1,1);
qm_comp = sdpvar(1,NT); % 
qm_turb = sdpvar(1,NT);

y1 = sdpvar(1,NT);% 线性化储气室压强约束
y2 = sdpvar(1,NT);
h1 = sdpvar(1,NT);% 线性化储热系统SOC
h2 = sdpvar(1,NT);

HM = 1e7; % big M
H_coll_s1 = sdpvar(1,NT); % collected heat by salt
H_coll_s2 = sdpvar(1,NT);
H_cons1 = sdpvar(1,NT);
H_cons2 = sdpvar(1,NT);
H_coll_sum = sdpvar(1,NT);
H_cons_sum = sdpvar(1,NT);
H_str = sdpvar(1,NT);   % 储热罐中存贮的热量
H_str0 = sdpvar(1,1);
Hg_CAES = sdpvar(1,NT); % 蓄热环节可供热负荷

for i = 1:N_ComCap  %% 线性化变量
    xd{i} = binvar(v+1,NT);
    delta{i} = sdpvar(v+1,NT);
end
for i = 1:N_OLTC
    rd{i} = binvar(n_OLTC,NT);
    h{i} = sdpvar(n_OLTC,NT);
end

%% Constraints
% % CAES
F_turb = [];    % 每级功率定义
F_comp = [];    % 每级功率定义
F_oper = [];    % 运行约束
F_power = [];   % 功率平衡约束
F_airstr = [];  % 储气罐压强动态约束
F_cold = [];
F_heat = [];
F_heatstr = [];

F_comp = [F_comp, Pcomp1 == 1/yita_comp(1)*k/(k-1)*Rg*qm_comp.*tao_comp_in1.*(y_comp1-1)];
F_comp = [F_comp, Pcomp2 == 1/yita_comp(2)*k/(k-1)*Rg*qm_comp.*tao_comp_in2.*(y_comp2-1)];
F_comp = [F_comp, Pcomp_min(1)*on_comp <= Pcomp1 <= Pcomp_max(1)*on_comp ];% 每级消耗的功率约束
F_comp = [F_comp, Pcomp_min(2)*on_comp <= Pcomp2 <= Pcomp_max(2)*on_comp ];
F_comp = [F_comp, Pcaes_d == Pcomp1 + Pcomp2];% 总功率定义

F_turb = [F_turb, Pturb1 == yita_turb(1)*k/(k-1)*Rg*qm_turb.*tao_turb_in1.*(1-y_turb1)];
F_turb = [F_turb, Pturb2 == yita_turb(2)*k/(k-1)*Rg*qm_turb.*tao_turb_in2.*(1-y_turb2)];
F_turb = [F_turb, Pturb_min(1)*on_turb <= Pturb1 <= Pturb_max(1)*on_turb];% 每级发出功率约束
F_turb = [F_turb, Pturb_min(2)*on_turb <= Pturb2 <= Pturb_max(2)*on_turb];
F_turb = [F_turb, Pcaes_g == Pturb1 + Pturb2];

F_oper = [F_oper, 0 <= on_comp + on_turb <= 1];%充放电不能同时进行
F_oper = [F_oper, qm_comp_min*on_comp <= qm_comp <= qm_comp_max*on_comp];%质量流量非否约束
F_oper = [F_oper, qm_turb_min*on_turb <= qm_turb <= qm_turb_max*on_turb];

F_airstr = [F_airstr, pr_st(1) == pr_st0 + 1/Vst * Rg * tao_str * 3600*(y1(1)- y2(1))];
F_airstr = [F_airstr, pr_st(1,2:NT) == pr_st(1,1:NT-1) + 1/Vst * Rg * tao_str * 3600*(y1(2:NT)- y2(2:NT))];
F_airstr = [F_airstr, pr_st(1,NT) == pr_st0]; % add 2016/06/15
F_airstr = [F_airstr, qm_comp_min*on_comp <= y1 <= qm_comp_max*on_comp];
F_airstr = [F_airstr, qm_turb_min*on_turb <= y2 <= qm_turb_max*on_turb];
F_airstr = [F_airstr, qm_comp_min*(1-on_comp) <= qm_comp - y1 <= qm_comp_max*(1-on_comp)];
F_airstr = [F_airstr, qm_turb_min*(1-on_turb) <= qm_turb - y2 <= qm_turb_max*(1-on_turb)];
F_airstr = [F_airstr, pr_st_min <= pr_st <= pr_st_max];
F_airstr = [F_airstr, pr_st_min <= pr_st0 <= pr_st_max]; % add 2016/06/15

F_cold = [F_cold, H_coll_s1 == cp_a*qm_comp.*(tao_cold_s_in1 - tao_cold_s_out1)];%% 回热系统热量约束% 冷却环节
F_cold = [F_cold, H_coll_s2 == cp_a*qm_comp.*(tao_cold_s_in2 - tao_cold_s_out2)];% 第1,2级导热油回收的热量
F_cold = [F_cold, H_coll_sum == H_coll_s1 + H_coll_s2]; % 只计及导热油收集的热量% 回收的总热量

F_heat = [F_heat, H_cons1 == cp_a*qm_turb.*(tao_heat_out1 - tao_heat_in1)];% 第1-2级消耗的热量
F_heat = [F_heat, H_cons2 == cp_a*qm_turb.*(tao_heat_out2 - tao_heat_in2)];
F_heat = [F_heat, H_cons_sum == H_cons1 + H_cons2]; % 加热环节消耗的总热量
% F_heat = [F_heat, sum(H_cons_sum) <= sum(H_coll_sum)]; % add 2016/06/15
% is needed？

F_heatstr = [F_heatstr, H_str(1,1) == H_str0 + h1(1,1) - h2(1,1) - Hg_CAES(1,1)];% 回热系统SOC% 高温储热罐中存贮的热量
F_heatstr = [F_heatstr, H_str(1,2:NT) == H_str(1,1:NT-1) + h1(1,2:NT) - h2(1,2:NT) - Hg_CAES(1,2:NT)];
% F_heatstr = [F_heatstr, H_str(1,1) == H_str0 + h1(1,1) - h2(1,1)];% 回热系统SOC% 高温储热罐中存贮的热量
% F_heatstr = [F_heatstr, H_str(1,2:NT) == H_str(1,1:NT-1) + h1(1,2:NT) - h2(1,2:NT)];
F_heatstr = [F_heatstr, -HM*on_comp <= h1 <= HM*on_comp];
F_heatstr = [F_heatstr, -HM*on_turb <= h2 <= HM*on_turb];
F_heatstr = [F_heatstr, -HM*(1-on_comp) <= H_coll_sum - h1 <= HM*(1-on_comp)];
F_heatstr = [F_heatstr, -HM*(1-on_turb) <= H_cons_sum - h2 <= HM*(1-on_turb)];
F_heatstr = [F_heatstr, H_str_min <= H_str <= H_str_max];% % Box 约束
F_heatstr = [F_heatstr, H_str_min <= H_str0 <= H_str_max];
F_heatstr = [F_heatstr, H_str(1,NT) == H_str0];  % add 2016/06/15 

%% PDN 约束
F_P = [];
F_Q = [];
F_U = [];
for t = 1:NT
    for i = 1:N_bus1 % 设置未放置发电机的母线节点的注入(P，Q)功率及a,b,c为 0
        if isempty(find(Ind_gen == i))
            F_P = [F_P, Pg(i,t) == 0];
%             F_Q = [F_Q, Qg(i,t) == 0];
        else
            Pg_max(i,:) = Wg(find(Ind_gen == i),:);
        end
        % 判断是否装有SVG
        if isempty(find(Ind_SVG ==i)) % 未安装SVG
            F_Q = [F_Q, Qc(i,t) == 0];
            else  % 安装SVG
            temp = find(Ind_SVG == i); % 索引
            F_Q = [F_Q, SVG(temp,3) <= Qc(i,t) <= SVG(temp,2)];
        end
    end
end

% 线路子线路潮流
for t = 1:NT
    for i = 1:N_line
    num_temp = size(find(Ind_subline(i,:) == 0),2);
    if num_temp == 1               
        F_P = [F_P, Psub(i,t) == P(Ind_subline(i,1),t)];
        F_Q = [F_Q, Qsub(i,t) == Q(Ind_subline(i,1),t)];
    elseif num_temp == 2
        F_P = [F_P, Psub(i,t) == 0];
        F_Q = [F_Q, Qsub(i,t) == 0];
    else
        F_P = [F_P, Psub(i,t) == P(Ind_subline(i,1),t) + P(Ind_subline(i,2),t)];
        F_Q = [F_Q, Qsub(i,t) == Q(Ind_subline(i,1),t) + Q(Ind_subline(i,2),t)];
    end
end
end

OTLC_count = zeros(1,NT);
ComCap_count = zeros(1,NT);
M = 1000;

F_P = [F_P, P(1,:) == Pgrid(1,:)]; %%
F_Q = [F_Q, Q(1,:) == 0]; %%
F_U = [F_U, U2(1,:) == 1.05^2];
Vs1 = 1.05^2 * ones(1,NT); 

for t = 1:NT  
    t
    for i = 1:N_line
        if line_j(i) == Ind_gen
%             F_P = [F_P, P(i,t) + Pg(line_j(i),t) + Pcaes_g(t)/1e3/Sb  == Psub(i,t) + Pd(line_j(i),t) + Pcaes_d(t)/1e3/Sb];
            F_P = [F_P, P(i,t) + Pg(line_j(i),t)  == Psub(i,t) + Pd(line_j(i),t)];
        else
            F_P = [F_P, P(i,t) + Pg(line_j(i),t) == Psub(i,t) + Pd(line_j(i),t)];
        end


    % 判断是否装有并联补偿电容
    if ~isempty(find(Ind_ComCap == line_j(i))) % 安装补偿电容
        ComCap_count(t) = ComCap_count(t) + 1;
%         F_Q = [F_Q, Q(i,t) + Qg(line_j(i),t) + 0.5*(U2(line_j(i),t)*ComCap(ComCap_count(t),3) + S(ComCap_count(t))*(2^0*delta{ComCap_count(t)}(1,t) + ...
%             2^1*delta{ComCap_count(t)}(2,t)+ 2^2*delta{ComCap_count(t)}(2,t)))+ Qc(line_j(i),t) - x(i)*I2(i,t) == Qsub(i,t) + Qd(line_j(i),t)];
%         F_Q = [F_Q, Q(i,t) + Qg(line_j(i),t) + 0.5*(U2(line_j(i),t)*ComCap(ComCap_count(t),3) + S(ComCap_count(t))*(2^0*delta{ComCap_count(t)}(1,t) + ...
%             2^1*delta{ComCap_count(t)}(2,t) + 2^2*delta{ComCap_count(t)}(2,t)))+ Qc(line_j(i),t)  == Qsub(i,t) + Qd(line_j(i),t)];
        F_Q = [F_Q, Q(i,t) + 0.5*(U2(line_j(i),t)*ComCap(ComCap_count(t),3) + S(ComCap_count(t))*(2^0*delta{ComCap_count(t)}(1,t) + ...
            2^1*delta{ComCap_count(t)}(2,t) + 2^2*delta{ComCap_count(t)}(2,t)))+ Qc(line_j(i),t)  == Qsub(i,t) + Qd(line_j(i),t)]; % 2016/06/15 去掉无功Qg
        for  m = 1:v+1  
             F_U = [F_U, U2(line_j(i),t) - M*(1-xd{ComCap_count(t)}(m,t)) <= delta{ComCap_count(t)}(m,t) <= U2(line_j(i),t) + M*(1-xd{ComCap_count(t)}(m,t))];
             F_U = [F_U, -M*xd{ComCap_count(t)}(m,t) <= delta{ComCap_count(t)}(m,t) <= M*xd{ComCap_count(t)}(m,t)];
        end
        F_U = [F_U, 0 <= 2^0*xd{ComCap_count(t)}(1,t)+ 2^1*xd{ComCap_count(t)}(2,t) + 2^2*xd{ComCap_count(t)}(3,t) <= (ComCap(ComCap_count(t),2) - ...
            ComCap(ComCap_count(t),3))/ComCap(ComCap_count(t),4)];
    else % 未安装补偿电容)
%         F_Q = [F_Q, Q(i,t) + Qg(line_j(i),t) + Qc(line_j(i),t) - x(i)*I2(i,t) == Qsub(i,t) + Qd(line_j(i),t)];
%         F_Q = [F_Q, Q(i,t) + Qg(line_j(i),t) + Qc(line_j(i),t) == Qsub(i,t) + Qd(line_j(i),t)];
        F_Q = [F_Q, Q(i,t) + Qc(line_j(i),t) == Qsub(i,t) + Qd(line_j(i),t)];% 2016/06/15 去掉无功Qg
    end

    if ~isempty(find(Ind_OLTC == i)) % 含OTLC的支路
        OTLC_count(t) = OTLC_count(t)+1;
        F_U = [F_U, sum(h{OTLC_count(t)}(:,t)./T_OLTC(:,t).^2,1) == U2(line_i(i),t)-(r(i)*P(i,t)+x(i)*Q(i,t))/Vs1(1,t)];
        for k = 1:n_OLTC
            F_U = [F_U, -M*(1-rd{OTLC_count(t)}(k,t)) + U2(line_j(i),t) <= h{OTLC_count(t)}(k,t) <= U2(line_j(i),t) + M*(1-rd{OTLC_count(t)}(k,t))];
            F_U = [F_U, -M*rd{OTLC_count(t)}(k,t) <= h{OTLC_count(t)}(k,t)<= M*rd{OTLC_count(t)}(k,t)];
        end
        F_U = [F_U,sum(rd{OTLC_count(t)}(:,t),1) == 1];
    else   % 不含OTLC的支路
        F_U = [F_U, U2(line_j(i),t)== U2(line_i(i),t)-(r(i)*P(i,t)+x(i)*Q(i,t))/Vs1(1,t)];
    end
    % 线路潮流约束
    F_P = [F_P, Pmin(i) <= P(i,t) <= Pmax(i)];
    end
  
    for i = 1:N_bus1
        F_P = [F_P, Pg_min(i,t) <= Pg(i,t) <= Pg_max(i,t)];
        F_U = [F_U, U2_min <= U2(i,t) <= U2_max];
    end
   
end
F_P = [F_P, Pgrid >= 0];
% F_Q = [F_Q, Qgrid >= 0];
% 热力站
F_H = [];
Count_HS = 0;
Count_Hd = 0;
for j = 1:N_bus2
    if ~isempty(find(Nd_HS == j)) % 若该节点设置供热机组
         Count_HS = Count_HS + 1; 
         F_H = [F_H, Hg_HP(j,:) == cp_w*m_HS(Count_HS)*(tao_NS(j,:) - tao_NR(j,:))];% add heat power generated by CAES
         F_H = [F_H, Hg_min <= Hg_HP(j,:) <= Hg_max];
    elseif ~isempty(find(Nd_Hd == j)) % 若该节点设置供热负荷
         Count_Hd = Count_Hd + 1; 
         F_H = [F_H, Hg_HP(j,:) == zeros(1,NT)];
         F_H = [F_H, cp_w*m_Hd(Count_Hd,:)*(tao_NS(j,:) - tao_NR(j,:)) == H_Hd(j,:)];
    else % 其他节点
         F_H = [F_H, Hg_HP(j,:) == zeros(1,NT)];
    end
         F_H = [F_H, tao_NS_min(j) <= tao_NS(j,:) <= tao_NS_max(j)];
         F_H = [F_H, tao_NR_min(j) <= tao_NR(j,:)<= tao_NR_max(j)];
end
F_H = [F_H, 0 <= Hg_CAES <= 0.2*H_str]; % 
% 供热网络
F_PH = [];
for i = 1:N_bus2
    temp1 = find(pipe_i == i); % 管道首节点为i
    temp2 = find(pipe_j == i); % 管道末节点为i
    S_pipe_F(i,1:length(temp1)) = temp1; % 首节点为i的管道集合
    S_pipe_T(i,1:length(temp2)) = temp2; % 末节点为i的管道集合
end

for i = 1:N_bus2
    % 供水管道温度节点混合
    num_temp1 = size(find(S_pipe_T(i,:) == 0),2);
    if num_temp1 == 1  % 以节点i末端的管道数为1
        b = S_pipe_T(i,1); % 管道编号
        F_PH = [F_PH, ms_pipe(b)* tao_PS_T(b,:) == ms_pipe(b) * tao_NS(i,:)];
        F_PH = [F_PH, tao_PR_F(b,:) == tao_NR(i,:)];
    elseif num_temp1 == 0  % 以节点i末端的管道数为2
        b = S_pipe_T(i,:); % 管道编号
        F_PH = [F_PH, ms_pipe(b(1))* tao_PS_T(b(1),:) + ms_pipe(b(2))* tao_PS_T(b(2),:) == (ms_pipe(b(1))+ms_pipe(b(2))) * tao_NS(i,:)];
 
        F_PH = [F_PH, tao_PR_F(b(1),:) == tao_NR(i,:)];
        F_PH = [F_PH, tao_PR_F(b(2),:) == tao_NR(i,:)];
    else  
        % % 以节点i末端的管道数为0  首节点
    end
    
    % 回水管道温度节点混合
    num_temp2 = size(find(S_pipe_F(i,:) == 0),2);
    if num_temp2 == 1  % 以节点i首端的管道数为1
        b = S_pipe_F(i,1); % 管道编号
        F_PH = [F_PH, mr_pipe(b)* tao_PR_T(b,:) == mr_pipe(b) * tao_NR(i,:)];
        
        F_PH = [F_PH, tao_PS_F(b,:) == tao_NS(i,:)];
    elseif num_temp2 == 0  % 以节点i首端的管道数为2
        b = S_pipe_F(i,:); % 管道编号
        F_PH = [F_PH, mr_pipe(b(1))* tao_PR_T(b(1),:) + mr_pipe(b(2))* tao_PR_T(b(2),:) == (mr_pipe(b(1))+mr_pipe(b(2))) * tao_NR(i,:)];
        
        F_PH = [F_PH, tao_PS_F(b(1),:) == tao_NS(i,:)];
        F_PH = [F_PH, tao_PS_F(b(2),:) == tao_NS(i,:)];
    else  
        % 以节点i末端的管道数为0  末节点
    end
end
for i = 1:N_pipe
    % 温度变化方程
    F_PH = [F_PH, tao_PS_T(i,:) == (tao_PS_F(i,:) - (tao_am-tao_K)).*exp(-lamada_pipe(i)*L_pipe(i)/(cp_w*ms_pipe(i))) + (tao_am-tao_K)];
    F_PH = [F_PH, tao_PR_T(i,:) == (tao_PR_F(i,:) - (tao_am-tao_K)).*exp(-lamada_pipe(i)*L_pipe(i)/(cp_w*mr_pipe(i))) + (tao_am-tao_K)];
end

%% 目标函数
Obj = 0;% 发电成本最小
% C_grid = 1500; % 0.15$/(kW.h) -> 1500$/(p.u)
C_grid = [0.05*ones(1,8) 0.10*ones(1,6) 0.08*ones(1,8) 0.05*ones(1,2)];
C_grid = C_grid * 1e4;
% C_wind = 2500; % 0.25$/(kW.h) -> 2500$/(p.u)
yite = 0.8;
Obj1 = C_grid*Pgrid'; % PDN 购电成本
Obj2 = C_grid*[(Hg_HP(CAES_ind,:))/1e3/Sb/yite]'; % DHN 购电成本
Obj = Obj1 + Obj2;
FP = [F_P,F_Q,F_U];
FH = [F_H,F_PH];
FHub = [F_comp,F_turb,F_power,F_oper,F_airstr,F_heat, F_cold ,F_heatstr];
F = [FP,FH,FHub];
opt= sdpsettings('solver','cplex');
sol = optimize(F,Obj,opt)
% checkset(F)
if sol.problem == 0
    disp('successed!')
    C_PDN = value(Obj1)
    C_DHN = value(Obj2)
    C_MEI = value(Obj)
    tte = 11;     % 绘图时刻
    tth = 11;
    nne = 2;    % 绘图电网节点编号
    lineno = 7; % 绘图电网线路编号
    pipeno = 8;
    
    figure(1)
    plot(1:NT,Pd0,'-b')
    hold on
    plot(1:NT,Qd0,'-.r')
    plot(1:NT,Sb*sum(Wg,1),'-g')
    plot(1:NT,H_hd0/1e3,'-.k')
    legend('P_d','Q_d','W_g')
    xlabel('Time (h)')
    ylabel('Power (MW/Mvar)')
    
%     figure(2)
%     plot(1:NT,Pd0,'-r')
%     hold on
%     plot(1:NT,Qd0,'-b')
%     xlabel('Time (h)')
%     ylabel('Power load (MW/Mvar)')
%     legend('P_d','Q_d')
    
%     figure(3)
%     plot(1:N_bus1,Pd_ratio,'-r')
%     hold on
%     plot(1:N_bus1,Qd_ratio,'-b')
%     xlabel('Bus No.')
%     ylabel('ratio')
%     legend('P_d','Q_d')
%     
%     figure(4)
%     plot(1:NT,Sb*Wg1,'-r')
%     hold on
%     plot(1:NT,Sb*Wg2,'-b')
%     xlabel('Time (h)')
%     ylabel('W_g(MW)')
%     legend('Wind Gen #1', 'Wind Gen #2-4')
    
    figure(5)
    plot(1:NT,Pd0,'-ko')
    hold on
    plot(1:NT,value(Pcaes_d)/1e3,'-.ks')
    plot(1:NT,Sb*value(Pg(2,:)),'-bo')
    plot(1:NT,Sb*value(Pg(7,:)),'-.bs')
    plot(1:NT,Sb*value(Pg(19,:)),'-ro')
    plot(1:NT,Sb*value(Pg(26,:)),'-.gs')
    plot(1:NT,Sb*value(Pgrid),'-go')
    plot(1:NT,value(Pcaes_g)/1e3,'-.rp')
    legend('P_d','P^{CAES}_d','W^g_1','W^g_2','W^g_3','W^g_4','\theta','P^{CAES}_g')
    xlabel('Time (h)')
    ylabel('Power (MW)')
    
%     figure(6)
%     plot(1:NT,-Sb*(value(Pg(2,:))-Wg(1,:)),'-.bo')
%     hold on
%     plot(1:NT,-Sb*(value(Pg(7,:))-Wg(2,:)),'-.rs')
%     plot(1:NT,-Sb*(value(Pg(19,:))-Wg(3,:)),'-bo')
%     plot(1:NT,-Sb*(value(Pg(26,:))-Wg(4,:)),'-rs')
%     legend('W^c_1','W^c_2','W^c_3','W^c_4')
%     xlabel('Time (h)')
%     ylabel('Wind  (MW)')
    
%     figure(7)
%     plot(1:N_line,Sb*value(P(:,tte)),'-ro')
%     hold on 
%     plot(1:N_line,Sb*value(Q(:,tte)),'-bs')
%     xlabel('Line No')
%     ylabel('Power flow (MW/Mvar)')
%     legend('P','Q')  
% 
%     figure(8)
%     plot(1:NT,Sb*value(P(lineno,:)),'-ro')
%     hold on 
%     plot(1:NT,Sb*value(Q(lineno,:)),'-bs')
%     xlabel('Time(h)')
% %     ylabel(['Power flow on line ', num2str(lineno),' (MW/Mvar) '])
%     ylabel('Power flow (MW/Mvar)')
%     legend('P','Q') 
%     
%     figure(9)
%     plot(1:N_bus1,sqrt(value(U2(:,tte))),'-ro')
%     xlabel('Bus No.')
%     ylabel('Voltage amplititude (p.u.)')
% %     ylabel(['Voltage amplititude at time period ',num2str(tte),' (p.u.)'])
% %    
%     figure(10)
%     plot(1:NT,sqrt(value(U2(nne,:))),'-ro')
%     xlabel('Time (h)')
% %     ylabel(['Voltage amplititude  at node ', num2str(nne), ' (p.u.)'])
%     ylabel('Voltage amplititude (p.u.)')
% 
% %     figure(11)
% %     plot(1:N_bus1,Sb*value(Qc(:,tte)),'-ro');
% %     xlabel('Bus No.')
% %     ylabel(['VSG at time period ', num2str(tte), ' (MVar)'])
% 
%     figure(11)
%     plot(1:NT,Sb*value(Qc(4,:)),'-ro');
%     xlabel('Time (h)')
% %     ylabel(['VSG output at Node 4 (MVar)'])
%     ylabel('VSG output (MVar)')
    
    
%     tij = zeros(N_OLTC,1);
%     for i = 1:N_OLTC
%         tij(i,1) = t_OLTC(find(value(rd{i}(:,tte)) == 1));
%     end
%     
% 
%     figure(12)
%     plot(1:N_OLTC,tij,'-ro')
%     xlabel('OLTC (No.)')
% %     ylabel('K at time period 1')    
%     ylabel('Value of OLTC Tap Changer K')  

    
%     CC = zeros(N_ComCap,1);
%     for i = 1:N_ComCap
%         CC(i,1) = (2^0*value(xd{i}(1,tte))+ 2^1*value(xd{i}(2,tte)) + 2^2*value(xd{i}(3,tte)))*S(i) + ComCap(i,3);
%     end
%     figure(13)
%     plot(1:N_ComCap,CC,'-ro')
%     xlabel('Compen (No.)')
%     ylabel('C ')   
%   
%     figure(14)
%     plot(1:N_pipe,value(tao_PS_F(:,tth)),'-ro')
%     hold on
%     plot(1:N_pipe,value(tao_PS_T(:,tth)),'-bs')
%     plot(1:N_pipe,value(tao_PR_F(:,tth)),'-.ro')
%     plot(1:N_pipe,value(tao_PR_T(:,tth)),'-.bs')
%     xlabel('Pipe No.')
% %     ylabel(['Temperature at time period ', num2str(tth), '(^{\circ}C)'])
%     ylabel('Pipe temperature (^{\circ}C)')
%     legend('PS_F','PS_T','PR_F','PR_T')
    
%     figure(15)
%     plot(1:N_bus2,value(tao_NS(:,tth)),'-ro')
%     hold on
%     plot(1:N_bus2,value(tao_NR(:,tth)),'-bs')
%     xlabel('Node No.')
% %     ylabel(['Temperature at time period ', num2str(tth), '(^{\circ}C)'])
%     ylabel(['Node temperature (^{\circ}C)'])
%     legend('\tau_{NS}','\tau_{NR}')
%     
%     figure(16)
%     plot(1:NT,value(Hg_HP(CAES_ind,:))/1e3,'-ro')
%     hold on 
% %     plot(0:NT,value(heatnode(:,2)),'-bs')
%     plot(1:NT,H_hd0/1e3)
%     plot(1:NT,value(Hg_CAES(1,:))/1e3,'-.bp')
%     xlabel('Time (h)')
%     ylabel('Heat(MW)')
%     legend('H^{HP}_g','H_d','H^{CAES}_g')
%     
%     figure(17)
%     plot(1:NT,value(tao_NS(4,:)),'-ro')
%     hold on
%     plot(1:NT,value(tao_NR(4,:)),'-bs')
%     xlabel('Time (h)')
%     ylabel(['Temperature of node 4 (^{\circ}C) '])
%     legend('\tau_{NS}','\tau_{NR}')
    
%     figure(18)
%     plot(0:NT,[0 value(on_comp)],'-bo')
%     hold on 
%     plot(0:NT,[0 value(on_turb)],'-.rs')
%     xlabel('Time (h)')
%     legend('compresser','turbine')
%     title('on/off of CAES')
%     xlim([0-0.1,NT+0.1])
% %     ylim([-0.1,1.1])

%     figure(19)
%     plot(1:NT,[value(Pcaes_d)]/1e3,'-bo')
%     hold on
%     plot(1:NT,[value(Pcaes_g)]/1e3,'-rs')
%     legend('P^{CAES}_d','P^{CAES}_g')
%     xlabel('Time (h)')
%     ylabel('P (WW)')
%     xlim([0-0.1,NT+0.1])

%     figure(20)
%     plot(1:NT,[value(pr_st)]/1e3,'-bo')
%     xlabel('Time (h)')
%     ylabel('Pr (Mpa)')
%     xlim([0-0.1,NT+0.1])
    figure(20)
    plot(1:NT, Sb*value(Pgrid),'-b')
    hold on
    plot(1:NT, Sb*Wg1,'-c')
    plot(1:NT, Sb*value(Pg(2,:)),'-g')
%     plot(1:NT, Sb*(Wg1-value(Pg(2,:))),'-r')
%     fill(zeros(1,NT), Sb*(Wg1-value(Pg(2,:))),'r')
    shadedplot(1:NT,Sb*(Wg1-value(Pg(2,:))),zeros(1,NT),'r');
    xlabel('Time (h)')
    ylabel('Power (MW)')
    legend('\theta','W^{g,u}_2','W^g_2','W^{cur}')

%     figure(21)
%     plot(1:NT, value(qm_comp),'-b')
%     hold on 
%     plot(1:NT, value(qm_turb),'-r')
%     legend('qm_{comp}','qm_{turb}')
%     xlim([0-0.1,NT+0.1])
% %     ylim([0-10,50])
%     xlabel('Time (h)')
%     ylabel('qm (kg/s)')
%     
%     figure(22)
%     x = 0:NT;
%     y1 = [0 value(Pcaes_d) - value(Pcaes_g)]/1e3;
%     y2 = [value(pr_st0),value(pr_st)]/1e3;
%     [AX,H1,H2] = plotyy(x,y1,x,y2,'bar','plot');
%     set(AX(1),'XColor','k','YColor','b');
%     set(AX(2),'XColor','k','YColor','r');
%     set(AX(1),'Xlim',[0,25]);
%     set(AX(2),'Xlim',[0,25]);
%     HH1=get(AX(1),'Ylabel');
%     set(HH1,'String','P^{CAES}(MW)');
%     set(HH1,'color','b');
%     HH2=get(AX(2),'Ylabel');
%     set(HH2,'String','Pr (MPa)');
%     set(HH2,'color','r');
%     set(H1,'LineStyle','-');
%     set(H2,'LineStyle','-');
%     set(H2,'color','r');
%     legend([H1,H2],{'Charge/Discharge Power(MW)';'State of Charge'});
%     xlabel('Time(h)');
%     
%     figure(23)
%     plot(1:NT, value(H_coll_sum)/1e3,'-.bo');
%     hold on
%     plot(1:NT, value(H_cons_sum)/1e3,'-.rs')
%     xlabel('Time (h)')
%     ylabel('Heat (MW.h)')
%     legend('H^{CAES}_g','H^{CAES}_d')
%   
%     figure(24)
%     x = 0:NT;
%     y1 = [0 value(H_coll_sum) - value(H_cons_sum)]/1e3;
%     y2 = [value(H_str0),value(H_str)]/1e3;
%     [AX,H1,H2] = plotyy(x,y1,x,y2,'bar','plot');
%     set(AX(1),'XColor','k','YColor','b');
%     set(AX(2),'XColor','k','YColor','r');
%     set(AX(1),'Xlim',[0,25]);
%     set(AX(2),'Xlim',[0,25]);
%     HH1=get(AX(1),'Ylabel');
%     set(HH1,'String','H^{CAES}(MW)');
%     set(HH1,'color','b');
%     HH2=get(AX(2),'Ylabel');
%     set(HH2,'String','H^{str} (MW.h)');
%     set(HH2,'color','r');
%     set(H1,'LineStyle','-');
%     set(H2,'LineStyle','-');
%     set(H2,'color','r');
%     legend([H1,H2],{'Charge/Discharge (MW)';'State of Charge'});
%     xlabel('Time (h)');
%     
%     figure(25)
%     plot(1:NT,C_grid/1e4,'-b')
%     xlabel('Time (h)')
%     ylabel('Price ($/(kW.h))')
%     ylim([0.04, 0.11])

else
    display('Hmm, something went wrong!');
    sol.info
    yalmiperror(sol.problem)
end

