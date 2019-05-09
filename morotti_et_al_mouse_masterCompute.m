%% Main file: morotti_et_al_mouse_masterCompute
% This model was built upon the code of the Soltis and Saucerman model
% of rabbit ventricular EC coupling.

% This file calls the ode solver and plots the results.

%%
close all;
clear all; 
clc;

%% Parameters for external modules

% ECC and CaM modules
freq = 1;                   % [Hz] CHANGE DEPENDING ON FREQUENCY
cycleLength = 1e3/freq;     % [ms]
CaMtotDyad = 418;           % [uM]
BtotDyad = 1.54/8.293e-4;   % [uM]
CaMKIItotDyad = 120;        % [uM]
CaNtotDyad = 3e-3/8.293e-4; % [uM]
PP1totDyad = 96.5;          % [uM]
CaMtotSL = 5.65;            % [uM]
BtotSL = 24.2;              % [uM]
CaMKIItotSL = 120*8.293e-4; % [uM]
CaNtotSL = 3e-3;            % [uM]
PP1totSL = 0.57;            % [uM]
CaMtotCyt = 5.65;           % [uM]
BtotCyt = 24.2;             % [uM]
CaMKIItotCyt = 120*8.293e-4;% [uM]
CaNtotCyt = 3e-3;           % [uM] 
PP1totCyt = 0.57;           % [uM]

% ADJUST CaMKII ACTIVITY LEVELS (expression = 'WT', 'OE', or 'KO')
expression = 'WT';
CKIIOE = 0; % Should be zero during 'WT' and 'KO' runs

if strcmp(expression,'OE') % OE
    CKIIOE = 1; % Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
    n_OE = 6;
    CaMKIItotDyad = 120*n_OE;        % [uM] 
    CaMKIItotSL = 120*8.293e-4*n_OE; % [uM]
    CaMKIItotCyt = 120*8.293e-4*n_OE;% [uM]
elseif strcmp(expression,'KO')
    CaMKIItotDyad = 0;          % [uM] 
    CaMKIItotSL = 0;            % [uM]
    CaMKIItotCyt = 0;           % [uM]
end

%plb_val=38; % RABBIT
plb_val=106; % MOUSE

% Parameters for CaMKII module
LCCtotDyad = 31.4*.9;       % [uM] - Total Dyadic [LCC] - (umol/l dyad)
LCCtotSL = 0.0846;          % [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
RyRtot = 382.6;             % [uM] - Total RyR (in Dyad)
PP1_dyad = 95.7;            % [uM] - Total dyadic [PP1]
PP1_SL = 0.57;              % [uM] - Total Subsarcolemmal [PP1]
PP2A_dyad = 95.76;          % [uM] - Total dyadic PP2A
OA = 0;                     % [uM] - PP1/PP2A inhibitor Okadaic Acid
PLBtot = plb_val;           % [uM] - Total [PLB] in cytosolic units

% Parameters for BAR module
Ligtot = 0                  % [uM] - SET LIGAND CONCENTRATION HERE 0-0.1
LCCtotBA = 0.025;           % [uM] - [umol/L cytosol]
RyRtotBA = 0.135;           % [uM] - [umol/L cytosol]
PLBtotBA = plb_val;         % [uM] - [umol/L cytosol]
TnItotBA = 70;              % [uM] - [umol/L cytosol]
IKstotBA = 0.025;           % [uM] - [umol/L cytosol]
ICFTRtotBA = 0.025;         % [uM] - [umol/L cytosol]
PP1_PLBtot = 0.89;          % [uM] - [umol/L cytosol]
IKurtotBA = 0.025;          % [uM] - [umol/L cytosol] MOUSE
PLMtotBA = 48;              % [uM] - [umol/L cytosol] MOUSE

% For Recovery from inactivation of LTCC
recoveryTime = 10; % initialize to smallest value

% Parameter varied in protocol simulation
variablePar = 20; % initilization

%% Collect all parameters and define mass matrix for BAR module
p = [cycleLength,recoveryTime,variablePar,CaMtotDyad,BtotDyad,CaMKIItotDyad,CaNtotDyad,PP1totDyad,...
    CaMtotSL,BtotSL,CaMKIItotSL,CaNtotSL,PP1totSL,...
    CaMtotCyt,BtotCyt,CaMKIItotCyt,CaNtotCyt,PP1totCyt...
    LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,LCCtotSL,PP1_SL...
    Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot,IKurtotBA,PLMtotBA...
    CKIIOE];

M = diag(ones(1,87+6+45+6+36));

%% Establish and define globals variables
global tStep tArray I_Ca_store I_to_store I_Na_store I_K1_store ibar_store 
global gates Jserca IKs_store Jleak ICFTR Incx
global I_ss_store dVm_store Ipca_store I_NaK_store I_Nabk_store I_kr_store
global I_kur1_store I_kur2_store
global Lmyo_store Fmyo_store
 
tStep = 1;
tArray = zeros(1,1e6);
I_Ca_store=zeros(1,1e6);
I_to_store=zeros(3,1e6);
I_Na_store = zeros(1,1e6);
I_K1_store = zeros(1,1e6);
ibar_store=zeros(1,1e6);
gates = zeros(2,1e6);
Jserca = zeros(1,1e6);
IKs_store = zeros(1,1e6);
Jleak = zeros(1e6,2);
ICFTR = zeros(1,1e6);
Incx = zeros(1,1e6);
I_kur1_store = zeros(1,1e6);
I_kur2_store = zeros(1,1e6);
I_ss_store = zeros(1,1e6);
dVm_store = zeros(1,1e6);
Ipca_store = zeros(1,1e6);
I_NaK_store = zeros(1,1e6);
I_Nabk_store = zeros(1,1e6);
I_kr_store = zeros(1,1e6);

Lmyo_store = zeros(1,1e6);
Fmyo_store = zeros(1,1e6);

%% Load initial conditions
% Isotonic contractions 

% Ligtot = 0 (steady-state @ 1 Hz)
load yfin_mouse_myofil_isoT_control_1Hz         % control

% Ligtot = 0.1 (120 s @ 1 Hz)
%load yfin_mouse_myofil_isoT_IBMX_1Hz_120s      % 100% cAMP @ all targets
%load yfin_mouse_myofil_isoT_ISOall_1Hz_120s    % -50% cAMP @ XBCa, XBcy, titin
%load yfin_mouse_myofil_isoT_ISOxbca_1Hz_120s   % -50% cAMP @ XBCa
%load yfin_mouse_myofil_isoT_ISOxbcy_1Hz_120s   % -50% cAMP @ XBcy
%load yfin_mouse_myofil_isoT_ISOtitin_1Hz_120s  % -50% cAMP @ titin

y0n = yfinal;

%% Run simulation
tic
tspan = [0 10e3]; % [ms]
options = odeset('RelTol',1e-5,'MaxStep',2);
[t,y] = ode15s(@morotti_et_al_mouse_masterODEfile,tspan,y0n,options,p);
toc

%% Save final conditions
yfinal = y(end,:)';

% Isotonic contractions 

% Ligtot = 0 (steady-state @ 1 Hz)
%save yfin_mouse_myofil_isoT_control_1Hz yfinal         % control

% Ligtot = 0.1 (120 s @ 1 Hz)
%save yfin_mouse_myofil_isoT_IBMX_1Hz_120s yfinal       % 100% cAMP  @ all targets
%save yfin_mouse_myofil_isoT_ISOall_1Hz_120s yfinal     % -50% cAMP @ XBCa, XBcy, titin
%save yfin_mouse_myofil_isoT_ISOxbca_1Hz_120s yfinal    % -50% cAMP @ XBCa
%save yfin_mouse_myofil_isoT_ISOxbcy_1Hz_120s yfinal    % -50% cAMP @ XBcy
%save yfin_mouse_myofil_isoT_ISOtitin_1Hz_120s yfinal   % -50% cAMP @ titin

%% Output variables
tArray = tArray(1:tStep);
Ica = I_Ca_store(1:tStep);
Ito = I_to_store(1,1:tStep);
Itof = I_to_store(2,1:tStep);
Itos = I_to_store(3,1:tStep);
INa = I_Na_store(1:tStep);
IK1 = I_K1_store(1:tStep);
s1 = gates(1,1:tStep);
k1 = gates(2,1:tStep);
Jserca = Jserca(1:tStep);
Iks = IKs_store(1:tStep);
Jleak = Jleak(1:tStep,:);
ICFTR = ICFTR(1:tStep);
Incx = Incx(1:tStep);
Ikur1 = I_kur1_store(1:tStep);
Ikur2 = I_kur2_store(1:tStep);
Iss = I_ss_store(1:tStep);
dVm = dVm_store(1:tStep);
Ipca = Ipca_store(1:tStep);
INaK = I_NaK_store(1:tStep);
INabk = I_Nabk_store(1:tStep);
Ikr = I_kr_store(1:tStep);

Lm = Lmyo_store(1:tStep);
Fm = Fmyo_store(1:tStep);

%% Plot results
ts = t/1e3; % time in s 
tss = tArray./1e3; % time in s

Cai = y(:,38);
Vm = y(:,39);
cAMPtot = y(:,87+6+45+6+11);

figure,set(gcf,'color','w')
subplot(2,2,1),hold on, plot(ts,Vm,'k'); ylabel('Em (mV)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(2,2,2),hold on; plot(ts,cAMPtot,'k'); ylabel('cAMPtot (uM)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(2,2,3),hold on, plot(ts,Cai,'k'); ylabel('[Ca]i (mM)');
xlabel('Time (s)'); set(gca,'box','off','tickdir','out','fontsize',12)
subplot(2,2,4),hold on, plot(tss,Lm,'k'); ylabel('Length (um)');
xlabel('Time (s)'); set(gca,'box','off','tickdir','out','fontsize',12)