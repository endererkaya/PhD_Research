%% Single Receptor Type Dynamic-Varying KD
%
% includes mobility - Diffusion TX-RX
% includes external molecule same with IM(same KD)
% noise on concentration at RX location
%
% Concentration >--mapped to(Sigmoid(KD))--> Receptor Occupation Ratio
% VARYING KD(DISSOCIATION CONSTANT)
%
% Detection with sum(ratio_occupied in one transmission period)~~receptor bound time
%
% Ender Erkaya
%
% April 2022
clearvars
%% Settings
Nbits   = 1000;
generated_bits = randi([0 1], 1, Nbits);

% Dissociation Constant of Receptor
KDinit  = 1e21;

% Diffusion Coefficients
DA  = 5e-9;
DRX = 0; %1e-11;
DTX = 0; %1e-11;

% Flow Coefficient
xflow = 0;
yflow = 0;
zflow = 0;

% Initial Receiver Positions
x0    = 5e-6;%um
y0    = 0;
z0    = 0;

% Symbol Interval 
Ts = 1e-2;%10ms
Q  = 1e7; % num of released molecules

% Discretization step
deltat = 1e-4; %10us
Tmax   = (Nbits+1)*Ts; %ms
Ns     = round(Ts/deltat);

% Effective Diffusion Coefficients
D1     = DA + DRX;
D2     = DRX + DTX;

% Initial receiver distance
r0     = sqrt(x0^2+y0^2+z0^2);
% arx    = 0.15*1e-6; % um
% Vr     = (4/3)*pi*arx^3;

% Function Handle
htau = @(tau,r) (4*pi*D1*tau).^(-3/2) .* exp(-(norm(r)^2)./(4*D1*tau));

%% Transmit Signal
t      = 0:deltat:Tmax;
H      = zeros(Nbits,length(t));

for s = 1:Nbits
    transmit_time = (s-1)*Ts;
    if s == 1
        r = r0;
    else 
        xr = x0 + sqrt(2*D2*transmit_time)*randn;
        yr = y0 + sqrt(2*D2*transmit_time)*randn;
        zr = z0 + sqrt(2*D2*transmit_time)*randn;
        r  = sqrt(xr^2+yr^2+zr^2);
    end
    tau = deltat:deltat:Ts*(Nbits+2-s);
    H(s,(s-1)*Ns+2:end) = htau(tau,r);
end

% Received Signal Strength
IMreceived = sum(diag(Q*generated_bits)*H,1); %IM concentration at point RX

%% Add Background Noise
% Add Noise(External Signal)
% noise = random('LogNormal', log(1e20), 1, size(t)); % noise concentration at point RX
noise = random('LogNormal', log(1e20), 0.25);
for k = 1:Nbits
    ntemp = random('LogNormal', log(noise(end)), 0.1);
    noise = [noise ntemp];
end
noise = repelem(noise, Ns);

% % Received RX Point Concentration Response
% yreceived = IMreceived + noise; %received concentration at point RX
yreceived = IMreceived(2:end) + noise;

%% Receiver  
% Received Ligand Receptor Response
% ADAPTIVE KD
FIR_filter = 0.8.^(0:19);
FIR_filter = FIR_filter/sum(FIR_filter); % normalize

KD        = zeros(1,Nbits);
yresponse = zeros(1,Nbits);
KD(1)        = KDinit;
yresponse(1) = mean(ligand_response(KDinit, yreceived(1:Ns)));
KDcurrent    = KDinit;
for k = 2:Nbits
    % ADAPTIVE KD
    temp = sum(yresponse(k-1:-1:k-min(k-1,length(FIR_filter))) .* FIR_filter(1:min(k-1,length(FIR_filter))));
    smoothed_occupied_ratio = temp;
    % Update KD
    KDend = inverse_ligand_response(KDcurrent,smoothed_occupied_ratio);
    KD(k) = KDend;    
    if k > 20
        KDcurrent = KDend
    end
    
    % Ligand Response
    ytemp = mean(ligand_response(KDcurrent, yreceived(Ns*(k-1)+1:Ns*k)));
    yresponse(k) = ytemp;
end

% Detection
detection_threshold = 0.45;
ydetect   = yresponse > detection_threshold;

%% Probability_Error
error      = xor(ydetect, logical(generated_bits));
prob_error = sum(error)/Nbits

%% Figuring 
figure
plot(t(2:end),IMreceived(2:end),'r','LineWidth',2);
hold on
plot(t(2:end),yreceived,'b','LineWidth',2);
hold on
plot(t(2:end-Ns),repelem(KD,Ns),'k','LineWidth',2);
grid on
xlabel('Time(ms)','Interpreter','latex');
ylabel('IM Received Signal','Interpreter','latex');

figure
stem(1:Nbits,yresponse,'ob','LineWidth',3)
hold on
plot(1:Nbits,xor(ydetect, logical(generated_bits)),'r')
grid on
xlabel('Number of Bits');
ylabel('Receptor Occupation')