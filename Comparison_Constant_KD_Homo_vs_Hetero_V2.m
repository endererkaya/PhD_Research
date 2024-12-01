% Received Signal Strength
generated_bits = randi([0 1], 1, Nbits);
temp_bits  = generated_bits;
temp_bits(generated_bits==1) = phi;
temp_bits(generated_bits==0) = phi2;

%% Add Background Noise
% Add Noise(External Signal)
% noise = random('LogNormal', log(1e20), 1, size(t)); % noise concentration at point RX
% noise = repelem(random('LogNormal', log(1e20), 1.2, [1 L+1]), Ns);
% sigma_noise  = log(1e5); %% OPEN IF RUN THE SCRIPT ONLY
% mu_noise     = 0;
% mean_noise   = exp(mu_noise + ((sigma_noise^2)/2));
% std_noise    = sqrt(exp(sigma_noise^2));

% noise = lognrnd(mu_noise, sigma_noise, 1, Nbits);
lognoise = mu_noise + sigma_noise * randn(1, Nbits);
noise = exp(lognoise);
% % Received RX Point Concentration Response
% yreceived = IMreceived + noise; %received concentration at point RX
% yreceived = IMreceived(2:end);
% ytemp = yreceived(Rtime:Ns:end-9*Ns);
% ysampled  =  exp(log(temp_bits) + lognoise);
ysampled = temp_bits + noise;

%% Qs
% tempQ = [Q/10 Q];
%% Single Homogeneous Receiver
% Set KD

yresponse1 = zeros(1,Nbits);
% ydetect1   = zeros(1,Nbits);
% detection_threshold1 = zeros(1,Nbits);

% Detection
for it = 1:Nbits
    yresponse1(it) = binornd(NR,ligand_response(KD, ysampled(it)));
    %     if it == 1
    %         detection_threshold1(it) = NR * 0.5;
    %         ydetect1(it) = yresponse1(it) > detection_threshold1(it);
    %     else
    %         if it == 2
    %             ISI_estimate = tempQ(ydetect1(it-1)+1) * h2;
    %         elseif it == 3
    %             ISI_estimate = tempQ(ydetect1(it-1)+1) * h2 + tempQ(ydetect1(it-2)+1) * h3;
    %         else
    %             ISI_estimate = tempQ(ydetect1(it-1)+1) * h2 + tempQ(ydetect1(it-2)+1) * h3 + ISI_constant;
    %         end
    %         ISI_estimate = 1;
    %         c1 = ISI_estimate + phi;%mu_noise
    %         c0 = ISI_estimate + phi2;%mu_noise
    %         detection_threshold1(it) = log((c1 + KD)/(c0 + KD))...
    %                                   /log((c1)/(c0));
    %         detection_threshold1(it) = find_homogeneous_detection_threshold(log(c0),log(c1),KD,mu_noise,sigma_noise,NR);
    %         ydetect1(it) = yresponse1(it) > detection_threshold1(it);
    %     end
end

ydetect1 = yresponse1 > detection_threshold1;

%% Probability_Error
prob_error1 = sum((xor(ydetect1, logical(generated_bits))))/Nbits;
C1 = confusionmat(logical(generated_bits),ydetect1);
%% Receiver Heterogeneous-2 KD

% Set KDs
% KD1 = KD;
% KD2 = KD * 100;
% tau  = (KD1 + KD2)/2;
% mulq = KD1 * KD2;

% Function Handle
% det_thresh = @(c1,c0) (log((c1 + KD1)/(c0 + KD1)) ...
%     + log((c1 + KD2)/(c0 + KD2))...
%     - log((c1*tau + mulq)/(c0*tau + mulq)))...
%     / (log(c1/c0) + log((c1 + tau)/(c0 + tau)) ...
%     - log((c1*tau + mulq)/(c0*tau + mulq)));

yresponse2 = zeros(1,Nbits);
% ydetect2   = zeros(1,Nbits);
% detection_threshold2 = zeros(1,Nbits);

% ISI_estimate2 = ISI_constant;
% c1 = ISI_estimate2 + phi;%%% mean_noise ???? %%%
% c0 = ISI_estimate2 + phi2; %%% ???? %%%
% detection_threshold2 = find_heterogeneous_detection_threshold(log(c0),log(c1),KD1,KD2,mu_noise,sigma_noise,NR);

% Received Ligand Receptor Response
% yresponse2 = NR * ligand_response2(KD1, KD2, 0.5, ysampled);

% KD1 = KD / 10 ;
% KD2 = KD * 10;
% Detection
for it = 1:Nbits
    yresponse2(it) = binornd(NR/2, ligand_response(KD1,ysampled(it))) + ...
        binornd(NR/2, ligand_response(KD2,ysampled(it)));
    %     if it == 1
    %         detection_threshold2(it) = NR * 0.5;
    %         ydetect2(it) = yresponse2(it) > detection_threshold2(it);
    %     else
    %         if it == 2
    %             ISI_estimate2 = tempQ(ydetect2(it-1)+1) * h2;
    %         elseif it == 3
    %             ISI_estimate2 = tempQ(ydetect2(it-1)+1) * h2 + tempQ(ydetect2(it-2)+1) * h3;
    %         else
    %             ISI_estimate2 = tempQ(ydetect2(it-1)+1) * h2 + tempQ(ydetect2(it-2)+1) * h3 + ISI_constant;
    %         end
    %         ISI_estimate2 = 1;
    %         c1 = ISI_estimate2 + phi;%%% mean_noise ???? %%%
    %         c0 = ISI_estimate2 + phi2; %%% ???? %%%
    %         detection_threshold2(it) = find_heterogeneous_detection_threshold(log(c0),log(c1),KD1,KD2,mu_noise,sigma_noise,NR);
    %         ydetect2(it) = yresponse2(it) > detection_threshold2(it);
    %     end
end

ydetect2 = yresponse2 > detection_threshold2;

%% Probability_Error
prob_error2 = sum((xor(ydetect2, logical(generated_bits))))/Nbits;
C2 = confusionmat(logical(generated_bits),ydetect2);

%% Receiver Homogeneous with Hill Coefficient

% yresponseH = zeros(1,Nbits);
% % ydetect1   = zeros(1,Nbits);
% % detection_threshold1 = zeros(1,Nbits);
% 
% % Detection
% for it = 1:Nbits
%     yresponseH(it) = binornd(NR,ligand_response(KDH, ysampled(it),hill));
% end
% ydetectH = yresponseH > detection_threshold_hill;
% 
% prob_errorH = sum((xor(ydetectH, logical(generated_bits))))/Nbits;
% CH = confusionmat(logical(generated_bits),ydetectH);

%% Receiver Heterogeneous-3 KD

% % Set KDs
% KD1 = KD / 1e2;
% KD2 = KD;
% KD3 = KD * 1e2;

% yresponse3 = zeros(1,Nbits);
% ydetect3   = zeros(1,Nbits);
% detection_threshold3 = zeros(1,Nbits);

% ISI_estimate3 = ISI_constant;
% c1 = ISI_estimate3 + phi;%%% mean_noise ???? %%%
% c0 = ISI_estimate3 + phi2; %%% ???? %%%
% detection_threshold3 = find_heterogeneous3_detection_threshold(log(c0),log(c1),KD,mu_noise,sigma_noise,NR);

% KD1 = KD / 1e2;
% KD2 = KD;
% KD3 = KD * 1e2;
%
% % Detection
% for it = 1:Nbits
%     yresponse3(it) = binornd(NR/3, ligand_response(KD1,ysampled(it))) + ...
%                      binornd(NR/3, ligand_response(KD2,ysampled(it))) + ...
%                      binornd(NR/3, ligand_response(KD3,ysampled(it)));
% %     if it == 1
% %         detection_threshold3(it) = NR * 0.5;
% %         ydetect3(it) = yresponse3(it) > detection_threshold3(it);
% %     else
% %         if it == 2
% %             ISI_estimate3 = tempQ(ydetect3(it-1)+1) * h2;
% %         elseif it == 3
% %             ISI_estimate3 = tempQ(ydetect3(it-1)+1) * h2 + tempQ(ydetect3(it-2)+1) * h3;
% %         else
% %             ISI_estimate3 = tempQ(ydetect3(it-1)+1) * h2 + tempQ(ydetect3(it-2)+1) * h3 + ISI_constant;
% %         end
% %         ydetect3(it) = yresponse3(it) > detection_threshold3(it);
% %     end
% end
%
% ydetect3 = yresponse3 > detection_threshold3;

%% Probability_Error
% prob_error3 = sum((xor(ydetect3, logical(generated_bits))))/Nbits;

% %% Receiver Heterogeneous-4 KD
%
% % Set KDs
% KD1 = KD / 1e3;
% KD2 = KD / 1e1;
% KD3 = KD * 1e1;
% KD4 = KD * 1e3;
%
% yresponse4 = zeros(1,Nbits);
% ydetect4   = zeros(1,Nbits);
% detection_threshold4 = zeros(1,Nbits);
%
% % Detection
% for it = 1:Nbits
%     yresponse4(it) = binornd(NR, ligand_response4(KD1,KD2,KD3,KD4,ysampled(it)));
%     if it == 1
%         detection_threshold = NR * 0.5;
%         ydetect4(it) = yresponse4(it) > detection_threshold;
%     else
%         if it == 2
%             ISI_estimate = Q * ydetect4(it-1) * h2;
%         elseif it == 3
%             ISI_estimate = Q * (ydetect4(it-1) * h2 + ydetect4(it-2) * h3);
%         else
%             ISI_estimate = Q * (ydetect4(it-1) * h2 + ydetect4(it-2) * h3) + ISI_constant;
%         end
% %         ISI_estimate = 1;
%         c1 = ISI_estimate + phi;%%% mean_noise ???? %%%
%         c0 = ISI_estimate; %%% ???? %%%
%         detection_threshold4(it) = find_heterogeneous4_detection_threshold(log(c0),log(c1),KD,mu_noise,sigma_noise,NR);
%         ydetect4(it) = yresponse4(it) > detection_threshold4(it);
%     end
% end
%
% %% Probability_Error
% prob_error4 = sum((xor(ydetect4, logical(generated_bits))))/Nbits;
%
% %% Figuring
% % figure
% % plot(t(2:end)*1e1,IMreceived(2:end),'r','LineWidth',2);
% % hold on
% % plot(t(2:end)*1e1,yreceived,'b','LineWidth',2);
% % grid on
% % xlabel('Time(ms)','Interpreter','latex');
% % ylabel('IM Received Signal','Interpreter','latex');
% %
% % figure
% % plot(t(2:end)*1e1,yresponse,'r','LineWidth',2);
% % hold on
% % stem((Ns/2:Ns:Ns/2+(L-1)*Ns)/1e3,yresponse_sum_stat,'ob','LineWidth',3)
% % grid on
% % xlabel('Time');
% % ylabel('Receptor Occupation')