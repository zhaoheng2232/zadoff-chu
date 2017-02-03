clc; 
% close all;
clear;

[simCfg,rfCfg,nodeCfg,algCfg,frameCfg] = init_config();
conops = 1;
if conops==1
    simCfg.topologyType = 'onCircle';'onCircle';'packedCircle';
    simCfg.areaRadiusM = 500;
else
    simCfg.topologyType = 'packedCircle';
    simCfg.areaRadiusM = 1000;
end
algCfg.maxSubactiveFrame = 5; inf;10;
algCfg.dtCts = 0.1;
algCfg.dtRts = 0.1;
algCfg.stTrainSlot = 0.1; %shutdown threshold in the training slot
algCfg.estChan = 1;
algCfg.compCFO = 1;

frameCfg.epochLen = 300;%250;%;1050;%500;
frameCfg.syncSeqLen = 251;%229;%1021;%479;
frameCfg.syncSeqCpLen = frameCfg.epochLen-frameCfg.syncSeqLen;
frameCfg.epch1LenUs = frameCfg.epochLen/rfCfg.bwHz*1e6;
frameCfg.epch2LenUs = frameCfg.epch1LenUs;
frameCfg.timeGapUs = 100;
frameCfg.framLenUs = 5000;
frameCfg.framLen= frameCfg.framLenUs*rfCfg.bwHz/1e6;
frameCfg.timeGap = frameCfg.timeGapUs/nodeCfg.fsHz*1e6;
frameCfg.soundingLen = 75;%80;% 100; %220;
frameCfg.ulPilotLen = 61;%61;% 79;%199;
frameCfg.ulPilotCpLen = frameCfg.soundingLen-frameCfg.ulPilotLen;


simCfg.runDetectAlg = 1;
frameCfg.numFrame = 10;
simCfg.maxSpeedMph = 2;
ptrn = {'r','g','b'};


tmp = load(['syncSeq' num2str(frameCfg.syncSeqLen) '.mat'],'syncSeqSet');
syncSeqSet = tmp.syncSeqSet;
tmp = load(['syncSeq' num2str(frameCfg.ulPilotLen) '.mat'],'syncSeqSet');
ulPilotSet = tmp.syncSeqSet;
clear tmp;

numLinkSet = 10;[10 20 30 40]; [2 3 4 5 8 10 12];
numMonteCarlo = 100;
matFileName = ['./fadingChanMatRtsLen' num2str(frameCfg.epochLen) 'SndLen' num2str(frameCfg.soundingLen) '_' num2str(frameCfg.numIniSlot) 'slots' num2str(simCfg.maxSpeedMph) 'Mph' num2str(rfCfg.maxCfoHz) 'Hz'];
load(matFileName);
% frameCfg.numIniSlot = 15;
chanSampTimeSec = cumsum([0 ones(1,frameCfg.epochLen) frameCfg.timeGap ... %epoch 1
    ones(1,frameCfg.epochLen) frameCfg.timeGap ... %epoch 2
    repmat([ones(1,frameCfg.soundingLen) frameCfg.timeGap ones(1,frameCfg.soundingLen) frameCfg.timeGap],1,frameCfg.numIniSlot), ... %Tx Rx BF iteration
    repmat([frameCfg.framLen frameCfg.timeGap ones(1,frameCfg.soundingLen) frameCfg.timeGap],1,frameCfg.numFrame)])/nodeCfg.fsHz;

for indxNumLink= 1 : length(numLinkSet)
    simCfg.numLink = numLinkSet(indxNumLink);
    rateTrace = zeros(5,frameCfg.numIniSlot+frameCfg.numFrame);
    numActiveLinkTrace = zeros(5,frameCfg.numIniSlot+frameCfg.numFrame+2);
    pwrTrace = zeros(1,frameCfg.numIniSlot+frameCfg.numFrame);
    avgNumActiveLink = zeros(2,frameCfg.numIniSlot+2);
    propDelayMat = zeros(simCfg.numLink);
    % loop for location realizations
    for indxMC = 1 : numMonteCarlo
        if mod(indxMC,2)==0
            fprintf('**Network User pairs: %d**, location number %d \n', simCfg.numLink, indxMC);
        end
        [txLocation, rxLocation] = gen_topology(simCfg); % generate topology of network user pairs
        linkInfo = gen_link_info(simCfg,rfCfg,nodeCfg,txLocation,rxLocation);
        
        %GENERATE PATHLOSS and FADING CHANNEL
        pathLossDbMat = zeros(simCfg.numLink, simCfg.numLink); %pathLossDbMat(i,j) is the path loss between Tx j and Rx i
        fadingChan = zeros(simCfg.numLink, simCfg.numLink, nodeCfg.numRx, nodeCfg.numTx, length(chanSampTimeSec));
        chanFlag = zeros(simCfg.numLink, simCfg.numLink);
        cc = 0; rndInd = randperm(100);
        for indxRx = 1 : simCfg.numLink
            for indxTx  = 1 : simCfg.numLink
                dist = norm(linkInfo(indxTx).txLoc - linkInfo(indxRx).rxLoc);
                propDelayMat(indxRx,indxTx) = dist/3e8*rfCfg.bwHz;
                pathLossDbMat(indxRx, indxTx) = calc_path_loss(rfCfg.anteGainDbi,rfCfg.pathloss1m,rfCfg.plExponent,rfCfg.noiseFloorDbm,dist);
                if pathLossDbMat(indxRx, indxTx)+nodeCfg.maxTxPowDbm>-5
                    chanFlag(indxRx,indxTx) = 1;
                    cc = mod(cc+1,100)+1;
                    ch = squeeze(fadingChanMat(rndInd(cc),:,:,:));
                    fadingChan(indxRx,indxTx,:,:,:) = ch(:,:,1:size(fadingChan,5));
                else
                    chanFlag(indxRx,indxTx) = 0;
                    fadingChan(indxRx,indxTx,:,:,:) = 0;
                end
            end
        end
        [sumRateRec,totTxPwrRec,numActiveLinkRec,sumRateRecIT] = optimize_iteration(nodeCfg,linkInfo,syncSeqSet,ulPilotSet(1:simCfg.numLink,:),pathLossDbMat,propDelayMat,fadingChan,chanFlag,frameCfg,algCfg,simCfg);
        rateTrace(1,:) = rateTrace(1,:)+sumRateRec/numMonteCarlo*rfCfg.bwHz/1e6;
        numActiveLinkTrace(1,:) = numActiveLinkTrace(1,:)+numActiveLinkRec(1,:)/numMonteCarlo;
        %-- test impact of channel estimation
%         for indxChanEst = 1:2
%             algCfg.estChan = indxChanEst-1;
%             [sumRateRec,totTxPwrRec,numActiveLinkRec] = optimize_iteration(nodeCfg,linkInfo,syncSeqSet,ulPilotSet(1:simCfg.numLink,:),pathLossDbMat,propDelayMat,fadingChan,chanFlag,frameCfg,algCfg,simCfg);
%             rateTrace(indxChanEst,:) = rateTrace(indxChanEst,:)+sumRateRec/numMonteCarlo*rfCfg.bwHz/1e6;
%         end
        
%         %-- test impact of maxSubActiveFrame
%         maxSubactiveFrameSet = [0 5 inf];
%         for indxSubActive = 1:3            
%             algCfg.maxSubactiveFrame = maxSubactiveFrameSet(indxSubActive);
%             [sumRateRec,totTxPwrRec,numActiveLinkRec] = optimize_iteration(nodeCfg,linkInfo,syncSeqSet,ulPilotSet(1:simCfg.numLink,:),pathLossDbMat,propDelayMat,fadingChan,chanFlag,frameCfg,algCfg,simCfg);
%             rateTrace(indxSubActive,:) = rateTrace(indxSubActive,:)+sumRateRec/numMonteCarlo*rfCfg.bwHz/1e6;
%             numActiveLinkTrace(indxSubActive,:) = numActiveLinkTrace(indxSubActive,:)+numActiveLinkRec(1,:)/numMonteCarlo;
%         end

%         %-- test impact of detection threshold
%         thresholdSet = [0 0.1 .2 .3 .5];
%         for indxThreshold = 1:4            
%             algCfg.dtRts = thresholdSet(indxThreshold);
%             algCfg.dtCts = thresholdSet(indxThreshold);
%             [sumRateRec,totTxPwrRec,numActiveLinkRec] = optimize_iteration(nodeCfg,linkInfo,syncSeqSet,ulPilotSet(1:simCfg.numLink,:),pathLossDbMat,propDelayMat,fadingChan,chanFlag,frameCfg,algCfg,simCfg);
%             rateTrace(indxThreshold,:) = rateTrace(indxThreshold,:)+sumRateRec/numMonteCarlo*rfCfg.bwHz/1e6;
%             numActiveLinkTrace(indxThreshold,:) = numActiveLinkTrace(indxThreshold,:)+numActiveLinkRec(1,:)/numMonteCarlo;
%         end

    end    
end
rateTrace = rateTrace/rfCfg.bwHz*20e6;
%%
% figure
% clrSet = {'r','g','b','m','c'};
% for n = 1:2
%     plot(rateTrace(n,:),clrSet{n},'LineWidth',2), hold on
% end
% ylabel('Aggregated Throughput (Mbps)')
% xlabel('Iteration #')
% grid on
% title([num2str(simCfg.numLink) '-link MANET, CONOPS #' num2str(conops) '; 10Mph']);
% legend('T_{rts,cts} = 0','T_{rts,cts} = 0.1','T_{rts,cts} = 0.2','T_{rts,cts} = 0.3','T_{rts,cts} = 0.5')

% 
% plot(rateTrace(1,:),'LineWidth',2), hold on
% ylabel('Aggregated Throughput (Mbps)')
% xlabel('Iteration #')
% grid on
% title([num2str(simCfg.numLink) '-link MANET; 10 Mph']);
% legend('Perfect Channel','Estimated Channel')
% 
% plot(numActiveLinkTrace(1,:),'LineWidth',2), hold on
% ylabel('Number of Active Links')
% xlabel('Iteration #')
% grid on
% title([num2str(simCfg.numLink) '-link MANET; 10 Mph']);
% legend('Perfect Channel','Estimated Channel')