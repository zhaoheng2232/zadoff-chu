function [sumRateRec,totTxPwrRec,numALinkRcd,rxRateRecIT] = optimize_iteration(nodeCfg,linkInfo,syncSeqSet,ulPilotSet,plDbMat,propDelayMat,fadingChan,chanFlag,frameCfg,algCfg,simCfg)
%OPTIMIZE_SLOTTED_ALOHA a joint physical and MAC optimization algorithm
% global linkInfo;

maxTxPowDbm = nodeCfg.maxTxPowDbm;
numLink = length(linkInfo);
epochLen = frameCfg.epochLen;
numALinkRcd(1:2,1) = numLink;

%------------- RTS
t = 1;
% if simCfg.runDetectAlg
[detectR,cfoEst] = init_sync_v2(fadingChan(:,:,:,:,t:t+epochLen-1),chanFlag,nodeCfg,frameCfg,algCfg,syncSeqSet,plDbMat,propDelayMat);
% else
%     detectR = ones(1,numLink);
%     cfoEst = [linkInfo(:).rxCfoHz]-[linkInfo(:).txCfoHz];
% end
aSetR = find(detectR);

%-------- CTS
t = t+epochLen+1;
[detectT,txVec] = init_train_tx_v2(fadingChan(:,:,:,:,t:t+epochLen-1),chanFlag,nodeCfg,frameCfg,algCfg,syncSeqSet,aSetR,plDbMat,propDelayMat);
for nn = 1:numLink
    linkInfo(nn).cfoEstHz = cfoEst(nn);
    linkInfo(nn).cfoTruHz = linkInfo(nn).rxCfoHz-linkInfo(nn).txCfoHz;
    if detectT(nn)
        linkInfo(nn).txPwrDbm = maxTxPowDbm;
        linkInfo(nn).txVec = txVec(:,nn);
    else
        linkInfo(nn).txPwrDbm = -inf;
        linkInfo(nn).txVec = zeros(nodeCfg.numTx,1);
    end
end
aSetT = find(detectT);
numALinkRcd(1,2) = length(aSetR);
numALinkRcd(2,2) = length(aSetT);

%----------- training of the IA beamforming vectors
t = t+epochLen+1; %+1 because there is a gap in between
if sum(detectT)==0 || sum(detectR)==0
    sumRateRec = 0; totTxPwrRec= 0; rxRateRecIT = 0;numALinkRcd(:,3:2+frameCfg.numIniSlot) = 0;
    return
    error('Number of Active link = 0 !!')
end
soundingLen = frameCfg.soundingLen;
[linkInfo, numALinkRcd(:,3:2+frameCfg.numIniSlot), totTxPwrRec, rxRateRec, rxRateRecIT, aSetT, aSetR] = iter_train_tx_rx(fadingChan(:,:,:,:,t:t+(2*soundingLen+2)*frameCfg.numIniSlot-1),chanFlag,nodeCfg,frameCfg,algCfg,linkInfo,ulPilotSet,aSetR,aSetT,plDbMat,propDelayMat);
t = t+(frameCfg.soundingLen+1)*frameCfg.numIniSlot*2+1;
[rxRateRec2,totTxPwrRec2,numALinkRcd2] = tx_data_frame(fadingChan(:,:,:,:,t:t+(soundingLen+3)*frameCfg.numFrame-1),chanFlag,nodeCfg,frameCfg,algCfg,linkInfo,ulPilotSet,aSetT,aSetR,plDbMat,propDelayMat);
sumRateRec = sum([rxRateRec rxRateRec2]);
totTxPwrRec = [totTxPwrRec totTxPwrRec2];
numALinkRcd = [numALinkRcd numALinkRcd2];
% run_data_frames
% covArea    = calc_coverage_area(simCfg,linkInfo);
% numActiveLink = sum(~isinf([linkInfo(:).txPwrDbm]));
