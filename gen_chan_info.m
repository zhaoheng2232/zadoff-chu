function [pathLossDbMat,fadingChan,propDelayMat,chanSampTimeSec] = gen_chan_info(simCfg,nodeCfg,frameCfg,rfCfg,linkInfo)

chanSampTimeSec = cumsum([0 ...
    ones(1,frameCfg.epochLen) frameCfg.timeGap ... % RTS
    ones(1,frameCfg.epochLen) frameCfg.timeGap ... % CTS
    repmat([ones(1,frameCfg.soundingLen) frameCfg.timeGap ones(1,frameCfg.soundingLen) frameCfg.timeGap],1,frameCfg.numIniSlot), ... %Tx Rx BF iteration
    repmat([frameCfg.framLen frameCfg.timeGap ones(1,frameCfg.soundingLen) frameCfg.timeGap],1,frameCfg.numFrame)])/nodeCfg.fsHz; % Data frame
pathLossDbMat = zeros(simCfg.numLink, simCfg.numLink); %pathLossDbMat(i,j) is the path loss between Tx j and Rx i
fadingChan = zeros(simCfg.numLink, simCfg.numLink, nodeCfg.numRx, nodeCfg.numTx, length(chanSampTimeSec));
propDelayMat = zeros(simCfg.numLink,simCfg.numLink);

for indxRx = 1 : simCfg.numLink
    for indxTx  = 1 : simCfg.numLink
        dist = norm(linkInfo(indxTx).txLoc - linkInfo(indxRx).rxLoc);
        propDelayMat(indxRx,indxTx) = dist/3e8*rfCfg.bwHz;
        pathLossDbMat(indxRx, indxTx) = 2*rfCfg.anteGainDbi + rfCfg.pathloss1m - 10*rfCfg.plExponent*log10(dist);
        % Generate fading channel
%         ch = squeeze(fadingChanMat(rndInd(cc),:,:,:));
%         fadingChan(indxRx,indxTx,:,:,:) = ch(:,:,1:size(fadingChan,5));
        
        fDopplerHz = (linkInfo(indxTx).txSpdMph-linkInfo(indxRx).rxSpdMph)*1.61/3.6/(3e8/rfCfg.fcHz);
        cfoHz = (linkInfo(indxRx).rxCfoHz-linkInfo(indxTx).txCfoHz);
        fadingChan(indxRx,indxTx,:,:,:) = gen_space_time_chan(chanSampTimeSec,fDopplerHz,nodeCfg.numTx,nodeCfg.numRx,cfoHz);
    end
end

end

