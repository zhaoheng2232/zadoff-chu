[simCfg,rfCfg,nodeCfg,algCfg,frameCfg] = init_config();
frameCfg.epochLen = 300;
frameCfg.soundingLen = 75;
for spd = [0 1 2 5 10 20]
    simCfg.maxSpeedMph = spd
    simCfg.numLink = 10;
    chanSampTimeSec = cumsum([0 ones(1,frameCfg.epochLen) frameCfg.timeGap ... %epoch 1
        ones(1,frameCfg.epochLen) frameCfg.timeGap ... %epoch 2
        repmat([ones(1,frameCfg.soundingLen) frameCfg.timeGap ones(1,frameCfg.soundingLen) frameCfg.timeGap],1,frameCfg.numIniSlot), ... %Tx Rx BF iteration
        repmat([frameCfg.framLen frameCfg.timeGap ones(1,frameCfg.soundingLen) frameCfg.timeGap],1,frameCfg.numFrame)])/nodeCfg.fsHz;
    
    rateTrace = zeros(1,frameCfg.numIniSlot+frameCfg.numFrame);
    pwrTrace = zeros(1,frameCfg.numIniSlot+frameCfg.numFrame);
    avgNumActiveLink = zeros(2,frameCfg.numIniSlot+2);
    propDelayMat = zeros(simCfg.numLink);
    % loop for location realizations
    fileName = ['./chan_mat/fadingChanMatRtsLen' num2str(frameCfg.epochLen) 'SndLen' num2str(frameCfg.soundingLen) '_' num2str(frameCfg.numIniSlot) 'slots' num2str(simCfg.maxSpeedMph) 'Mph' num2str(rfCfg.maxCfoHz) 'Hz']
    [txLocation, rxLocation] = gen_topology(simCfg); % generate topology of network user pairs
    linkInfo = gen_link_info(simCfg,rfCfg,nodeCfg,txLocation,rxLocation);
    
    %GENERATE PATHLOSS and FADING CHANNEL
    pathLossDbMat = zeros(simCfg.numLink, simCfg.numLink); %pathLossDbMat(i,j) is the path loss between Tx j and Rx i
    fadingChanMat = zeros(simCfg.numLink*simCfg.numLink, nodeCfg.numRx, nodeCfg.numTx, length(chanSampTimeSec));
    nn = 1;
    for indxRx = 1 : simCfg.numLink
        for indxTx  = 1 : simCfg.numLink
            dist = norm(linkInfo(indxTx).txLoc - linkInfo(indxRx).rxLoc);
            fDopplerHz = (linkInfo(indxTx).txSpdMph-linkInfo(indxRx).rxSpdMph)*1.61/3.6/(3e8/rfCfg.fcHz);
            cfoHz = (linkInfo(indxRx).rxCfoHz-linkInfo(indxTx).txCfoHz);
            ch = gen_space_time_chan(chanSampTimeSec,fDopplerHz,nodeCfg.numTx,nodeCfg.numRx,cfoHz);
            fadingChanMat(nn,:,:,:) = ch;
            nn = nn+1;
        end
    end
    save(fileName,'fadingChanMat');
end