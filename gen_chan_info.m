function [pathLossDbMat,fadingChan,propDelayMat] = gen_chan_info(simCfg,nodeCfg,frameCfg,rfCfg,linkInfo)

numChanSamp = (frameCfg.epochLen + frameCfg.timeGap) + ...  % Epoch 1
    (frameCfg.epochLen + frameCfg.timeGap) + ...    % Epoch 2
    (frameCfg.soundingLen + frameCfg.timeGap)*2*frameCfg.numIniSlot + ...   % Beamforming training
    (frameCfg.framLen + frameCfg.timeGap + frameCfg.soundingLen + frameCfg.timeGap)*frameCfg.numFrame;  % Data Frame
    

pathLossDbMat = zeros(simCfg.numLink, simCfg.numLink); %pathLossDbMat(i,j) is the path loss between Tx j and Rx i
fadingChan = zeros(simCfg.numLink, simCfg.numLink, nodeCfg.numRx, nodeCfg.numTx, numChanSamp);
propDelayMat = zeros(simCfg.numLink,simCfg.numLink);

switch simCfg.fadingType
    case 'Rayleigh'
        chan = comm.MIMOChannel(...
                'SampleRate', nodeCfg.fsHz,...
                'FadingDistribution', simCfg.fadingType,...
                'MaximumDopplerShift', 0,...
                'SpatialCorrelation', false,...
                'NumTransmitAntennas', nodeCfg.numTx,...
                'NumReceiveAntennas', nodeCfg.numRx,...
                'PathGainsOutputPort', true);
    case 'Rician'   
        chan = comm.MIMOChannel(...
                'SampleRate', nodeCfg.fsHz,...
                'FadingDistribution', simCfg.fadingType,...
                'KFactor', simCfg.fadingKFactor,...
                'MaximumDopplerShift', 0,...
                'SpatialCorrelation', false,...
                'NumTransmitAntennas', nodeCfg.numTx,...
                'NumReceiveAntennas', nodeCfg.numRx,...
                'PathGainsOutputPort', true);
    otherwise
        disp('Unknown Channel Fading Type');
end


for indxRx = 1 : simCfg.numLink
    for indxTx  = 1 : simCfg.numLink
        dist = norm(linkInfo(indxTx).txLoc - linkInfo(indxRx).rxLoc);
        propDelayMat(indxRx,indxTx) = dist/3e8*rfCfg.bwHz;
        pathLossDbMat(indxRx, indxTx) = 2*rfCfg.anteGainDbi + rfCfg.pathloss1m - 10*rfCfg.plExponent*log10(dist);
        fDopplerHz = (abs(linkInfo(indxTx).txSpdMph-linkInfo(indxRx).rxSpdMph))*1.61/3.6/(3e8/rfCfg.fcHz);
        cfoHz = (linkInfo(indxRx).rxCfoHz-linkInfo(indxTx).txCfoHz);
        
%         reset(chan);
        release(chan);
        chan.MaximumDopplerShift = fDopplerHz;
        [~,ch] = step(chan,randi([0 1], numChanSamp, nodeCfg.numTx));
        ch = squeeze(ch);
        ch = permute(ch, [3 2 1]);
        
        cfoSeq = exp(1i*2*pi*cfoHz/nodeCfg.fsHz*(0:numChanSamp-1));
        cfoSeq = repmat(cfoSeq, [nodeCfg.numRx 1 nodeCfg.numTx]);
        cfoSeq = permute(cfoSeq, [1 3 2]);
        
        fadingChan(indxRx,indxTx,:,:,:) = ch.*cfoSeq;
%         fadingChan(indxRx,indxTx,:,:,:) = gen_space_time_chan(chanSampTimeSec,fDopplerHz,nodeCfg.numTx,nodeCfg.numRx,cfoHz);
    end
end

end

