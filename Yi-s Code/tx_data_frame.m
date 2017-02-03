function [rxRateRec,totTxPwrRec,numALinkRcd] = tx_data_frame(fadingChan,chanFlag,nodeCfg,frameCfg,algCfg,linkInfo,ulPilotSet,aSetT,aSetR,plDbMat,propDelayMat)

numDataFrame= frameCfg.numFrame;
numRxAnte = nodeCfg.numRx;
numTxAnte = nodeCfg.numTx;
numLink = size(fadingChan,1);
ulPilotLen = frameCfg.ulPilotLen;
ulPilotCpLen = frameCfg.ulPilotCpLen;

soundingLen = frameCfg.soundingLen;
decodeSnrThrsDb = algCfg.snrThrsDb;

L = 2^(ceil(log2(ulPilotLen))+3);
soundSeq(aSetR,:) = [ulPilotSet(aSetR,end-ulPilotCpLen+1:end) ulPilotSet(aSetR,:)];

maxTxPowLin = db2pow(nodeCfg.maxTxPowDbm);
pathlossLin = db2pow(plDbMat);

stChan = zeros(numRxAnte,soundingLen);
rxRateRec = zeros(numLink,numDataFrame);
totTxPwrRec = zeros(1,numDataFrame);
numALinkRcd = zeros(2,numDataFrame);
propLen = 10;

t = 0;
for iterIndx = 1:numDataFrame
    for indxRx = aSetR    
        chanSig = sqrt(db2pow(linkInfo(indxRx).txPwrDbm)*pathlossLin(indxRx,indxRx))*squeeze(fadingChan(indxRx,indxRx,:,:,t+1))*linkInfo(indxRx).txVec;
        chanIntf = zeros(numRxAnte,numLink);
        for indxTx = aSetT
            if chanFlag(indxRx,indxTx)==0
                continue
            end
            Hij = squeeze(fadingChan(indxRx,indxTx,:,:,t+1));
            txVec = sqrt(db2pow(linkInfo(indxTx).txPwrDbm)*pathlossLin(indxRx,indxTx))*linkInfo(indxTx).txVec;
            if indxRx ~= indxTx  %&& linkInfo(indxTx).txPow > 0
                chanIntf(:,indxTx) = Hij*txVec;
            end
        end
        Q = chanIntf*chanIntf'+eye(numTxAnte);
        if all(chanSig)==0
            linkInfo(indxRx).rxVec = zeros(numRxAnte,1);
            ppsnrDb = -inf;
        else
            rxVec = Q\chanSig;
            linkInfo(indxRx).rxVec = rxVec/norm(rxVec);
            ppsnr = abs(rxVec'*chanSig)^2/real(rxVec'*Q*rxVec);
            ppsnrDb = pow2db(ppsnr);
        end
        linkInfo(indxRx).ppsnrDb = ppsnrDb;
        feasInd = find(ppsnrDb >= decodeSnrThrsDb);
        if isempty(feasInd)
            linkInfo(indxRx).rate = 0;
        else
            linkInfo(indxRx).rate = algCfg.rateSet(feasInd(end));
        end
        [linkInfo(indxRx).txPwrDeltaDb linkInfo(indxRx).rate] = calc_pwr_delta(ppsnrDb,algCfg.snrThrsDb,linkInfo(indxRx).txPwrDbm,nodeCfg.maxTxPowDbm,algCfg.dataTxPwrAdjSet,algCfg.rateSet);
        linkInfo(indxRx).txPwrDbm = linkInfo(indxRx).txPwrDbm+linkInfo(indxRx).txPwrDeltaDb;
    end
    rxRateRec(:,iterIndx)  = [linkInfo(:).rate];
    totTxPwrRec(iterIndx) = pow2db(sum(db2pow([linkInfo(:).txPwrDbm])));
    %     % Tx-> Rx, the receiver can estimate the rate and feedback to the Tx the power adjustment
    %     for indxRx = aSet
    %         % re_measure the channel matrices
    %         nn = 1;
    %         chanAll = zeros(numRxAnte,numActiveLink);
    %         for indxLink = 1:numActiveLink
    %             indxTx = aSet(indxLink);
    %             if chanFlag(indxRx,indxTx)==0
    %                 continue
    %             end
    %             if indxRx == indxTx  %&& linkInfo(indxTx).txPow > 0
    %                 Hkk    = squeeze(fadingChan(indxRx,indxTx,:,:,t));
    %                 chanSig = sqrt(maxTxPowLin*pathlossLin(indxRx,indxTx))*Hkk*linkInfo(indxTx).txVec;
    %                 chanAll(:,indxLink) = chanSig;
    %                 sigLinkIndx = indxLink;
    %             elseif ~isinf(linkInfo(indxTx).txPwrDbm)
    %                 Hij   = squeeze(fadingChan(indxRx,indxTx,:,:,t));
    %                 chanTmp  = sqrt(maxTxPowLin*pathlossLin(indxRx,indxTx))*Hij*linkInfo(indxTx).txVec;
    %                 chanAll(:,indxLink) = chanTmp;
    %                 chanIntf(:,nn) = chanTmp;   nn = nn+1;
    %             end
    %         end
    %         pilotActive = ulPilotSet(aSet,:); % pilot PN sequences
    %         Y = chanAll*pilotActive + sqrt(1/2)*(randn(numTxAnte,ulPilotLen)+1i*randn(numTxAnte,ulPilotLen));
    %         % R = Y*Y'/size(Y,2);
    %         % joint channel estimation
    %         estChanTmp     = Y*pilotActive' / (pilotActive*pilotActive');
    %         estChan        = estChanTmp(:,sigLinkIndx);
    %         % W_mmse = (R\estChan)';
    %         estIntfChan = estChanTmp(:,setdiff(1:numActiveLink,sigLinkIndx));
    %         Qest = estIntfChan*estIntfChan'+eye(numTxAnte);
    %         W_mmse = (Qest\estChan)';
    %         %sig_pilot    = ulPilotSet(indxTx,:);
    % %         Q = estIntfChan*estIntfChan'+eye(numTxAnte);
    %         Q = chanIntf*chanIntf'+eye(numTxAnte);
    %         ppsnr = abs(W_mmse*chanSig)^2/real(W_mmse*Q*W_mmse');
    %
    % %         total_cov = chanSig*chanSig'+chanIntf*chanIntf'+eye(numRxAnte);
    % %         W_mmse_orig   = chanSig' / total_cov;
    % %         W_mmse        = W_mmse_orig ./ (chanSig'*W_mmse_orig');
    % %         MSE      = 1-W_mmse_orig*chanSig;
    % %         ppsnr0    = 1/real(MSE)-1;
    %         ppsnrDb = pow2db(ppsnr);
    %         feasInd = find(ppsnrDb >= decodeSnrThrsDb);
    % %         if ppsnrDb < algCfg.snrThrsDb(linkInfo(indxRx).rate==algCfg.rateSet)-1
    % %             linkInfo(indxRx).per = 1;
    % %         else
    % %             linkInfo(indxRx).per = 0;
    % %         end
    %         if isempty(feasInd)
    %             linkInfo(indxRx).rate   = 0;
    %         else
    %             linkInfo(indxRx).rate = algCfg.rateSet(feasInd(end));
    %         end
    % %         linkInfo(indxRx).rate = log2(1+ppsnr);
    % %         linkInfo(indxRx).txPwrDeltaDb = calc_pwr_delta(ppsnrDb,decodeSnrThrsDb,linkInfo(indxRx).txPwrDbm,nodeCfg.maxTxPowDbm,nodeCfg.ctrlPwrDbGrid2);
    % %     if sum([linkInfo(indxRx).txPwrDeltaDb])
    % %         [linkInfo(:).txPwrDeltaDb]
    % %     end
    %         linkInfo(indxRx).rxVec = W_mmse/norm(W_mmse);
    %     end
    %     %disp([linkInfo(:).txPwrDeltaDb])
    %     % receiver vector finish
    %     rxRateRec(:,numIniSlot+iterIndx)  = [linkInfo(:).rate];
    %     totTxPwrRec(1,numIniSlot+iterIndx) = pow2db(sum(db2pow([linkInfo(:).txPwrDbm])));
    
    %% Rx-> Tx, the Tx fine tune the Tx beamforming vector
    t = t+1;
    aSetT1 = aSetT;
    for indxTx = aSetT
        chanIntf = zeros(numTxAnte,numLink);
        rxVec = sqrt(maxTxPowLin*pathlossLin(indxTx,indxTx))*linkInfo(indxTx).rxVec;
        Hii = squeeze(fadingChan(indxTx,indxTx,:,:,t+1));
        chanSig = Hii.'*conj(rxVec);
        Y = 1/sqrt(2).*(randn(numTxAnte,soundingLen+propLen)+1i*randn(numTxAnte,soundingLen+propLen));
        for indxRx = aSetR
            if chanFlag(indxRx,indxTx)==0
                continue
            end            
            dly = round(propDelayMat(indxRx,indxTx));
            Hij = squeeze(fadingChan(indxRx,indxTx,:,:,t+1:t+soundingLen));
            rxVec = sqrt(maxTxPowLin*pathlossLin(indxRx,indxTx))*linkInfo(indxRx).rxVec;
            for tt = 1:soundingLen
                stChan(:,tt) = Hij(:,:,tt).'*conj(rxVec);
            end
            Y(:,dly+1:dly+soundingLen) = Y(:,dly+1:dly+soundingLen)+stChan.*repmat(soundSeq(indxRx,:),numTxAnte,1);
            if indxRx ~= indxTx  %&& linkInfo(indxTx).txPow > 0
                chanIntf(:,indxRx) = stChan(:,1);
            end
        end
        % the Tx doesn't need to know the CFO
        d = ulPilotCpLen;
        zc1 = ulPilotSet(indxTx,:);
        Y1 = Y(:,d+1:d+ulPilotLen).*(ones(numRxAnte,1)*conj(zc1));
        Rx = Y1*Y1';    [U,S,V] = svd(Rx);    T = U*diag(diag(S).^(-1/2))*V';        %Y1 = T*Y1;
        [q1,Yf] = est_q(T,Y1,L);
        [mxQ,indx] = max(q1);
        if mxQ>algCfg.ctTrainSlot
            if algCfg.estChan
                estChan = Yf(:,indx);
                txVec = Rx\estChan;
            else
                Q = chanIntf*chanIntf'+eye(numTxAnte);
                txVec = Q\chanSig;
            end
            if all(chanSig)==0
                error('false detection!!')
            end
            linkInfo(indxTx).txVec = conj(txVec)/norm(txVec);  % there was no conj()
        else
            linkInfo(indxTx).txVec = zeros(numTxAnte,1);
            aSetT1 = setdiff(aSetT1,indxTx);
            linkInfo(indxTx).txPwrDbm = -inf;
        end
    end
    t = t+soundingLen+1;
    aSetT = aSetT1;
    numALinkRcd(1,iterIndx) = length(aSetR);
    numALinkRcd(2,iterIndx) = length(aSetT);
    
%     for indxTx = aSetT
%         % estimate channel
%         chanAll     = zeros(numTxAnte,numLink);
%         for indxRx = aSetR;
%             if indxRx == indxTx
%                 Hkk    = squeeze(fadingChan(indxRx,indxTx,:,:,t));
%                 chanSig = sqrt(maxTxPowLin/numTxAnte.*pathlossLin(indxRx,indxTx)).*linkInfo(indxRx).rxVec*Hkk; % true channel
%                 chanAll(:,indxRx) = chanSig.';
%                 sigLinkIndx       = indxLink;
%                 % use transpose since the sequence is sent from Rx to Tx
%             else
%                 Hij   = squeeze(fadingChan(indxRx,indxTx,:,:,t));
%                 chanIntf1 = sqrt(maxTxPowLin/numTxAnte.*pathlossLin(indxRx,indxTx))*linkInfo(indxRx).rxVec*Hij;
%                 chanAll(:,indxLink) = chanIntf1.';
%             end
%         end
%         pilotActive = ulPilotSet(aSet,:); % pilot PN sequences
%         Y = chanAll*pilotActive + sqrt(1/2)*(randn(numTxAnte,ulPilotLen)+1i*randn(numTxAnte,ulPilotLen));
%         % joint channel estimation
%         estChanTmp     = Y*pilotActive' / (pilotActive*pilotActive');
%         estChan        = estChanTmp(:,sigLinkIndx);
%         estIntfChan = estChanTmp(:,setdiff(1:numActiveLink,sigLinkIndx));
%         %sig_pilot    = ulPilotSet(indxTx,:);
%         estQ        = (estIntfChan*estIntfChan').';
%         txVec = (estQ+1e-6*eye(numTxAnte))\estChan;
%         txVec = txVec/norm(txVec);
%         %[txVec shutdown] = optimize_tx_bf_pow(estQ,estChan,simCfg,indxTx);
%         %         if iterIndx == 1
%         %             linkInfo(indxTx).timeStamp  = 0;
%         %         end
%         %         if linkInfo(indxTx).txPwrDbm>=nodeCfg.maxTxPowDbm && linkInfo(indxTx).txPwrDeltaDb>0
%         %             linkInfo(indxTx).timeStamp  = linkInfo(indxTx).timeStamp + 1;
%         %         else
%         %             linkInfo(indxTx).timeStamp  = 0;
%         %         end
%         %         if linkInfo(indxTx).timeStamp > maxWaitTime
%         %             linkInfo(indxTx).txPwrDbm = -inf;
%         %         end
%         linkInfo(indxTx).txPwrDbm =  max(min(nodeCfg.maxTxPowDbm,linkInfo(indxTx).txPwrDbm+linkInfo(indxTx).txPwrDeltaDb),nodeCfg.minTxPowDbm);
%         linkInfo(indxTx).txVec = db2mag(linkInfo(indxTx).txPwrDbm-nodeCfg.maxTxPowDbm)*txVec;
%     end
%     for indxTx = 1:numLink
%         if ~ismember(indxTx,aSet)
%             continue
%         end        
%         if iterIndx == 1
%             linkInfo(indxTx).timeStamp  = 0;
%         end
%         if linkInfo(indxTx).txPwrDbm>=nodeCfg.maxTxPowDbm && linkInfo(indxTx).txPwrDeltaDb>0 %linkInfo(indxRx).rate==0
%             linkInfo(indxTx).timeStamp  = linkInfo(indxTx).timeStamp + 1;
%         else
%             linkInfo(indxTx).timeStamp  = 0;
%         end
%         if linkInfo(indxTx).timeStamp > maxWaitTime(indxTx)
%             linkInfo(indxTx).txPwrDbm = -inf;
%             aSet = setdiff(aSet,indxTx);
%         end
% %         linkInfo(indxTx).txPwrDbm =  max(min(nodeCfg.maxTxPowDbm,linkInfo(indxTx).txPwrDbm+linkInfo(indxTx).txPwrDeltaDb),nodeCfg.minTxPowDbm);
%     end
%     numActiveLink = numel(aSet);
end
end

function [q,Yf] = est_q(T,Y,Lf)
L = size(Y,2);
if nargin<3
    Lf = 2^(ceil(log2(L*2-1)))*1;
end
Yf = fft(Y,Lf,2);
Yf1 = T*Yf;
p = sum(Yf1.*conj(Yf1),1);
q = p/L;
end