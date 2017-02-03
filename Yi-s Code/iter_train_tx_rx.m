function [linkInfo, numALinkRcd, totTxPwrRec, rxRateRec, rxRateRecIT, aSetT, aSetR] =  iter_train_tx_rx(fadingChan,chanFlag,nodeCfg,frameCfg,algCfg,linkInfo,ulPilotSet,aSetR,aSetT,plDbMat,propDelayMat)
% EPOCH 3:Tx <-> Rx spatial fine tuning

numRxAnte = nodeCfg.numRx;
numTxAnte = nodeCfg.numTx;
numLink = size(fadingChan,1);
numIniSlot = frameCfg.numIniSlot;
ulPilotLen = frameCfg.ulPilotLen;
ulPilotCpLen = frameCfg.ulPilotCpLen;
soundingLen = frameCfg.soundingLen;
decodeSnrThrsDb = algCfg.snrThrsDb;

L = 2^(ceil(log2(ulPilotLen))+3);
soundSeq = [ulPilotSet(1:numLink,end-ulPilotCpLen+1:end) ulPilotSet(1:numLink,:)];

maxTxPowLin = db2pow(nodeCfg.maxTxPowDbm);
pathlossLin = db2pow(plDbMat);

rxRateRec = zeros(numLink,numIniSlot);
rxRateRecIT = zeros(numLink,numIniSlot);
totTxPwrRec = zeros(1,numIniSlot);
maxWaitTime = ceil(algCfg.maxSubactiveFrame*rand(1,numLink)+1);
% maxWaitTime = algCfg.maxSubactiveFrame*ones(1,numLink);
stChan = zeros(numRxAnte,soundingLen);
timeStampRx = zeros(1,numLink);
numALinkRcd = zeros(2,numIniSlot);
propLen = 50;

t = 0;
for iterIndx = 1:numIniSlot
    %------------------- Tx-> Rx, the receiver can estimate the rate and feedback to the Tx the power adjustment
    for indxRx = aSetR
        chanSig = sqrt(db2pow(linkInfo(indxRx).txPwrDbm)*pathlossLin(indxRx,indxRx))*squeeze(fadingChan(indxRx,indxRx,:,:,t+1))*linkInfo(indxRx).txVec;
        chanIntf = zeros(numRxAnte,numLink);
        Y = 1/sqrt(2).*(randn(numRxAnte,soundingLen+propLen)+1i*randn(numRxAnte,soundingLen+propLen));
        for indxTx = aSetT
            if chanFlag(indxRx,indxTx)==0
                continue
            end
            dly = round(propDelayMat(indxRx,indxTx));
            Hij = squeeze(fadingChan(indxRx,indxTx,:,:,t+1:t+soundingLen));
            txVec = sqrt(db2pow(linkInfo(indxTx).txPwrDbm)*pathlossLin(indxRx,indxTx))*linkInfo(indxTx).txVec;
            for tt = 1:soundingLen
                stChan(:,tt) = Hij(:,:,tt)*txVec;
            end
            Y(:,dly+1:dly+soundingLen) = Y(:,dly+1:dly+soundingLen)+stChan.*repmat(soundSeq(indxTx,:),numRxAnte,1);
            if indxRx ~= indxTx  %&& linkInfo(indxTx).txPow > 0
                chanIntf(:,indxTx) = stChan(:,1);
            end
        end
        if algCfg.compCFO == 1
            Y = Y.*repmat(exp(-1j*2*pi*linkInfo(indxRx).cfoEstHz/nodeCfg.fsHz*(1:soundingLen+propLen)),numRxAnte,1);
        end
        d = ulPilotCpLen;
        zc1 = ulPilotSet(indxRx,:);
        Y1 = Y(:,d+1:d+ulPilotLen).*(ones(numRxAnte,1)*conj(zc1));
        Rx = Y1*Y1';    [U,S,V] = svd(Rx);    T = U*diag(diag(S).^(-1/2))*V';        %Y1 = T*Y1;
        [q1,Yf] = est_q(T,Y1,L);
        [mxQ,indx] = max(q1);
        estChan = Yf(:,indx);
%         estChan = mean(Y1,2);
        
        ppsnrEst = mxQ/(1-mxQ);
        ppsnrEstDb = pow2db(ppsnrEst);
        Q = chanIntf*chanIntf'+eye(numTxAnte);
        if algCfg.estChan
            rxVec = Rx\estChan;
            ppsnr = abs(rxVec'*chanSig)^2/real(rxVec'*Q*rxVec);
        else
            if all(chanSig)==0
                mxQ = 0;
                rxVec = 0;
                ppsnrEstDb = -inf;
            else
                rxVec = Q\chanSig;
                ppsnr = abs(rxVec'*chanSig)^2/real(rxVec'*Q*rxVec);
                ppsnrEstDb = pow2db(ppsnr);
                mxQ = ppsnr/(1+ppsnr);
            end
        end
        if mxQ<algCfg.stTrainSlot
            linkInfo(indxRx).rate = 0;
            timeStampRx(indxRx) = inf;
        else
            rxRateRecIT(indxRx,iterIndx) = log2(1+ppsnr);
            if mxQ>algCfg.ctTrainSlot
                ppsnrDb = pow2db(ppsnr);
                linkInfo(indxRx).ppsnrDb = ppsnrDb;
                % linkInfo(indxRx).rate = log2(1+ppsnr);
                feasInd = find(ppsnrEstDb >= decodeSnrThrsDb);
                if isempty(feasInd)
                    linkInfo(indxRx).rate = 0;
                else
                    linkInfo(indxRx).rate = algCfg.rateSet(feasInd(end));
                end
                timeStampRx(indxRx) = 0;
            else
                linkInfo(indxRx).rate = 0;
                timeStampRx(indxRx) = timeStampRx(indxRx)+1;
            end
        end
        if timeStampRx(indxRx)<= maxWaitTime
            linkInfo(indxRx).rxVec = rxVec/norm(rxVec); % W_mmse; %W_mmse/norm(W_mmse);
        else
            linkInfo(indxRx).rxVec = zeros(numRxAnte,1);
        end
        linkInfo(indxRx).ppsnrEstDb = ppsnrEstDb;
    end
    aSetR = setdiff(aSetR,find(timeStampRx>=maxWaitTime));
    rxRateRec(:,iterIndx)  = [linkInfo(:).rate];
    totTxPwrRec(iterIndx) = pow2db(sum(db2pow([linkInfo(:).txPwrDbm])));
    
    %---------- Rx-> Tx, the Tx fine tune the Tx beamforming vector
    t = t+soundingLen+1; %+1 because there is a gap in between
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
            if algCfg.compCFO == 1
                Y(:,dly+1:dly+soundingLen) = Y(:,dly+1:dly+soundingLen)+stChan.*repmat(soundSeq(indxRx,:),numTxAnte,1).*repmat(exp(-1j*2*pi*linkInfo(indxRx).cfoEstHz/nodeCfg.fsHz*(1:soundingLen)),numRxAnte,1);
            else
                Y(:,dly+1:dly+soundingLen) = Y(:,dly+1:dly+soundingLen)+stChan.*repmat(soundSeq(indxRx,:),numTxAnte,1);
            end
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
            linkInfo(indxTx).txVec = conj(txVec)/norm(txVec); % there was no conj()
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
end
for indxLink = aSetR
    [linkInfo(indxLink).txPwrDeltaDb linkInfo(indxLink).rate] = calc_pwr_delta(linkInfo(indxLink).ppsnrEstDb,algCfg.snrThrsDb,linkInfo(indxLink).txPwrDbm,nodeCfg.maxTxPowDbm,algCfg.trainTxPwrAdjSet,algCfg.rateSet);
    linkInfo(indxLink).txPwrDbm = linkInfo(indxLink).txPwrDbm+linkInfo(indxLink).txPwrDeltaDb;
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