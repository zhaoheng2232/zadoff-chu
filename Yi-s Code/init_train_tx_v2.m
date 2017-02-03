function [detect,txVec] = init_train_tx_v2(fadingChan,chanFlag,nodeCfg,frameCfg,algCfg,syncSeqSet,aSetR,plDbMat,propDelayMat)
% in Epoch 2, the Rx nodes send Z-C sequence, the Tx figure out the Tx eigen-beamforming vector
% the difference between v2 from the previous version (init_sync) is that here
% we do use beta distribution based detection algorithm

syncSeqLen = frameCfg.syncSeqLen;
epochLen = frameCfg.epochLen;
syncSeqCpLen = frameCfg.syncSeqCpLen;
numRxAnte = nodeCfg.numRx;
numTxAnte = nodeCfg.numTx;
numLink = size(fadingChan,1);
L = 2^(ceil(log2(syncSeqLen*2-1)))*2;
winT = round(-L/40):round(L/40);

detect = ones(1,numLink);
maxTxPowLin = db2pow(nodeCfg.maxTxPowDbm);
pathlossLin = db2pow(plDbMat);
propLen = 50;
% generate the sequence for epoch 2. it is a numRxAnte x epochLen matrix for each node
syncSeqEpoch2 = zeros(numRxAnte,epochLen,numLink);
for indxTx = 1:numLink
    for indxAnt = 1:numRxAnte
        seqTmp = circshift(syncSeqSet(indxTx,:),[0 round(syncSeqLen/numRxAnte*(indxAnt-1)/indxTx)]);
        syncSeqEpoch2(indxAnt,:,indxTx) = [seqTmp(1,end-syncSeqCpLen+1:end) seqTmp];
    end
end

txVec = zeros(numTxAnte,numLink);
for indxTx = 1:numLink
    Y = 0/sqrt(2).*(randn(numTxAnte,epochLen+propLen)+1i*randn(numTxAnte,epochLen+propLen));
    for indxRx = aSetR
        if chanFlag(indxRx,indxTx)==0
            continue
        end
        for indxAnte = 1:numRxAnte
            dly = round(propDelayMat(indxRx,indxTx));
            Hij = squeeze(fadingChan(indxRx,indxTx,indxAnte,:,:));
            h = sqrt(maxTxPowLin/numRxAnte*pathlossLin(indxRx,indxTx))*Hij;        
            Y(:,dly+1:dly+epochLen) = Y(:,dly+1:dly+epochLen)+h.*repmat(syncSeqEpoch2(indxAnte,:,indxRx),numRxAnte,1);
        end
    end
    d = syncSeqCpLen;
    zc1 = syncSeqSet(indxTx,:);
    Y1 = Y(:,d+1:d+syncSeqLen).*(ones(numRxAnte,1)*conj(zc1));
    fY = fft(Y1,L,2);
    Rx = Y1*Y1';    [U,S,V] = svd(Rx);    T = U*diag(diag(S).^(-1/2))*V';
    Y1 = T*Y1;
    q1 = est_q(Y1,L);
    indx = zeros(1,numRxAnte);
    [mx,indx(1)] = max(q1);
    if mx<algCfg.dtCts
        detect(indxTx) = false;
        continue
    end
    if algCfg.eigIniTx == 0
        tmp = randn(numTxAnte,1)+1j*randn(numTxAnte,1);
        txVec(:,indxTx) = tmp/norm(tmp);
        continue
    end
    for nn = 2:numRxAnte
        win2 = mod(indx(1)+L/numRxAnte*(nn-1)+winT,L)+1;
        [~,indxTmp] = max(q1(win2));
        indx(nn) = win2(indxTmp);
    end
    estChan = fY(:,indx).';
    [~,~,V]  = svd(estChan);
    txVec(:,indxTx) = V(:,1);
%     truChan = squeeze(fadingChan(indxTx,indxTx,:,:,1)); for n = 1:4, abs(dot(estChan(:,n),truChan(:,n)))/norm(estChan(:,n))/norm(truChan(:,n)),end
end
end

function [q,Yf] = est_q(Y,Lf)
L = size(Y,2);
if nargin<2
    Lf = 2^(ceil(log2(L*2-1)))*1;
end
Yf = fft(Y,Lf,2);
p = sum(Yf.*conj(Yf),1);
q = p/L;
end