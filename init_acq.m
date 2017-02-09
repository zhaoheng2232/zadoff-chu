function [detected,cfoEstHz,tauEst] = init_acq(nodeCfg,frameCfg,algCfg,rfCfg,pathLossDbMat,fadingChan,propDelayMat)

syncSeqLen = frameCfg.syncSeqLen;
syncSeqRateSet = frameCfg.syncSeqRateSet;
syncSeqCpLen = frameCfg.syncSeqCpLen;
epochLen = frameCfg.epochLen;
fsHz = nodeCfg.fsHz;
upSampRate = nodeCfg.upSampRate;

numTxAnte = nodeCfg.numTx;
numRxAnte = nodeCfg.numRx;
[numLink,~,~,~,~,numMultipath] = size(fadingChan);

pathlossLin = db2pow(pathLossDbMat);
maxTxPowLin = db2pow(nodeCfg.maxTxPowDbm);
noiseFloorLin = db2pow(rfCfg.noiseFloorDbm);

propLen = 50; delayOffset = 10;
cfoEstHz = zeros(1,numLink);
tauEst = zeros(1,numLink);
detected = ones(1,numLink);
for indxRx = 1:numLink
    %--- generate Rx samples
    Y = 1/sqrt(2).*(randn(numRxAnte,epochLen+propLen)+1i*randn(numRxAnte,epochLen+propLen));
    for indxTx = 1:numLink
        zc1 = gen_zc_seq(syncSeqRateSet(indxTx),syncSeqLen);
        zc2 = conj(zc1);
        seqTmp = [zc1(end-syncSeqCpLen+1:end) zc1 zc2(end-syncSeqCpLen+1:end) zc2];
        for indxMp = 1:numMultipath
            dly = round(propDelayMat(indxRx,indxTx,indxMp))+delayOffset;
            Hij = squeeze(fadingChan(indxRx,indxTx,:,:,:,indxMp));
            h = sqrt(maxTxPowLin/numTxAnte/noiseFloorLin*pathlossLin(indxRx,indxTx))*squeeze(Hij(:,1,:));   % Noise normalized to 0 dBm
            Y(:,dly+1:dly+epochLen) = Y(:,dly+1:dly+epochLen)+h.*repmat(seqTmp,numRxAnte,1);
        end
    end
    
    %--- estimate zeta
    start = syncSeqCpLen; %+8 because the pulse shaping introduce a 8-sample delay
    zc1 = gen_zc_seq(syncSeqRateSet(indxRx),syncSeqLen);
	zc2 = conj(zc1);
    Y1 = Y(:,start+1:start+syncSeqLen).*(ones(numRxAnte,1)*conj(zc1));
    Rx = Y1*Y1';    [U,S,V] = svd(Rx);    T = U*diag(diag(S).^(-1/2))*V';
    Y1 = T*Y1;
    q1 = est_q(Y1);    
    [mxQ1,loc1] = max(q1);
    if mxQ1<algCfg.dtRts
        detected(indxRx) = false;
        continue
    end    
    zeta1 = fine_est_peak(q1,syncSeqLen,loc1,1); 
    
    Y2 = Y(:,start+syncSeqLen+syncSeqCpLen+1:start+syncSeqLen*2+syncSeqCpLen).*(ones(numRxAnte,1)*conj(zc2));
    Rx = Y2*Y2';    [U,S,V] = svd(Rx);    T = U*diag(diag(S).^(-1/2))*V';
    Y2 = T*Y2;
    q2 = est_q(Y2);
    [mxQ2,loc2] = max(q2);
    if mxQ2<algCfg.dtRts
        detected(indxRx) = false;
        continue
    end
    zeta2 = fine_est_peak(q2,syncSeqLen,loc2,1); 
    fHat = mean([zeta1 zeta2]);
    cfoEstHz(indxRx) = fHat*fsHz/upSampRate;
    tauEst(indxRx) = (zeta2-zeta1)*syncSeqLen/syncSeqRateSet(indxRx)/2;
end
% fprintf('---------------\n');

end

function [q,Yf] = est_q(Y,Lf)
L = size(Y,2);
if nargin<2
    Lf = 2^(ceil(log2(L*2-1)))*1;
end
Yf = fft(Y,Lf,2);
% if nargin<3
p = sum(Yf.*conj(Yf),1);
% else
%     iR = inv(R);
%     p = zeros(1,Lf);
%     for n = 1:Lf
%         p(n) = real(Yf(:,n)'*iR*Yf(:,n));
%     end
% end

q = p/L;
end

function x = fine_est_peak(q,L,indx,flag)
% flag indicate if apply fine channel estimation
% output: x is the fine estimation of the frequency
Lf = length(q);
if nargin<3
    [~,indx] = max(q);
end
% indx = indx+indx3;
if indx <= Lf/2
    zeta = (indx-1)/Lf;
else
    zeta = (indx-1)/Lf-1;
end
if nargin==4 && flag==0
    x = zeta;
    return
end
a = ifft(q);
x0 = zeta;
while 1
    e = a(2:L).*exp(-1j*2*pi*x0*(1:L-1));
    f0 = a(1)+2*real(sum(e));
    df = 4*pi*imag(dot(1:L-1,e));
    v = sign(df);
    t = 1e-3; bt = 0.3; alf = 0.3;
    while 1
        x = x0+t*v;
        e = a(2:L).*exp(-1j*2*pi*x*(1:L-1));
        f = a(1)+2*real(sum(e));
        if f > f0+alf*t*df*v
            x0 = x;
            break
        end
        t = bt*t;
    end
    if abs(t*v)<1e-5
        break
    end
end
% plot(-.5:1/Lf:.5-1/Lf,fftshift(q)), hold on, plot(x,f,'o','MarkerFaceColor','b')
% -- newton iteratio wouldn't work
% for ii = 1:7
%     e = a(2:L).*exp(-1j*2*pi*x*(1:L-1));
%     f = a(1)+2*real(sum(e));
%     d1 = 4*pi*imag(dot(1:L-1,e));
%     d2 = -8*pi^2*real(dot((1:L-1).^2,e));
%     dx = d1/d2;
%     x = x-0.5*dx;
%     if abs(dx)<1e-6
%         break
%     end
% end
end


