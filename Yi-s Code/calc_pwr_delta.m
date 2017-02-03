function [txPwrDeltaDb,rate] = calc_pwr_delta(ppsnrDb,decodeSnrThrsDb,txPwrDbm,maxTxPowDbm,deltaDbGrid,rateVec)
% upMax = max(deltaDbGrid);
% dnMax = min(deltaDbGrid);
if length(rateVec)~=length(decodeSnrThrsDb)
    error('length(rateVec)~=length(decodeSnrThrsDb)!!');
end
if ppsnrDb > max(decodeSnrThrsDb)
    txPwrDeltaDb = deltaDbGrid(find(deltaDbGrid>max(decodeSnrThrsDb)-ppsnrDb,1));
    rate = rateVec(end);
    %txPwrDeltaDb = max(ceil(max(decodeSnrThrsDb)-ppsnrDb),min(deltaDbGrid));
elseif ppsnrDb < min(decodeSnrThrsDb)
    txPwrDeltaDb = deltaDbGrid(find(deltaDbGrid>min(decodeSnrThrsDb)-ppsnrDb,1));
    if isempty(txPwrDeltaDb)
        txPwrDeltaDb = deltaDbGrid(end);
    end
    rate = 0;
    %txPwrDeltaDb = min(ceil(min(decodeSnrThrsDb)-ppsnrDb),max(deltaDbGrid));  
else
    indx = find(ppsnrDb < decodeSnrThrsDb,1);
    rate = rateVec(indx-1);
    up = decodeSnrThrsDb(indx)-ppsnrDb;
%     dn = decodeSnrThrsDb(indx-1)-ppsnrDb;
    if up < maxTxPowDbm-txPwrDbm
        txPwrDeltaDb = deltaDbGrid(find(deltaDbGrid>decodeSnrThrsDb(indx)-ppsnrDb,1));
        if isempty(txPwrDeltaDb)
            txPwrDeltaDb = deltaDbGrid(end);
        end
%         txPwrDeltaDb =  min(ceil(up),max(deltaDbGrid));
    else
        txPwrDeltaDb = deltaDbGrid(find(deltaDbGrid>decodeSnrThrsDb(indx-1)-ppsnrDb,1));        
%         txPwrDeltaDb =  max(ceil(dn),min(deltaDbGrid));
    end
end
