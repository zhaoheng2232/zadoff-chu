function plDb =  calc_path_loss(anteGain,pathloss1m,plExponent,noiseFloorDbm,dist)
% pl =  calc_path_loss(anteGain,pathloss1m,plExponent,noiseFloorDbm,dist)
% calculate path loss normalized by noise floor, in other words, if the Tx
% power is 0dBm, then the Rx SNR will be plDb

plDb = 2*anteGain + pathloss1m - 10*plExponent*log10(dist) - noiseFloorDbm;
