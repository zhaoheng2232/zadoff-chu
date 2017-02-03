function [ simCfg,rfCfg,nodeCfg,algCfg,frameCfg ] = init_config()

% ---------- Configure simulation parameters
simCfg.conopsId = 1;
simCfg.topologyType = 'packedCircle';
simCfg.numLink = 3;
simCfg.areaRadiusM = 2900;
simCfg.maxTxRxDistM = 500;
simCfg.minTxRxDistM = 500;
simCfg.maxSpeedMph = 0;
simCfg.fadingType = 'Rayleigh';
simCfg.numMultipath = 1;
simCfg.runDetectAlg = 1;
simCfg.maxTxRange = 1000;%2900;%10^((nodeCfg.maxTxPowDbm+2*2 + rfCfg.pathloss1m - rfCfg.noiseFloorDbm)/(10*rfCfg.plExponent));

% ---------- Configure RF parameters
rfCfg.fcHz = 2.4e9;
rfCfg.bwHz = 1.25e6;
rfCfg.maxCfoHz = rfCfg.fcHz*5e-6; %5ppm
rfCfg.plExponent = 3.8;
rfCfg.pathloss1m = 20*log10( 3e8/rfCfg.fcHz/4/pi/1 ); %1m loss
rfCfg.nfDb = 3;
rfCfg.anteGainDbi = 2;
rfCfg.noiseFloorDbm = -113 + rfCfg.nfDb+10*log10(rfCfg.bwHz); %dBm

% ---------- Configure node parameters
nodeCfg.numTx = 4;
nodeCfg.numRx = 4;
nodeCfg.maxTxPowDbm = 40;
nodeCfg.minTxPowDbm = -20;
nodeCfg.ctrlPwrDbGrid1 = -9:3:0;
nodeCfg.ctrlPwrDbGrid2 = -2:1;
nodeCfg.fsHz = rfCfg.bwHz;
nodeCfg.upSampRate = 1;
if nodeCfg.upSampRate==1
    nodeCfg.srrcf = [0 0 1 0 0];
else
    nodeCfg.srrcf = srrc(nodeCfg.upSampRate,0.5);
end

% ---------- Configure algorithm parameters
algCfg.maxSubactiveFrame = 6; inf;10;
algCfg.dtCts = 0.2;
algCfg.dtRts = 0.1;
algCfg.eigIniTx = 1;
algCfg.stTrainSlot = 0.1; %shutdown threshold in the training slot
% algCfg.snrThrsDb = [7.4 11.6 18.8];  % QPSK 1/2; 16QAM, 1/2; 64QAM, 2/3
% algCfg.rateSet = [1 2 4];
algCfg.snrThrsDb = [4.1 7.4 9.5 11.6 15.2 18.8 20.7 22.9];  % QPSK 1/2; 16QAM, 1/2; 64QAM, 2/3
algCfg.rateSet = [0.5 1 1.5 2 3 4 4.5 5];
algCfg.dataTxPwrAdjSet = -1:2;
algCfg.trainTxPwrAdjSet = -9:3:0;
algCfg.estChan = 0;
algCfg.minSnrDb = min(algCfg.snrThrsDb);
minSnr = db2pow(min(algCfg.snrThrsDb));
algCfg.ctTrainSlot = minSnr/(1+minSnr); %countdown threshold in the training slot

% ---------- Configure Frame parameters
frameCfg.numFrame = 10;
frameCfg.numIniSlot = 10;
frameCfg.epochLen = 550;%250;%;1050;%500;
frameCfg.syncSeqLen = 509;%229;%1021;%479;
frameCfg.syncSeqCpLen = frameCfg.epochLen-frameCfg.syncSeqLen;
frameCfg.epch1LenUs = frameCfg.epochLen/rfCfg.bwHz*1e6;
frameCfg.epch2LenUs = frameCfg.epch1LenUs;
frameCfg.timeGapUs = 100;
frameCfg.framLenUs = 5000;
frameCfg.framLen= frameCfg.framLenUs*rfCfg.bwHz/1e6;
frameCfg.timeGap = frameCfg.timeGapUs/nodeCfg.fsHz*1e6;
frameCfg.soundingLen = 150;%80;% 100; %220;
frameCfg.ulPilotLen = 127;%61;% 79;%199;
frameCfg.ulPilotCpLen = frameCfg.soundingLen-frameCfg.ulPilotLen;

end

