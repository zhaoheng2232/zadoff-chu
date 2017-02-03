function linkInfo = gen_link_info(simCfg,rfCfg,nodeCfg,txLocation,rxLocation)

txSpdMph = rand(1,simCfg.numLink)*simCfg.maxSpeedMph;
rxSpdMph = rand(1,simCfg.numLink)*simCfg.maxSpeedMph;
txCfoHz = (rand(1,simCfg.numLink)*2-1)*rfCfg.maxCfoHz;
rxCfoHz = (rand(1,simCfg.numLink)*2-1)*rfCfg.maxCfoHz;
linkInfo = struct('txLoc',num2cell(txLocation,1),'rxLoc',num2cell(rxLocation,1),'txSpdMph',num2cell(txSpdMph),'rxSpdMph',...
    num2cell(rxSpdMph),'txCfoHz',num2cell(txCfoHz),'rxCfoHz',num2cell(rxCfoHz),'txPwrDbm',num2cell(-inf(1,simCfg.numLink)),...
    'rate',num2cell(zeros(1,simCfg.numLink)),'ppsnrDb',num2cell(-inf(1,simCfg.numLink)),'txVec',num2cell(zeros(nodeCfg.numTx,simCfg.numLink),1),'rxVec',num2cell(zeros(nodeCfg.numRx,simCfg.numLink),1));
