clc
% close all
clear all

%% Configuration
[simCfg,rfCfg,nodeCfg,algCfg,frameCfg] = init_config();

numMonteCarlo = 1000;
cfoErrDump = zeros(simCfg.numLink,numMonteCarlo);
for indxMC = 1 : numMonteCarlo
    fprintf('**Monte Carlo Simulation Index: %d**\n', indxMC);
    %% Topology and Link Infomation Realization 
    linkInfo = gen_link_info(simCfg,rfCfg,nodeCfg);

    %% Pathloss and Fading Channel Generation
    [pathLossDbMat,fadingChan,propDelayMat] = gen_chan_info(simCfg,nodeCfg,frameCfg,rfCfg,linkInfo);

    %% Optimization Itearations
    % ------------------------------RTS----------------------------------------
    t = 1;
    [detectR,cfoEst] = init_acq(nodeCfg,frameCfg,algCfg,rfCfg,pathLossDbMat,fadingChan(:,:,:,:,t:t+frameCfg.epochLen-1),propDelayMat);
    aSetR = find(detectR);

    for nn = 1:simCfg.numLink
        linkInfo(nn).cfoEstHz = cfoEst(nn);
        linkInfo(nn).cfoTruHz = linkInfo(nn).rxCfoHz-linkInfo(nn).txCfoHz;
        cfoErrDump(nn,indxMC) = linkInfo(nn).cfoEstHz - linkInfo(nn).cfoTruHz;
    end
    % ------------------------------CTS----------------------------------------
    % t = t + frameCfg.epochLen + frameCfg.timeGap;
end