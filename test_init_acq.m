clc
% close all
% clear all

%% Configuration
[simCfg,rfCfg,nodeCfg,algCfg,frameCfg] = init_config();

%% Topology and Link Infomation Realization 
linkInfo = gen_link_info(simCfg,rfCfg,nodeCfg);

%% Pathloss and Fading Channel Generation
[pathLossDbMat,fadingChan,propDelayMat,chanSampTimeSec] = gen_chan_info(simCfg,nodeCfg,frameCfg,rfCfg,linkInfo);
