function [txLocation, rxLocation] = gen_topology(simCfg)
% link_info = gen_topology( sim_ctrl, link_info )
% generates topology of network user pairs based on sim_ctrl and stores the
% information in link_info

txLocation = zeros(2, simCfg.numLink);
rxLocation = zeros(2, simCfg.numLink);

switch simCfg.topologyType
    case 'packedCircle'
        for indxLink = 1:simCfg.numLink
            %rx distributed randomly within coverage circle
            random_phase = rand*2*pi;
            dist = 1-abs(rand+rand-1);
            rxLocation(:,indxLink) = dist*[ simCfg.areaRadiusM*cos(random_phase) simCfg.areaRadiusM*sin(random_phase) ].';
            
            %tx distributed randomly within max_tx_range around rx
            random_phase = rand*2*pi;
           % dist = 1-abs(rand+rand-1); 
            txRxDist = rand * (simCfg.maxTxRxDistM - simCfg.minTxRxDistM) + simCfg.minTxRxDistM;
            relativeTxPos = txRxDist *[ cos(random_phase) sin(random_phase) ];
            txLocation(:,indxLink) = rxLocation(:,indxLink) + relativeTxPos.';

            %ensures that tx is also located within the same coverage circle as rx
            while (norm(txLocation(:, indxLink)) > simCfg.areaRadiusM )
                random_phase = rand*2*pi;
             %   r = 1-abs(rand+rand-1);
                txRxDist = rand * (simCfg.maxTxRxDistM - simCfg.minTxRxDistM) + simCfg.minTxRxDistM;
                relativeTxPos = txRxDist * [cos(random_phase) sin(random_phase) ];
                txLocation(:,indxLink) = rxLocation(:,indxLink) + relativeTxPos.';
            end
        end
        
    case 'onCircle'
        txLocation = [cos(0:2*pi/simCfg.numLink:2*pi-2*pi/simCfg.numLink);sin(0:2*pi/simCfg.numLink:2*pi-2*pi/simCfg.numLink)]*simCfg.areaRadiusM;
    otherwise
        disp('Unknown topology');
end