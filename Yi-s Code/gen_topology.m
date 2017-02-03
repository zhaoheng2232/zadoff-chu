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

    case 2 %convoy operation
        for indxLink = 1:simCfg.numLink
            %rx distributed linearly with spacing convoy_gap
            rxLocation(:,indxLink) = [((mod(indxLink-1,8)+1)-4.5)*sim_ctrl.convoy_gap 0].';
            
            %tx distributed randomly within max_tx_range around rx
            random_phase = rand*2*pi;
            relativeTxPos = rand*[ sim_ctrl.max_tx_range*cos(random_phase) sim_ctrl.max_tx_range*sin(random_phase) ];
            txLocation(:,indxLink) = rxLocation(:,indxLink) + relativeTxPos.';
        end

    case 3 %UAV
        for indxLink = 1:simCfg.numLink
            %2 rx located uav_gap distance from each other
            rxLocation(:,indxLink) = [(mod(indxLink-1,2)-0.5)*sim_ctrl.uav_gap 0].';
            
            %tx distributed randomly within max_tx_range around rx
            random_phase = rand*2*pi;
            relativeTxPos = rand*[ sim_ctrl.max_tx_range*cos(random_phase) sim_ctrl.max_tx_range*sin(random_phase) ];
            txLocation(:,indxLink) = rxLocation(:,indxLink) + relativeTxPos.';
        end
    case 4 %fixed distance btw Tx and Rx
        for indxLink = 1 : simCfg.numLink
            random_phase = rand*2*pi;
            dist = 1-abs(rand+rand-1);
            rxLocation(:,indxLink) = dist*[ simCfg.areaRadiusM*cos(random_phase) simCfg.areaRadiusM*sin(random_phase) ].';           
            
            % tx_position
            random_phase = rand*2*pi;
            relativeTxPosition = sim_ctrl.tx_rx_radius .* [cos(random_phase) sin(random_phase) ];
            txLocation(:, indxLink) = rxLocation(:,indxLink) + relativeTxPosition.';
            while (norm(txLocation(:,indxLink)) > simCfg.areaRadiusM)
                random_phase = rand*2*pi;
                relativeTxPosition = sim_ctrl.tx_rx_radius .* [cos(random_phase) sin(random_phase) ];
                txLocation(:, indxLink) = rxLocation(:,indxLink) + relativeTxPosition.';
            end
        end
    otherwise
end

%for debug purposes
if 0*simCfg.debugPlot
    plot(simCfg.areaRadiusM*sin(0:pi/180:2*pi),simCfg.areaRadiusM*cos(0:pi/180:2*pi),'b'), hold on
    axis equal
    for indxLink = 1:simCfg.numLink
        plot(txLocation(1, indxLink),txLocation(2, indxLink),'rp','MarkerFaceColor','r','MarkerSize',10)
        plot(rxLocation(1, indxLink),rxLocation(2, indxLink),'go','MarkerFaceColor','g','MarkerSize',10)
        plot_arrow(txLocation(:, indxLink),rxLocation(:, indxLink),100)
    end
%     if (sim_ctrl.scenario_ID == 1 ||sim_ctrl.scenario_ID == 4 )
%         ring_num = 1000;
%         ring_tx = zeros(2,ring_num);
%         TX_radius = sim_ctrl.max_tx_range;
%         for ring_idx = 1:1000
%             phase_res = ring_idx/1000*2*pi;
%             ring_tx(:,ring_idx) = [TX_radius*cos(phase_res) TX_radius*sin(phase_res)];
%         end
% 
%         hold on;
%         plot(txLocation(1,:)/1e3, txLocation(2,:)/1e3, 'x', 'LineWidth',2.0, 'MarkerSize',8.0);
%         hold on;
%         axis([-1 1 -1 1]*10);
%         plot(rxLocation(1,:)/1e3, rxLocation(2,:)/1e3, 'ro', 'LineWidth',2.0, 'MarkerSize',8.0);
%         legend('TX User', 'RX User');
%         for indxLink = 1:size(rxLocation,2)
%             plot_arrow( txLocation(1,indxLink)/1e3,txLocation(2,indxLink)/1e3,rxLocation(1,indxLink)/1e3,rxLocation(2,indxLink)/1e3,'headwidth',0.00001,'headheight',0.00001);
%         end
%         for indxLink = 1:size(rxLocation,2)
%             plot(ring_tx(1,:)/1e3 + rxLocation(1,indxLink)/1e3, ring_tx(2,:)/1e3 + rxLocation(2,indxLink)/1e3, 'LineWidth',1.0, 'MarkerSize',6.0);
%         end
% 
%         cov_radius = simCfg.areaRadiusM;
%         for ring_idx = 1:1000
%             phase_res = ring_idx/1000*2*pi;
%             ring_tx(:,ring_idx) = [cov_radius*cos(phase_res) cov_radius*sin(phase_res)];
%         end
%         plot(ring_tx(1,:)/1e3, ring_tx(2,:)/1e3, 'LineWidth',1.0, 'MarkerSize',6.0,'Color',[0 0 0]);
% 
%         xlabel('Unit: Km');
%         drawnow;
%     elseif (sim_ctrl.scenario_ID == 2 || sim_ctrl.scenario_ID == 3)
%         ring_num = 1000;
%         ring_tx = zeros(2,ring_num);
%         TX_radius = sim_ctrl.max_tx_range;
%         for ring_idx = 1:1000
%             phase_res = ring_idx/1000*2*pi;
%             ring_tx(:,ring_idx) = [TX_radius*cos(phase_res) TX_radius*sin(phase_res)];
%         end
%         hold on;
%         plot(txLocation(1,:)/1e3, txLocation(2,:)/1e3, 'x', 'LineWidth',2.0, 'MarkerSize',8.0);
%         hold on;
%         axis([-1 1 -1 1]*12);
%         plot(rxLocation(1,:)/1e3, rxLocation(2,:)/1e3, 'ro', 'LineWidth',2.0, 'MarkerSize',8.0);
%         legend('TX User', 'RX User');
%         for indxLink = 1:size(rxLocation,2)
%             plot_arrow( txLocation(1,indxLink)/1e3,txLocation(2,indxLink)/1e3,rxLocation(1,indxLink)/1e3,rxLocation(2,indxLink)/1e3,'headwidth',0.00001,'headheight',0.00001);
%         end
%         for indxLink = 1:size(rxLocation,2)
%             plot(ring_tx(1,:)/1e3 + rxLocation(1,indxLink)/1e3, ring_tx(2,:)/1e3 + rxLocation(2,indxLink)/1e3, 'LineWidth',1.0, 'MarkerSize',6.0);
%         end
% 
%         xlabel('Unit: Km');
%         drawnow;            
%     end
end