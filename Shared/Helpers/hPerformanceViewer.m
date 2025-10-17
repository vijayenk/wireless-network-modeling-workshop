classdef hPerformanceViewer < handle
    %helperPerformanceViewer Performance metrics viewer (throughput, latency, and
    %packet loss ratio)
    %
    %   helperPerformanceViewer(NODES,SIMULATIONTIME) plots the packet loss
    %   ratio, throughput, and average application packet latency for the
    %   Bluetooth and WLAN nodes specified.
    %
    %   NODES is a normal array or a cell array of objects of type
    %   bluetoothNode, bluetoothLENode, wlanNode, or helperCoexNode.
    %
    %   SIMULATIONTIME is a finite positive scalar indicating the simulation
    %   time in seconds.
    %
    %   helperPerformanceViewer methods:
    %
    %   plotNetworkStats        - Plots the performance metrics (throughput,
    %                             average packet latency, and packet loss ratio)
    %   throughput              - Returns transmission throughput of all the nodes
    %   averageReceiveLatency   - Returns average receive latency of all the nodes
    %   packetLossRatio         - Returns packet loss ratio of all the nodes

    %   Copyright 2023-2025 The MathWorks, Inc.

    properties (SetAccess=private)
        %WLANNodes List of WLAN nodes for visualization
        WLANNodes = wlanNode.empty()

        %BluetoothBREDRNodes List of Bluetooth BR/EDR nodes for visualization
        BluetoothBREDRNodes

        %BluetoothLENodes List of Bluetooth LE nodes for visualization
        BluetoothLENodes

        %CoexNodes List of coexistence nodes for visualization
        CoexNodes

        %AveragePacketLatency This property specifies the average latency computed
        %at each node in seconds. The first column represents node ID, second
        %column represents the average latency value, third column represents the
        %node type, and the fourth column represents the device type.
        AveragePacketLatency

        %Throughput This property specifies the throughput computed at each node in
        %seconds. The first column represents node ID, second column represents the
        %throughput value, third column represents the node type, and the fourth
        %column represents the device type.
        Throughput

        %PacketLossRatio This property specifies the packet loss ratio computed at
        %each node in seconds. The first column represents node ID, second column
        %represents the PLR value, third column represents the node type, and the
        %fourth column represents the device type.
        PacketLossRatio

        %PacketLatencyAll This property specifies the latency value of all the
        %received packet in seconds. Node IDs indicating bluetoothLENode is not
        %populated.
        PacketLatencyAll
    end

    properties (Access=private)
        %pPacketLatency This property specifies the latency computed at each node
        %in seconds. The first column represents node ID, second column represents
        %the latency values of the received packets, third column represents the
        %node type, and the fourth column represents the device type. Node IDs
        %indicating bluetoothLENode is not populated.
        pPacketLatency

        %pNumPackets This property specifies the number of application packets
        %received at each node. The first column represents node ID, second column
        %represents the number of packets received, third column represents the
        %node type, and the fourth column represents the device type.
        pNumPackets

        %pSimulationTime Simulation time in seconds
        pSimulationTime = 0

        %pNodeTypeList This property is a vector specifying the node IDs and the
        %type of the nodes. The first column represents the node ID, the second
        %column represents the node type, and the third column represents the
        %device type. This property is used for identifying the nodes and devices
        %for reference.
        pNodeTypeList

        %pNodes This property is a cell array containing all the nodes added in the
        %same order as the node IDs stored in pNodeTypeList.
        pNodes
    end

    %% Constructor
    methods
        function obj = hPerformanceViewer(nodes,simulationTime)

            % Validate the simulation time
            validateattributes(simulationTime,{'numeric'},{'nonempty','scalar','positive','finite'},mfilename,"simulationTime");

            % Return if no nodes are available
            if isempty(nodes)
                return;
            end

            % Validate the nodes
            if iscell(nodes)
                for idx = 1:numel(nodes)
                    validateattributes(nodes{idx},["bluetoothLENode","bluetoothNode","wlanNode","helperCoexNode"],{'scalar'},mfilename,"nodes");
                end
            else
                validateattributes(nodes(1),["bluetoothLENode","bluetoothNode","wlanNode","helperCoexNode"],{'scalar'},mfilename,"nodes");
                nodes = num2cell(nodes);
            end

            % Use weak-references for cross-linking handle objects
            objWeakRef = matlab.lang.WeakReference(obj);

            % Store the nodes
            [numBluetoothBREDRNodes,numBluetoothLENodes,numWLANNodes,numCoexNodes] = deal(0);
            numNodes = numel(nodes);
            for idx = 1:numNodes
                if isa(nodes{idx},"helperCoexNode")
                    numCoexNodes = numCoexNodes+1;
                    obj.CoexNodes{numCoexNodes} = nodes{idx};
                    addlistener(nodes{idx},"AppDataReceived",@(srcNode,eventData) objWeakRef.Handle.calculatePacketLatency(srcNode,eventData));
                elseif isa(nodes{idx},"bluetoothNode")
                    numBluetoothBREDRNodes = numBluetoothBREDRNodes+1;
                    obj.BluetoothBREDRNodes{numBluetoothBREDRNodes} = nodes{idx};
                    addlistener(nodes{idx},"AppDataReceived",@(srcNode,eventData) objWeakRef.Handle.calculatePacketLatency(srcNode,eventData));
                elseif isa(nodes{idx},"bluetoothLENode")
                    numBluetoothLENodes = numBluetoothLENodes+1;
                    obj.BluetoothLENodes{numBluetoothLENodes} = nodes{idx};
                    if isa(nodes{idx},"helperBluetoothLE6GHzNode")
                        addlistener(nodes{idx},"AppDataReceived",@(srcNode,eventData) objWeakRef.Handle.calculatePacketLatency(srcNode,eventData));
                    end
                else
                    numWLANNodes = numWLANNodes+1;
                    obj.WLANNodes(numWLANNodes) = nodes{idx};
                    addlistener(nodes{idx},"AppDataReceived",@(srcNode,eventData) objWeakRef.Handle.calculatePacketLatency(srcNode,eventData));
                end
            end

            % Initialize the throughput and latency
            obj.pNodeTypeList = zeros(numNodes,3);
            for idx = 1:numBluetoothBREDRNodes
                obj.pNodeTypeList(idx,1) = obj.BluetoothBREDRNodes{idx}.ID;
                obj.pNodeTypeList(idx,2:3) = 4;
                obj.pNodes{idx} = obj.BluetoothBREDRNodes{idx};
            end
            for idx = 1:numBluetoothLENodes
                obj.pNodeTypeList(numBluetoothBREDRNodes+idx,1) = obj.BluetoothLENodes{idx}.ID;
                obj.pNodeTypeList(numBluetoothBREDRNodes+idx,2:3) = 3;
                obj.pNodes{numBluetoothBREDRNodes+idx} = obj.BluetoothLENodes{idx};
            end
            for idx = 1:numWLANNodes
                obj.pNodeTypeList(numBluetoothBREDRNodes+numBluetoothLENodes+idx,1) = obj.WLANNodes(idx).ID;
                obj.pNodeTypeList(numBluetoothBREDRNodes+numBluetoothLENodes+idx,2:3) = 1;
                obj.pNodes{numBluetoothBREDRNodes+numBluetoothLENodes+idx} = obj.WLANNodes(idx);
            end

            % Attach ID and packet type for coexistence nodes
            countDevices = 1;
            for idx = 1:numCoexNodes
                obj.pNodeTypeList(numBluetoothBREDRNodes+numWLANNodes+numBluetoothLENodes+countDevices,2) = 5;
                if ~isempty(obj.CoexNodes{idx}.BluetoothDevice.DeviceName)
                    obj.pNodeTypeList(numBluetoothBREDRNodes+numWLANNodes+numBluetoothLENodes+countDevices,1) = obj.CoexNodes{idx}.ID;
                    obj.pNodeTypeList(numBluetoothBREDRNodes+numWLANNodes+numBluetoothLENodes+countDevices,3) = 4;
                    obj.pNodes{numBluetoothBREDRNodes+numWLANNodes+numBluetoothLENodes+countDevices} = obj.CoexNodes{idx};
                    countDevices = countDevices+1;
                end
                if ~isempty(obj.CoexNodes{idx}.BluetoothLEDevice.DeviceName)
                    obj.pNodeTypeList(numBluetoothBREDRNodes+numWLANNodes+numBluetoothLENodes+countDevices,1) = obj.CoexNodes{idx}.ID;
                    obj.pNodeTypeList(numBluetoothBREDRNodes+numWLANNodes+numBluetoothLENodes+countDevices,3) = 3;
                    obj.pNodes{numBluetoothBREDRNodes+numWLANNodes+numBluetoothLENodes+countDevices} = obj.CoexNodes{idx};
                    countDevices = countDevices+1;
                end
                if ~isempty(obj.CoexNodes{idx}.WLANDevice.DeviceName)
                    obj.pNodeTypeList(numBluetoothBREDRNodes+numWLANNodes+numBluetoothLENodes+countDevices,1) = obj.CoexNodes{idx}.ID;
                    obj.pNodeTypeList(numBluetoothBREDRNodes+numWLANNodes+numBluetoothLENodes+countDevices,3) = 1;
                    obj.pNodes{numBluetoothBREDRNodes+numWLANNodes+numBluetoothLENodes+countDevices} = obj.CoexNodes{idx};
                    countDevices = countDevices+1;
                end
            end
            numDevices = size(obj.pNodeTypeList,1);
            obj.PacketLatencyAll = cell(numDevices,2);
            obj.PacketLatencyAll(:,1) = num2cell(obj.pNodeTypeList(:,1));
            obj.pPacketLatency = [obj.pNodeTypeList(:,1) zeros(numDevices,1) obj.pNodeTypeList(:,2) zeros(numDevices,1)];
            [obj.AveragePacketLatency, obj.Throughput, obj.PacketLossRatio, obj.pNumPackets] = deal(obj.pPacketLatency);

            % Schedule action to perform at the end of simulation
            obj.pSimulationTime = simulationTime;
            networkSimulator = wirelessNetworkSimulator.getInstance;
            if numel(obj.BluetoothBREDRNodes)>0
                scheduleAction(networkSimulator,@(varargin) objWeakRef.Handle.calculateStats(obj.BluetoothBREDRNodes),[],obj.pSimulationTime);
            end
            if numel(obj.BluetoothLENodes)>0
                scheduleAction(networkSimulator,@(varargin) objWeakRef.Handle.calculateStats(obj.BluetoothLENodes),[],obj.pSimulationTime);
            end
            if numel(obj.WLANNodes)>0
                scheduleAction(networkSimulator,@(varargin) objWeakRef.Handle.calculateStats(num2cell(obj.WLANNodes)),[],obj.pSimulationTime);
            end
            if numel(obj.CoexNodes)>0
                scheduleAction(networkSimulator,@(varargin) objWeakRef.Handle.calculateStats(obj.CoexNodes),[],obj.pSimulationTime);
            end
        end

        function plotNetworkStats(obj)
            %plotNetworkStats Plots the performance metrics (throughput, average
            %receive latency, and packet loss ratio)
            %
            %   plotNetworkStats(OBJ) plots throughput (in Mbps), average packet
            %   latency (in seconds), and packet loss ratio at each node.
            %
            %   OBJ is an object of type helperPerformanceViewer.

            % Calculate and plot the graph for WLAN Nodes
            if numel(obj.WLANNodes)>0
                [nodeNames,plr,throughput,latency] = getMetrics(obj,num2cell(obj.WLANNodes));
                figHandle = plotBarGraph(obj,nodeNames,plr,throughput,latency,"Mbps",[0 ones(1,numel(obj.WLANNodes))],"Node Name");
                figHandle.Tag = "WLAN Performance";
                sgtitle("Performance of WLAN Nodes");
            end

            % Calculate and plot the graph for Bluetooth Nodes
            if numel(obj.BluetoothBREDRNodes)>0 || numel(obj.BluetoothLENodes)>0
                bluetoothNodes = [obj.BluetoothBREDRNodes obj.BluetoothLENodes];
                [nodeNames,plr,throughput,latency] = getMetrics(obj,bluetoothNodes);
                figHandle = plotBarGraph(obj,nodeNames,plr,throughput,latency,"Mbps",[0 ones(1,numel(bluetoothNodes))],"Node Name");
                figHandle.Tag = "Bluetooth Performance";
                sgtitle("Performance of Bluetooth Nodes");
            end

            % Calculate and plot the graph for coexistence Nodes
            if numel(obj.CoexNodes)>0
                nodeNames = []; plr = []; throughput = []; latency = [];
                numDevices = zeros(numel(obj.CoexNodes)+1,1);
                for idx = 1:numel(obj.CoexNodes)
                    [stackNames,stackPLR,stackThroughput,stackLatency] = getMetrics(obj,obj.CoexNodes(idx));
                    nodeNames = [nodeNames stackNames]; %#ok<*AGROW>
                    plr = [plr stackPLR];
                    throughput = [throughput stackThroughput];
                    latency = [latency stackLatency];
                    numDevices(idx+1) = sum(obj.CoexNodes{idx}.NumDevicesType);
                end
                figHandle = plotBarGraph(obj,nodeNames,plr,throughput,latency,"Mbps",numDevices,"Device Name");
                figHandle.Tag = "Coexistence node Performance";
                sgtitle("Performance of Coexistence Nodes");
            end
        end

        function tput = throughput(obj, nodeID, varargin)
            %throughput Returns transmission throughput of given nodes.
            %
            %   TPUT = throughput(OBJ, NODEID) returns transmission throughput for the
            %   given node, NODEID.
            %
            %   TPUT is a vector of size 1-by-N where N is the number of values in the
            %   NODEID. If the NODEID specified, indicates a coexistence node of type
            %   helperCoexNode, TPUT, is returned as a cell array of size 1-by-N where
            %   N is the number of values in the NODEID. Each cell is a size 1-by-Nd
            %   where Nd is the number of devices in the coexistence node.
            %
            %   OBJ is an object of type helperPerformanceViewer.
            %
            %   NODEID is a scalar or a vector of integers indicating valid node
            %   identifier(s).
            %
            %   TPUT = throughput(..., LinkID=Value) returns transmission throughput
            %   for the specified link, in the given node, NODEID. LinkID indicates a
            %   valid link identifier in an MLD node specified by NODEID. If NODEID is
            %   a vector of values, the LinkID must be a vector of the same size,
            %   indicating valid link identifiers for the corresponding NODEID values.
            %   This is applicable when the specified node ID corresponds to a WLAN
            %   node.

            tput = calculateMetrics(obj,"Throughput",nodeID,varargin);

        end

        function plr = packetLossRatio(obj, nodeID, varargin)
            %packetLossRatio Returns packet loss ratio of given nodes.
            %
            %   PLR = packetLossRatio(OBJ, NODEID) returns packet loss ratio for the
            %   given node, NODEID.
            %
            %   PLR is a vector of size 1-by-N where N is the number of values in the
            %   NODEID. If the NODEID specified, indicates a coexistence node of type
            %   helperCoexNode, PLR, is returned as a cell array of size 1-by-N where N
            %   is the number of values in the NODEID. Each cell is a size 1-by-Nd
            %   where Nd is the number of devices in the coexistence node.
            %
            %   OBJ is an object of type helperPerformanceViewer.
            %
            %   NODEID is a scalar or a vector of integers indicating valid node
            %   identifier(s).
            %
            %   PLR = packetLossRatio(..., LinkID=Value) returns packet loss ratio for
            %   the specified link, in the given node, NODEID. LinkID indicates a valid
            %   link identifier in an MLD node specified by NODEID. If NODEID is a
            %   vector of values, the LinkID must be a vector of the same size,
            %   indicating valid link identifiers for the corresponding NODEID values.
            %   This is applicable when the specified node ID corresponds to a WLAN
            %   node.

            plr = calculateMetrics(obj,"PLR",nodeID,varargin);
        end

        function avgLatency = averageReceiveLatency(obj, nodeID, latencyMetricVal)
            %averageReceiveLatency Returns average receive application latency of given
            %nodes.
            %
            %   AVGLATENCY = averageReceiveLatency(OBJ, NODEID) returns average
            %   application receive latency for the given node, NODEID.
            %
            %   AVGLATENCY is a vector or cell array of size 1-by-N where N is the
            %   number of values in the NODEID. If the NODEID specified, indicates a
            %   coexistence node of type helperCoexNode, AVGLATENCY, is returned as a
            %   cell array of size 1-by-N where N is the number of values in the
            %   NODEID. Each cell is a size 1-by-Nd where Nd is the number of devices
            %   in the coexistence node.
            %
            %   OBJ is an object of type helperPerformanceViewer.
            %
            %   NODEID is a scalar or a vector of integers indicating valid node
            %   identifier(s).

            latencyMetric = "Average";
            narginchk(2,3)
            if nargin>2 
                possibleLatencyMetric = ["Average";"P95";"P99"];
                latencyMetric = validatestring(latencyMetricVal,possibleLatencyMetric,mfilename,"LatencyMetric");
            end

            avgLatency = calculateMetrics(obj,"Latency",nodeID,[],latencyMetric);
        end
    end

    %% Standard Based Methods
    methods (Access=private)
        function calculatePacketLatency(obj,srcNode,eventData)
            %calculatePacketLatency Calculate average packet latency

            notificationData = eventData.Data;
            % Find the node index in the stored aggregate packet latency
            nodeIdx = find(srcNode.ID==obj.pNodeTypeList(:,1));

            % Number of indices found will be >1 for coexistence node. Hence find the
            % index in the list based on the device
            if numel(nodeIdx)>1
                if isfield(notificationData,"PacketType")
                    nodeIdx = nodeIdx(obj.pNodeTypeList(nodeIdx,3)==notificationData.PacketType);
                else
                    nodeIdx = nodeIdx(obj.pNodeTypeList(nodeIdx,3)==4);
                end
            end

            latencyVal = notificationData.CurrentTime-notificationData.PacketGenerationTime;

            % Store the list of packet latencies of each packet
            packetLatencyList = obj.PacketLatencyAll{nodeIdx,2};
            packetLatencyList = [packetLatencyList latencyVal];
            obj.PacketLatencyAll{nodeIdx,2} = packetLatencyList;

            % Aggregate the packet latency
            obj.pPacketLatency(nodeIdx,2) = obj.pPacketLatency(nodeIdx,2)+latencyVal;

            % Increment the number of App packets received
            obj.pNumPackets(nodeIdx,2) = obj.pNumPackets(nodeIdx,2)+1;

            % Calculate the average packet latency
            obj.AveragePacketLatency(nodeIdx,2) = obj.pPacketLatency(nodeIdx,2)/obj.pNumPackets(nodeIdx,2);
        end

        function calculateStats(obj,nodes)
            %calculateStats Calculate all the statistics of the nodes added to the
            %performance viewer at the end of simulation

            % Initialize
            numNodes = numel(nodes);

            % Calculate the performance of all the nodes
            for idx = 1:numNodes
                % Get the index at the packet latency for the node
                nodeIdx = find(nodes{idx}.ID==obj.pNodeTypeList(:,1));

                % Get the node name and statistics of the node
                statsNode = statistics(nodes{idx});
                statsDevice = statsNode;

                % Get the device types of the node
                deviceList = obj.pNodeTypeList(nodeIdx,3);
                numDevices = numel(deviceList);

                % Initialize the metrics
                packetLossRatio = zeros(1,numDevices);
                throughput = zeros(1,numDevices);
                latency = zeros(1,numDevices);
                count = 1;

                % Calculate the metric based on the node and device type
                for packetTypeIdx = 1:numel(deviceList)
                    deviceType = deviceList(packetTypeIdx);
                    if deviceType==1
                        if isa(nodes{idx},"helperCoexNode")
                            statsDevice = statsNode.WLANDevice;
                        end
                        techName = "WLAN";
                    elseif deviceType==4
                        if isa(nodes{idx},"helperCoexNode")
                            statsDevice = statsNode.BluetoothDevice;
                        end
                        techName = "BREDR";
                    elseif deviceType==3
                        if isa(nodes{idx},"helperCoexNode")
                            statsDevice = statsNode.BluetoothLEDevice;
                        else
                            obj.AveragePacketLatency(nodeIdx,2) = statsNode.App.AveragePacketLatency;
                        end
                        techName = "LE";
                    end
                    packetLossRatio(count) = calculatePerformanceMetric(obj,"PLR",techName,statsDevice,obj.pSimulationTime,nodeIdx(packetTypeIdx));
                    throughput(count) = calculatePerformanceMetric(obj,"Throughput",techName,statsDevice,obj.pSimulationTime,nodeIdx(packetTypeIdx));
                    latency(count) = calculatePerformanceMetric(obj,"Latency",techName,statsDevice,obj.pSimulationTime,nodeIdx(packetTypeIdx),"Average");

                    count = count+1;
                end

                % Store the PLR and throughput at the respective node indices
                obj.PacketLossRatio(nodeIdx,2) = packetLossRatio;
                obj.Throughput(nodeIdx,2) = throughput;
                obj.AveragePacketLatency(nodeIdx,2) = latency;
            end
        end

        function metricVal = calculateMetrics(obj,metricName,nodeID,linkIDNameValue,latencyMetric)
            %calculateMetrics Calculate the performance metrics for the specified
            %performance metric and node ID

            % Validate the inputs
            narginchk(2,5)

            % Initialize
            metricVal = zeros(1,numel(nodeID));
            if nargin~=5
                latencyMetric = "Average";
            end

            % For loop for all the nodes specified
            for idx = 1:numel(nodeID)

                % Get the ID of the node
                idVal = nodeID(idx);

                % Find the index of the node in the performance viewer object
                nodeIdx = find(idVal==obj.pNodeTypeList(:,1));

                % Get the type of the node
                type = [obj.pNodeTypeList(nodeIdx,2)];
                if numel(type)>1
                    type = 5;
                end

                if isempty(type)
                    error(['Node with ID=' num2str(idVal) ' not found in the list of nodes added to this object.']);
                end

                % Get the statistics of the node specified
                switch type
                    case 1 % WLAN node
                        nodeIdxWLAN = find(idVal == [obj.WLANNodes(:).ID]);
                        statsNode = statistics(obj.WLANNodes(nodeIdxWLAN),"all");
                    case 3 % Bluetooth LE node
                        ids = zeros(1,numel(obj.BluetoothLENodes));
                        for btIdx = 1:numel(obj.BluetoothLENodes)
                            ids(btIdx) = obj.BluetoothLENodes{btIdx}.ID;
                        end
                        nodeIdxBluetooth = idVal == ids;
                        statsNode = statistics(obj.BluetoothLENodes{nodeIdxBluetooth});
                    case 4 % Bluetooth BR/EDR node
                        ids = zeros(1,numel(obj.BluetoothBREDRNodes));
                        for btIdx = 1:numel(obj.BluetoothBREDRNodes)
                            ids(btIdx) = obj.BluetoothBREDRNodes{btIdx}.ID;
                        end
                        nodeIdxBluetooth = idVal == ids;
                        statsNode = statistics(obj.BluetoothBREDRNodes{nodeIdxBluetooth});
                    case 5 % Coexistence node
                        ids = zeros(1,numel(obj.CoexNodes));
                        for coexIdx = 1:numel(obj.CoexNodes)
                            ids(coexIdx) = obj.CoexNodes{coexIdx}.ID;
                        end
                        nodeIdxCoex = idVal == ids;
                        statsNode = statistics(obj.CoexNodes{nodeIdxCoex});
                end

                % Calculate the performance metric for the node specified based on the
                % retrieved statistics
                switch type
                    case 1
                        % Calculate link level performance metric
                        if nargin>3 && ~isempty(linkIDNameValue)
                            assert((nargin-3==1) && strcmpi(linkIDNameValue{1},"LinkID"), 'Invalid N/V pair');
                            linkID = linkIDNameValue{2};
                            assert((numel(linkID)==numel(nodeID)), 'Number of elements in LinkID must be the same as the number of elements in NodeID');
                            if ~isa(obj.WLANNodes(nodeIdxWLAN).DeviceConfig, "wlanMultilinkDeviceConfig") || numel(obj.WLANNodes(nodeIdxWLAN).DeviceConfig.LinkConfig) < linkID(idx)
                                error('Specified node ID must belong to an MLD node and the link ID must be less than or equal to number of links in the MLD');
                            end
                            if strcmp(metricName,"Throughput")
                                metricVal(idx) = ((statsNode.MAC(1).Link(linkID(idx)).TransmittedPayloadBytes)*8*1e-6)/obj.pSimulationTime;
                            else % PLR
                                metricVal(idx) = ((statsNode.MAC(1).Link(linkID(idx)).RetransmittedDataFrames)/statsNode.MAC(1).Link(linkID(idx)).TransmittedDataFrames);
                            end
                        else % Calculate node level performance metric
                            metricVal(idx) = calculatePerformanceMetric(obj,metricName,"WLAN",statsNode,obj.pSimulationTime,nodeIdx,latencyMetric);
                        end
                    case 3
                        if ~isa(obj.pNodes{nodeIdx},"helperBluetoothLE6GHzNode")
                            if ~strcmp(latencyMetric,"Average")
                                error("Bluetooth LE nodes support only average latency")
                            end
                        end
                        metricVal(idx) = calculatePerformanceMetric(obj,metricName,"LE",statsNode,obj.pSimulationTime,nodeIdx,latencyMetric);
                    case 4
                        metricVal(idx) = calculatePerformanceMetric(obj,metricName,"BREDR",statsNode,obj.pSimulationTime,nodeIdx,latencyMetric);
                    case 5
                        deviceList = obj.pNodeTypeList(nodeIdx,3);
                        metricValCoexDevices = zeros(1,numel(deviceList));
                        if ~iscell(metricVal)
                            metricVal = num2cell(metricVal);
                        end
                        for packetTypeIdx = 1:numel(deviceList)
                            deviceType = deviceList(packetTypeIdx);
                            if isfield(statsNode,"WLANDevice") && deviceType==1
                                statsDevice = statsNode.WLANDevice;
                                techName = "WLAN";
                            elseif isfield(statsNode,"BluetoothDevice") && deviceType==4
                                statsDevice = statsNode.BluetoothDevice;
                                techName = "BREDR";
                            elseif isfield(statsNode,"BluetoothLEDevice") && deviceType==3
                                statsDevice = statsNode.BluetoothLEDevice;
                                techName = "LE";
                            end
                            metricValCoexDevices(packetTypeIdx) = calculatePerformanceMetric(obj,metricName,techName,statsDevice,obj.pSimulationTime,nodeIdx(packetTypeIdx),latencyMetric);
                        end
                        metricVal{idx} = metricValCoexDevices;
                end
            end
        end

        function metricVal = calculatePerformanceMetric(obj,metricName,techName,stats,simulationTime,nodeIdx,latencyMetric)
            %calculatePerformanceMetric Calculate and return the metric value for
            %specified performance metric and technology

            node = obj.pNodes{nodeIdx};
            isaLENode = isa(node,"bluetoothLENode") && ~isa(node,"helperCoexNode");
            if ~isaLENode
                if strcmp(metricName,"PLR")
                    metricVal = calculatePLR(obj,techName,stats,nodeIdx);
                elseif strcmp(metricName,"Throughput")
                    metricVal = calculateThroughput(obj,techName,stats,simulationTime,nodeIdx);
                else
                    if strcmp(latencyMetric,"Average")
                        metricVal = [obj.AveragePacketLatency(nodeIdx,2)];
                    else
                        latencyVal = obj.PacketLatencyAll{nodeIdx,2};
                        if ~isempty(latencyVal)
                            if strcmp(latencyMetric,"P99")
                                metricVal = prctile(latencyVal,99);
                            elseif strcmp(latencyMetric,"P95")
                                metricVal = prctile(latencyVal,95);
                            end
                        else
                            metricVal = 0;
                        end
                    end
                end
            else
                nodeConnectedIDs = obj.pNodes{nodeIdx}.ConnectedNodeIDs;
                for idx = 1:numel(nodeConnectedIDs)
                    connectedID = nodeConnectedIDs(idx);
                    connectedNodeIdx = find(connectedID==obj.pNodeTypeList(:,1));
                    if strcmp(metricName,"PLR")
                        metricVal(idx) = calculatePLR(obj,techName,stats,nodeIdx,connectedNodeIdx);
                    elseif strcmp(metricName,"Throughput")
                        metricVal(idx) = calculateThroughput(obj,techName,stats,simulationTime,nodeIdx,connectedNodeIdx)/1e3;
                    else
                        if isa(node,"helperBluetoothLE6GHzNode")
                            if strcmp(latencyMetric,"Average")
                                metricVal = [obj.AveragePacketLatency(nodeIdx,2)];
                            else
                                latencyVal = obj.PacketLatencyAll{nodeIdx,2};
                                if ~isempty(latencyVal)
                                    if strcmp(latencyMetric,"P99")
                                        metricVal = prctile(latencyVal,99);
                                    elseif strcmp(latencyMetric,"P95")
                                        metricVal = prctile(latencyVal,95);
                                    end
                                else
                                    metricVal = 0;
                                end
                            end
                            break;
                        else
                            metricVal(idx) = kpi(obj.pNodes{connectedNodeIdx},obj.pNodes{nodeIdx},"latency",Layer="App");
                        end
                    end
                end
                metricVal = mean(metricVal);
            end
        end

        function packetLossRatio = calculatePLR(obj,techName,stats,nodeIdx,connectedNodeIdx)
            %calculatePLR Calculate the packet loss ratio based on the technology

            switch techName
                case "WLAN"
                    packetLossRatio = sum(([stats.MAC.RetransmittedDataFrames])/[stats.MAC.TransmittedDataFrames]);
                case "BREDR"
                    transmittedDataPackets = 0; retransmittedDataPackets = 0;
                    for peripheralIdx = 1:numel(stats.Baseband.ConnectionStats)
                        statsConn = stats.Baseband.ConnectionStats(peripheralIdx);
                        transmittedDataPackets = transmittedDataPackets+statsConn.TransmittedACLPackets+ ...
                            statsConn.TransmittedDVPackets;
                        retransmittedDataPackets = retransmittedDataPackets+statsConn.RetransmittedACLPackets+ ...
                            statsConn.DVWithRetransmittedData;
                    end
                    packetLossRatio = retransmittedDataPackets/(transmittedDataPackets+retransmittedDataPackets);
                case "LE"
                    if nargin==4
                        for peripheralIdx = 1:numel(stats.LL)
                            statsConn = stats.LL(peripheralIdx);
                            transmittedDataPackets = statsConn.TransmittedDataPackets;
                            retransmittedDataPackets = statsConn.RetransmittedDataPackets;
                        end
                        packetLossRatio = retransmittedDataPackets/(transmittedDataPackets+retransmittedDataPackets);
                    else
                        packetLossRatio = kpi(obj.pNodes{nodeIdx},obj.pNodes{connectedNodeIdx},"PLR",Layer="LL");
                    end

            end
        end

        function throughput = calculateThroughput(obj,techName,stats,simulationTime,nodeIdx,connectedNodeIdx)
            %calculateThroughput Calculate the throughput based on the technology

            switch techName
                case "WLAN"
                    throughput = (sum([stats.MAC.TransmittedPayloadBytes])*8*1e-6)/simulationTime;
                case "BREDR"
                    transmittedPayloadBytes = 0;
                    for peripheralIdx = 1:numel(stats.Baseband.ConnectionStats)
                        statsConn = stats.Baseband.ConnectionStats(peripheralIdx);
                        transmittedPayloadBytes = transmittedPayloadBytes+statsConn.TransmittedDataBytes;
                    end
                    throughput = transmittedPayloadBytes*8*1e-6/simulationTime;
                case "LE"
                    if nargin==5
                        for peripheralIdx = 1:numel(stats.LL)
                            statsConn = stats.LL(peripheralIdx);
                            transmittedPayloadBytes = statsConn.TransmittedPayloadBytes;
                        end
                        throughput = transmittedPayloadBytes*8*1e-6/simulationTime;
                    else
                        throughput = kpi(obj.pNodes{nodeIdx},obj.pNodes{connectedNodeIdx},"throughput",Layer="LL");
                    end
            end
        end
    end

    %% Figure Based Methods
    methods (Access=private)
        function [nodeNames,packetLossRatio,throughput,averageAppPacketLatency] = getMetrics(obj,nodes)
            %getMetrics Retrieve the performance metrics for the nodes added to the
            %performance viewer for the plot

            % Initialize
            numNodes = numel(nodes);
            packetLossRatio = zeros(1,numNodes);
            throughput = zeros(1,numNodes);
            averageAppPacketLatency = zeros(1,numNodes);
            count = 1;
            nodeNames = "";

            % Get the performance of all the nodes
            for idx = 1:numNodes
                % Get the index at the packet latency for the node
                nodeIdx = find(nodes{idx}.ID==obj.pNodeTypeList(:,1));
                if ~isa(nodes{idx},"helperCoexNode")
                    nodeNames(idx) = nodes{idx}.Name;
                end

                deviceList = obj.pNodeTypeList(nodeIdx,3);
                for packetTypeIdx = 1:numel(deviceList)
                    if isa(nodes{idx},"helperCoexNode")
                        switch deviceList(packetTypeIdx)
                            case 1
                                techName = "WLAN";
                            case 3
                                techName = "LE";
                            case 4
                                techName = "BREDR";
                        end
                        if contains(nodes{idx}.Name,techName)
                            nodeNames(count) = nodes{idx}.Name;
                        else
                            nodeNames(count) = strjoin([nodes{idx}.Name techName]);
                        end
                    end

                    % Retrieve the performance metrics for the specified nodes
                    averageAppPacketLatency(count) = obj.AveragePacketLatency(nodeIdx(packetTypeIdx),2);
                    packetLossRatio(count) = obj.PacketLossRatio(nodeIdx(packetTypeIdx),2);
                    throughput(count) = obj.Throughput(nodeIdx(packetTypeIdx),2);
                    count = count+1;
                end
            end
        end

        function figHandle = plotBarGraph(~,nodeNames,plr,throughput,latency,throughputUnits,numDevices,xlabelVal)
            %plotBarGraph Plots the bar graph

            % Get screen resolution and create figure
            figHandle = figure;
            matlab.graphics.internal.themes.figureUseDesktopTheme(figHandle);
            numNodes = numel(nodeNames);
            if numNodes>15
                xTickAngleVal = 90;
            else
                xTickAngleVal = 20;
            end

            % Calculate how the bars need to be grouped. Same node device bars will be
            % spaced nearer
            xValue = zeros(sum(numDevices),1);
            xValueCount = 0;
            for idx = 2:numel(numDevices)
                xValue(xValueCount+1:xValueCount+numDevices(idx),:) = (1:numDevices(idx))+sum(numDevices(1:idx-1));
                xValueCount = xValueCount+numDevices(idx);
                numDevices(idx) = numDevices(idx)+0.5;
            end

            % Plot the throughput
            s1 = subplot(3,1,1);
            bar(s1,xValue,throughput,"BarWidth",0.95);
            title("Throughput at Each Node");
            xlabel(xlabelVal);
            xticks(xValue)
            xtickangle(xTickAngleVal);
            xticklabels(nodeNames);
            ylabel('Throughput ('+ string(throughputUnits) + ')');

            % Plot the Packet Loss Ratio
            s2 = subplot(3,1,2);
            bar(s2,xValue,plr,"BarWidth",0.95);
            title("Packet Loss at Each Node");
            xlabel(xlabelVal);
            xticks(xValue)
            xtickangle(xTickAngleVal);
            xticklabels(nodeNames);
            ylabel("Packet Loss Ratio");

            % Plot the Average Application Packet Latency
            s3 = subplot(3,1,3);
            bar(s3,xValue,latency,"BarWidth",0.95);
            title("Average Packet Latency at Each Node");
            xlabel(xlabelVal);
            xticks(xValue)
            xtickangle(xTickAngleVal);
            xticklabels(nodeNames);
            ylabel("Latency (s)");
        end
    end
end
