classdef hSLSTGaxMultiFrequencySystemChannel < handle
%hSLSTGaxMultiFrequencySystemChannel Create a system channel object
%
% CHAN = hSLSTGaxMultiFrequencySystemChannel(NODES) returns a system
% channel object for an array of wlanNode and bluetoothLENode objects,
% NODES. This assumes all nodes can transmit and receive and channels
% between the nodes are reciprocal.
%
% CHAN = hSLSTGaxMultiFrequencySystemChannel(NODES,PROTLINKCHAN) uses the
% prototype channel for a link PROTLINKCHAN to create channels between all
% nodes. PROTLINKCHAN is a wlanTGaxChannel, wlanTGacChannel or
% wlanTGnChannel object.
%
% CHAN =
% hSLSTGaxMultiFrequencySystemChannel(...,ShadowFadingStandardDeviation=val)
% sets the shadow fading standard deviation in dB. The default is 0 dB (no
% shadow fading).
%
% CHAN = hSLSTGaxMultiFrequencySystemChannel(...,PathLossModel=model)
% sets the path loss model to 'free-space', 'residential', or 'enterprise'.
% The default is 'free-space').
%
%   hSLSTGaxMultiFrequencySystemChannel properties:
%
%   Channels   - Array of system channels; one channel per frequency.
%   ChannelFcn - Function handle to channel for all nodes and links in
%                the simulation
%
%   hSLSTGaxMultiFrequencySystemChannel methods:
%
%   getChannel - Returns the channel object for a signal and receiver
%   getLink    - Returns the link structure for a signal and receiver

%   Copyright 2022-2025 The MathWorks, Inc.

    properties
        %Channels Array of system channels; one channel per frequency
        %   The class of Channels depends on PHYAbstractionMethod specified
        %   in NODES:
        %    "none" - array of hSLSTGaxSystemChannel objects
        %    otherwise - array of hSLSTGaxAbtractSystemChannel objects
        Channels;
        %ChannelFcn Function handle to channel for all nodes and links in
        %the simulation
        ChannelFcn;
    end

    properties (Access=private)
        UseFullPHY = false;
        ChannelInfo;
        DeviceInfo;
        LinkActiveWithinChannel; % NumLinks-by-NumChannels logical matrix indicating if a link is active within a channel
        LinkNodeID; % NumLinks-by-2 matrix of the sorted Node IDs in each link
    end

    methods
        function obj = hSLSTGaxMultiFrequencySystemChannel(nodes,varargin)
            % CHAN = hSLSTGaxMultiFrequencySystemChannel(NODES) returns a system
            % channel object for an array of wlanNode and bluetoothLENode objects,
            % NODES. This assumes all nodes can transmit and receive and channels
            % between the nodes are reciprocal.
            %
            % CHAN = hSLSTGaxMultiFrequencySystemChannel(NODES,PROTLINKCHAN)
            % uses the prototype channel for a link PROTLINKCHAN to create
            % channels between all nodes.
            %
            % CHAN =
            % hSLSTGaxMultiFrequencySystemChannel(...,ShadowFadingStandardDeviation=val)
            % sets the shadow fading standard deviation in dB. The default
            % is 0 dB (no shadow fading).

            % Set the base channel properties - this is used for all links Note
            % the number of transmit and receive antennas, and channel
            % bandwidth are sets from the node configuration.
            prototypeChannel = wlanTGaxChannel(...
                'DelayProfile','Model-D',...
                'TransmitReceiveDistance',15,...
                'ChannelFiltering',false,...
                'OutputDataType','single', ...
                'EnvironmentalSpeed',0); % Stationary
            validate = true; % Validate channels by default
            nvpairStart = 1;
            if nargin>1
                if isa(varargin{1}, 'wlan.internal.ChannelBase')
                    prototypeChannel = varargin{1};
                    nvpairStart = 2;
                    if nargin>2 && islogical(varargin{2})
                        % hSLSTGaxMultiFrequencySystemChannel(nodes,chan,(true|false), ....)
                        validate = varargin{2};
                        nvpairStart = 3;
                    end
                elseif islogical(varargin{1})
                    % hSLSTGaxMultiFrequencySystemChannel(nodes,(true|false), ...)
                    validate = varargin{1};
                    nvpairStart = 2;
                end
            end

            % Validate the nodes
            if iscell(nodes)
                for idx = 1:numel(nodes)
                    validateattributes(nodes{idx},["bluetoothLENode","wlanNode"],{'scalar'},mfilename,"nodes");
                end
            else
                validateattributes(nodes(1),["bluetoothLENode","wlanNode"],{'scalar'},mfilename,"nodes");
                nodes = num2cell(nodes);
            end

            % Segregate the nodes to WLAN and Bluetooth nodes
            wlanNodes = nodes(cellfun(@(x) isa(x,"wlanNode"),nodes));
            numWLANNodes = numel(wlanNodes);
            leNodes = nodes(cellfun(@(x) isa(x,"bluetoothLENode"),nodes));
            numLENodes = numel(leNodes);

            % Select type of channel model depending on PHY abstraction used. If
            % Bluetooth LE nodes exist, always full PHY.
            obj.UseFullPHY = numLENodes>0 || wlanNodes{1}.PHYAbstractionMethod == "none";

            % Create vectors of parameters for all devices in all nodes
            freqAllDevices = [];                 % Center frequency
            bwAllDevices = [];                   % Channel bandwidth (Hz)
            nodeNumAllDevices = [];              % Node number (1, 2, ..., numNodes)
            nodeIDsAllDevices = [];              % Node identifiers
            deviceNumAllDevices = [];            % Number of device in each nodes (1, ..., numDevices)
            numAntsAllDevices = [];              % Number of antennas
            interferenceModelingAllDevices = []; % Interference modeling
            % Get PHY parameters of WLAN
            numWLANDevicesTotal = 0;
            for n = 1:numWLANNodes
                if numLENodes>0
                    assert((wlanNodes{n}.PHYAbstractionMethod=="none") == obj.UseFullPHY, ...
                        'When LE nodes are added to the simulation all WLAN nodes must be configured to use no PHY abstraction.')
                else
                    assert((wlanNodes{n}.PHYAbstractionMethod=="none") == obj.UseFullPHY, ...
                        'All nodes must be configured to use the same PHY abstraction method')
                end
                devCfg = getDeviceConfig(wlanNodes{n});
                for d = 1:numel(devCfg)
                    numWLANDevicesTotal = numWLANDevicesTotal+1;
                    freqAllDevices = [freqAllDevices wlanChannelFrequency(devCfg(d).BandAndChannel(2),devCfg(d).BandAndChannel(1))]; %#ok<*AGROW>
                    bwAllDevices = [bwAllDevices devCfg(d).ChannelBandwidth];
                    nodeNumAllDevices = [nodeNumAllDevices n];
                    nodeIDsAllDevices = [nodeIDsAllDevices wlanNodes{n}.ID];
                    deviceNumAllDevices = [deviceNumAllDevices d];
                    numAntsAllDevices = [numAntsAllDevices devCfg(d).NumTransmitAntennas];
                    interferenceModelingAllDevices = [interferenceModelingAllDevices devCfg(d).InterferenceModeling];
                end
            end
            % Get PHY parameters of LE
            for n = 1:numLENodes
                freqAllDevices = [freqAllDevices leNodes{n}.ReceiveFrequency];
                bwAllDevices = [bwAllDevices leNodes{n}.ReceiveBandwidth];
                nodeNumAllDevices = [nodeNumAllDevices n+numWLANNodes];
                nodeIDsAllDevices = [nodeIDsAllDevices leNodes{n}.ID];
                deviceNumAllDevices = [deviceNumAllDevices 1];
                numAntsAllDevices = [numAntsAllDevices 1];
                interferenceModelingAllDevices = [interferenceModelingAllDevices leNodes{n}.InterferenceModeling];
            end
            obj.DeviceInfo = table;
            obj.DeviceInfo.NodeID = nodeIDsAllDevices';
            obj.DeviceInfo.DeviceNumber = deviceNumAllDevices';
            obj.DeviceInfo.FrequencyStartStop =  permute(freqAllDevices'+permute([-0.5 0.5],[3 1 2]).*bwAllDevices',[1 3 2]);

            % Get all combinations of devices, i.e. all possible links.
            % Each row contains two indices, one for each device making up
            % a channel.
            linkDeviceIndices = nchoosek(1:numel(freqAllDevices),2);

            % Remove combinations of devices which are on the same node
            isSameNode = (nodeNumAllDevices(linkDeviceIndices(:,1))==nodeNumAllDevices(linkDeviceIndices(:,2)))';
            linkDeviceIndices = linkDeviceIndices(~isSameNode,:);

            % Get the channel bandwidth of all devices, for all combinations
            bwCombinations = bwAllDevices(linkDeviceIndices);

            % Get the frequencies of all devices combinations
            freqCombinations = freqAllDevices(linkDeviceIndices);

            % Calculate the start and end frequencies of the channels
            % required for each device, for all combinations
            rxStartFreq = freqCombinations-bwCombinations/2;
            rxEndFreq = freqCombinations+bwCombinations/2;

            % Get device interference modeling, for all combinations
            % All devices must use the same interference modeling
            if numel(unique(interferenceModelingAllDevices))>1
                error('All devices and links must use the same interference mode')
            end

            % Determine which links to model. For each pair of devices
            % generate a packet from the first device to the second to
            % determine where there is a frequency overlap and it should be
            % modeled (this can be co-channel, overlapping, or adjacent).
            packet = wirelessnetwork.internal.wirelessPacket;
            packet.Abstraction = true; % Set to avoid validation
            packet.Metadata.OversamplingFactor = 1.125; % Set to avoid validation
            modelChannel = false(height(linkDeviceIndices),1);
            for i = 1:height(linkDeviceIndices)
                packet.CenterFrequency = freqCombinations(i,1);
                packet.Bandwidth = bwCombinations(i,1);
                rxNodeNum = nodeNumAllDevices(linkDeviceIndices(i,2));
                deviceNum = deviceNumAllDevices(linkDeviceIndices(i,2));
                node = nodes{rxNodeNum};
                if isa(node,"wlanNode")
                    packet.Type = 1; % WLAN
                    % Create channels between all nodes which overlap in
                    % frequency, therefore ignore co-channel interference.
                    interferenceFidelity = node.InterferenceFidelity;
                    node.InterferenceFidelity = max(node.InterferenceFidelity,1); % Force interference model to not be co-channel
                    modelChannel(i) = wlan.internal.sls.isFrequencyOverlapping(node,packet,deviceNum);
                    node.InterferenceFidelity = interferenceFidelity; % Restore
                else
                    packet.Type = 3; % Bluetooth LE
                    % Initialize
                    modelChannel(i) = false;
                    txStartFrequency = packet.CenterFrequency - packet.Bandwidth/2;
                    txEndFrequency = packet.CenterFrequency + packet.Bandwidth/2;

                    % Get the max interference offset based on the interference modeling
                    if strcmp(node.InterferenceModeling,node.InterferenceModeling_Values(1))
                        offset = 0;
                    else
                        offset = node.MaxInterferenceOffset;
                    end

                    % Invoke channel if signal lies in the 2.4 GHz band. The 2.4 GHz
                    % band starts at 2.4 GHz and ends at 2.4835 GHz.
                    if (txStartFrequency >= node.BluetoothLEStartBand-offset && txStartFrequency <= node.BluetoothLEEndBand+offset) || ...
                            (txEndFrequency >= node.BluetoothLEStartBand-offset && txEndFrequency <= node.BluetoothLEEndBand+offset)
                        % Model the channel
                        modelChannel(i) = true;
                    end
                end
            end

            linkDeviceIndToModel = linkDeviceIndices(modelChannel,:); % Device index for each channel
            linkNodeIDs = nodeIDsAllDevices(linkDeviceIndToModel);
            deviceNumToModel = deviceNumAllDevices(linkDeviceIndToModel);
            deviceFreqsToModel = freqAllDevices(linkDeviceIndToModel);
            deviceBWToModel = bwAllDevices(linkDeviceIndToModel);

            % Store the node IDs within this link sorted
            obj.LinkNodeID = sort(linkNodeIDs,2);

            % The combined span of the bandwidth of two devices is used to
            % model the channel between them. Therefore set the channel
            % bandwidth to this total span.
            channelStartFreq = min(rxStartFreq(modelChannel,:),[],2);
            channelStopFreq = max(rxEndFreq(modelChannel,:),[],2);
            overallChannelBandwidth = channelStopFreq-channelStartFreq;

            % Calculate the center frequencies of each channel, given the
            % span of the required bandwidth.
            centerFrequencies = (channelStartFreq+channelStopFreq)./2;

            % Get the index of the available channel bandwidth which the
            % overall channel bandwidth will fit in. If it exceeds the
            % maximum channel bandwidth string, then use the largest. Note
            % max() returns the first element in a row if the values are
            % the same - which is what we need.
            availableBWs = [20 40 80 160 320 Inf]*1e6;
            availableCBWs = ["CBW20" "CBW40" "CBW80" "CBW160" "CBW320" "CBW320"];
            [~,bwIndex] = max(availableBWs>=overallChannelBandwidth,[],2);
            channelBandwidths = availableCBWs(bwIndex);

            % Combine channels with the same center frequencies
            uniqueCenterFrequencies = unique(centerFrequencies);
            numChannels = numel(uniqueCenterFrequencies);

            % Create a channel manager for each unique channel center
            % frequency
            if obj.UseFullPHY
                obj.Channels = hSLSTGaxSystemChannel.empty(1,0);
            else
                obj.Channels = hSLSTGaxAbstractSystemChannel.empty(1,0);
            end
            obj.LinkActiveWithinChannel = false(height(linkNodeIDs),numChannels);
            for i = 1:numChannels
                % Create a channel between all devices operating within the
                % same channel. This includes adjacent/overlapping channels
                % if enabled.

                % Set the bandwidth and frequency of the channel. For the
                % channel, get the maximum bandwidth of all links using
                % that channel.
                isFreq = centerFrequencies==uniqueCenterFrequencies(i);
                [maxBW,maxBWIdx] = max(overallChannelBandwidth(isFreq));
                linkCBWsWithinThisChannel = channelBandwidths(isFreq);
                prototypeChannel.SampleRate = maxBW; % Used to create the path filters used in abstract PHY
                prototypeChannel.ChannelBandwidth = linkCBWsWithinThisChannel(maxBWIdx);
                prototypeChannel.CarrierFrequency = uniqueCenterFrequencies(i);

                % Get the node IDs of the pairs of nodes in this channel.
                linkNodeIDsThisChannel = linkNodeIDs(isFreq,:);

                % Indicate the link is active within this channel
                obj.LinkActiveWithinChannel(isFreq,i) = true;

                % Capture channel info for validation and plotting
                channelInfo = struct;
                channelInfo.LinkNodeIDs = linkNodeIDsThisChannel;
                channelInfo.LinkDevice = deviceNumToModel(isFreq,:);
                channelInfo.LinkDeviceFrequencies = deviceFreqsToModel(isFreq,:);
                channelInfo.LinkDeviceIndices = linkDeviceIndToModel(isFreq,:);
                channelInfo.linkDeviceBW = deviceBWToModel(isFreq,:);
                channelInfo.CarrierFrequency = prototypeChannel.CarrierFrequency;
                channelInfo.ChannelBandwidth = prototypeChannel.ChannelBandwidth;
                channelInfo.SampleRate = prototypeChannel.SampleRate;
                obj.ChannelInfo = [obj.ChannelInfo channelInfo];

                % Create numAntennasChannel, which contiains the number of
                % antennas for each device used in this channel.
                devicesInChannel = unique(linkDeviceIndToModel(isFreq,:));
                numAntennasChannel = numAntsAllDevices(devicesInChannel);

                if obj.UseFullPHY
                    obj.Channels(i) = hSLSTGaxSystemChannel(prototypeChannel,numAntennasChannel,'ChannelIndicesLUT',linkNodeIDsThisChannel,varargin{nvpairStart:end});
                else
                    obj.Channels(i) = hSLSTGaxAbstractSystemChannel(prototypeChannel,numAntennasChannel,'ChannelIndicesLUT',linkNodeIDsThisChannel,varargin{nvpairStart:end});
                end
            end

            if validate
                validateChannels(obj);
            end

            % Function handle to return impaired signal
            obj.ChannelFcn = @(rxInfo,signal)impairSignal(obj,signal,rxInfo);
        end
    end

    methods
        function channel = getChannel(obj,sig,rxInfo)
            % Returns the system channel object for a link between
            % signal SIG and receiver info RXINFO

            % Get all channels between the transmitter and receiver - there
            % may be multiple devices connected between a transmitter and
            if sig.TransmitterID>rxInfo.ID
                % LinkNodeID is stored sorted, so reverse IDs if required
                idsToTest = [rxInfo.ID sig.TransmitterID];
            else
                idsToTest = [sig.TransmitterID rxInfo.ID];
            end
            linkMatch = all(obj.LinkNodeID==idsToTest,2);
            candidateIndices = any(obj.LinkActiveWithinChannel(linkMatch,:),1);
            if ~any(candidateIndices)
                error('hSLSTGaxMultiFrequencySystemChannel:NoChannelExists','Channel does not exist between Node ID %d and Node ID %d.',sig.TransmitterID,rxInfo.ID)
            end

            % Select the channel with the closest center frequency to the
            % signal as the desired channel
            candidateChannels = [obj.Channels(candidateIndices)];
            [~,idx] = min(abs(sig.CenterFrequency-[candidateChannels.CenterFrequency]));
            channel = candidateChannels(idx);
        end

        function link = getLink(obj,sig,rxInfo)
            % Returns the link structure between signal SIG the receiver
            % info RXINFO

            % Get system channel
            channel = getChannel(obj,sig,rxInfo);

            % Extract channel information
            link = getLink(channel,sig.TransmitterID,rxInfo.ID);
        end
    end

    methods (Access=private)
        function [sig,varargout] = impairSignal(obj,sig,rxInfo)
            %impairSignal Apply path loss, log-normal shadow fading, and
            %frequency selective fading to the packet and update relevant
            %fields of the output data.

            % Get system channel between devices
            channel = getChannel(obj,sig,rxInfo);

            % Model path loss
            sig = pathLoss(obj,channel,sig,rxInfo);

            % Model shadow fading
            sig = shadowFading(obj,channel,sig,rxInfo);

            % Model frequency-selective fading
            if obj.UseFullPHY
                nargs = nargout;
                [sig,varargout{1:nargs-1}] = applyChannelToSignalStructure(channel,sig,rxInfo);
            else
                sig.Metadata.Channel = getChannelStatistics(channel,sig,rxInfo);
            end
        end

        function sig = pathLoss(obj, channel, sig, rxInfo)
            %pathLoss Apply path loss to the packet and update relevant fields
            %of the output data.

            pl = getPathLoss(channel,sig,rxInfo); % dB

            % Apply path loss on the power of the packet
            sig.Power = sig.Power - pl;

            if obj.UseFullPHY
                % Scale signal by path loss
                sig.Data = sig.Data*db2mag(-pl);
            end
        end

        function sig = shadowFading(obj, channel, sig, rxInfo)
            %shadowFading Apply log-normal shadow fading to the packet and
            %update relevant fields of the output data.

            txIdx = sig.TransmitterID;
            rxIdx = rxInfo.ID;
            l = getShadowFading(channel,txIdx,rxIdx); % dB
            % Apply shadow fading to the power of the packet
            sig.Power = sig.Power + l; % power in dBm

            if obj.UseFullPHY
                % Scale signal by shadow fading
                sig.Data = sig.Data*db2mag(l);
            end
        end

        function validateChannels(obj)
            for ic = 1:numel(obj.ChannelInfo)
                channelInfo = obj.ChannelInfo(ic);
                % Validate that multiple devices in a node do not share the
                % same channel with another node. A transmission from one
                % node cannot be simultaneously received by multiple
                % devices on a single node.
                uniqueNodeLinks = unique(channelInfo.LinkNodeIDs,"rows");
                if height(uniqueNodeLinks)~=height(channelInfo.LinkNodeIDs)
                    % This condition is met if multiple device links exist
                    % between nodes in this one channel

                    for iu = 1:height(uniqueNodeLinks)
                        % Loop over the pair of nodes in a link until we
                        % find one with multiple devices
                        match = all(uniqueNodeLinks(iu,:)==channelInfo.LinkNodeIDs,2);
                        if nnz(match)>1
                            error('hSLSTGaxMultiFrequencySystemChannel:MultipleDevicesShareChannel', ...
                                'Multiple devices or links in a node cannot receive a transmission from another node. Set BandAndChannel of devices or links in Node ID %d and Node ID %d so that the operating channels do not overlap within a node.', ...
                                uniqueNodeLinks(iu,1),uniqueNodeLinks(iu,2));
                        end
                    end
                end
            end
        end
    end
end

function devCfg = getDeviceConfig(node)
    %getDeviceConfig Returns the object holding MAC/PHY configuration
    %
    %   DEVCFG = getDeviceConfig(NODE) returns the object that holds the
    %   MAC/PHY configuration.
    %
    %   DEVCFG is an object of type wlanDeviceConfig if the input is a non-MLD
    %   node and it is an object of type wlanLinkConfig if the input is an MLD
    %   node.
    %
    %   NODE is an object of type wlanNode.

    if isa(node.DeviceConfig, "wlanMultilinkDeviceConfig")
        devCfg = node.DeviceConfig.LinkConfig;
        if strcmp(node.DeviceConfig.Mode, "STA") && isprop(node.DeviceConfig,'EnhancedMultilinkMode') && strcmp(node.DeviceConfig.EnhancedMultilinkMode, "EMLSR")
            % EMLSR case. For each device link set the number of antennas
            % to the sum for all devices. This will generate a channel for
            % each link (band) using the maximum number of antennas
            % possible. When we send data though the channel we will only
            % use the number of antennas appropriate at that time.
            totalNumAntennas = sum([devCfg.NumTransmitAntennas]);
            for i = 1:numel(devCfg)
                devCfg(i).NumTransmitAntennas = totalNumAntennas;
            end
        end
    else
        devCfg = node.DeviceConfig;
    end
end
