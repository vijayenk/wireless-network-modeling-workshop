classdef helperCaptureIQSamples < handle
    %helperCaptureIQSamples Captures the IQ samples of a Bluetooth, WLAN, 5G
    %NR, or coexistence node
    %
    %   helperCaptureIQSamples(NODE,varargin) captures the IQ samples for
    %   Bluetooth, WLAN, 5G NR, or coexistence node.
    %
    %   helperCaptureIQSamples(NODE,Name=Value) captures the IQ samples for
    %   Bluetooth, WLAN, 5G NR, or coexistence node, with the specified
    %   property Name set to the specified Value. You can specify additional
    %   name-value arguments in any order as (Name1=Value1, ..., NameN=ValueN).
    %
    %   NODE is the wireless node for which the IQ samples needs to be
    %   captured. This is a scalar object of type bluetoothNode,
    %   bluetoothLENode, wlanNode, nrGNB, nrUE or helperCoexNode.
    %
    %   helperCaptureIQSamples properties (Configurable as NV Pair):
    %
    %   FileName           - MAT file name to store IQ samples
    %   CapturePeriodicity - Periodicity of IQ sample capture
    %   CenterFrequency    - Center frequency to capture in Hz
    %   Bandwidth          - Bandwidth around the center frequency to capture in Hz

    %   Copyright 2023 The MathWorks, Inc.

    properties (SetAccess=private) % Configurable as NV pair
        %FileName MAT file name to store IQ samples
        %   Name of the MAT file in which IQ samples will be stored. The IQ samples
        %   will be stored in a variable "IQSamples". Specify this property as a
        %   character vector or string scalar. If the file name is not specified,
        %   the file name is by the format: <NodeName_NodeID_yyyyMMdd_HHmmss>.
        FileName {mustBeTextScalar} = ""

        %CapturePeriodicity Periodicity of IQ sample capture in seconds
        %   Periodicity at which the IQ samples must be retrieved and captured from
        %   the interference buffer in the PHY receiver. If this is not specified
        %   the default value is 0.125s. If the periodicity is less than simulation
        %   time, samples will be captured at the end of simulation. A high capture
        %   periodicity value will force the interference buffer to retain the
        %   packets in it for a duration of periodicity after the packet end time.
        %   This will result in a high memory usage by the interference buffer and
        %   might affect the performance of the simulation. Specify this property
        %   as a positive finite numeric scalar. The default value is 0.125 s.
        CapturePeriodicity (1,1) {mustBeNumeric,mustBePositive,mustBeFinite} = 0.125

        %CenterFrequency Center frequency to capture in Hz
        %   Center frequency to capture IQ samples. Specify this property as
        %   numeric nonnegative value in Hz. If the center frequency is not
        %   provided, for Bluetooth basic rate/enhanced data rate (BR/EDR) and low
        %   energy (LE) nodes, the default value is 2.441 MHz. For WLAN and 5G
        %   nodes (gNB and UE), the default value is the node's frequency of
        %   operation.
        CenterFrequency (1,:) {mustBeNumeric,mustBeReal,mustBeNonnegative,mustBeFinite} = 2.441e6

        %Bandwidth Bandwidth around the center frequency to capture in Hz
        %   Bandwidth around the center frequency to capture IQ samples. Specify
        %   this property as a numeric positive finite value. If the bandwidth is
        %   not provided, for Bluetooth BR/EDR and LE nodes, the default value is
        %   82 MHz which corresponds to default center frequency of 2.441e6. For
        %   WLAN node, the default value is the node's channel bandwidth of
        %   operation. For gNB and UE nodes, you must specify the bandwidth of
        %   channel on which the node is operating.
        Bandwidth (1,:) {mustBeNumeric,mustBeReal,mustBePositive,mustBeFinite} = 82e6
    end

    properties (SetAccess=private)
        %IQSamples Captured IQ Signal, stored as timeseries
        IQSamples
    end

    properties (Constant,Hidden)
        % Bluetooth LE channel center frequencies in Hz
        LECenterFrequencies = [2404:2:2424 2428:2:2478 2402 2426 2480]*1e6
    end

    properties (Access=private)
        %pData Data information of the IQ samples. The size of pData is [MxNxP],
        %where M represent the number of samples, N represents the number of
        %streams, and P represents the number of devices. Incase of Bluetooth, N
        %and P are always 1.
        pData

        %pTimestamps Timestamp for the IQ samples. The size of pTimestamps is
        %[Mx1xP], where M represent the number of samples and P represents the
        %number of devices. Incase of Bluetooth, P is always 1.
        pTimestamps

        %pStartTime Start time to capture IQ samples based on the current time
        pStartTime = 0

        %pIsInitialized Flag to know if the object is initialized
        pIsInitialized = false

        %pSize Size of the IQ sample data and timestamp
        pSize = 0

        %pExpectedSampleRate Expected sample rate of the IQ signals
        pExpectedSampleRate = 20e6

        %pInterferenceBuffer Interference buffers of coexistence node
        pInterferenceBuffer

        %pDeviceNames Device names of the node
        pDeviceNames = string.empty
    end

    methods
        function obj = helperCaptureIQSamples(node,varargin)
            % Constructor

            % Validate the specified node
            validateattributes(node,{'bluetoothNode','bluetoothLENode','wlanNode','nrGNB', 'nrUE','helperCoexNode'},"scalar",1);

            % Name-value pair check
            coder.internal.errorIf(mod(nargin-1,2)==1,"MATLAB:system:invalidPVPairs");

            % Default file name
            nodeName = strrep(node.Name," ","_");
            obj.FileName = strjoin([nodeName "_" node.ID "_" char(datetime("now","Format","yyyyMMdd_HHmmss"))],"");

            if isa(node,"helperCoexNode")
                obj.CenterFrequency = [];
                obj.Bandwidth = [];
                countDevice = 1;
                obj.pExpectedSampleRate = [];
                if ~isempty(node.WLANDevice.DeviceName)
                    devCfg = getWLANDeviceConfig(obj, node.WLANDevice);
                    obj.CenterFrequency = [obj.CenterFrequency node.WLANDevice.ReceiveFrequency];
                    obj.Bandwidth = [obj.Bandwidth devCfg.ChannelBandwidth];
                    obj.pExpectedSampleRate = [obj.pExpectedSampleRate devCfg.ChannelBandwidth];
                    obj.pDeviceNames(countDevice) = strjoin([nodeName "WLAN" num2str(node.WLANDevice.ReceiveFrequency/1e6) "MHz"],"_");
                    countDevice = countDevice+1;
                end
                if ~isempty(node.BluetoothDevice.DeviceName)
                    obj.CenterFrequency = [obj.CenterFrequency 2441e6];
                    obj.Bandwidth = [obj.Bandwidth 82e6];
                    obj.pExpectedSampleRate = [obj.pExpectedSampleRate 20e6];
                    obj.pDeviceNames(countDevice) = strjoin([nodeName "BREDR"],"_");
                    countDevice = countDevice+1;
                end
                if ~isempty(node.BluetoothLEDevice.DeviceName)
                    obj.CenterFrequency = [obj.CenterFrequency 2441e6];
                    obj.Bandwidth = [obj.Bandwidth 82e6];
                    obj.pExpectedSampleRate = [obj.pExpectedSampleRate 20e6];
                    obj.pDeviceNames(countDevice) = strjoin([nodeName "LE"],"_");
                end
                if isempty(obj.CenterFrequency)
                    error("Add a device to the coexistence node to capture IQ samples");
                end
            elseif isa(node,"wlanNode")
                obj.CenterFrequency = node.ReceiveFrequency;
                devCfg = getWLANDeviceConfig(obj, node);
                obj.Bandwidth = [devCfg.ChannelBandwidth];
                for idx=1:numel(obj.CenterFrequency)
                    obj.pDeviceNames(idx) = strjoin([nodeName num2str(node.ReceiveFrequency(idx)/1e6) "MHz"],"_");
                end
            elseif isa(node,"nrGNB")
                obj.CenterFrequency = node.ULCarrierFrequency;
                obj.pDeviceNames = nodeName;
            elseif isa(node,"nrUE")
                obj.CenterFrequency = node.DLCarrierFrequency;
                obj.pDeviceNames = nodeName;
            else
                obj.CenterFrequency = 2441e6;
                obj.Bandwidth = 82e6;
                obj.pDeviceNames = nodeName;
            end

            isBandwidthSupplied = false;
            % Assign name-value arguments
            for idx = 1:2:nargin-1
                if isa(node,"helperCoexNode") && (strcmp(varargin{idx},"CenterFrequency")||strcmp(varargin{idx},"Bandwidth"))
                    error(strjoin(["NV Pair" varargin{idx} "is not supported for helperCoexNode."]))
                end
                if strcmp(varargin{idx}, "Bandwidth")
                    isBandwidthSupplied = true;
                end
                obj.(char(varargin{idx})) = varargin{idx+1};
            end
            if (isa(node,"nrUE") || isa(node,"nrGNB")) && isBandwidthSupplied==false
                error("Bandwidth must be supplied for IQ sample capture on 5G nodes (nrGNB and nrUE)");
            end

            % Use weak-references for cross-linking handle objects
            objWeakRef = matlab.lang.WeakReference(obj);

            networkSimulator = wirelessNetworkSimulator.getInstance;
            % Schedule action to initialize the packet capture
            captureIQFcn = @(varargin) objWeakRef.Handle.captureSignals(node,networkSimulator,[]);
            scheduleAction(networkSimulator,captureIQFcn,[],0);

            % Validate and set Center frequency and bandwidth
            if isa(node,"wlanNode") % WLAN validation
                validateWLANNode(obj,node)
            elseif isa(node,"bluetoothNode") || isa(node,"bluetoothLENode") % LE and BREDR validation
                validateBluetoothNode(obj);
            elseif isa(node,"nrGNB") || isa(node,"nrUE") % For 5G nodes
                obj.pExpectedSampleRate = node.PhyEntity.RxBuffer.SampleRate;
                validate5GNode(obj);
            end
        end

        function value = get.IQSamples(obj)
            % Return the IQ samples as a time series
            value = timeseries.empty;
            for idx = 1:numel(obj.pData)
                value(idx) = timeseries(obj.pData{idx}(1:obj.pSize(idx),:),obj.pTimestamps{idx}(1:obj.pSize(idx),:),Name=obj.pDeviceNames(idx));
            end
        end
    end

    methods (Access=private)
        function captureSignals(obj,node,networkSimulator,duration)
            %captureSignal Captures the IQ signals at periodic intervals

            currentTime = networkSimulator.CurrentTime;
            endTime = networkSimulator.EndTime;
            if currentTime==0
                if ~obj.pIsInitialized
                    % Validate and update the capture periodicity
                    duration = obj.CapturePeriodicity;
                    if duration>endTime
                        duration = endTime;
                    end
                    obj.CapturePeriodicity = duration;

                    % Update the interference buffer to store the signal atleast for duration
                    % of periodicity after the signal end time
                    initializeIQData(obj,node,endTime,duration)

                    % Use weak-references for cross-linking handle objects
                    objWeakRef = matlab.lang.WeakReference(obj);

                    % Schedule action at network simulator to capture signals at a particular
                    % periodicity
                    captureIQFcn = @(varargin) objWeakRef.Handle.captureSignals(node,networkSimulator,duration);
                    scheduleAction(networkSimulator,captureIQFcn,[],0,duration);

                    % Schedule action at the network simulator to capture signal at the end of
                    % simulation if periodicity doesn't include end time
                    if round(mod(endTime,duration),9)~=0
                        captureIQFcn = @(varargin) objWeakRef.Handle.captureSignals(node,networkSimulator,round(mod(endTime,duration),9));
                        scheduleAction(networkSimulator,captureIQFcn,[],endTime);
                    end
                    obj.pIsInitialized = true;
                end
            else
                % Capture IQ samples
                if isa(node,"helperCoexNode")
                    % Capture IQ samples for Coex node
                    captureIQSamples(obj,obj.pInterferenceBuffer,duration,1);
                elseif isa(node,"wlanNode")
                    % Capture IQ samples for WLAN node. For WLAN node, number of Tx antennas is
                    % passed because the Tx and Rx antenna count is assumed to be same
                    devCfg = getWLANDeviceConfig(obj, node);
                    captureIQSamples(obj,[node.PHYRx.Interference],duration,[devCfg.NumTransmitAntennas]);
                elseif isa(node,"nrGNB") || isa(node,"nrUE")
                    % Capture IQ samples for 5G node
                    captureIQSamples(obj,node.PhyEntity.RxBuffer,duration,node.NumReceiveAntennas);
                else
                    % Capture IQ samples for Bluetooth node
                    captureIQSamples(obj,node.PHYReceiver.Interference,duration,1);
                end

                % At the end of simulation save the IQ samples as MAT file
                if currentTime==endTime
                    saveIQSamples(obj);
                end
            end
        end

        function initializeIQData(obj,node,endTime,duration)
            %initializeIQData Initialize the IQ data information for the node

            if isa(node,"helperCoexNode")
                if ~isempty(node.WLANDevice.DeviceName)
                    obj.pInterferenceBuffer = [obj.pInterferenceBuffer [node.WLANDevice.PHYRx.Interference]];
                end
                if ~isempty(node.BluetoothDevice.DeviceName)
                    obj.pInterferenceBuffer = [obj.pInterferenceBuffer node.BluetoothDevice.PHYReceiver.Interference];
                end
                if ~isempty(node.BluetoothLEDevice.DeviceName)
                    obj.pInterferenceBuffer = [obj.pInterferenceBuffer node.BluetoothLEDevice.PHYReceiver.Interference];
                end
                interferenceBuffer = obj.pInterferenceBuffer;
                numRxAnts = ones(numel(obj.CenterFrequency),1);
                obj.pSize = zeros(numel(obj.CenterFrequency),1);
                obj.pStartTime = zeros(numel(obj.CenterFrequency),1);
            elseif isa(node,"wlanNode") % WLAN node
                interferenceBuffer = [node.PHYRx.Interference];
                devCfg = getWLANDeviceConfig(obj, node);
                numRxAnts = [devCfg.NumTransmitAntennas]; % Tx and Rx antenna count is same
                obj.pSize = zeros(numel(node.ReceiveFrequency),1);
                obj.pStartTime = zeros(numel(node.ReceiveFrequency),1);
            elseif isa(node,"nrGNB") || isa(node,"nrUE") % 5G node
                interferenceBuffer = node.PhyEntity.RxBuffer;
                numRxAnts = node.NumReceiveAntennas;
            else % Bluetooth node
                interferenceBuffer = node.PHYReceiver.Interference;
                numRxAnts = 1;
            end
            for idx = 1:numel(obj.pExpectedSampleRate) % For each device initialize the IQ data information
                interferenceBuffer(idx).BufferCleanupTime = duration;
                obj.pData{idx}(:,:) = zeros(obj.pExpectedSampleRate(idx)*endTime,numRxAnts(idx));
                obj.pTimestamps{idx}(:,:) = zeros(obj.pExpectedSampleRate(idx)*endTime,1);
            end
        end

        function captureIQSamples(obj,interferenceBuffer,duration,numRxAntennas)
            %captureIQSignals Capture the IQ samples from the signals stored in the
            %interference buffer

            for idx=1:numel(obj.CenterFrequency) % For each device
                startTime = obj.pStartTime(idx);
                endTime = round(startTime+duration,9);

                % Combine the signals and return the resultant waveform
                [waveform,~,waveformSampleRate] = resultantWaveform(interferenceBuffer(idx),startTime,endTime, ...
                    CenterFrequency=obj.CenterFrequency(idx),Bandwidth=obj.Bandwidth(idx));
                expectedSampleRate = obj.pExpectedSampleRate(idx);

                if ~isempty(waveform)
                    % Resample the signal for the expected sample rate
                    if waveformSampleRate~=expectedSampleRate
                        [L,M] = rat(expectedSampleRate/waveformSampleRate);
                        waveform = resample(waveform,L,M);
                    end
                    % Calculate the associated timestamps. Rounding to 15 for better accurate
                    % timestamp
                    timeStep = round(1/expectedSampleRate,15);
                    timestamps = startTime:timeStep:endTime;
                    timestamps = timestamps(1:size(waveform,1))';

                    % Store the IQ data
                    storeIQData(obj,waveform,timestamps,startTime,idx);
                else
                    % Fill empty waveform for the duration
                    fillEmptyWaveform(obj,expectedSampleRate,startTime,duration,numRxAntennas(idx),idx);
                end
                % Update the start time for the next capture
                obj.pStartTime(idx) = obj.pStartTime(idx)+duration;
            end
        end

        function fillEmptyWaveform(obj,expectedSampleRate,startTime,duration,numRxAntennas,idx)
            %fillEmptyWaveform Fills empty waveform for the specified duration

            timestamps = (startTime:1/expectedSampleRate:startTime+duration)';
            timestamps = timestamps(1:round(duration*expectedSampleRate));
            if ~isempty(timestamps)
                storeIQData(obj,zeros(size(timestamps,1),numRxAntennas),timestamps,startTime,idx);
            end
        end

        function storeIQData(obj,waveform,timestamps,startTime,idx)
            %storeIQData Store the IQ waveform and corresponding timestamps

            % Check if there is a time overlap
            isOverlap = checkOverlap(obj,obj.pSize(idx),obj.pTimestamps{idx},startTime);

            % Calculate the range in with the IQ data and timestamps needs to be stored
            rangeIQ = obj.pSize(idx)+isOverlap:obj.pSize(idx)+size(timestamps,1)-~isOverlap;

            obj.pData{idx}(rangeIQ,:) = waveform;
            obj.pTimestamps{idx}(rangeIQ,:) = timestamps;
            obj.pSize(idx) = rangeIQ(end);
        end

        function isOverlap = checkOverlap(~,size,timestamp,startTime)
            %checkOverlap Check if there is an overlap of the previous last timestamp
            %with the first current timestamp.

            % Since previous end time and the current start time can be same, the start
            % index for the range of the waveform must be calculated accordingly.

            if size>0 && timestamp(size,1)==startTime
                isOverlap = 0;
            else
                isOverlap = 1;
            end
        end

        function saveIQSamples(obj)
            %saveIQSamples Save the captured IQ samples as a MAT file

            iqSamples = obj.IQSamples;
            save(strjoin([obj.FileName ".mat"],""),"iqSamples","-v7.3","-nocompression");
        end

        function validateWLANNode(obj,node)
            %validateWLANNode Validate the WLAN node

            devCfg = getWLANDeviceConfig(obj, node);
            validateattributes(obj.CenterFrequency,"double",{'vector','numel',numel(devCfg)},mfilename,"CenterFrequency");
            validateattributes(obj.Bandwidth,"double",{'vector','numel',numel(devCfg)},mfilename,"Bandwidth");
            if node.MACFrameAbstraction==1 || ~strcmpi(node.PHYAbstractionMethod,"none")
                error("IQ sample capture is not supported with abstraction.");
            end
            % If the number of center frequency is more than 1, order it as same order
            % of PHY Rx in WLAN node
            if numel(obj.CenterFrequency)>1
                centerFreqBand = round(obj.CenterFrequency/1e6);
                centerFreqBandNode = round(node.ReceiveFrequency/1e6);
                freqOrder = zeros(numel(centerFreqBandNode),1);
                for idx=1:numel(centerFreqBandNode)
                    freqIdx = find(centerFreqBand==centerFreqBandNode(idx));
                    if numel(freqIdx)==0
                        error(strjoin(["The specified center frequency" string(obj.CenterFrequency(idx)) "does not correspond to any device in the node."]));
                    elseif numel(freqIdx)>1
                        error("More than 1 specified center frequencies are in the same frequency band.");
                    end
                    freqOrder(idx) = freqIdx;
                end
                obj.CenterFrequency = obj.CenterFrequency(freqOrder);
                obj.Bandwidth = obj.Bandwidth(freqOrder);
            end
            obj.pExpectedSampleRate = obj.Bandwidth;
        end

        function validateBluetoothNode(obj)
            %validateBluetoothNode Validate the Bluetooth node

            % For Bluetooth BR/EDR and LE node Center frequency and bandwidth must be
            % scalar
            validateattributes(obj.CenterFrequency,"double",{'scalar'},mfilename,"CenterFrequency");
            validateattributes(obj.Bandwidth,"double",{'scalar'},mfilename,"Bandwidth");
        end

        function validate5GNode(obj)
            %validate5GNode Validate the 5G node

            % Center frequency and bandwidth must be scalar
            validateattributes(obj.CenterFrequency,"double",{'scalar'},mfilename,"CenterFrequency");
            validateattributes(obj.Bandwidth,"double",{'scalar'},mfilename,"Bandwidth");
        end

        function deviceConfig = getWLANDeviceConfig(~, node)
            %getWLANDeviceConfig Returns the object holding WLAN MAC/PHY configuration

            devCfg = node.DeviceConfig;
            if isa(devCfg, "wlanDeviceConfig")
                deviceConfig = devCfg;
            else % wlanMultilinkDeviceConfig
                deviceConfig = devCfg.LinkConfig;
            end
        end
    end
end