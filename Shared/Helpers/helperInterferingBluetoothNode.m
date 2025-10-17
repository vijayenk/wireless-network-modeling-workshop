classdef helperInterferingBluetoothNode < wirelessnetwork.internal.wirelessNode
    %helperInterferingBluetoothNode Create an object for a Bluetooth LE
    %node that generates Bluetooth LE signals that act as interference
    %   INTLENODE = helperInterferingBluetoothNode creates an object for
    %   Bluetooth LE node that generates a signal periodically to introduce
    %   Bluetooth LE interference. The node does not do protocol modeling
    %   and only transmits PHY signals (IQ samples) based on the
    %   configuration.
    %
    %   INTLENODE = helperInterferingWLANNode(Name=Value) creates an object
    %   for Bluetooth LE node with the specified property Name set to the
    %   specified Value. You can specify additional name-value pair
    %   arguments in any order as (Name1=Value1, ..., NameN=ValueN).
    %
    %   helperInterferingWLANNode properties:
    %
    %   ID                - Node identifier
    %   Position          - Node position
    %   Name              - Node namee
    %   SignalPeriodicity - Periodicity of signals in seconds
    %   TransmitterPower  - Signal transmission power in dBm
    %   PHYMode           - Physical layer (PHY) mode
    %   AccessAddress     - Unique connection address of the nod
    %   DataLength        - Length of data to be transmitted in bytes
    %   UsedChannels      - List of used (good) data channels
    %   Algorithm         - Type of channel selection algorithm

    %   Copyright 2025 The MathWorks, Inc.

    properties
        %SignalPeriodicity Periodicity of signals in seconds
        %   Specify the signal periodicity as a positive double value in
        %   seconds. This property specifies the interval between the start
        %   of two successive signals. The default value is 2e-3 seconds.
        SignalPeriodicity (1,1) {mustBePositive} = 2e-3

        %TransmitterPower Signal transmission power in dBm
        %   Specify the transmit power as a scalar double value. It
        %   specifies the signal transmission power in dBm. The default
        %   value is 20 dBm.
        TransmitterPower (1,1) {mustBeNumeric, mustBeFinite, ...
            mustBeLessThanOrEqual(TransmitterPower,100)} = 20

        %PHYMode Physical layer (PHY) mode
        %   Specify the PHY mode as "LE1M","LE2M","LE500K" or "LE125K".
        %   This property indicates the physical layer mode used for
        %   generating the signals. The default value is "LE1M".
        PHYMode (1,1) {mustBeMember(PHYMode,["LE1M";"LE2M";"LE125K";"LE500K"])} = "LE1M"

        %AccessAddress Unique connection address of the node
        %   Specify the access address as an 8-element character vector or
        %   a string scalar denoting a 4-octet hexadecimal value. This
        %   property specifies a unique 32-bit address for the link layer
        %   connection which is used for generating a pcket. The default
        %   value is "5DA44270". You can easily create a valid access
        %   address using ble.internal.generateAccessAddress.
        AccessAddress = "5DA44270"

        %DataLength Length of data to be transmitted in bytes
        %   Specify the length of data to be transmitted in bytes as a
        %   scalar integer value less than or equal to 251. It specifies
        %   the transmitted packet data length in bytes. This will be used
        %   to form the link layer (LL) packet data unit (PDU) which will
        %   then be used to generate the waveform. The default is 100
        %   bytes.
        DataLength (1,1) {mustBePositive,mustBeInteger,mustBeLessThanOrEqual(DataLength,251)} = 100

        %UsedChannels List of used (good) data channels
        %   Specify the list of used data channels as an integer vector with
        %   element values in the range [0, 36]. This value specifies the
        %   indices of the assigned data channels. To ensure that at least two
        %   channels are set as used (good) channels, specify the vector length
        %   greater than 1. This property indicates the set of good channels
        %   used by the connection. The default value is a row vector
        %   containing all the channel indices [0:36].
        UsedChannels = 0:36

        %Algorithm Type of channel selection algorithm
        %   Specify the channel selection algorithm as 1 or 2 representing
        %   "Algorithm #1" or "Algorithm #2", respectively. This property
        %   specifies the channel selection algorithm used for hopping of this
        %   connection. The default value is 1.
        Algorithm (1,1) {mustBeMember(Algorithm,[1;2])} = 1
    end

    properties (SetAccess = private)
        %LastCurrentTime Current simulation time
        LastCurrentTime = 0;

        %TransmitBuffer Buffer contains the data to be transmitted from the
        %node
        TransmitBuffer

        %TransmittedSignals Total number of transmitted signals
        TransmittedSignals = 0

        %Waveform IQ samples of the generated waveform
        Waveform
    end

    properties (Constant,Hidden)
        %Bandwidth Bandwidth of the Bluetooth LE signal
        Bandwidth = 2e6
    end

    properties (Access = private)
        %NextSignalTimer Timer to track the start of next signal
        pNextSignalTimer = 0

        %pTransmitBuffer Buffer contains the data to be transmitted from the
        %node
        pTransmitBuffer

        %pMetadata Structure containing the metadata of the signal
        pMetadata

        %pIsInitialized Flag to check whether the node is initialized
        pIsInitialized = false

        %SampleRate Sample rate of the signal in Hz
        SampleRate = 20e6

        %AccessAddressBin Access address in binary
        pAccessAddressBin

        %pChannelSelection Channel selection object for Bluetooth LE
        pChannelSelection = bleChannelSelection

        %pCFList List of LE center frequencies
        pCFList
    end

    events (Hidden)
        %PacketTransmissionStarted is triggered when the node starts
        %transmitting a packet. PacketTransmissionStarted passes the event
        %notification along with a structure to the registered callback. The
        %structure has fields in pTransmitBuffer property along woth PacketDuration
        %and TransmittedPower property
        PacketTransmissionStarted
    end

    methods
        % Constructor
        function obj = helperInterferingBluetoothNode(varargin)
            % Name-value pairs
            for idx = 1:2:nargin-1
                obj.(varargin{idx}) = varargin{idx+1};
            end

            % Initialize the WLAN signal
            obj.pMetadata = struct('NumSamples', 0, ...
                'Duration', 0);

            % Initialize the transmission buffer with empty signal
            obj.pTransmitBuffer = wirelessnetwork.internal.wirelessPacket;

            % Initialize the list of center frequencies
            obj.pCFList = ble.internal.networkUtilities.ChannelFrequencies;

            % Initialize the channel selection algorithm
            obj.pChannelSelection = bleChannelSelection( ...
                        "Algorithm",obj.Algorithm,...
                        "AccessAddress",obj.AccessAddress,...
                        "UsedChannels",obj.UsedChannels);
        end

        % Auto-completion
        function v = set(obj, prop)
            v = obj.([prop, 'Values']);
        end

        % Set used channels
        function set.UsedChannels(obj, value)
            validateattributes(value, {'numeric'}, {'vector', 'nonnegative', ...
                'integer', '<=', 36}, obj.FileName, "UsedChannels");
            value = unique(value);
            % At least two channels must be set as used (good) channels as
            % specified in Bluetooth Core Specification v5.3 | Vol 6, Part B,
            % Section 4.5.8.1
            if numel(value) < 2
                error(message("bluetooth:bleLL:InvalidUsedChannels"));
            end
            obj.UsedChannels = value;
        end

        % Set access address
        function set.AccessAddress(obj, value)
            ble.internal.validateHex(value, 8, "AccessAddress");
            obj.AccessAddress = string(value);
        end

        function nextInvokeTime = run(obj, currentTime)
            %run Runs the WLAN node
            %
            %   NEXTINVOKETIME = run(OBJ, CURRENTTIME) runs the
            %   functionality of WLAN node at the current time,
            %   CURRENTTIME, and returns the time, NEXTINVOKETIME, to run
            %   the node again.
            %
            %   NEXTINVOKETIME is the time instant (in seconds) at which
            %   the node runs again.
            %
            %   OBJ is an object of type helperInterferingWLANNode.
            %
            %   CURRENTTIME is the current simulation time in seconds.

            % Initialize the node when the node is run for the first time
            if ~obj.pIsInitialized
                init(obj);
            end

            % Update the simulation time (in seconds)
            elapsedTime = currentTime - obj.LastCurrentTime;
            obj.LastCurrentTime = currentTime;

            % Update the timer
            obj.pNextSignalTimer = round(obj.pNextSignalTimer - elapsedTime, 9);

            % Reset the transmission buffer
            obj.TransmitBuffer = [];

            % Update the transmission buffer with the waveform and its
            % metadata
            if obj.pNextSignalTimer <= 0
                % Update the transmission buffer with the signal
                % information
                obj.TransmitBuffer = obj.pTransmitBuffer;
                obj.TransmitBuffer.StartTime = obj.LastCurrentTime;

                % Choose the channel index and generate the waveform
                channelIndex = obj.pChannelSelection();
                obj.TransmitBuffer.CenterFrequency = obj.pCFList(channelIndex+1);
                generateWaveform(obj,channelIndex);
    
                % Apply Tx power on the waveform
                scale = 10.^((-30 + obj.TransmitterPower)/20);
                obj.Waveform = obj.Waveform * scale;

                % Add the waveform information
                obj.TransmitBuffer.Metadata.NumSamples = numel(obj.Waveform);
                obj.TransmitBuffer.Data = obj.Waveform;
                obj.TransmitBuffer.Duration = round(obj.TransmitBuffer.Metadata.NumSamples/obj.SampleRate, 9);

                % Validate the waveform duration
                if obj.pTransmitBuffer.Duration > obj.SignalPeriodicity
                    error('Transmission duration of signal must be less than signal periodicity');
                end

                % Update the position and velocity of the node
                obj.TransmitBuffer.TransmitterPosition = obj.Position;
                obj.TransmitBuffer.TransmitterVelocity = obj.Velocity;

                % Update the timer for next signal transmission
                obj.pNextSignalTimer = obj.SignalPeriodicity;

                % Update the number of signals transmitted
                obj.TransmittedSignals = obj.TransmittedSignals + 1;

                % Notify the packet transmission started event
                txBuffer = obj.TransmitBuffer;
                txBuffer.PacketDuration = obj.pTransmitBuffer.Duration;
                txBuffer.TransmittedPower = obj.pTransmitBuffer.Power;
                triggerEvent(obj,"PacketTransmissionStarted",txBuffer);
            end
            nextInvokeTime = round(obj.pNextSignalTimer + obj.LastCurrentTime, 9);
        end

        function txPackets = pullTransmittedData(obj)
            txPackets = obj.TransmitBuffer;
            obj.TransmitBuffer = [];
        end

        function pushReceivedData(~, ~)
            % Do nothing
        end

        function [flag, rxInfo] = isPacketRelevant(~, ~)
            %isPacketRelevant Return flag to indicate whether channel
            %has to be applied on incoming signal

            %The WLAN node is only for transmission and not for reception
            %hence the flag is set as false and rxInfo is set empty
            flag = false;
            rxInfo = [];
        end

        function updateChannelMap(obj,newUsedChannelList)
            %updateChannelList Provide updated channel list to Bluetooth LE node
            %
            %   STATUS = updateChannelList(OBJ, NEWUSEDCHANNELSLIST) updates
            %   the channel map by providing a new list of used channels,
            %   NEWUSEDCHANNELSLIST, to the channel selection algorithm.
            %
            %   OBJ is an object of type helperInterferingBluetoothLENode.
            %
            %   NEWUSEDCHANNELSLIST is the list of good (used) channels,
            %   specified as an integer vector with element values in the range
            %   [0, 36].

            obj.UsedChannels = newUsedChannelList;

            obj.pChannelSelection.UsedChannels = unique(newUsedChannelList);
        end
    end

    methods (Access = private)
        function init(obj)
            %init Initializes the waveform and transmit buffer

            % Convert access address to binary
            obj.pAccessAddressBin = int2bit(hex2dec(obj.AccessAddress),32,0);

            % Update the transmission buffer with the signal information
            obj.pTransmitBuffer.Type = 3;
            obj.pTransmitBuffer.SampleRate = obj.SampleRate;
            obj.pTransmitBuffer.NumTransmitAntennas = 1;
            obj.pTransmitBuffer.Bandwidth = 2e6;
            obj.pTransmitBuffer.TransmitterID = obj.ID;
            obj.pTransmitBuffer.Power = obj.TransmitterPower;

            % Set initialization flag
            obj.pIsInitialized = true;
        end

        function generateWaveform(obj,channelIndex)
            % generateWaveform generates the Bluetooth LE waveform using the
            % features of Bluetooth Toolbox(TM)

            % Create random data
            data = randi([0 1], obj.DataLength, 1);

            % Generate the LL PDU
            llPDU = bleLLDataChannelPDU(bleLLDataChannelPDUConfig(LLID="Data (start fragment/complete)"),data);

            % Generate Bluetooth LE waveform
            bleWaveform = bleWaveformGenerator(single(llPDU),"Mode",obj.PHYMode,...
                "ChannelIndex",channelIndex,"SamplesPerSymbol",...
                getSamplesPerSymbol(obj),"AccessAddress",obj.pAccessAddressBin);
            obj.Waveform = bleWaveform;
        end

        function triggerEvent(obj, eventName, eventData)
            %triggerEvent Trigger the event to notify all the listeners

            if event.hasListener(obj, eventName)
                eventData.NodeName = obj.Name;
                eventData.NodeID = obj.ID;
                eventData.CurrentTime = obj.LastCurrentTime;
                eventDataObj = wirelessnetwork.internal.nodeEventData;
                eventDataObj.Data = eventData;
                notify(obj, eventName, eventDataObj);
            end
        end

        % Get samples per symbol
        function samplesPerSymbol = getSamplesPerSymbol(obj)
            %getSamplesPerSymbol Return the samples per symbol based on the PHY
            %mode and sample rate
    
            % PHY mode
            switch obj.PHYMode
                case {"LE1M","LE500K","LE125K"}
                    % 1Msps (one mega symbols per second)
                    symbolRate = 1e6;
                case "LE2M"
                    % 2Msps (two mega symbols per second)
                    symbolRate = 2e6;
            end
            % Multiply with sps (samples per symbol)
            samplesPerSymbol = obj.SampleRate/symbolRate;
        end
    end
end