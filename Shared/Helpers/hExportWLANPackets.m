classdef hExportWLANPackets < handle
%hExportWLANPackets Export transmitted and received MAC frames to PCAP or
%PCAPNG file
%
%   OBJ = hExportWLANPackets(WLANNODES) exports the transmitted and
%   received MAC frames to a PCAP file format, by using the default value
%   of FILEEXTENSION.
%
%   OBJ = hExportWLANPackets(WLANNODES, FILEEXTENSION) exports the
%   transmitted and received MAC frames, captured at WLANNODES, to a file
%   with PCAP (Packet Capture) or PCAPNG (Packet Capture Next Generation)
%   extension, FILEEXTENSION. The node name, node ID, and timestamp of
%   creation are included in the file name. OBJ captures the MAC frames by
%   using the "MPDUGenerated" and "MPDUDecoded" events of the objects in
%   WLANNODES. When the events are triggered, this function invokes the
%   callback function, packetWriterCallback, that writes MAC frames to a
%   PCAP or PCAPNG file.
%
%   OBJ = hExportWLANPackets(WLANNODES, FILEEXTENSION, ENABLERADIOTAP)
%   additionally adds the radiotap information in the captured MAC frames,
%   and export it to the PCAP or PCAPNG file format, if ENABLERADIOTAP is
%   true. The second and third arguments of this syntax can be interchanged
%   and specified independently of each other.
%
%   WLANNODES is an array or cell array of objects of type <a
%   href="matlab:help('wlanNode')">wlanNode</a>.
%
%   FILEEXTENSION is a file extension of a PCAP or PCAPNG file that stores
%   the MAC frames. Specify this value as "pcap" or "pcapng". The default
%   value is "pcap".
%
%   ENABLERADIOTAP is a logical scalar that enables or disables adding
%   radiotap information to the captured MAC frames. The default value is
%   "false". These are the radiotap fields formed for the transmitted
%   frames: Channel, Flags, dBm TX power. These are the radiotap fields
%   formed for the received frames: Channel, Flags, dBm TX power, HE, MCS,
%   Rate, and VHT.
%
%   hExportWLANPackets properties (read-only):
%
%   PCAPObjList             - Array of objects of type wlanPCAPWriter
%   WLANNodes               - Configured WLAN Nodes in the network

%   Copyright 2022-2024 The MathWorks, Inc.

    properties(GetAccess = public, SetAccess = private)
        %PCAPObjList Array of objects of type wlanPCAPWriter
        %   This property is an array of objects of type <a
        %   href="matlab:help('wlanPCAPWriter')">wlanPCAPWriter</a>
        PCAPObjList;

        %WLANNodes Configured WLAN Nodes in the network
        %   Specify the WLAN nodes configured in the network as an array of
        %   objects of type <a href="matlab:help('wlanNode')">wlanNode</a>.
        WLANNodes;

        %RadiotapFields Structure of radiotap fields
        %   This property stores the logical values to enable different
        %   radiotap fields.
        RadiotapFields;

        %EnableRadiotap Enable radiotap information
        %   This property is a logical scalar. If set to true, radiotap
        %   information is written into the PCAP/PCAPNG file for the
        %   captured MAC frames. The default is false.
        EnableRadiotap = false;
    end

    properties(Constant, Hidden)
        % Frame configuration objects used for calculating data rate
        HTConfigObject = wlanHTConfig;
        VHTConfigObject = wlanVHTConfig;
        HESUConfigObject = wlanHESUConfig;
    end

    methods
        function obj = hExportWLANPackets(wlanNodes, varargin)

            % Check for number of input arguments
            narginchk(1,3);

            % Set default values for non-mandatory inputs
            fileExtension = "pcap";
            enableRadiotap = false;

            % Check for non-mandatory inputs
            switch nargin
                % Two inputs provided:
                % wlanNodes and fileExtension (OR)
                % wlanNodes and enableRadiotap
                case 2
                    if isstring(varargin{1}) || ischar(varargin{1})
                        fileExtension = setFileExtension(obj,varargin{1});
                    else
                        validateattributes(varargin{1}, "logical", {'scalar'});
                        enableRadiotap = varargin{1};
                    end
                    % Three inputs provided:
                    % wlanNodes, fileExtension, enableRadiotap (OR)
                    % wlanNodes, enableRadiotap, fileExtension
                case 3
                    if isstring(varargin{1}) || ischar(varargin{1})
                        fileExtension = setFileExtension(obj,varargin{1});
                        validateattributes(varargin{2}, "logical", {'scalar'});
                        enableRadiotap = varargin{2};
                    else
                        validateattributes(varargin{1}, "logical", {'scalar'});
                        enableRadiotap = varargin{1};
                        fileExtension = setFileExtension(obj,varargin{2});
                    end
            end

            % Validate WLAN node inputs
            if iscell(wlanNodes)
                for idx = 1:numel(wlanNodes)
                    validateattributes(wlanNodes{idx}, "wlanNode", {'scalar'});
                end
            else
                validateattributes(wlanNodes(1), "wlanNode", {'scalar'});
                wlanNodes = num2cell(wlanNodes);
            end

            obj.WLANNodes = wlanNodes;
            obj.EnableRadiotap = enableRadiotap;
            obj.PCAPObjList = wlanPCAPWriter.empty(0,numel(wlanNodes));
            fileName = strings(1,numel(wlanNodes));
            for nodeIDx = 1:numel(wlanNodes)
                if wlanNodes{nodeIDx}.MACFrameAbstraction
                    warning("hExportWLANPackets:AbstractedMACNotValid","MACFrameAbstraction must be set to false for capturing MAC frames.");
                else
                    wlanNodes{nodeIDx}.IncludeRxVector = true;

                    % PCAP file name is by the format:
                    % <NodeName_NodeID_yyyyMMdd_HHmmss>
                    fileName(nodeIDx) = strjoin([wlanNodes{nodeIDx}.Name "_" wlanNodes{nodeIDx}.ID "_" char(datetime('now','Format','yyyyMMdd_HHmmss'))],"");

                    % Create a WLAN PCAP file writer object with the specified
                    % file name and extension by using the wlanPCAPWriter
                    % object.
                    obj.PCAPObjList(nodeIDx) = wlanPCAPWriter('FileName', fileName(nodeIDx), 'FileExtension',fileExtension, 'RadiotapPresent', obj.EnableRadiotap);
                    addlistener(wlanNodes{nodeIDx}, 'MPDUDecoded', ...
                        @(src, eventData)packetWriterCallback(obj,obj.PCAPObjList(nodeIDx), eventData));
                    addlistener(wlanNodes{nodeIDx}, 'MPDUGenerated', ...
                        @(src, eventData)packetWriterCallback(obj,obj.PCAPObjList(nodeIDx), eventData));
                end
            end
            resetRadiotapFields(obj);
        end
    end

    methods(Access = private)
        function packetWriterCallback(obj, pcapObj, eventData)
            %packetWriterCallback Callback function that writes the packets
            %into PCAP/PCAPNG format
            %
            %   packetWriterCallback(PCAPOBJ, NODEOBJ, EVENTDATA) is a
            %   callback function that writes the packets into PCAP/PCAPNG
            %   format when 'MPDUDecoded' or 'MPDUGenerated' is triggered
            %   from node. It also forms the radiotap bytes for the frames
            %   received at each node.
            %
            %   PCAPOBJ is an object of type wlanPCAPWriter.
            %
            %   EVENTDATA is a structure containing the following fields:
            %       Data      - Structure containing the following fields
            %                   when EventName is MPDUGenerated:
            %                   DeviceID - Scalar representing device
            %                              identifier
            %                   CurrentTime - Scalar representing current
            %                   simulation time in seconds
            %                   MPDU - For full MAC, it is a cell array of
            %                          MPDU(s) where each element is a
            %                          vector containing MPDU bytes in
            %                          decimal format. In case of
            %                          abstracted MAC frames, it is a
            %                          structure containing MAC frame
            %                          information
            %                   Frequency - Scalar representing receiving
            %                   frequency in Hz If EventName is
            %                   MPDUDecoded, the following two additional
            %                   fields are present:
            %                   FCSFail - Flag indicating whether frame
            %                             check sequence(FCS) failed. In
            %                             case of multiple MPDUs, it is a
            %                             vector with values for each MPDU
            %                   RxVector - Structure containing parameters
            %                              provided by PHY upon receipt of
            %                              a valid PHY header
            %
            %       Source    - An object of type wlanNode representing the
            %                   node from which the event is triggered
            %       EventName - Name of the event

            if strcmp(eventData.EventName, "MPDUDecoded") || strcmp(eventData.EventName, "MPDUGenerated")
                macFrame = eventData.Data.MPDU;
                if strcmp(eventData.EventName, "MPDUGenerated")
                    macFrame = macFrame';
                end

                % Calculate the timestamp of the packet in microseconds
                timestamp = round(eventData.Data.CurrentTime*1e6);
                if obj.EnableRadiotap
                    radiotapBytes = formRadiotap(obj,eventData);
                    for frameIDx = 1:numel(macFrame)
                        % For received frames, radiotapBytes is 2D array
                        % containing the radiotap bytes for each subframe.
                        % The radiotap flags byte (has FCS bit) must be
                        % filled accordingly
                        if strcmp(eventData.EventName, "MPDUDecoded")
                            write(pcapObj, macFrame{frameIDx},timestamp,'Radiotap',radiotapBytes(frameIDx,:));
                        else % FCS information is not available for transmitted frames
                            write(pcapObj, macFrame{frameIDx},timestamp,'Radiotap',radiotapBytes);
                        end
                    end
                else
                    for frameIDx = 1:numel(macFrame)
                        write(pcapObj, macFrame{frameIDx},timestamp);
                    end
                end
            end
        end

        function radiotapBytes = formRadiotap(obj, eventData, varargin)
            %formRadiotap Returns the radiotap bytes
            %
            %   formRadiotap(EVENTDATA) is a function that forms and
            %   returns the radiotap bytes.
            %
            %   formRadiotap(EVENTDATA, CAPTURETIME) additionally fills
            %   the timestamp bytes in radiotap.
            %
            %   EVENTDATA is a structure containing event data, source
            %   node object, and event name.
            %
            %   CAPTURETIME is the current time captured from event data.
            %
            %   The following radiotap fields are formed using context
            %   captured from EVENTDATA.
            %   For received frames (captured via MPDUDecoded event):
            %   Channel, Flags, HE, MCS, Rate, VHT, dBm TX power,
            %   timestamp, and LSIG (disabled by default).
            %   For transmitted frames (captured via MPDUGenerated event):
            %   Channel, Flags, dBm TX power, timestamp

            headerRevision = 0;
            headerPad = 0;
            radiotapBytes = [headerRevision headerPad 0 0]; % header length is 0 initially

            formRadiotapBytes = (strcmp(eventData.EventName, "MPDUDecoded") && ~eventData.Data.PHYDecodeFail) || ...
                strcmp(eventData.EventName, "MPDUGenerated");

            if formRadiotapBytes
                % Reset radiotap flags to default values for each captured
                % frame
                resetRadiotapFields(obj);

                % Include frame timestamp if capture time is provided as an
                % input
                if ~isempty(varargin)
                    obj.RadiotapFields.FrameTimestamp = 1;
                    captureTime = varargin{1};
                end

                if strcmp(eventData.EventName, "MPDUDecoded")
                    packetInfo = eventData.Data.RxVector;
                    switch packetInfo.PPDUFormat
                        case 1 % NonHT, MCS information is not captured for NonHT
                            obj.RadiotapFields.MCSInformation = 0;
                        case 3 % VHT
                            obj.RadiotapFields.VHTInformation = 1;
                        case {4, 5} % HE-SU or HE-EXT-SU
                            obj.RadiotapFields.HEInformation = 1;
                        case 7 % HE-TB
                            obj.RadiotapFields.LSIG = 1;
                    end

                    % Do not fill rate byte if it exceeds 255
                    dataRate = formRadiotapRateByte(obj,packetInfo);
                    dataRate = round(dataRate);
                    if dataRate > 255
                        obj.RadiotapFields.Rate = 0;
                    end
                else % MPDUGenerated
                    % Rate and MCS fields are not formed for transmitted
                    % frames, so disable the respective radiotap fields
                    obj.RadiotapFields.Rate = 0;
                    obj.RadiotapFields.MCSInformation = 0;
                end

                % Get the logical values as an array
                fieldVector = cellfun(@(x)(obj.RadiotapFields.(x)),fieldnames(obj.RadiotapFields));
                % Add zeros to fill up the unsupported bits
                fieldVector = [fieldVector(1:16);0;fieldVector(17:24); 0;fieldVector(25:30)].';
                % Create the hexadecimal word for radiotap flags
                radiotapFlagsWord = dec2hex(flip((bit2int(fieldVector',4,false))'))';
                radiotapFlagsWord = reshape(radiotapFlagsWord, 2, []).';
                % Create the byte representation to pass in the radiotap vector
                radiotapFlagBytes = flip(hex2dec(radiotapFlagsWord).');

                % Append to radiotap
                radiotapBytes = [radiotapBytes radiotapFlagBytes];

                if obj.RadiotapFields.Flags
                    % Form Flags byte, only 'FCS at end' is true by default,
                    % decimal 16;
                    radiotapBytes = [radiotapBytes 16];
                    if strcmp(eventData.EventName, "MPDUDecoded")
                        flagsIndex = numel(radiotapBytes);
                    end
                end

                if obj.RadiotapFields.Rate
                    % Fill previously computed Rate byte
                    radiotapBytes = [radiotapBytes dataRate];
                end

                if obj.RadiotapFields.Channel
                    % Check for 2 bytes alignment
                    aligned= mod(numel(radiotapBytes), 2);
                    if aligned~=0
                        numZeroPad = calculateNumberOfZeroPads(obj,radiotapBytes,2);
                        radiotapBytes = [radiotapBytes zeros(1,numZeroPad)];
                    end
                    channelBytes = formRadiotapChannelBytes(obj,eventData);
                    radiotapBytes = [radiotapBytes channelBytes];
                end

                if obj.RadiotapFields.dBmTxPower
                    % dBm Tx Power
                    if strcmp(eventData.EventName, "MPDUDecoded")
                        txPower = round(packetInfo.PerUserInfo.TxPower,0);
                    else
                        devCfg = getDeviceConfig(obj, eventData);
                        txPower = round(devCfg.TransmitPower,0);
                    end
                    radiotapBytes = [radiotapBytes txPower];
                end

                if obj.RadiotapFields.MCSInformation
                    % Form MCS information bytes
                    mcsBytes = formRadiotapMCSBytes(obj,packetInfo);
                    radiotapBytes = [radiotapBytes mcsBytes];
                end

                if obj.RadiotapFields.VHTInformation
                    % Check for 2 bytes alignment
                    aligned= mod(numel(radiotapBytes), 2);
                    if aligned~=0
                        numZeroPad = calculateNumberOfZeroPads(obj,radiotapBytes,2);
                        radiotapBytes = [radiotapBytes zeros(1,numZeroPad)];
                    end

                    % Form VHT information bytes
                    VHTInformationBytes = formRadiotapVHTBytes(obj,packetInfo);
                    radiotapBytes = [radiotapBytes VHTInformationBytes];
                end

                if obj.RadiotapFields.FrameTimestamp
                    % Check for 8 bytes alignment
                    aligned = mod(numel(radiotapBytes), 8);
                    if aligned~=0
                        numZeroPad = calculateNumberOfZeroPads(obj,radiotapBytes,8);
                        radiotapBytes = [radiotapBytes zeros(1,numZeroPad)];
                    end

                    % Form timestamp bytes
                    timestampBytes = formRadiotapTimestampBytes(obj,captureTime);
                    % Append to radiotap
                    radiotapBytes = [radiotapBytes timestampBytes];
                end

                if obj.RadiotapFields.HEInformation
                    % Check for 2 bytes alignment
                    aligned = mod(numel(radiotapBytes), 2);
                    if aligned~=0
                        numZeroPad = calculateNumberOfZeroPads(obj,radiotapBytes,2);
                        radiotapBytes = [radiotapBytes zeros(1,numZeroPad)];
                    end

                    % Form HE radiotap bytes
                    HEInformationBytes = formRadiotapHEBytes(obj,packetInfo);
                    radiotapBytes = [radiotapBytes HEInformationBytes];
                end

                if obj.RadiotapFields.LSIG
                    % Check for 2 bytes alignment
                    aligned = mod(numel(radiotapBytes), 2);
                    if aligned~=0
                        numZeroPad = calculateNumberOfZeroPads(obj,radiotapBytes,2);
                        radiotapBytes = [radiotapBytes zeros(1,numZeroPad)];
                    end

                    % Form LSIG bytes
                    LSIGBytes = formRadiotapLSIGBytes(obj,packetInfo);
                    radiotapBytes = [radiotapBytes LSIGBytes];
                end

                % Final header length
                if numel(radiotapBytes)>255
                    headerLength = [255 mod(numel(radiotapBytes), 255)];
                else
                    headerLength = [numel(radiotapBytes) 0];
                end

                % Header Pad value
                radiotapBytes(1,3:4) = headerLength;

                % Fill FCS at end flag for each subframe
                if strcmp(eventData.EventName, "MPDUDecoded")
                    numSubframes = numel(eventData.Data.FCSFail);
                    fcsBytes = eventData.Data.FCSFail;
                    outputRadiotap = repmat(radiotapBytes, numSubframes, 1);
                    for i = 1:numSubframes
                        if fcsBytes(i)
                            outputRadiotap(i,flagsIndex) = 80;
                        end
                    end
                    radiotapBytes =  outputRadiotap;
                end
            end
        end

        function numZeroPad = calculateNumberOfZeroPads(~,radiotapBytes, byteAlignment)
            %calculateNumberOfZeroPads Returns the number of zeros to be
            %padded to radiotap bytes for a given byte alignment

            counter = 1;
            while(numel(radiotapBytes)>(byteAlignment*counter))
                counter = counter+1;
            end
            numZeroPad = mod(byteAlignment*counter,numel(radiotapBytes));
        end

        function resetRadiotapFields(obj)
            %resetRadiotapFields Reset the radiotap flags to default logical
            %values

            obj.RadiotapFields = struct('TSFT', 0, ... % 802.11 Time Synchronization Function timer
                'Flags', 1, ... % Properties of transmitted and received frames
                'Rate', 1, ... % TX/RX data rate
                'Channel', 1, ... % Tx/Rx frequency in MHz, followed by flags
                'FHSS', 0, ... % The hop set and pattern for frequency-hopping radios
                'dBmAntennaSignal', 0, ... % RF signal power at the antenna(decibels from a 1 milliwatt reference)
                'dBmAntennaNoise', 0, ... % RF noise power at the antenna(decibels from a 1 milliwatt reference)
                'LockQuality', 0, ... % Quality of Barker code lock, unitless
                'TxAttenuation', 0, ... % Transmit power expressed as unitless distance from max power set at factory calibration
                'dBTxAttenuation', 0, ... % Transmit power expressed as decibel distance from max power set at factory calibration
                'dBmTxPower', 1, ... % Transmit power expressed as dBm (decibels from a 1 milliwatt reference)
                'Antenna', 0, ... % Unitless indication of the Rx/Tx antenna for this packet
                'dBAntennaSignal', 0, ... % RF signal power at the antenna, decibel difference from an arbitrary, fixed reference
                'dBAntennaNoise', 0, ... % RF noise power at the antenna, decibel difference from an arbitrary, fixed reference
                'RxFlags', 0, ... % Properties of received frames
                'TxFlags', 0, ... % Properties of transmitted frames
                'DataRetries', 0, ... % Reserved field
                'ChannelPlus', 0, ... % Reserved field
                'MCSInformation', 1, ... % MCS rate index as in IEEE_802.11n-2009
                'AMPDUStatus', 0, ... % Frame was received as part of an a-MPDU
                'VHTInformation', 0, ... % Information for VHT frames
                'FrameTimestamp', 0, ... % Frame timestamp
                'HEInformation', 0, ... % Frame was received or transmitted using the HE PHY
                'HEMUInformation', 0, ... % PPDUs of HE_MU type that wasn’t already captured in the regular HE field.
                'ZeroLengthPSDU', 0, ... % no PSDU captured for this PPDU
                'LSIG', 0, ... % L-SIG contents, if known
                'TLVs', 0, ... % type-length-value item list
                'RadiotapNSnext', 0, ... % Reserved field
                'VendorNSnext', 0, ... % Reserved field
                'Ext', 0); % Reserved field
        end

        function dataRate = formRadiotapRateByte(obj,packetInfo)
            %formRadiotapRateByte Returns the data rate for captured frame
            %
            % PACKETINFO is a structure containing parameters provided by PHY 
            % upon receipt of a valid PHY header.

            mcs = packetInfo.PerUserInfo.MCS;
            numSTS = packetInfo.PerUserInfo.NumSpaceTimeStreams;
            cbw = strcat('CBW', int2str(packetInfo.ChannelBandwidth));
            switch packetInfo.PPDUFormat
                case 1 % NonHT
                    switch mcs
                        case 0 % 6 Mbps
                            dataRate = 6;
                        case 1 % 9 Mbps
                            dataRate = 9;
                        case 2 % 12 Mbps
                            dataRate = 12;
                        case 3 % 18 Mbps
                            dataRate = 18;
                        case 4 % 24 Mbps
                            dataRate = 24;
                        case 5 % 36 Mbps
                            dataRate = 36;
                        case 6  % 48 Mbps
                            dataRate = 48;
                        otherwise % 54 Mbps
                            dataRate = 54;
                    end
                    dataRate = dataRate*1e3/500; % in units of 500 Kbps

                case 2 % HT-mixed
                    htConfig = obj.HTConfigObject;
                    htConfig.ChannelBandwidth = cbw;
                    htConfig.AggregatedMPDU = packetInfo.AggregatedMPDU;
                    htConfig.MCS = mcs;
                    htConfig.NumTransmitAntennas = numSTS; % Set as NumSTS for validation
                    htConfig.NumSpaceTimeStreams = numSTS;

                    r = wlan.internal.getRateTable(htConfig);
                    symbolTime = 4; % in microseconds
                    dataRate = (r.NDBPS/symbolTime)*1e3/500; % in units of 500 Kbps

                case 3 % VHT
                    vhtConfig = obj.VHTConfigObject;
                    vhtConfig.ChannelBandwidth = cbw;
                    vhtConfig.MCS = mcs;
                    vhtConfig.NumTransmitAntennas = numSTS;
                    vhtConfig.NumSpaceTimeStreams = numSTS;

                    r = wlan.internal.getRateTable(vhtConfig);
                    symbolTime = 4; % in microseconds
                    dataRate = (r.NDBPS(1)/symbolTime)*1e3/500;

                case {4, 5} % HE-SU, HE-EXT-SU
                    heSUConfig = obj.HESUConfigObject;
                    heSUConfig.ChannelBandwidth = cbw;
                    heSUConfig.MCS = mcs;
                    heSUConfig.ExtendedRange = (packetInfo.PPDUFormat == 5);
                    heSUConfig.NumTransmitAntennas = numSTS;
                    heSUConfig.NumSpaceTimeStreams = numSTS;
                    r = wlan.internal.heRateDependentParameters(ruInfo(heSUConfig).RUSizes,mcs,numSTS,heSUConfig.DCM);
                    symbolTime = 16; % in microseconds
                    dataRate = (r.NDBPS/symbolTime)*1e3/500;

                otherwise % EHT-SU
                    ehtSUConfig = wlanEHTMUConfig(cbw);
                    ehtSUConfig.User{1}.MCS = mcs;
                    ehtSUConfig.NumTransmitAntennas = numSTS;
                    [~,r] = wlan.internal.ehtCodingParameters(ehtSUConfig);
                    symbolTime = 16; % in microseconds
                    dataRate = (r.NDBPS/symbolTime)*1e3/500;
            end
        end

        function channelBytes = formRadiotapChannelBytes(obj,eventData)
            %formRadiotapChannelBytes Return the channel bytes for captured
            %frame
            %
            %   EVENTDATA is a structure containing event data, source
            %   node object, and event name.

            % u16 frequency, u16 flags
            channelBytes = zeros(1,4);
            frequencyInMHz = round(eventData.Data.Frequency/1e6);
            frequencyHexBits = dec2hex(frequencyInMHz);
            % Scale to 4 bits if required
            if numel(frequencyHexBits) < 4
                frequencyHexBits = [repmat('0', 1, ...
                    4 - numel(frequencyHexBits)), frequencyHexBits];
            end
            channelBytes(1:2) = [hex2dec(frequencyHexBits(3:4)) hex2dec(frequencyHexBits(1:2))];

            devCfg = getDeviceConfig(obj, eventData);
            band = devCfg.BandAndChannel(1);

            if band == 5
                channelFlagsHexBits = '0140';
            elseif band == 2.4
                channelFlagsHexBits = '00c0';
            else % 6 GHz band
                channelFlagsHexBits = '0040';
            end
            channelBytes(3:4) = [hex2dec(channelFlagsHexBits(3:4)) hex2dec(channelFlagsHexBits(1:2))];
        end

        function mcsBytes = formRadiotapMCSBytes(~, packetInfo)
            %formRadiotapMCSBytes Returns the MCS bytes for captured frame
            %
            % PACKETINFO is a structure containing parameters provided by PHY 
            % upon receipt of a valid PHY header.

            % MCS information, bandwidth, MCS index present, and check for
            % HT-format enabled
            if packetInfo.PPDUFormat==2 % HT-Mixed
                cbwknown = 1;
                switch(packetInfo.ChannelBandwidth)
                    case 20
                        chanBW = 0;
                    case 40
                        chanBW = 1;
                    otherwise
                        chanBW = 0; % 20L and 20U sidebands are not implemented, default 0
                end
                HTFormatknown = 1;
            else
                cbwknown = 0;
                HTFormatknown = 0;
            end

            mcsKnown = [0 0 0 0 HTFormatknown 0 1 cbwknown];
            mcsKnown_dec = bit2int(mcsKnown',8);
            if cbwknown == 1
                mcsBytes = [mcsKnown_dec chanBW packetInfo.PerUserInfo.MCS];
            else
                mcsBytes = [mcsKnown_dec 0 packetInfo.PerUserInfo.MCS];
            end
        end

        function VHTInformationBytes = formRadiotapVHTBytes(~,packetInfo)
            %VHTInformationBytes Return the radiotap VHT bytes for captured
            %frame
            %
            % PACKETINFO is a structure containing parameters provided by PHY 
            % upon receipt of a valid PHY header.


            VHTInformationBytes = zeros(1,12);
            % Form VHT known byte, only bandwidth known
            VHTKnownWord = ['00';'40'];
            VHTKnownBytes = flip(hex2dec(VHTKnownWord).');
            VHTInformationBytes(1:2) = VHTKnownBytes;

            % Form bandwidth byte
            switch(packetInfo.ChannelBandwidth)
                case 20
                    cbw = 0;
                case 40
                    cbw = 1;
                case 80
                    cbw = 2;
                case 160
                    cbw = 3;
            end

            VHTInformationBytes(4) = cbw;
            % Form mcs_nss byte for 1 user Number of spatial streams is
            % same as number of space time streams if STBC is not in use
            assert(packetInfo.PerUserInfo.MCS<=9,'Max MCS is 9');
            assert(packetInfo.PerUserInfo.NumSpaceTimeStreams<=8,'Max NSS is 8');
            VHTWord_mcsnss = int2str(10*packetInfo.PerUserInfo.MCS+packetInfo.PerUserInfo.NumSpaceTimeStreams);
            VHTBytes_mcsnss = hex2dec(VHTWord_mcsnss);
            VHTInformationBytes(5:8) = [VHTBytes_mcsnss 0 0 0];
        end

        function timestampBytes = formRadiotapTimestampBytes(~,captureTime)
            %formRadiotapTimestampBytes Returns the radiotap timestamp
            %bytes for captured frame

            % Hex Timestamp in nanoseconds
            timestamp_hexbits = dec2hex(round(captureTime*1e9, 0));
            timestampWord = strings(1,8);
            num = numel(timestamp_hexbits);
            % Form the timestamp word
            idx = 1;
            count = 1;
            while num>0
                if num==1
                    timestampWord(count) = timestamp_hexbits(idx);
                else
                    timestampWord(count) = timestamp_hexbits(idx:idx+1);
                end

                idx = idx+2;
                num = num -2;
                count = count+1;
            end
            x = find(cellfun(@isempty,timestampWord),1);
            timestampWord(1:x-1) = timestampWord(x-1:-1:1);
            % Form timestamp bytes
            timestampBytes = hex2dec(timestampWord);
            accuracy = [1 0]; % 1
            position = 2; % End of frame
            unit = 2; % Nanoseconds
            % Pos value (4 bits) followed by unit value(4 bits)
            unit_pos = bin2dec(strcat(dec2bin(position,4), dec2bin(unit,4)));
            flags = 2; % Accuracy known
            timestampBytes = [timestampBytes accuracy unit_pos flags];
        end

        function HEInformationBytes = formRadiotapHEBytes(~,packetInfo)
            %HEInformationBytes Returns the radiotap HE bytes for the
            %captured frame
            %
            % PACKETINFO is a structure containing parameters provided by PHY 
            % upon receipt of a valid PHY header.

            % Form HE information bytes. Each of data1, data2, data3,
            % data4, data5, and data6 is 2 bytes. These bytes are formed as
            % per the radiotap definition.
            HEInformationBytes = zeros(1, 12);
            % Form data1
            switch packetInfo.PPDUFormat
                case 4 % HE_SU(4)
                    format = 0; % 0=HE_SU, 1=HE_EXT_SU, 2=HE_MU, 3=HE_TRIG (radiotap)
                case 5 % HE_EXT_SU(5)
                    format = 1;
                case 6 % HE_MU(6)
                    format = 2;
                case 7 % HE_TB(7)
                    format = 3;
                otherwise
                    format = 0; % Show HE_SU for all other frames
            end
            formatBits = dec2bin(format);
            % Add ones for data bandwidth, data MCS known, BSS color
            if numel(formatBits)==1
                data1Bits = [0 1 0 0, 0 0 0 0, 0 0 1 0, 0 1 0 str2double(formatBits(1))];
            else
                data1Bits = [0 1 0 0, 0 0 0 0, 0 0 1 0, 0 1 str2double(formatBits(1)) str2double(formatBits(2))];
            end

            a = bit2int(data1Bits(9:16)',8);
            b = bit2int(data1Bits(1:8)',8);
            data1Bytes = [a b];
            HEInformationBytes(1:2) = data1Bytes;
            % Form data2
            data2Bytes = [64 0]; % Hardcoded to include only TxOP known flag
            HEInformationBytes(3:4) = data2Bytes;
            % Form data3
            mcs_hex = dec2hex(packetInfo.PerUserInfo.MCS);
            bss_hex = dec2hex(packetInfo.BSSColor);
            data3Word = ['0' mcs_hex '0' bss_hex];
            data3Word = reshape(data3Word, 2, []).';
            data3Bytes = flip(hex2dec(data3Word).');
            HEInformationBytes(5:6) = data3Bytes;
            % Form data5, data4 not necessary
            switch(packetInfo.ChannelBandwidth)
                case 20
                    cbw = 0;
                case 40
                    cbw = 1;
                case 80
                    cbw = 2;
                case 160
                    cbw = 3;
            end
            cbw_hex = dec2hex(cbw);
            data5Word = ['0' '0' '0' cbw_hex];
            data5Word = reshape(data5Word, 2, []).';
            data5Bytes = flip(hex2dec(data5Word).');
            HEInformationBytes(9:10) = data5Bytes;
            % Form data6
            txOP_hex = dec2hex(packetInfo.TXOPDuration);
            numSTS_hex = dec2hex(packetInfo.PerUserInfo.NumSpaceTimeStreams);
            data6Word = [txOP_hex '0' numSTS_hex];
            if numel(data6Word)<4
                data6Word = ['0' data6Word];
            end

            data6Word = reshape(data6Word, 2, []).';
            data6Bytes = flip(hex2dec(data6Word).');
            HEInformationBytes(11:12) = data6Bytes;
        end

        function LSIGBytes = formRadiotapLSIGBytes(~,packetInfo)
            %formRadiotapLSIGBytes Returns the LSIG radiotap bytes for the
            %captured frame
            %
            % PACKETINFO is a structure containing parameters provided by PHY 
            % upon receipt of a valid PHY header.

            % Form L-SIG DURATION bytes, u16 data1, data2
            LSIGBytes = zeros(1, 4);
            LSIGdata1Word = ['00';'02'];
            LSIGdata1Bytes = flip(hex2dec(LSIGdata1Word).');
            LSIGBytes(1:2) = LSIGdata1Bytes;
            LSIGlength_hex = dec2hex(packetInfo.LSIGLength);
            switch numel(LSIGlength_hex)
                case 1
                    LSIGdata2Word = ['0' '0' LSIGlength_hex '0'];
                case 2
                    LSIGdata2Word = ['0' LSIGlength_hex '0'];
                case 3
                    LSIGdata2Word = [LSIGlength_hex '0'];
            end
            LSIGdata2Word = reshape(LSIGdata2Word, 2, []).';
            LSIGdata2Bytes = flip(hex2dec(LSIGdata2Word).');
            LSIGBytes(3:4) = LSIGdata2Bytes;
        end

        function fileExtension = setFileExtension(~,inputExtension)
            if strcmp(inputExtension, "pcapng") || strcmp(inputExtension, "pcap")
                fileExtension = inputExtension;
            else
                error("File extension must be ""pcap"" or ""pcapng""");
            end
        end

        function deviceConfig = getDeviceConfig(~, eventData)
        %getDeviceConfig Returns the object holding MAC/PHY configuration

            devCfg = eventData.Source.DeviceConfig;
            capturedPacketdeviceID = eventData.Data.DeviceID;
            if isa(devCfg, "wlanDeviceConfig")
                deviceConfig = devCfg(capturedPacketdeviceID);
            else % wlanMultilinkDeviceConfig
                deviceConfig = devCfg.LinkConfig(capturedPacketdeviceID);
            end
        end
    end
end
