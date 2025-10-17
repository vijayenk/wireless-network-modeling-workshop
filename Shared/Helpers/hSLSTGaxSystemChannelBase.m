classdef hSLSTGaxSystemChannelBase < handle
%hSLSTGaxSystemChannelBase Create a hSLSTGaxSystemChannelBase object
%
%   CM = hSLSTGaxSystemChannelBase(CHAN,NUMANTENNAS) returns a channel
%   manager object for the specified channel configuration object CHAN and 
%   an array containing the number of antennas per node NUMANTENNAS. This 
%   assumes all nodes can transmit and receive and channels between the 
%   nodes are reciprocal.
%
%   CHAN is a wlanTGaxChannel, wlanTGacChannel or wlanTGnChannel object. 
%   The ChannelFiltering property must be set to false and the NumSamples
%   property must be set. The channel configuration is assumed to be the
%   same between all nodes.
%
%   hSLSTGaxSystemChannelBase properties:
%
%   Links           - Array of structures containing the channel for
%                     each link with the following fields:
%                     Channel - Fading channel object. The channel
%                               sample rate is set to the lowest
%                               possible to generate samples given
%                               Doppler requirements of the channel.
%                               You can edit channel properties, in
%                               which case call the reset() method of
%                               hSLSTGaxAbstractChannel to update the
%                               channel sampling rate.
%                     Node1   - Identifier of the first node in the link
%                     Node2   - Identifier of the second node in the
%                              link
%                     SampleRate - Sample rate in Hz used to generate
%                                  path gains.
%                     ShadowFading - Shadow fading in dB.
%   CenterFrequency  - Center frequency of all channels in Hertz
%   ShadowFadingStandardDeviation - Shadow fading standard deviation
%                                   in dB
%   PathLossModel    - Path loss model
%   PathLossModelFcn - Custom path loss model function handle
%
%   hSLSTGaxSystemChannelBase methods:
%
%   reset           - Reset channel models
%   initialize      - Reset channels and calculate first path gains
%   getLink         - Return the link between nodes
%   getChannel      - Return the channel object between nodes
%   getPathDelays   - Calculate channel path delays between nodes
%   getPathFilters  - Calculate the channel path filters between nodes
%   getPathGains    - Calculate the channel path gain between nodes
%   getShadowFading - Calculate the shadow fading between nodes
%   getPathLoss     - Calculate the path loss between nodes

%   Copyright 2022-2025 The MathWorks, Inc.

   properties (Access=private)
        Channels;      % Cell array of WLAN channel model objects
        PathGains;     % Cell array of path gains for generated links
        PathFilters;   % Cell array of path gains for generated links
        PathFilterSampleRate; % Vector of sample rate for each set of path filters
        PathDelays;    % Cell array of path delays for all links
        PathTimes;     % Cell array of path times for generated links
        LastPathTime;  % Vector of last path time generated for each link
        SampleTimeOffset; % Vector of last sample time generated for each link
        PathTimeOffset;   % Vector of path time offsets for each link
    end

    properties (Access=protected)
        NumLinks;
    end

    properties
        %Links Structure array containing the channel for each link
        %   Structure with fields:
        %     Channel      - Channel model object
        %     Node1        - Identifier of first node in link
        %     Node2        - Identifier of second node in link
        %     SampleRate   - Nominal sample rate of the link in Hz 
        %     ShadowFading - Shadow fading loss for the link in dB
        Links;
        %CenterFrequency Channel operating frequency in Hz
        CenterFrequency;
        %ShadowFadingStandardDeviation Shadow fading standard deviation in dB
        %    Set to 0 for no shadow fading
        ShadowFadingStandardDeviation = 0;
        %PathLossModel Path loss model
        %    Set to 'free-space', 'enterprise', 'residential', or 'custom'.
        %    'enterprise' and  'residential' use the formulas given in IEEE
        %    802.11-14/0980r16. 'custom' uses the path loss function
        %    handled PathLossModelFcn. The default is 'free-space'.
        PathLossModel = 'free-space';
        % PathLossModelFcn Custom path loss model function handle
        %    Function handle which returns path loss in dB:
        %    PathLossModelFcn(sig,rxInfo). sig is a structure containing
        %    information about the transmitted packet and transmitter.
        %    rxInfo is a structure containing information about the
        %    receiver.
        PathLossModelFcn;
        % ChannelIndicesLUT Lookup table of node indices for each channel
        %    Matrix of NumLinks-by-2 storing the node indices for each
        %    link. This is created when the object is constructed but
        %    may be overwritten.
        ChannelIndicesLUT; % Array used to map link index
    end

    properties (Constant, Hidden=true)
        PacketIterationSimTime = 10e-3; % Simulate 10 ms of channel at a time
        LightSpeed = physconst('lightspeed');
    end

    properties (Access=private)
        % Vector of unique node identifiers in order that they were
        % specified when creating the channel
        UniqueNodeIndices;
    end

    methods
        function obj = hSLSTGaxSystemChannelBase(chan,numAntennas,varargin)
            % CM = hSLSTGaxSystemChannelBase(CHAN,NUMANTENNAS) returns a
            % channel manager object for the specified channel
            % configuration object CHAN and an array containing the number
            % of antennas per node NUMANTENNAS. This assumes all nodes can
            % transmit and receive and channels between the nodes are
            % reciprocal.
            %
            % CM = hSLSTGaxSystemChannelBase(...,ShadowFadingStandardDeviation=val)
            % sets the shadow fading standard deviation in dB. The default
            % is 0 dB (no shadow fading).

            % Generate downlink and uplink channels between transmitters
            % and receivers assuming they may have different numbers of
            % antennas.
            numNodes = numel(numAntennas);
            if numNodes<2
                error('More than 1 node is required to create a channel')
            end

            % Generate channels between all transmitters and receivers.
            % * Assume a transmitter and receiver are the same node,
            %   therefore do not generate a channel between a node and
            %   itself.
            % * Assume the channel is reciprocal between a transmitter
            %   and receiver, but antenna configuration can be different.
            %
            % Create a mapping function to map requested channel between a
            % pair of nodes (tx,rx) to a channel. Each row contains the
            % transmitter index and receiver index of a channel. The
            % transmitter index will always be lower than the receiver
            % index.
            obj.ChannelIndicesLUT = nchoosek(1:numNodes,2);

            % Set name value pair properties now to allow
            % ChannelIndicesLUT to be overwritten to allow for node IDs to
            % be used.
            for i = 1:2:nargin-3
                obj.(varargin{i}) = varargin{i+1};
            end

            % Allocate arrays to store channels generated
            obj.NumLinks = height(obj.ChannelIndicesLUT);
            obj.Links = repmat(struct('Channel',[],'Node1',[],'Node2',[],'SampleRate',chan.SampleRate,'ShadowFading',[]),1,obj.NumLinks);
            obj.Channels = cell(1,obj.NumLinks);
            obj.PathDelays = cell(obj.NumLinks,1);
            obj.PathFilters = cell(obj.NumLinks,1);
            obj.PathFilterSampleRate = zeros(obj.NumLinks,1);
            obj.PathGains = cell(obj.NumLinks,1);
            obj.PathTimes = cell(obj.NumLinks,1);
            obj.CenterFrequency = chan.CarrierFrequency;

            obj.UniqueNodeIndices = unique(obj.ChannelIndicesLUT,'stable');

            % Generate channels
            release(chan);
            for i = 1:obj.NumLinks
                node1ID = obj.ChannelIndicesLUT(i,1);
                node2ID = obj.ChannelIndicesLUT(i,2);

                obj.Channels{i} = clone(chan);
                obj.Channels{i}.ChannelFiltering = false; % Disable channel filtering as will be done externally
                obj.Channels{i}.NumTransmitAntennas = numAntennas(obj.UniqueNodeIndices==node1ID);
                obj.Channels{i}.NumReceiveAntennas = numAntennas(obj.UniqueNodeIndices==node2ID);

                % Create public array of channels user can control
                obj.Links(i).Node1 = node1ID;
                obj.Links(i).Node2 = node2ID;
                obj.Links(i).Channel = obj.Channels{i};

                % Log-normal shadow fading in dB
                obj.Links(i).ShadowFading = obj.ShadowFadingStandardDeviation*randn;
            end

            % Reset channels and calculate Doppler dependent parameters
            reset(obj);
        end

        function reset(obj,varargin)
            % reset(OBJ) reset all channel models reset(OBJ,TXIDX,RXIDX)
            % channel model specified by transmit and receive index

            if nargin>1
                txIdx = varargin{1};
                rxIdx = varargin{2};
                linksToSet = obj.sub2linkInd(txIdx,rxIdx);
            else
                linksToSet = 1:obj.NumLinks;
            end

            % Set sample rate and number of samples for each channel at
            % lowest rate given Doppler frequency
            for i = linksToSet
                % Reset
                obj.PathGains{i} = [];
                obj.PathTimes{i} = [];
                obj.PathDelays{i} = [];
                obj.PathFilters{i} = [];
                obj.PathFilterSampleRate(i) = 0;
                obj.SampleTimeOffset(i) = 0;
                obj.LastPathTime(i) = -1;
                obj.PathTimeOffset(i) = 0;

                if isempty(obj.Channels) || isempty(obj.Channels{i})
                    % No channel exists
                    continue
                end

                % Log-normal shadow fading in dB
                if obj.ShadowFadingStandardDeviation>0
                    obj.Links(i).ShadowFading = obj.ShadowFadingStandardDeviation*randn;
                end

                release(obj.Channels{i});
                obj.Channels{i}.SampleRate = obj.Links(i).SampleRate; % Reset sample rate for regenerating path filters as potentially reduced by code.
                obj.Channels{i}.ChannelFiltering = false; % Disable channel filtering as will be done externally

                % Set channel sampling rate required for Doppler component
                wavelength = 3e8/obj.Channels{i}.CarrierFrequency;
                fdoppler = (obj.Channels{i}.EnvironmentalSpeed*(5/18))/wavelength; % Cut-off frequency (Hz), change km/h to m/s
                normalizationFactor = 1/300;
                fc = fdoppler/normalizationFactor; % Channel sampling frequency
                if fc>0
                    fc = fc+1e-6; % Add to avoid numeric issues comparing oversampled and input sample rates when filtering
                end
                interpolationFactor = 1/40; % Interpolation required for Fluorescent effect

                % Set sample rate of channel to lowest possible to generate
                % path gains based on Doppler frequency
                obj.Channels{i}.SampleRate = max(fc/interpolationFactor,1e-3); % minimum very low sample rate (TODO handle 0 speed better)

                % Calculate how many path gain samples to generate so that
                % it will have at least as many required for one packet
                % step time
                numPathGainSamples = max(ceil(obj.PacketIterationSimTime*obj.Channels{i}.SampleRate),2); % At least 2 samples required
                while (numPathGainSamples-1)/obj.Channels{i}.SampleRate < obj.PacketIterationSimTime
                    % As the first sample is time 0, make sure we have
                    % enough samples to capture the simulation time
                    numPathGainSamples = numPathGainSamples+1;
                end
                obj.Channels{i}.NumSamples = numPathGainSamples;
            end
        end

        function l = getLink(obj,varargin)
            % L = getLink(OBJ,TXIDX,RXIDX) returns the link structure
            % between node index TXIDX and RXIDX.
            %
            % L = getLink(OBJ,IDX) returns the link for link index
            % LINKIDX.

            l = obj.Links(linkIndex(obj,varargin{:}));
        end

        function chan = getChannel(obj,varargin)
            % CHAN = getChannel(OBJ,TXIDX,RXIDX) returns the channel model
            % object between node index TXIDX and RXIDX.
            %
            % CHAN = getChannel(OBJ,LINKIDX) returns the channel object for
            % the link index LINKIDX.

            chan = obj.Channels{linkIndex(obj,varargin{:})};
        end

        function pd = getPathDelays(obj,varargin)
            % PD = getPathDelays(OBJ,TXIDX,RXIDX) calculates channel path
            % delays between node index TXIDX and RXIDX.
            %
            % PD = getPathDelays(OBJ,LINKIDX) calculates the path delays
            % for the link index LINKIDX.

            idx = linkIndex(obj,varargin{:});

            if isempty(obj.PathDelays{idx})
                % Get path delays and filters if they do not already exist
                chanInfo = info(obj.Channels{idx});
                obj.PathDelays{idx} = chanInfo.PathDelays;
            end
            pd = obj.PathDelays{idx};
        end

        function [pf,pd] = getPathFilters(obj,txIdx,rxIdx,varargin)
            % PF = getPathFilters(OBJ,TXIDX,RXIDX) calculates channel path
            % filters between node index TXIDX and RXIDX for the sample
            % rate set in the base channel model.
            %
            % PF = getPathFilters(...,SR) specifies the sample rate in Hz.
            %
            % [PF,PD] = getPathFilters(...) additionally returns the path
            % delays.

            idx = linkIndex(obj,txIdx,rxIdx);
            if nargin>3
                sr = varargin{1};
            else
                sr = obj.Links(idx).SampleRate;
            end
            pd = getPathDelays(obj,idx);
            if isempty(obj.PathFilters{idx}) || obj.PathFilterSampleRate(idx)~=sr
                % Get path filters if they do not already exist or the
                % sample rate has changed
                obj.PathFilters{idx} = wireless.internal.L2SM.channelFilterCoefficients(pd,sr);
                obj.PathFilterSampleRate(idx) = sr;
            end
            pf = obj.PathFilters{idx};
        end

        function [pg,st] = getPathGains(obj,txIdx,rxIdx,numSamples,varargin)
            % [PG,ST] = getPathGains(OBJ,TXIDX,RXIDX,NUMSAMPLES) calculates
            % channel path gains and associated sample time between node
            % index TXIDX and RXIDX for NUMSAMPLES since the last call.
            % Linear interpolation is used to generate path gains at the
            % required sample rate.
            %
            % [PG,ST] =
            % getPathGains(OBJ,TXIDX,RXIDX,NUMSAMPLES,STARTTIME,METHOD)
            % specifies the STARTTIME of path gains.
            %
            % METHOD controls the interpolation method: 0 - returns closest
            % path gain for each sample time 2 - interpolate path gains
            % over sample times (linear)
            %
            % [PG,ST] = getPathGains(OBJ,TXIDX,RXIDX,2,[STARTTIME
            % ENDTIME],1) calculates channel path gains and associated
            % sample times. Returns all path gains within the sample
            % period, and one either side to allow for interpolation given
            % STARTIME and ENDTIME.

            interpMethod = 2; % linear

            % Extract channel information
            [idx,switched] = obj.sub2linkInd(txIdx,rxIdx);

            fs = obj.Links(idx).SampleRate;
            samplesSimTime = numSamples/fs;

            lastPathTime = obj.LastPathTime(idx);
            pathTimeOffset = obj.PathTimeOffset(idx);
            if nargin>4
                if nargin>5 && varargin{2} == 1 % method = start stop (1)
                    % getPathGains(OBJ,TXIDX,RXIDX,NUMSAMPLES,TIME,1)
                    assert(numSamples==2)
                    assert(numel(varargin{1})==2)
                    waveformStartTime = varargin{1}(1);
                    waveformEndTime = varargin{1}(2);
                    interpMethod = 1; % Start stop
                else
                    % getPathGains(OBJ,TXIDX,RXIDX,NUMSAMPLES,STARTTIME,[METHOD])
                    sampleTimeOffset = varargin{1};
                    if nargin>5
                        interpMethod = varargin{2};
                        assert(any(interpMethod==[0 2]),'Expected interpolation method to be linear or closest')
                    end
                end
            else
                % getPathGains(OBJ,TXIDX,RXIDX,NUMSAMPLES) Continue from
                % where we left off
                sampleTimeOffset = obj.SampleTimeOffset(idx);
            end

            if interpMethod == 1 % start end
                sampleTimes = ([waveformStartTime waveformEndTime]);
            else
                % Update the times for each sample to generate
                sampleTimes = sampleTimeOffset+((0):(1/fs):(samplesSimTime-(1/fs)));
            end
            obj.SampleTimeOffset(idx) = sampleTimes(end)+(1/fs);

            while sampleTimes(end)>lastPathTime
                % Generate new path gains if required
                chan = obj.Channels{idx};
                pgFc = chan();

                % Simulation time for each sample
                pathTimes = cast(pathTimeOffset+((0):(1/chan.SampleRate):((chan.NumSamples-1)/chan.SampleRate)),chan.OutputDataType);
                lastPathTime = pathTimes(end);
                pathTimeOffset = lastPathTime+(1/chan.SampleRate);

                if ~isempty(obj.PathTimes{idx})
                    % Keep path gains at times which are still needed,
                    % discard rest. Note we need 1 less than the minimum
                    % sample time to allow interpolation
                    keepIdx = find((obj.PathTimes{idx}>=sampleTimes(1))',1);
                    keepIdx = max(keepIdx-1,1);
                    if isempty(keepIdx)
                        % Keep at least one path time to make sure we will
                        % always have one that is one before the
                        % lastPathTime
                        keepIdx = numel(obj.PathTimes{idx});
                    end
                    obj.PathTimes{idx} = [obj.PathTimes{idx}(keepIdx:end); pathTimes.'];
                    obj.PathGains{idx} = [obj.PathGains{idx}(keepIdx:end,:,:,:); pgFc];
                else
                    % First time
                    obj.PathTimes{idx} = pathTimes.';
                    obj.PathGains{idx} = pgFc;
                end
                obj.LastPathTime(idx) = lastPathTime;
                obj.PathTimeOffset(idx) = pathTimeOffset;
            end

            % Channel path gains are generated (and stored) for one
            % direction. Therefore if the reciprocal channel is required
            % switch the transmit and receive antenna dimension.
            pgUse = obj.PathGains{idx};
            if switched
                pgUse = permute(pgUse,[1 2 4 3]);
            end

            switch interpMethod
                case 0 % closest
                    % Return path closest to time requested
                    [~,closestInd] = min(abs(obj.PathTimes{idx}-sampleTimes));
                    pg = pgUse(closestInd,:,:,:);
                    st = obj.PathTimes{idx}(closestInd)';
                case 1
                    % Return path gains which in packet duration, and one
                    % before and after to allow for interpolation
                    firstIdx = find(obj.PathTimes{idx}<=waveformStartTime,1,'last');
                    lastIdx = find(obj.PathTimes{idx}>=waveformEndTime,1,'first');
                    st = obj.PathTimes{idx}(firstIdx:lastIdx)';
                    pg = pgUse(firstIdx:lastIdx,:,:,:);
                otherwise % linear
                    % Interpolate path gains over sample times
                    st = sampleTimes';
                    pg = interp1(obj.PathTimes{idx},pgUse,sampleTimes');
            end
        end

        function pg = frequencyShiftPathGains(obj,pg,txCenterFrequency,txIdx,rxIdx)
            % Frequency shift path gains so they are centered at the
            % transmission center frequency if the channel center frequency
            % is different than the transmission center frequency
            frequencyOffset = txCenterFrequency-obj.CenterFrequency;
            if abs(frequencyOffset)>0
                pd = getPathDelays(obj,txIdx,rxIdx);
                pg = pg.*exp(1i*2*pi*frequencyOffset.*pd);
            end
        end

        function l = getShadowFading(obj,txIdx,rxIdx)
            % L = getShadowFading(OBJ,TXIDX,RXIDX) returns shadow fading in
            % dB between node index TXIDX and RXIDX.

            % Extract channel information
            idx = obj.sub2linkInd(txIdx,rxIdx);

            l = obj.Links(idx).ShadowFading;
        end

        function pl = getPathLoss(obj,sig,rxInfo)
            %pathLoss Calculates path loss based on the signal and receiver
            %information

            d = norm(sig.TransmitterPosition - rxInfo.Position);

            switch obj.PathLossModel
                case 'free-space'
                    pl = freeSpacePathLoss(obj, d);
                case 'residential'
                    pl = tgaxResidentialPathLoss(obj, d);
                case 'enterprise'
                    pl = tgaxEnterprisePathLoss(obj, d);
                case 'custom'
                    pl = obj.PathLossModelFcn(sig,rxInfo);
            end
        end
    end

    methods (Access=protected)
        function [idx,switched] = sub2linkInd(obj,txIdx,rxIdx)
            % Returns the link index given the transmit and receive node
            % indices
            idx = all(obj.ChannelIndicesLUT == [txIdx rxIdx],2);
            switched = false;
            if ~any(idx)
                % Try the switched indices if no match
                idx = all(obj.ChannelIndicesLUT == [rxIdx txIdx],2);
                switched = true;
            end
            % Check that the channel has been created and not more than one
            % channel exists between the pair. A channel does not exist
            % between a node and itself.
            assert(nnz(idx)<=1,'More than one channel exists between node #%d and #%d at this frequency.',txIdx,rxIdx)
            if ~any(idx) || isempty(obj.Links(idx).Channel)
                error('hSLSTGaxSystemChannelBase:NoChannelExists','Channel does not exist between Node ID %d and  Node ID %d.',txIdx,rxIdx)
            end
        end

        function idx = linkIndex(obj,varargin)
            %linkIndex returns the link index given either the
            %link index or transmitter and receiver node index.
            if nargin==2
                idx = varargin{1};
            else
                % Extract channel information
                txIdx = varargin{1};
                rxIdx = varargin{2};
                idx = obj.sub2linkInd(txIdx,rxIdx);
            end
        end

        function pathGains = extractRequiredPathGains(~,pathGains,numTxAnts,numRxAnts)
            % PATHGAINS = extractRequiredPathGains(OBJ,PATHGAINS,NUMTXANTS,NUMRXANTS)
            % extracts the path gains required given the number of transmit
            % and receive antennas used in a link from the maximum number
            % of path gains.

            [~,~,maxNumTxAnts,maxNumRxAnts] = size(pathGains);
            if maxNumTxAnts~=numTxAnts || maxNumRxAnts~=numRxAnts
                txAntIdx = 1:numTxAnts;
                rxAntIdx = 1:numRxAnts;
                pathGains = pathGains(:,:,txAntIdx,rxAntIdx);
            end
        end
    end

    methods (Access=private)
        function pl = tgaxEnterprisePathLoss(obj, d)
            %tgaxEnterprisePathLoss Apply distance-based path loss on the
            %packet and update relevant fields of the output data. As per
            %IEEE 802.11-14/0980r16

            d = max(d,1);
            W = 0; % Number of walls penetrated
            % Enterprise
            dBP = 10; % breakpoint distance
            pl = 40.052 + 20*log10((obj.CenterFrequency/1e9)/2.4) + 20*log10(min(d,dBP)) + (d>dBP) * 35*log10(d/dBP) + 7*W;
        end

        function pl = tgaxResidentialPathLoss(obj, d)
            %tgaxResidentialPathLoss Apply distance-based path loss on the
            %packet and update relevant fields of the output data. As per
            %IEEE 802.11-14/0980r16

            d = max(d,1);
            F = 0; % Number of floors penetrated
            W = 0; % Number of walls penetrated
            % Residential
            dBP = 5;
            pl = 40.052 + 20*log10((obj.CenterFrequency/1e9)/2.4) + 20*log10(min(d,dBP)) + (d>dBP) * 35*log10(d/dBP) + 18.3*F^((F+2)/(F+1)-0.46) + 5*W;
        end

        function pl = freeSpacePathLoss(obj, d)
            %freeSpacePathLoss Apply free space path loss on the packet and
            %update relevant fields of the output data

            % Calculate free space path loss (in dB)
            pl = fspl(d, obj.LightSpeed/obj.CenterFrequency);
        end
    end
end