classdef hSLSTGaxAbstractSystemChannel < hSLSTGaxSystemChannelBase
%hSLSTGaxAbstractSystemChannel Create a channel manager object for
%abstracted PHY
%
%   CM = hSLSTGaxAbstractSystemChannel(CHAN,NUMANTENNAS) returns a channel
%   manager object for the specified channel configuration object CHAN and
%   an array containing the number of antennas per node NUMANTENNAS. CHAN
%   is a wlanTGaxChannel, wlanTGacChannel or wlanTGnChannel object. The
%   ChannelFiltering property must be set to false and the NumSamples
%   property must be set. The channel configuration is assumed to be the
%   same between all nodes. This helper object assumes all nodes can
%   transmit and receive, and channels between the nodes are reciprocal.
%
%   hSLSTGaxAbstractSystemChannel properties:
%
%   Links - Array of structures containing the channel for each link
%
%   hSLSTGaxAbstractSystemChannel methods:
%
%   getChannelStatistics - returns the channel stats between a pair of
%                          nodes

%   Copyright 2022-2025 The MathWorks, Inc.

    methods
        function obj = hSLSTGaxAbstractSystemChannel(varargin)
            % CM = hSLSTGaxAbstractSystemChannel(CHAN,NUMANTENNAS) returns
            % a channel manager object for the specified channel
            % configuration object CHAN and an array containing the number
            % of antennas per node NUMANTENNAS. This assumes all nodes can
            % transmit and receive and channels between the nodes are
            % reciprocal.
            obj = obj@hSLSTGaxSystemChannelBase(varargin{:})
        end

        function s = getChannelStatistics(obj,varargin)
            % S = getChannelStatistics(CM,TXIDX,RXIDX,SIMTIME) returns a
            % structure containing channel statisticschannel between a pair
            % of nodes with transmitter index TXIDX and receiver index
            % RXIDX at the given time in seconds SIMTIME.
            %
            % Once a channel realization is created using the above method,
            % subsequent calls to the method with the same TXIDX, RXIDX
            % pair will return the same channel matrix. A new channel
            % manager object must be created to obtain a new realization.
            %
            % S = getChannelStatistics(CM,SIG,RXINFO) returns a structure
            % containing channel statistics at the given time in seconds
            % given signal structure SIG and receiver info structure
            % RXINFO.

            sigStructPresent = isstruct(varargin{1});
            if ~sigStructPresent
                % Tx node index (and sim time) passed
                txIdx = varargin{1};
                rxIdx = varargin{2};
                simTime = varargin{3};
            else
                % Structure containing signal passed
                sig = varargin{1};
                rxInfo = varargin{2};
                txIdx = sig.TransmitterID;
                rxIdx = rxInfo.ID;
                numTxAnts = sig.NumTransmitAntennas;
                numRxAnts = rxInfo.NumReceiveAntennas;
                sr = sig.SampleRate;

                % Return channel statistics closest to midpoint of packet
                simTime = sig.StartTime+sig.Duration/2; % seconds
            end

            % Evolving channel, get path gain for desired simulation time
            numSamples = 1; % Get a single path gain for each packet
            interpMethod = 0; % 0 = closest, 2 = linear
            [pathGains,sampleTimes] = getPathGains(obj,txIdx,rxIdx,numSamples,simTime,interpMethod);

            if sigStructPresent
                % Extract only the tx/rx antennas needed if that
                % information is provided. To support EMLSR, the pathgains
                % are generated for the maximum possible number of
                % antennas.
                pathGains = extractRequiredPathGains(obj,pathGains,numTxAnts,numRxAnts);

                % Frequency shift path gains so they are centered at the
                % transmission center frequency if the channel center
                % frequency is different than the transmission center
                % frequency
                pathGains = frequencyShiftPathGains(obj,pathGains,sig.CenterFrequency,sig.TransmitterID,rxInfo.ID);

                [pathFilters,pathDelays] = getPathFilters(obj,txIdx,rxIdx,sr); % Sample rate can change so get appropriate path filter
            else
                [pathFilters,pathDelays] = getPathFilters(obj,txIdx,rxIdx); % Same filters and delays for all channels
            end

            s = struct;
            s.PathGains = pathGains;
            s.PathFilters = pathFilters;
            s.PathDelays = pathDelays;
            s.SampleTimes = sampleTimes;
        end
    end
end
