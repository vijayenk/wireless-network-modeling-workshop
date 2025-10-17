classdef hSLSTGaxSystemChannel < hSLSTGaxSystemChannelBase
%hSLSTGaxSystemChannel Create a channel manager object for full PHY
%
% CM = hSLSTGaxSystemChannel(CHAN,NUMANTSPERNODE) returns a channel
% manager object for the specified channel configuration object CHAN
% and an array containing the number of antennas per node
% NUMANTSPERNODE. This assumes all nodes can transmit and receive and
% channels between the nodes are reciprocal.
%
% CHAN is a wlanTGaxChannel, wlanTGacChannel or wlanTGnChannel object.
% The ChannelFiltering property must be set to false and the NumSamples
% property must be set. The channel configuration is assumed to be the
% same between all nodes.
%
%   hSLSTGaxSystemChannel properties:
%
%   Links           - Array of structures containing the channel for
%                     each link
%   CenterFrequency - Center frequency of all channels in Hertz
%
%   hSLSTGaxSystemChannel methods:
%
%   applyChannelToWaveform        - Filters a waveform through a
%                                   channel between two nodes
%   applyChannelToSignalStructure - Filters a waveform in a signal
%                                   structure through a channel between
%                                   two nodes

%   Copyright 2022-2025 The MathWorks, Inc.

    properties (Constant)
        ChannelFiltersStruct = struct('LinkIndex',zeros(0,1),'NumTxAnts',zeros(0,1),'SampleRate',zeros(0,1),'ChannelFilter',zeros(0,1),'Length',zeros(0,1),'Delay',zeros(0,1));
    end

    properties (Access=private)
        % Structure array containing channel filters with the fields:
        %   LinkIndex - Index of the link between a pair of nodes. Assumes
        %               each pair of nodes has a unique set of path delays.
        %   NumTxAnts - Number of transmit antennas
        %   SampleRate    - Sample rate in Hz of the waveform to filter
        %   ChannelFilter - comm.ChannelFilter object
        %   Length        - Filter length
        %   Delay         - Filter delay in samples
        ChannelFilters;
    end

    methods
        function obj = hSLSTGaxSystemChannel(varargin)
            % CM = hSLSTGaxSystemChannel(CHAN,NUMANTENNAS) returns a
            % channel object for the specified channel configuration object
            % CHAN and an array containing the number of antennas per node
            % NUMANTENNAS. This assumes all nodes can transmit and receive
            % and channels between the nodes are reciprocal.
            obj = obj@hSLSTGaxSystemChannelBase(varargin{:})
        end

        function reset(obj)
            % Clear channel filters
            obj.ChannelFilters = repmat(obj.ChannelFiltersStruct,1,0);
            reset@hSLSTGaxSystemChannelBase(obj);
        end

        function [sig,pg,chanInfo] = applyChannelToSignalStructure(obj,sig,rxInfo)
            % SIG = applyChannelToSignalStructure(OBJ,SIG,RXIFNO) filters
            % the waveform in a signal structure SIG through the channel
            % between two nodes. The receiver is specified by the structure
            % RXINFO.

            [numSamples,numTxAnts] = size(sig.Data);
            [chanFilt,filterLen,filterDelay] = getChannelFilter(obj,sig.TransmitterID,rxInfo.ID,sig.SampleRate,numTxAnts);
            if nargout>2
                chanInfo = info(chanFilt);
            end

            % Trailing zeros will be added to data to allow for channel delay
            numPadSamples = filterLen-1;

            % Get path gains for all samples of input data
            numSamplesToSim = numSamples+numPadSamples;
            simTime = sig.StartTime; % seconds
            pg = getPathGains(obj,sig.TransmitterID,rxInfo.ID,numSamplesToSim,simTime);

            % Extract only the tx/rx antennas needed if that information is
            % provided. To support EMLSR, the pathgains are generated for
            % the maximum possible number of antennas. Extract the
            % appropriate path gains for channel filtering.
            pg = extractRequiredPathGains(obj,pg,sig.NumTransmitAntennas,rxInfo.NumReceiveAntennas);

            % Frequency shift path gains so they are centered at the
            % transmission center frequency if the channel center frequency
            % is different than the transmission center frequency
            pg = frequencyShiftPathGains(obj,pg,sig.CenterFrequency,sig.TransmitterID,rxInfo.ID);

            % Add trailing zeros and pad for antenna selection as required
            dataPad = [sig.Data; zeros(numPadSamples,numTxAnts)];

            % Reset filter as we assume one packet filtered at a time and
            % we are jumping ahead in time and we don't want any internal
            % state
            reset(chanFilt);

            % Filter waveform
            filteredData = chanFilt(dataPad,pg);

            % Remove implementation delay
            sig.Data = filteredData(filterDelay+1:end,:);
            pg = pg(filterDelay+1:end,:,:,:);

            % Add trailing transient to packet duration in seconds
            numTransientSamples = filterLen-1-filterDelay;
            sig.Duration = sig.Duration+(numTransientSamples/chanFilt.SampleRate);
        end

        function [y,pg] = applyChannelToWaveform(obj,x,fs,txIdx,rxIdx,varargin)
            % Y = applyChannelToWaveform(OBJ,X,FS,TXIDX,RXIDX,[TIMEOFFSET])
            % filters the waveform X at sample rate FS through the channel
            % between node index TXIDX and RXIDX. TIMEOFFSET is optional
            % and specifies the time of the first sample to pass through
            % the channel in seconds. If a time offset is specified the
            % channel filter is reset as time is assumed to progress beyond
            % the filter group delay.

            numTx = size(x,2);
            chanFilt = getChannelFilter(obj,txIdx,rxIdx,fs,numTx);

            if nargin>4
                % If time offset provided reset filter as we assume one
                % packet filtered at a time and we are jumping ahead in
                % time and we don't want any internal state.
                reset(chanFilt);
            end
            % Get path gains for the required number of samples.
            numSamplesToSim = size(x,1);
            pg = getPathGains(obj,txIdx,rxIdx,numSamplesToSim,varargin{:});
            %  Filter waveform
            y = chanFilt(x,pg);
        end
    end

    methods (Access=private)
        function [chanFilt,filterLen,filterDelay] = getChannelFilter(obj,varargin)
            % CHANFILT = getChannelFilter(OBJ,TXIDX,RXIDX,FS,NUMTXANTS)
            % returns the channel filter between node index TXIDX and RXIDX
            % for sample rate FS Hz and number of transmit antennas
            % NUMTXANTS.
            %
            % CHANFILT = getChannelFilter(OBJ,IDX,FS)  returns the chnanel
            % filter for link index IDX with sampel rate FS Hz.

            if nargin==3
                % CHANFILT = getChannelFilter(OBJ,IDX,FS)
                idx = varargin{1};
                fs = varargin{2};
                numTxAnts = 0;
            else
                % CHANFILT = getChannelFilter(OBJ,TXIDX,RXIDX,FS,NUMTXANTS)
                idx = linkIndex(obj,varargin{1:2});
                fs = varargin{3};
                numTxAnts = varargin{4};
            end

            % Create a unique channel filter between  each pair of nodes
            % for the specified number of transmit antennas and sample
            % rate. The number of transmit antennas and sample rate cannot
            % change in comm.ChannelFilter once created, hence why we
            % create a new one. Before creating a new channel filter check
            % if one is already created
            if isempty(obj.ChannelFilters)
                % No channel filters, create new one
                [chanFilt,filterLen,filterDelay] = createChannelFilter(obj,idx,fs,numTxAnts);
            else
                existingFilterIdx = all([obj.ChannelFilters.LinkIndex]==idx,1) & [obj.ChannelFilters.SampleRate]==fs & [obj.ChannelFilters.NumTxAnts]==numTxAnts;
                if any(existingFilterIdx)
                    % Filter exists, use it
                    chanFilt = obj.ChannelFilters(existingFilterIdx).ChannelFilter;
                    filterLen = obj.ChannelFilters(existingFilterIdx).Length;
                    filterDelay = obj.ChannelFilters(existingFilterIdx).Delay;
                else
                    % Filter does not exist, create new one
                    [chanFilt,filterLen,filterDelay] = createChannelFilter(obj,idx,fs,numTxAnts);
                end
            end
        end

        function [chanFilt,filtLen,filtDelay] = createChannelFilter(obj,idx,fs,numTxAnts)
            % Create new channel filter and store in structure array. Force
            % NormalizeChannelOutputs to false as the signal strength
            % calculations in the simulator is not normalized.
            chanFilt = comm.ChannelFilter('PathDelays',getPathDelays(obj,idx),'SampleRate',fs,'NormalizeChannelOutputs',false);
            cinfo = info(chanFilt);
            filtLen = size(cinfo.ChannelFilterCoefficients,2);
            filtDelay = cinfo.ChannelFilterDelay;
            newChannelFilter = struct('LinkIndex',idx,'NumTxAnts',numTxAnts,'SampleRate',fs,'ChannelFilter',chanFilt,'Length',filtLen,'Delay',filtDelay);
            obj.ChannelFilters = cat(2,obj.ChannelFilters,newChannelFilter);
        end
    end
end
