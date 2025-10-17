classdef h38901Channel < handle
%h38901Channel TR 38.901 channel model
%
%   h38901Channel properties:
%
%   Scenario           - Deployment scenario ("UMi", "UMa", "RMa", "InH",
%                        "InF-SL", "InF-DL", "InF-SH", "InF-DH", "InF-HH")
%                        (default "UMa")
%   InterSiteDistance  - Intersite distance in meters (default 500)
%   Wrapping           - Geographical distance-based wrapping (true, false)
%                        (default true)
%   SpatialConsistency - Spatial consistency ("None", "Static", 
%                        "ProcedureA", "ProcedureB") (default "None")
%   UpdateDistance     - Update distance for spatially consistent mobility
%   Seed               - Random number generator (RNG) seed (default 0)
%   OfficeType         - Office type for InH scenario ("Mixed","Open")
%                        (default "Mixed")
%   ScenarioExtents    - Location and size of the scenario
%                        (default [])
%   HallSize           - Dimensions of hall for InF scenarios
%                        (default [120 60 10])
%   ClutterSize        - Clutter size for InF scenarios (default 2)
%   ClutterDensity     - Clutter density for InF scenarios (default 0.6)
%   ClutterHeight      - Clutter height for InF scenarios (default 6)
%   AbsoluteTOA        - Absolute time of arrival (false, true) 
%                        (default false)
%   LOSProbability     - LOS probability (default [])
%
%   h38901Channel object functions:
%
%   h38901Channel               - Create channel model
%   connectNodes                - Connect simulator nodes to channel
%   channelFunction             - Simulator custom channel function
%   createChannelLink           - Create a single channel link
%   applyBeamforming            - Apply beamforming (TXRU virtualization)
%   createAutoCorrMatrices      - Create a set of autocorrelation matrices
%   uniformAutoCorrRVs          - Generate uniformly distributed spatially 
%                                 correlated random variables
%   wrappingOffsets             - Distance offsets for wrap-around 
%                                 calculations
%   sitePolygon                 - Vertices of site boundary polygon
%   spatiallyConsistentMobility - Apply spatially consistent mobility

%   Copyright 2022-2025 The MathWorks, Inc.

    % =====================================================================
    % public interface

    properties (Access=public)

        % Deployment scenario, used to determine properties of the channel
        % links ("UMi", "UMa", "RMa", "InH", "InF-SL", "InF-DL", "InF-SH",
        % "InF-DH", "InF-HH") (default "UMa")
        Scenario = "UMa";

        % Intersite distance in meters. Only required for UMi, UMa or RMa
        % scenarios (default 500)
        InterSiteDistance = 500;

        % Enable wrap around calculations, as defined in Rec. ITU-R
        % M.2101-0 Attachment 2 to Annex 1. Only required for UMi, UMa or
        % RMa scenarios (default true)
        Wrapping = true;

        % Spatial consistency. Set to "None" (or false) to apply no spatial
        % consistency procedure. Set to "Static" (or true) to apply 
        % TR 38.901 Section 7.6.3.1 "Spatial consistency procedure". Set to
        % "ProcedureA" or "ProcedureB" to apply Procedure A or Procedure B
        % from TR 38.901 Secton 7.6.3.2 "Spatially-consistent UT/BS
        % mobility modelling" (default "None")
        SpatialConsistency = "None";

        % The change in BS-to-UE distance, in meters, that will trigger the
        % spatial consistency update procedure when SpatialConsistency is
        % set to "ProcedureA" or "ProcedureB". TR 38.901 Section 7.6.3.2 
        % states that this distance "should be limited within 1 meter" 
        % (default 1.0)
        UpdateDistance = 1.0;

        % Random number generator seed (default 0)
        Seed = 0;

        % Office type for InH scenario ("Mixed", "Open") (default "Mixed")
        OfficeType = "Mixed";

        % Location and size of the scenario, a four-element vector of the
        % form [left bottom width height]. The elements are defined as 
        % follows:
        %   left -   The X coordinate of the left edge of the scenario in 
        %            meters
        %   bottom - The Y coordinate of the bottom edge of the scenario in 
        %            meters
        %   width  - The width of the scenario in meters, that is, the 
        %            right edge of the scenario is left + width
        %   height - The height of the scenario in meters, that is, the 
        %            top edge of the scenario is bottom + height
        % Use empty ([]) to automatically calculate the value. For UMi, UMa
        % and RMa scenarios, the value is calculated assuming that each
        % nrGNB node lies at the center of a hexagonal cell with the size
        % given by the InterSiteDistance property. For InF scenarios, the
        % value is calculated from the HallSize property assuming that the
        % hall is centered on (0,0). For InH scenarios, the value is
        % calculated from the locations of the nodes attached to the
        % simulator. (default [])
        ScenarioExtents = [];

        % Dimensions of hall for InF scenarios, a 3-by-1 vector [L W H]
        % where L is the hall length, W is the hall width and H is the hall
        % height in meters (default [120 60 10]). Not required for UMi,
        % UMa, RMa or InH scenarios.
        HallSize = [120 60 10];

        % Clutter size in meters for InF scenarios (default 2)
        ClutterSize = 2;

        % Clutter density for InF scenarios (0...1) (default 0.6)
        ClutterDensity = 0.6;

        % Clutter height in meters for InF scenarios (default 6)
        ClutterHeight = 6;

        % Absolute time of arrival as defined in TR 38.901 Section 7.6.9.
        % Only applicable for InF scenarios (false, true) (default false)
        AbsoluteTOA = false;

        % LOS probability. A scalar between 0 and 1 giving the probability
        % that a randomly-dropped UE will be in line of sight (LOS)
        % condition. If empty, the LOS probability is determined using
        % TR 38.901 Table 7.4.2-1 (default [])
        LOSProbability = [];

    end

    methods

        function set.Scenario(obj,val)

            persistent sv;
            if (isempty(sv))
                sv = nrPathLossConfig.Scenario_Values;
            end
            obj.Scenario = validatestring(val,sv,'h38901Channel','Scenario');

        end

    end

    methods (Access=public)

        function channel = h38901Channel(varargin)
        % Create channel model

            % Set properties from name-value arguments
            setProperties(channel,varargin{:});

            % Set up path loss configuration
            channel.thePathLossConfig = nrPathLossConfig(Scenario=channel.Scenario);

            % Create map of attachments between UEs and BSs
            channel.theUEtoBSMap = dictionary([],[]);

            % Create map between IDs and UEs
            channel.theUEMap = dictionary([],struct());

            % Create map between IDs and BSs
            channel.theBSMap = dictionary([],struct());

            % Set up empty cache of site positions
            channel.theSitePositions = [];

            % Create map between link ID pairs and channels
            c = newChannel();
            channel.theLinkToChannelMap = dictionary([],c(false));

        end

        function connectNodes(channel,varargin)
        % connectNodes(CHANNEL,SLS) obtains the connections between
        % BSs and UEs by querying the wireless network simulator SLS.
        %
        % connectNodes(CHANNEL,SLS,CHCFG) additionally specifies
        % configuration properties for the channel links to be created.
        % CHCFG is a structure with the following fields:
        %   Site    
        %     - A row vector specifying the 1-based site index for each 
        %       gNB.
        %   Sector
        %     - A row vector specifying the 1-based sector index for each
        %       gNB.
        %   TXRUVirtualization 
        %     - A row vector of structures specifying the antenna
        %       virtualization parameters for each gNB. The
        %       parameterization is according to TR 36.897 Section 5.2.2
        %       TXRU virtualization model option-1B. Each structure has the
        %       following fields:
        %         K    - Vertical weight vector length
        %         Tilt - Tilting angle in degrees
        %         L    - Horizontal weight vector length
        %         Pan  - Panning angle in degrees
        %       The default value is struct(K=1,Tilt=0,L=1,Pan=0).
        %   TransmitArrayOrientation 
        %     - A matrix where the columns specify the transmit antenna
        %       orientations for each gNB. See
        %       nrCDLChannel/TransmitArrayOrientation. If this field is
        %       present then the Sector field is ignored, as the Sector
        %       field is only used to determine the bearing angle alpha of
        %       the array orientation.
        %   n_fl    
        %     - A row vector specifying the 1-based floor number for each
        %       UE (1 = ground floor). Only required for UMi, UMa or RMa
        %       scenarios. If absent, floor numbers will be determined
        %       automatically from UE node heights according to TR 36.873
        %       Table 6-1.
        %   d_2D_in 
        %     - A row vector specifying the 2-D indoor distance each UE in
        %       meters. Only required for UMi, UMa or RMa scenarios. If
        %       absent, the 2-D indoor distances default to 0 for all UEs
        %       i.e. all UEs are outdoor.
        %   ReceiveArrayOrientation 
        %     - A matrix where the columns specify the receive antenna
        %       orientations for each UE. See
        %       nrCDLChannel/ReceiveArrayOrientation.
        % The order of the gNBs and UEs in the columns must be the order
        % in which the nrGNB and nrUE nodes appear in the SLS.Nodes array.
        % Alternatively if the TXRUVirtualization, n_fl, d_2D_in, 
        % TransmitArrayOrientation or ReceiveArrayOrientation field has a
        % single column, it applies to all gNBs or all UEs as appropriate.

            % Determine the number of input arguments that are not
            % name-value arguments and the position of the first name (if
            % present)
            firstnvname = find(cellfun(@(x)(ischar(x) || isstring(x)),varargin),1,'first');
            if (isempty(firstnvname))
                ninarg = nargin;
            else
                % -1 for the first NV name, +1 for obj
                ninarg = firstnvname - 1 + 1;
            end

            % Get gNBs and UEs
            sls = varargin{1};
            toMat = @(x)[x{:}];
            nodes = sls.Nodes;
            nodeTypes = cellfun(@class,nodes,UniformOutput=false);
            nodesOfType = @(t)toMat(nodes(strcmp(nodeTypes,t)));
            gNBs = nodesOfType('nrGNB');
            UEs = nodesOfType('nrUE');
            if (~isempty(UEs))
                ueIDs = [UEs.ID];
                UEs = cellfun(@(x)UEs(ueIndices(ueIDs,x)),{gNBs.UENodeIDs},UniformOutput=false);
            else
                ueIDs = [];
            end

            % Get channel configuration structure
            if (ninarg==3)
                chCfg = varargin{2};
                chCfg.UEIDs = ueIDs;
            else % ninarg==2
                chCfg = [];
            end

            % Set properties from name-value arguments; this allows the
            % hidden properties InterfererHasSmallScale, PhaseLOS_d_3D
            % and InterfererSameLinkEnd to be controlled
            setProperties(channel,varargin{firstnvname:end});

            % Record list of gNBs, connections between gNBs and UEs, and
            % other scenario information that can be established from node
            % names or the channel configuration structure
            recordNodes(channel,gNBs,UEs,chCfg);

            function ind = ueIndices(ueIDs,x)
                if (isempty(x))
                    ind = [];
                else
                    ind = any(ueIDs==x.',1);
                end
            end

        end
        
        function packet = channelFunction(obj,rxinfo,packet)
        % RXPACKET = CHANNEL.channelFunction(RXINFO,TXPACKET) is the custom
        % channel function CUSTOMCHANNELFCN described in
        % wirelessNetworkSimulator/addChannelModel. Call
        % wirelessNetworkSimulator/addChannelModel and pass a handle to
        % CHANNEL.channelFunction to connect the channel model to the
        % simulator.
            
            % Check if connectNodes method has been called
            if (~(size(obj.theRecordedSitePositions,2)==3))
                error('nr5g:h38901Channel:NoConnectNodes','Call the connectNodes method to connect simulator nodes to the channel before executing the simulation.');
            end

            % -------------------------------------------------------------
            % TR 38.901 Section 7.5 Steps 2 - 10
            % Get channel for the current link
            [ch,bsID,ueID,linkind] = getSLSChannel(obj,rxinfo,packet);
            % -------------------------------------------------------------

            % If the channel is empty, signifying a BS-to-BS or UE-to-UE
            % link when InterfererSameLinkEnd=false, return an empty packet
            % that will be dropped by the simulator
            if (isempty(ch))
                packet = [];
                return;
            end

            % Apply spatially consistent mobility if required
            if (isSpatiallyConsistentMobility(obj))
                ch = packetMobilitySLS(obj,ch,packet,rxinfo);
                obj.theLinkToChannelMap(linkind) = ch;
            end

            % If the small scale channel is a CDL channel
            if (~isstruct(ch.SmallScale))

                % Configure channel according to packet StartTime and
                % Duration
                ch.SmallScale.InitialTime = packet.StartTime;
                ch.SmallScale.NumTimeSamples = ceil(packet.Duration * ch.SmallScale.SampleRate);

                % ---------------------------------------------------------
                % TR 38.901 Section 7.5 Step 11
                % Execute the channel
                [pathGains,sampleTimes] = ch.SmallScale();
                pathDelays = ch.PathDelays;
                pathFilters = ch.PathFilters;

                % Apply d_3D-related term of Eq 7.5-29
                if (obj.PhaseLOS_d_3D)
                    channelInfo = info(ch.SmallScale);
                    los = strcmpi(channelInfo.ClusterTypes,'LOS');
                    pathGains = applyPhaseLOS_d_3D(ch,pathGains,los);
                end

                % Apply beamforming to the channel output, this allows for
                % TXRU virtualization
                pathGains = h38901Channel.applyBeamforming(ch.SmallScale,pathGains,ch.TXRUVirtualization);
                % ---------------------------------------------------------

                % Ensure that path gains and sample times span at least one
                % slot
                [pathGains,sampleTimes] = spanSlot(obj,bsID,ueID,pathGains,sampleTimes);

            else % LOS ray only
                
                pathGains = h38901Channel.noFastFadingPathGains(ch);
                if (ch.LargeScale.TransmitAndReceiveSwapped())
                    pathGains = permute(pathGains,[1 2 4 3]);
                end
                if (obj.PhaseLOS_d_3D)
                    pathGains = applyPhaseLOS_d_3D(ch,pathGains,true);
                end
                ch.SmallScale.TransmitAndReceiveSwapped = ch.LargeScale.TransmitAndReceiveSwapped();
                pathGains = h38901Channel.applyBeamforming(ch.SmallScale,pathGains,ch.TXRUVirtualization);
                pathDelays = ch.PathDelays;
                pathFilters = ch.PathFilters;
                sampleTimes = 0;

            end

            % For full PHY, apply small scale channel to packet data
            if (~packet.Abstraction)
                if (~isequal(pathFilters,1))
                    % Channel filtering is required as path filters are not
                    % a unit scalar
                    packet.Data = channelFiltering(ch,packet.Data,pathGains,sampleTimes);
                    T = size(packet.Data,1);
                    packet.Duration = T / ch.SmallScale.SampleRate;
                else
                    % No channel filtering is required, channel is a matrix
                    % between transmit and receive antennas
                    H = permute(pathGains,[3 4 1 2]);
                    packet.Data = packet.Data * H;
                end
            end

            % -------------------------------------------------------------
            % TR 38.901 Section 7.5 Step 12
            % Update packet power with large scale channel effects
            PLdB = ch.LargeScale.execute(packet.TransmitterPosition,rxinfo.Position,packet.CenterFrequency);
            packet.Power = packet.Power - PLdB;

            % For full PHY, apply large scale channel to packet data
            if (~packet.Abstraction)
                packet.Data = packet.Data * db2mag(-PLdB);
            end
            % -------------------------------------------------------------

            % Update the channel metadata in the packet
            packet.Metadata.Channel.PathGains = pathGains;
            packet.Metadata.Channel.PathDelays = pathDelays;
            packet.Metadata.Channel.PathFilters = pathFilters;
            packet.Metadata.Channel.SampleTimes = sampleTimes;

        end

    end

    methods (Static, Access=public)

        function [channel,chinfo] = createChannelLink(chcfg)
        % [CHANNEL,CHINFO] = createChannelLink(CHCFG) creates a single 
        % channel link for channel link configuration structure CHCFG.

            % If no inputs are provided, return an empty channel and empty
            % channel information (this helps with building arrays having
            % the appropriate structure fields)
            if (nargin==0)
                channel = newChannel();
                chinfo = newChannelInfo();
                return;
            end

            % Provide default values for configuration parameters
            chcfg = setDefaults(chcfg,channelLinkDefaults(chcfg));

            % Validate ScenarioExtents
            if (isempty(chcfg.ScenarioExtents))
                error('nr5g:h38901Channel:EmptyExtents','ScenarioExtents must not be empty.');
            end

            % Validate that the UE position is inside the system boundary
            % (required for correlated LSPs and/or spatial consistency)
            if (~isnan(chcfg.InterSiteDistance) || spatialConsistency(chcfg))
                checkUEPosition(chcfg.ScenarioExtents,chcfg.UEPosition);
            end

            % Prepare the TR 38.901 Section 7.5 fast fading model
            [channel,chinfo,chcfg,siteuers,SCRVs] = fastFadingChannelModel(chcfg);

            % Section 7.6.9 "Absolute time of arrival" for InF scenarios
            if (startsWith(chcfg.Scenario,"InF") && chcfg.AbsoluteTOA)
                [channel,chinfo] = absoluteTOA(chcfg,channel,chinfo,siteuers,SCRVs);
            end

            % Store the channel link configuration in the channel output            
            channel.ChannelConfiguration = chcfg;

        end

        function pathGains = applyBeamforming(cdl,pathGains,virt)
        % PATHGAINS = applyBeamforming(CDL,PATHGAINS,VIRT) applies
        % beamforming (TXRU virtualization) to the set of path gains
        % PATHGAINS for nrCDLChannel object CDL and TXRU virtualization
        % parameters VIRT, as defined in TR 36.897 Section 5.2.2 TXRU
        % virtualization model option-1B. VIRT is a structure containing
        % the following fields:
        % K    - Vertical weight vector length
        % Tilt - Tilting angle in degrees
        % L    - Horizontal weight vector length
        % Pan  - Panning angle in degrees

            persistent c;
            if (isempty(c))
                c = physconst('LightSpeed');
            end

            % TR 36.897 Section 5.2.2 TXRU virtualization model option-1B
            K = virt.K;
            L = virt.L;
            if (K > 1 || L > 1)

                k = 1:K;
                lambda = c / cdl.CarrierFrequency;
                d_V = cdl.TransmitAntennaArray.ElementSpacing(1) * lambda;
                l = 1:L;
                d_H = cdl.TransmitAntennaArray.ElementSpacing(2) * lambda;

                theta_etilt = virt.Tilt;
                W = 1 / sqrt(K) * exp(-1i * 2 * pi / lambda * (k-1) * d_V * cosd(theta_etilt));
                theta = virt.Pan;
                V = 1 / sqrt(L) * exp(-1i * 2 * pi / lambda * (l-1) * d_H * sind(theta));

                % 'Ct' is the number of antennas in the channel (that is,
                % the non-virtual antenna count). 'Nt' is the number of
                % antennas that the link is aware of (that is, the virtual
                % antenna count)
                if (cdl.TransmitAndReceiveSwapped)
                    pathGains = permute(pathGains,[1 2 4 3]);
                end
                [Ncs,Np,Ct,Nr] = size(pathGains);
                Nt = Ct / (K*L);

                % Apply virtualization to the path gains
                X = kron(eye(Nt),kron(V,W));
                pathGains = permute(reshape(X * reshape(permute(pathGains,[3 4 1 2]),Ct,[]),[Nt Nr Ncs Np]),[3 4 1 2]);
                if (cdl.TransmitAndReceiveSwapped)
                    pathGains = permute(pathGains,[1 2 4 3]);
                end

            end

        end

        function [autoCorrMatrices,firstCoord] = createAutoCorrMatrices(rs,minpos,maxpos,distances)
            % [AUTOCORRMATRICES,FIRSTCOORD] =
            % createAutoCorrMatrices(RS,MINPOS,MAXPOS,DISTANCES) returns a
            % 3-D array containing autocorrelation matrices,
            % AUTOCORRMATRICES, and the coordinate pair of the first
            % autocorrelaton matrix element, FIRSTCOORD. Each matrix (i.e.
            % plane) of AUTOCORRMATRICES is the autocorrelation matrix for
            % an element of DISTANCES, a vector of correlation distances.
            % RS is a RandomStream object used to generate the normal
            % random variables prior to spatial filtering. MINPOS and
            % MAXPOS are the coordinate pairs of the lower-left and
            % upper-right corners of the rectangular region for which the
            % autocorrelation matrices are defined.

            [autoCorrMatrices,firstCoord] = createAutoCorrMatricesLocal(rs,minpos,maxpos,distances);
        end

        function rvs = uniformAutoCorrRVs(autoCorrMatrix,firstCoord,pos)
            % RVS = uniformAutoCorrRVs(AUTOCORRMATRIX,FIRSTCOORD,POS)
            % creates uniformly distributed spatially correlated random
            % variables RVS given an autocorrelaton matrix, AUTOCORRMATRIX,
            % the coordinate pair of the first autocorrelaton matrix
            % element, FIRSTCOORD, and a N-by-2 matrix of UE coordinate POS
            % (each row is the [X Y] coordinate of a UE).

            % Get normally distributed spatially correlated RVs
            rvs = normalAutoCorrRVs(autoCorrMatrix,firstCoord,pos);
            
            % Transform to uniform distribution U(0,1) using the
            % probability integral transform
            rvs = probabilityIntegralTransform(rvs);
        end

        function offsets = wrappingOffsets(ISD,numCellSites,numSectors)
        % OFFSETS = wrappingOffsets(ISD,numCellSites,numSectors) returns distance
        % offsets OFFSETS for wrap-around calculations, according to Rec. ITU-R
        % M.2101-0 Attachment 2 to Annex 1, for a specified intersite distance ISD,
        % number of cell sites numCellSites and number of sector per cell site
        % numSectors.

            % Rec. ITU-R M.2101-0 Attachment 2 to Annex 1, for 19 cell
            % sites modified for the difference in the orientation of the
            % rings of 6 and 12 sites around the central site between the
            % ITU-R document and TR 38.901. Similar equations are derived
            % for 3 and 7 cell sites.

            if numCellSites == 3 && numSectors == 3
                offsets = [[0 -0.5 0.5 -0.5 0.5 -1 1]*sqrt(3)*ISD; [0 -1.5 1.5 1.5 -1.5 0 0]*ISD].';
            elseif numCellSites == 7
                offsets = [[0 -0.5 0.5 -1 1 -1.5 1.5]*sqrt(3)*ISD; [0 -2.5 2.5 2 -2 -0.5 0.5]*ISD].';
            elseif numCellSites == 19
                offsets = [[0 -1 1 -1.5 1.5 -2.5 2.5]*sqrt(3)*ISD; [0 -4 4 3.5 -3.5 -0.5 0.5]*ISD].';
            else
                % Throw an error if the number of cell sites is neither 3 (with trisectorization), 
                % nor 7, nor 19 when toroidal wrap-around modeling is enabled
                error('nr5g:h38901Channel:InvalidCellConfiguration','Toroidal wrap-around modeling does not support (%d) cell sites with (%d) sectors. Use 3 cell sites with trisectorization, or 7 or 19 cell sites.', numCellSites,numSectors);
            end

        end

        function [x,y] = sitePolygon(ISD)
        % [X,Y] = sitePolygon(ISD) returns vertices of the polygon that
        % forms the boundary of a site with specified intersite distance
        % ISD.

            hexang = 0:60:360;
            [x,y] = pol2cart(deg2rad(hexang),ISD/sqrt(3));

        end

        % Calculate path gains corresponding to LOS channel ray
        function pathgains = noFastFadingPathGains(varargin)
        % PATHGAINS = noFastFadingPathGains(CH) returns the path gains
        % array PATHGAINS for the LOS ray of specified channel
        % configuration structure CH. PATHGAINS is of size
        % 1-by-1-by-Nt-by-Nr where Nt is the number of transmit antennas
        % and Nr is the number of receive antennas.

            pathgains = noFastFadingPathGainsLocal(varargin{:});

        end

        % Apply spatially consistent mobility to a channel link
        function varargout = spatiallyConsistentMobility(channel,cfg,t,tx,rx)
        % [CHANNEL,UPDATED,LOS_SOFT] = spatiallyConsistentMobility(CHANNEL,CFG,T,TX,RX)
        % applies spatially consistent mobility to input CHANNEL according
        % to configuration structure CFG, current simulation time T,
        % transmitter information structure TX, and receiver information
        % structure RX. If the output UPDATED is true, then the output
        % CHANNEL is the updated channel after a mobility update. If the
        % output UPDATED is false, then the output CHANNEL is equal to the
        % input, and no mobility update has been applied. The output 
        % LOS_SOFT provides the soft LOS state.
        %
        % CHANNEL is a channel link structure, as returned by the 
        % createChannelLink object function. 
        %
        % CFG is a structure containing the following fields:
        % SpatialConsistency - ("ProcedureA", "ProcedureB")
        % UpdateDistance     - The change in BS-to-UE distance, in meters, 
        %                      that will trigger the spatial consistency 
        %                      update procedure. TR 38.901 Section 7.6.3.2 
        %                      states that this distance "should be limited
        %                      within 1 meter" (default 1.0)
        %
        % TX and RX are structures containing the following fields:
        % Position          - Position of node, specified as a real-valued
        %                     vector in Cartesian coordinates [x y z] in 
        %                     meters
        % Velocity          - Velocity (v) of node in the x-, y-, and 
        %                     z-directions, specified as a real-valued 
        %                     vector of the form [vx vy vz] in meters per 
        %                     second
        % RotationVelocity  - Rotation velocity (ω) of the node in the 
        %                     direction of positive bearing angle α,
        %                     downtilt angle β, and slant angle γ,
        %                     specified as a real-valued vector of the form
        %                     [ωα; ωβ; ωγ] in rotations per minute (RPM).
        %                     These values are used to update the
        %                     orientation of the antenna array of the node

            if (isempty(channel.Mobility.manager))
                cfg = setDefaults(cfg,struct(UpdateDistance=1.0));
                if (~isSpatiallyConsistentMobility(cfg))
                    error('nr5g:h38901Channel:InvalidSpatialConsistency',"SpatialConsistency must be 'ProcedureA' or 'ProcedureB' when calling h38901Channel.spatiallyConsistentMobility");
                end
                if (cfg.UpdateDistance < 0)
                    error('nr5g:h38901Channel:InvalidUpdateDistance','UpdateDistance must be nonnegative.');
                end
                channel.Mobility.manager = mobilityManager(cfg);
            end
            [varargout{1:nargout}] = channel.Mobility.manager.update(channel,t,tx,rx);

        end

    end

    % =====================================================================
    % private

    properties (SetAccess=private,Hidden)

        InterfererHasSmallScale = false;
        PhaseLOS_d_3D = true;
        InterfererSameLinkEnd = false;

    end

    properties (Access=private)

        theLinkToChannelMap;
        thePathLossConfig;
        theUEtoBSMap;
        theUEMap;
        theBSMap;
        theRecordedSitePositions;
        theRecordedUEPositions;
        theScenarioInfo = struct(MaxBSID=NaN,Wrapping=false,SpatialConsistency=false,InterSiteDistance=NaN,NodeSiz=[1 1 1],MaxLinkID=NaN);
        theSitePositions;

    end

end

%% ========================================================================
%  local functions related to wirelessNetworkSimulator
%  ========================================================================

function [channel,bsID,ueID,linkind] = getSLSChannel(obj,rxinfo,packet)

    % Get BS and UE IDs for this link and establish if the link is uplink
    % (that is, the UE is the transmitter)
    linkID = [packet.TransmitterID rxinfo.ID];
    [bsID,ueID,isUplink] = getBSandUE(obj,linkID);

    % If either the BS or UE ID is undefined, signifying a BS-to-BS or
    % UE-to-UE link, return an empty channel if InterfererSameLinkEnd=false
    if (~obj.InterfererSameLinkEnd && (isnan(bsID) || isnan(ueID)))
        channel = [];
        linkind = [];
        return;
    end

    % When InterfererSameLinkEnd=true and the link is BS-to-BS and
    % corresponds to two sectors for the same site, return an empty channel
    if (isnan(ueID))
        nodesubs = cat(1,obj.theBSMap(linkID).NodeSubs);
        sitesubs = nodesubs(:,1);
        if (sitesubs(1)==sitesubs(2))
            channel = [];
            linkind = [];
            return;
        end
    end

    % ---------------------------------------------------------------------
    % TR 38.901 Section 7.5 Steps 2 - 10
    % Get or create the channel for the appropriate link direction
    if (isnan(isUplink))
        % BS-to-BS or UE-to-UE link, do not allow reciprocity as nodes can
        % have different values for properties such as TXRUVirtualization
        % and d_2D_in
        thisLinkID = linkID;
        otherLinkID = linkID;
    elseif (isUplink)
        thisLinkID = [ueID bsID];
        otherLinkID = [bsID ueID];
    else
        thisLinkID = [bsID ueID];
        otherLinkID = [ueID bsID];
    end
    [channel,isUplink,linkind] = getSLSChannelLink(obj,thisLinkID,otherLinkID,linkID,rxinfo,packet,isUplink);

    % ---------------------------------------------------------------------

    % Now that a channel is selected, ensure that the channel is set for
    % the correct link direction. Note that channels are always created in
    % the downlink direction, so uplink links must always be configured as
    % reciprocal links (that is, with transmit and receive swapped). For
    % BS-to-BS or UE-to-UE links reciprocity is not used as nodes can have
    % different values for properties such as TXRUVirtualization and
    % d_2D_in
    if (~isnan(isUplink) && xor(isUplink,channel.LargeScale.TransmitAndReceiveSwapped()))
        channel.LargeScale.swapTransmitAndReceive();
        if (~isempty(channel.SmallScale) && ~isstruct(channel.SmallScale))
            swapTransmitAndReceive(channel.SmallScale);
        end
    end

end

% Get or create the channel for the specified link direction
function [channel,isUplink,linkindout] = getSLSChannelLink(obj,thisLinkID,otherLinkID,linkID,rxinfo,packet,isUplink)

    % If a channel exists in this link direction
    linkind = linkIndicesForID(obj,thisLinkID);
    if (isKey(obj.theLinkToChannelMap,linkind))

        % Use it
        channel = obj.theLinkToChannelMap(linkind);
        linkindout = linkind;

    else % a channel does not exist in this link direction

        % If a channel exists in the other link direction
        linkind = linkIndicesForID(obj,otherLinkID);
        if (isKey(obj.theLinkToChannelMap,linkind))

            % If the center frequencies of that channel and this packet
            % match
            ch = obj.theLinkToChannelMap(linkind);
            sameFrequency = isequal(ch.CenterFrequency,packet.CenterFrequency);
            chAnts = [ch.NumTransmitAntennas ch.NumReceiveAntennas];
            nodeAnts = [rxinfo.NumReceiveAntennas packet.NumTransmitAntennas];
            sameAnts = all(chAnts == nodeAnts);
            if (sameFrequency && sameAnts)

                % TDD and the same antenna count, the channel can be
                % re-used for this link direction - use it
                channel = ch;
                linkindout = linkind;

            else

                % FDD and/or a different antenna count, the channel cannot
                % be re-used for this link direction - create a new channel
                channel = [];

            end

        else

            % A channel does not exist for either link direction - create a
            % new channel for this link direction
            channel = [];

        end

    end

    % ---------------------------------------------------------------------
    % TR 38.901 Section 7.5 Steps 2 - 10
    if (isempty(channel))
        [channel,isUplink,linkindout] = createSLSChannelLink(obj,linkID,packet,isUplink);
    end
    % ---------------------------------------------------------------------

end

% Create a channel link for the SLS
function [channel,isUplink,linkind] = createSLSChannelLink(obj,linkID,packet,isUplink)

    % Establish if this link needs a CDL channel (that is, between two
    % nodes that are attached or InterfererHasSmallScale is true).
    % Establish the BS and UE IDs
    [attached,bsID,ueID] = isAttached(obj,linkID);
    hasSmallScale = true;
    if (obj.InterfererHasSmallScale)
        fastFading = true;
    else
        fastFading = attached;
    end

    % Prepare vectors describing node counts and subscripts, and get the
    % number of transmit and receive antennas
    if (isnan(isUplink))
        % BS-to-BS or UE-to-UE link
        if (~isnan(bsID))
            % BS-to-BS
            theMap = obj.theBSMap;
        else
            % UE-to-UE
            theMap = obj.theUEMap;
        end
        % 'f' and 'r' are the forward and reverse nodes, i.e. adopt the
        % convention that the communication path linkID(1) -> linkID(2) is
        % the forward link
        f = theMap(linkID(1));
        r = theMap(linkID(2));
        % Create 'nodeSiz' and 'nodeSubs', allowing for UEs in the forward
        % link and BSs (including sectors) in the reverse link
        nodeSiz = obj.theScenarioInfo.NodeSiz;
        [nodeSiz,nodeSubs] = interferingNodeVectors(bsID,f,r,nodeSiz);
        % Get number of transmit and receive antennas
        Nt = f.Node.NumTransmitAntennas;
        Nr = r.Node.NumReceiveAntennas;
    else
        % BS-to-UE or UE-to-BS link
        BS = obj.theBSMap(bsID);
        UE = obj.theUEMap(ueID);
        % Create 'nodeSiz'
        nodeSiz = obj.theScenarioInfo.NodeSiz;
        % Create 'nodeSubs'
        nodeSubs = [BS.NodeSubs(1:2) UE.NodeSubs(3)];
        % Get number of transmit and receive antennas
        if (isUplink)
            % Note that channels are in the downlink direction and operate
            % as reciprocal links for the uplink direction. Therefore, the
            % antenna counts below are for the downlink
            Nt = BS.Node.NumReceiveAntennas;
            Nr = UE.Node.NumTransmitAntennas;
        else % downlink
            Nt = BS.Node.NumTransmitAntennas;
            Nr = UE.Node.NumReceiveAntennas;
        end
    end

    % Update the scenario in the path loss configuration, in case the
    % object property has changed    
    obj.thePathLossConfig.Scenario = obj.Scenario;

    % ---------------------------------------------------------------------
    % TR 38.901 Section 7.5 Steps 2 - 10
    % Create channel from low-level parameters
    ISD = obj.theScenarioInfo.InterSiteDistance;
    chcfg = struct();
    chcfg.Seed = obj.Seed;
    chcfg.Scenario = string(obj.Scenario);
    chcfg.InterSiteDistance = ISD;
    % If site positions have changed since the last call to
    % createSLSChannelLink, or if this is the first call
    sitePositions = obj.theRecordedSitePositions;
    if (~isequal(obj.theSitePositions,sitePositions))
        % Store those positions and pass them to
        % h38901Channel.createChannelLink in order to reset
        % spatially consistent random variables (SCRVs)
        chcfg.SitePositions = sitePositions;
        obj.theSitePositions = sitePositions;
    else
        % Otherwise pass empty site positions
        chcfg.SitePositions = [];
    end
    chcfg.HasSmallScale = hasSmallScale;
    chcfg.FastFading = fastFading;
    chcfg.NodeSubs = nodeSubs;
    chcfg.NodeSiz = nodeSiz;
    if obj.InterfererSameLinkEnd
        chcfg.NumCellSites = obj.theScenarioInfo.NodeSiz(1);
    end
    chcfg.NumTransmitAntennas = Nt;
    chcfg.NumReceiveAntennas = Nr;
    if (isnan(isUplink))
        % BS-to-BS or UE-to-UE link
        chcfg.BSPosition = f.Node.Position;
        chcfg.UEPosition = r.Node.Position;
        chcfg.SampleRate = f.OFDMInfo.SampleRate;
        if (~isnan(bsID))
            % BS-to-BS link, treat receiving BS like an outdoor UE
            chcfg.TXRUVirtualization = f.TXRUVirtualization;
            chcfg.TransmitArrayOrientation = f.TransmitArrayOrientation;
            chcfg.n_fl = 0;
            chcfg.d_2D_in = 0;
            chcfg.ReceiveArrayOrientation = [];
        else
            % UE-to-UE link, treat transmitting UE like a BS with no TXRU
            % virtualization
            chcfg.TXRUVirtualization = struct(K=1,Tilt=0,L=1,Pan=0);
            chcfg.TransmitArrayOrientation = [];
            chcfg.n_fl = r.n_fl;
            chcfg.d_2D_in = r.d_2D_in;
            chcfg.ReceiveArrayOrientation = r.ReceiveArrayOrientation;
        end
    else
        % BS-to-UE or UE-to-BS link
        chcfg.SampleRate = BS.OFDMInfo.SampleRate;
        chcfg.TXRUVirtualization = BS.TXRUVirtualization;
        chcfg.TransmitArrayOrientation = BS.TransmitArrayOrientation;
        chcfg.BSPosition = BS.Node.Position;
        chcfg.UEPosition = UE.Node.Position;
        chcfg.n_fl = UE.n_fl;
        chcfg.d_2D_in = UE.d_2D_in;
        chcfg.ReceiveArrayOrientation = UE.ReceiveArrayOrientation;
    end
    chcfg.CenterFrequency = packet.CenterFrequency;
    chcfg.Wrapping = obj.theScenarioInfo.Wrapping;
    chcfg.SpatialConsistency = obj.theScenarioInfo.SpatialConsistency;
    chcfg.OfficeType = obj.OfficeType;
    if (~isempty(obj.ScenarioExtents))
        extents = obj.ScenarioExtents;
    else
        if (startsWith(chcfg.Scenario,"InF"))
            % Defaulted inside h38901Channel.createChannelLink from
            % chCfg.HallSize
            extents = [];
        else
            if (any(chcfg.Scenario==["UMi" "UMa" "RMa"]))
                % Get polygons that are the boundaries for each site
                [sitex,sitey] = h38901Channel.sitePolygon(ISD);
                % Get bounding box of the union of the site polygons
                sysx = sitePositions(:,1) + sitex;
                sysy = sitePositions(:,2) + sitey;
                minpos = [min(sysx, [], 'all'), min(sysy, [], 'all')];
                maxpos = [max(sysx, [], 'all'), max(sysy, [], 'all')];
            else % InH
                allnodepos = [sitePositions; obj.theRecordedUEPositions];
                % Note: the bounding box around the node positions is
                % extended by 1 meter on all sides here to avoid precision
                % issues when determining if nodes lie inside the scenario
                % extents
                maxpos = max(allnodepos(:,1:2),[],1) + 1;
                minpos = min(allnodepos(:,1:2),[],1) - 1;
            end
            extents = [minpos maxpos-minpos];
        end
    end
    chcfg.ScenarioExtents = extents;
    chcfg.HallSize = obj.HallSize;
    chcfg.ClutterSize = obj.ClutterSize;
    chcfg.ClutterDensity = obj.ClutterDensity;
    chcfg.ClutterHeight = obj.ClutterHeight;
    chcfg.AbsoluteTOA = obj.AbsoluteTOA;
    chcfg.LOSProbability = obj.LOSProbability;
    channel = h38901Channel.createChannelLink(chcfg);
    channel.LargeScale = channel.LargeScale(obj.thePathLossConfig);
    % ---------------------------------------------------------------------

    % For full PHY, construct channel filters if required (that is, if path
    % delays are non-zero, due to small scale fading and/or absolute TOA)
    if (~packet.Abstraction && ~isequal(channel.PathDelays,0))
        f = comm.ChannelFilter();
        f.SampleRate = channel.SmallScale.SampleRate;
        f.PathDelays = channel.PathDelays;
        f.NormalizeChannelOutputs = false;
        r = clone(f);
        fInfo = info(f);
        channel.PathFilters = fInfo.ChannelFilterCoefficients;
        channel.ChannelFilter = f;
        channel.ChannelFilterReciprocal = r;
    end

    % Record the channel
    linkind = linkIndicesForID(obj,linkID);
    obj.theLinkToChannelMap(linkind) = channel;

    if (isnan(isUplink))
        % BS-to-BS or UE-to-UE link, treat as downlink i.e. reciprocity
        % will not be used
        isUplink = false;
    end

end

function [nodeSiz,nodeSubs] = interferingNodeVectors(bsID,f,r,nodeSiz)

    oldSiz = nodeSiz;

    % 'nUEasSite' is the number of extra sites needed to account for all
    % UEs treated as sites
    nUEasSite = nodeSiz(3);

    % 'nBSasUE' is the number of extra UEs needed to account for all BS
    % sites and sectors treated as UEs
    nBSasUE = prod(nodeSiz(1:2));

    % Update site and UE entries in 'nodeSiz'
    nodeSiz([1 3]) = nodeSiz([1 3]) + [nUEasSite nBSasUE];

    % Create 'nodeSubs'
    if (~isnan(bsID))

        % BS-to-BS link

        % 'thisBSasUE' combines 'r' BS site and sector indices into a UE
        % index
        i = r.NodeSubs(1);
        j = r.NodeSubs(2);
        thisBSasUE = (i-1)*oldSiz(2) + j;

        % 'nodeSubs' consists of site and sector index of the 'f' BS, and
        % index of the 'r' BS treated as a UE
        nodeSubs = [f.NodeSubs(1:2) oldSiz(3)+thisBSasUE];

    else

        % UE-to-UE link

        % 'thisUEasSite' uses 'f' UE index as a site index
        thisUEasSite = f.NodeSubs(3);

        % 'nodeSubs' consists of index of the 'f' UE treated as a site,
        % sector set to 1, and the index of the 'r' UE
        nodeSubs = [oldSiz(1)+thisUEasSite 1 r.NodeSubs(3)];

    end

end

% Apply channel filtering to signal 'x' using path gains 'pg' and sample
% times 't'
function y = channelFiltering(ch,x,pg,t)

    if (~ch.SmallScale.TransmitAndReceiveSwapped)
        f = ch.ChannelFilter;
    else
        f = ch.ChannelFilterReciprocal;
    end
    fInfo = info(f);
    pathFilters = fInfo.ChannelFilterCoefficients;
    Nh = size(pathFilters,2);
    Nt = size(x,2);
    x = [x; zeros([Nh Nt])];
    insize = size(x);
    if (~isstruct(ch.SmallScale))
        sampleDensity = ch.SmallScale.SampleDensity;
    else
        sampleDensity = 1;
    end
    t = t - t(1);
    outputtype = class(x);
    % If the channel filter is not locked (i.e. hasn't yet been used)
    if (~isLocked(f))
        % Determine the best filter option for the current input
        sizg = size(pg,1:4);
        T = insize(1);
        p = wireless.internal.channelmodels.getFilterPolicy(outputtype,[sizg T]);
        % If the option is to filter waveform sections where each section
        % corresponds to a first dimension element of 'pg', configure the
        % channel filter for efficient execution of this case
        f.OptimizeScalarGainFilter = (p.FilterOption==2);
    end
    filterOption = 1 + f.OptimizeScalarGainFilter;
    y = wireless.internal.channelmodels.smartChannelFiltering(x,f,f.SampleRate,pg,insize,sampleDensity,t,outputtype,filterOption);

end

% Get BS and UE IDs for this link and establish if the link is uplink
function [bsID,ueID,isUplink] = getBSandUE(obj,linkID)

    maxBSID = obj.theScenarioInfo.MaxBSID;
    txIsUE = isKey(obj.theUEtoBSMap,linkID(1));
    if (txIsUE)
        ueID = linkID(1);
        bsID = linkID(2);
        isUplink = true;
    else
        ueID = linkID(2);
        bsID = linkID(1);
        isUplink = false;
        if (ueID <= maxBSID)
            ueID = NaN;
            isUplink = NaN;
            return;
        end
    end
    if (bsID > maxBSID)
        bsID = NaN;
        isUplink = NaN;
    end

end

% Establish if a link is between two nodes that are attached, and also
% establish the BS and UE IDs
function [c,bsID,ueID] = isAttached(obj,linkID)

    [bsID,ueID] = getBSandUE(obj,linkID);
    if (isnan(bsID) || isnan(ueID))
        % BS-to-BS or UE-to-UE link, not attached
        c = false;
    else
        c = (obj.theUEtoBSMap(ueID)==bsID);
    end

end

% Get linear index (used as a dictionary hash) from a link ID pair
function ind = linkIndicesForID(obj,linkID)

    ind = (linkID(1)-1)*obj.theScenarioInfo.MaxLinkID + linkID(2);

end

% Record list of gNBs, connections between gNBs and UEs, and other scenario
% information that can be established from node names or channel
% configuration structure
function recordNodes(obj,gNBs,UEs,chCfg)

    noNodeWarning = @(s)warning(['nr5g:h38901Channel:No' s 'Nodes'],'No nr%s nodes have been added to the wirelessNetworkSimulator. Call h38901Channel/connectNodes after all nodes have been added to the wirelessNetworkSimulator.',s);
    if (~isempty(gNBs))
        maxBSID = max([gNBs.ID]);
        obj.theScenarioInfo.MaxBSID = maxBSID;
        obj.theScenarioInfo.NodeSiz(1) = maxBSID;
        obj.theScenarioInfo.MaxLinkID = maxBSID;
    else
        maxBSID = [];
        noNodeWarning('GNB');
    end

    if (~isempty(UEs))
        minUEID = min([cat(2,UEs{:}).ID]);
        maxUEID = max([cat(2,UEs{:}).ID]);
        if (~isempty(maxUEID))
            obj.theScenarioInfo.NodeSiz(3) = maxUEID;
        end
    else
        minUEID = [];
        noNodeWarning('UE');
    end

    if (~isempty(maxBSID) && ~isempty(minUEID))
        if (minUEID < maxBSID)
            error('nr5g:h38901Channel:UENodeFirst','The minimum UE ID (%d) is less than maximum gNB ID (%d). Create all gNB nodes before creating any UE nodes, to ensure that all UE IDs are greater than all gNB IDs.',minUEID,maxBSID);
        end
    end

    obj.theScenarioInfo.SpatialConsistency = obj.SpatialConsistency;

    chCfgSpecified = ~isempty(chCfg);

    if (chCfgSpecified)
        if (isfield(chCfg,"TransmitArrayOrientation"))
            [~,~,ic] = unique(chCfg.TransmitArrayOrientation.','rows');
            chCfg.Sector = ic.';
        end
        if (~any(string(obj.Scenario)==["UMi" "UMa" "RMa"]))
            for field = ["d_2D_in" "n_fl"]
                if (isfield(chCfg,field))
                    chCfg = rmfield(chCfg,field);
                end
            end
        end
    end

    firstUE = true;
    for i = 1:numel(gNBs)

        bs = gNBs(i);
        site = getFieldColumnOrDefault(chCfg,"Site",i,bs.ID);
        sector = getFieldColumnOrDefault(chCfg,"Sector",i,1);

        nvps = getNVPs(bs.Name);
        if (size(nvps,2)==2)
            subs = getNodeIndices(nvps);
        else
            subs = [site sector];
        end
        ofdmInfo = nrOFDMInfo(bs.NumResourceBlocks,bs.SubcarrierSpacing/1e3);

        virt = getFieldColumnOrDefault(chCfg,"TXRUVirtualization",i,struct(K=1,Tilt=0,L=1,Pan=0));
        txorientation = getFieldColumnOrDefault(chCfg,"TransmitArrayOrientation",i,[]);
        obj.theBSMap(bs.ID) = struct(Node=bs,NodeSubs=subs,OFDMInfo=ofdmInfo,TXRUVirtualization=virt,TransmitArrayOrientation=txorientation);

        for j = 1:numel(UEs{i})

            ue = UEs{i}(j);

            nvps = getNVPs(ue.Name);
            if (isCellular(obj.Scenario) && size(nvps,2)==2)
                [subs,siz] = getNodeIndices(nvps);
                if (firstUE)
                    obj.theScenarioInfo.InterSiteDistance = str2double(getV(nvps,"ISD"));
                    obj.theScenarioInfo.Wrapping = str2double(getV(nvps,"Wrapping"));
                    spatialConsistency = getV(nvps,"SpatialConsistency");
                    d = str2double(spatialConsistency);
                    if (~isnan(d))
                        spatialConsistency = d;
                    end
                    obj.theScenarioInfo.SpatialConsistency = spatialConsistency;
                    obj.theScenarioInfo.NodeSiz = siz;
                else
                    obj.theScenarioInfo.NodeSiz = max([obj.theScenarioInfo.NodeSiz; subs],[],1);
                end
                n_fl = str2double(getV(nvps,"n_fl"));
                d_2D_in = str2double(getV(nvps,"d_2D_in"));
                rxorientation = [];
            else
                if (chCfgSpecified)
                    subs = [site sector j];
                else
                    subs = [site sector ue.ID];
                end
                if (chCfgSpecified && firstUE)
                    obj.theScenarioInfo.InterSiteDistance = obj.InterSiteDistance;
                    obj.theScenarioInfo.Wrapping = obj.Wrapping;
                    siz = [max(chCfg.Site) max(chCfg.Sector) max(cellfun(@numel,UEs))];
                    obj.theScenarioInfo.NodeSiz = siz;
                end
                obj.theScenarioInfo.NodeSiz = max([obj.theScenarioInfo.NodeSiz; subs],[],1);
                if (chCfgSpecified)
                    default_n_fl = round((ue.Position(3)-1.5)/3 + 1);
                    el = find(chCfg.UEIDs==ue.ID,1,'first');
                else
                    default_n_fl = 0;
                    el = [];
                end
                n_fl = getFieldColumnOrDefault(chCfg,"n_fl",el,default_n_fl);
                d_2D_in = getFieldColumnOrDefault(chCfg,"d_2D_in",el,0);
                rxorientation = getFieldColumnOrDefault(chCfg,"ReceiveArrayOrientation",el,[]);
            end

            firstUE = false;

            obj.theUEMap(ue.ID) = struct(Node=ue,NodeSubs=subs,OFDMInfo=ofdmInfo,n_fl=n_fl,d_2D_in=d_2D_in,ReceiveArrayOrientation=rxorientation);
            obj.theScenarioInfo.MaxLinkID = max([obj.theScenarioInfo.MaxLinkID ue.ID]);

            obj.theUEtoBSMap(ue.ID) = bs.ID;

        end

    end

    obj.theRecordedSitePositions = getNodePositions(obj.theBSMap);
    obj.theRecordedUEPositions = getNodePositions(obj.theUEMap);

end

function p = getNodePositions(m)

    v = values(m);
    if (~isempty(v))
        p = cat(1,cat(1,v.Node).Position);
    else
        p = zeros(0,3);
    end

end

function v = getFieldColumnOrDefault(s,f,i,d)

    if (isfield(s,f))
        v = s.(f);
        if (size(v,2) > 1)
            v = v(:,i);
        end
    else
        v = d;
    end

end

% Get name-value pairs from a node Name property
function nvp = getNVPs(name)

    nvp = arrayfun(@(x)strsplit(x,'='),strsplit(name,',').','UniformOutput',false);
    nvp = cat(1,nvp{:});

end

% Get value (from the output of getNVPs) corresponding to a name 
function v = getV(nvps,n)

    v = nvps(nvps(:,1)==n,2);

end

% Get indices of site / sector / UE from the name-value pairs extracted
% from the node name
function [subs,siz] = getNodeIndices(nvps)

    names = ["Site" "Sector"];
    if (any(nvps(:,1)=="UE"))
        names = [names "UE"];
    end
    y = arrayfun(@(x)strsplit(x,'/'),arrayfun(@(n)getV(nvps,n),names).','UniformOutput',false);
    y = cat(1,y{:});

    subs = str2double(y(:,1)).';
    siz = str2double(y(:,2)).';

end

% Calculation of BS antenna array panel dimensions based on number of
% CSI-RS ports, using the same logic as for the PMI panel dimensions in
% the SLS
function [Ng,M,N] = panelDimensions(numCSIRSPorts)

    Ng = 1;
    if (numCSIRSPorts > 2)
        panelConfigs = [4 2 1; 8 2 2; 8 4 1; 12 3 2; 12 6 1; 16 4 2; ...
            16 8 1; 24 4 3; 24 6 2; 24 12 1; 32 4 4; 32 8 2; 32 16 1];
        idx = find(panelConfigs(:,1)==numCSIRSPorts,1,'first');
        M = panelConfigs(idx,2);
        N = panelConfigs(idx,3);
    else
        M = 1;
        N = 1;
    end

end

function ch = packetMobilitySLS(obj,ch,packet,rxinfo)

    % Set up mobility configuration structure from object properties
    cfg.SpatialConsistency = obj.SpatialConsistency;
    cfg.UpdateDistance = obj.UpdateDistance;

    % Get input time, equal to packet start time
    t = packet.StartTime;

    % Set up transmitter info structure from transmitted packet in SLS
    tx = struct();
    tx.Position = packet.TransmitterPosition;
    tx.Velocity = packet.TransmitterVelocity;

    % Set up receiver info structure from receiver info in SLS
    rx = struct();
    rx.Position = rxinfo.Position;
    rx.Velocity = rxinfo.Velocity;

    % SLS does not provide rotation information. The rotation velocity is
    % set to zero. The current direction of travel is stored. The array
    % orientations will be updated after the mobility update in order to
    % keep the antenna array orientations fixed relative to the direction
    % of travel
    tx.RotationVelocity = [0; 0; 0];
    rx.RotationVelocity = [0; 0; 0];
    old_direction = ch.Mobility.UTDirectionOfTravel;

    % Perform mobility update
    if (isempty(ch.Mobility.manager))
        ch.Mobility.manager = mobilityManager(cfg);
    end
    chSpatialConsistency = ch.ChannelConfiguration.SpatialConsistency;
    if (~strcmpi(chSpatialConsistency,cfg.SpatialConsistency))
        error('nr5g:h38901Channel:ChangedSpatialConsistencySLS','Spatial consistency configured when initially dropping UEs, h38901Scenario.SpatialConsistency (%s), must match SpatialConsistency configured in h38901Channel (%s).',chSpatialConsistency,cfg.SpatialConsistency);
    end
    [ch,updated] = ch.Mobility.manager.update(ch,t,tx,rx);

    % Update antenna array orientations in order to keep them fixed
    % relative to the direction of travel
    if (updated)
        ch = updateArrayOrientations(ch,old_direction);
    end

end

function ch = updateArrayOrientations(ch,old_direction)

    ss = ch.SmallScale;

    delta_direction = ch.Mobility.UTDirectionOfTravel - old_direction;

    delta_direction_rx = delta_direction(:,1);
    ss.ReceiveArrayOrientation = updateArrayOrientation(ss.ReceiveArrayOrientation,delta_direction_rx);

    delta_direction_tx = delta_direction(:,2);
    ss.TransmitArrayOrientation = updateArrayOrientation(ss.TransmitArrayOrientation,delta_direction_tx);

    ch.SmallScale = ss;

end

function ao = updateArrayOrientation(ao,delta_azel)

    % The bearing angle of the array is updated by the change in azimuth
    % angle of the direction of travel
    delta_phi = delta_azel(1);
    ao(1) = ao(1) + delta_phi;

    % The downtilt angle of the array is updated by the change in elevation
    % angle of the direction of travel
    delta_theta = delta_azel(2);
    ao(2) = ao(2) + delta_theta;

end

%% ========================================================================
%  local functions independent of wirelessNetworkSimulator
%  ========================================================================

% Create a new channel structure
function c = newChannel()

    c = struct();
    c.CenterFrequency = NaN;
    c.NumTransmitAntennas = NaN;
    c.NumReceiveAntennas = NaN;
    c.LargeScale = [];
    c.SmallScale = [];
    c.TXRUVirtualization = struct();
    c.PathFilters = [];
    c.NodeSubs = [];
    c.NodeSiz = [];
    c.PathDelays = [];
    c.ChannelFilter = [];
    c.ChannelFilterReciprocal = [];
    c.ChannelConfiguration = [];
    c.Mobility = newMobility();

end

% Create a new channel information structure
function c = newChannelInfo()

    c = struct();
    c.O2I = NaN;
    c.pLOS = NaN;
    c.LOS = NaN;
    c.h_E = NaN;
    c.SF = NaN;
    c.K = NaN;
    c.d_3D = NaN;

end

% Create new mobility structure
function m = newMobility()

    m = struct();
    m.tau_prime = [];
    m.tau = [];
    m.DS = [];
    m.r_tau = [];
    m.zeta = [];
    m.K = [];
    m.indoor = [];
    m.Deltatau = [];
    m.UTDirectionOfTravel = [];
    m.manager = [];

end

% Set object properties from name-value arguments
function setProperties(obj,varargin)

    for i = 1:2:(nargin-1)
        n = varargin{i};
        v = varargin{i+1};
        obj.(n) = v;
    end

end

% Wrapper for large scale channel which handles reciprocity, caching for
% mobility and toroidal wrap-around
function obj = largeScaleWrapper(wrapping,ISD,numCellSites,numSectors,in)

    cache.txpos = [NaN NaN NaN];
    cache.rxpos = [NaN NaN NaN];
    cache.out = NaN;
    swapped = false;
    if (wrapping)
        offsets = h38901Channel.wrappingOffsets(ISD,numCellSites,numSectors);
    else
        offsets = [0 0];
    end

    fn = in;

    obj.execute = @execute;
    obj.swapTransmitAndReceive = @swapTransmitAndReceive;
    obj.TransmitAndReceiveSwapped = @getSwapped;

    function swapTransmitAndReceive()

        swapped = ~swapped;

    end

    function s = getSwapped()

        s = swapped;

    end

    function out = execute(txpos,rxpos,varargin)

        if (swapped)
            temp = rxpos;
            rxpos = txpos;
            txpos = temp;
        end

        if (any(txpos~=cache.txpos) || any(rxpos~=cache.rxpos))
            d = vecnorm((txpos(1:2) + offsets) - rxpos(1:2),2,2);
            [~,idx] = min(d);
            cache.out = fn(txpos + [offsets(idx,:) 0],rxpos,varargin{:});
            cache.txpos = txpos;
            cache.rxpos = rxpos;
        end

        out = cache.out;

    end

end

% Select a value or call a function based on the scenario
function x = scenarioSwitch(s,umi,uma,rma,indoorh,indoorf)

    if (s=="UMi")
        x = select(umi);
    elseif (s=="UMa")
        x = select(uma);
    elseif (s=="RMa")
        x = select(rma);
    elseif (startsWith(s,"InH"))
        indoorh = desuffix(s,indoorh,["M","O"]);
        x = select(indoorh);
    elseif (startsWith(s,"InF"))
        indoorf = desuffix(s,indoorf,["SL","DL","SH","DH","HH"]);
        x = select(indoorf);
    end

    function x = desuffix(s,x,suffixes)
        if (~isscalar(x))
            x = x{suffixes==extractAfter(s,'-')};
        end
    end

    function x = select(x)
        if (isa(x,'function_handle'))
            x = x();
        end
    end

end

% Outdoor-to-indoor penetration loss in dB, TR 38.901 Section 7.4.3
function PL = outdoorToIndoor(rs,s,d_2D_in,fc,sc)

    % For InH and InF, although TR 38.901 Section 7.4.3.1 says "only the
    % high-loss model is applicable to InF", it is assumed here that O2I is
    % not applicable for these scenarios
    hasO2I = isCellular(s);

    spatialConsistency = ~isempty(sc);
    indoor = (d_2D_in~=0);

    if (~hasO2I || ~indoor)

        PL = 0;

    else

        sub6GHz = (fc < 6e9);
        if (sub6GHz)

            % TR 38.901 Table 7.4.3-3 "single-frequency simulations <6 GHz"
            % NOTE: d_2D_in will already be established during scenario
            % building and passed in via the CHCFG input to the
            % connectNodes method. sigma_SF is handled in the
            % getShadowFading function in this file
            PL_tw = 20;
            PL_in = 0.5 * d_2D_in;
            sigma_P = 0;

        else

            % Low loss ratio from TR 38.901 Table 7.8-1 for UMi and UMa.
            % For RMa, TR 38.901 Section 7.4.3.1 says "only the low-loss
            % model is applicable to RMa"
            lowlossRatio = scenarioSwitch(s,0.5,0.5,1.0,NaN,NaN);
            if (spatialConsistency)
                % TR 38.901 Section 7.6.3.3 "building type is determined
                % using a spatially consistent uniform random variable"
                uepos = sc.UEPosition;
                rv = h38901Channel.uniformAutoCorrRVs(sc.O2I.AutoCorrMatrices(:,:,1),sc.O2I.FirstCoord,uepos(1:2));
            else
                rv = rs.rand;
            end
            lowloss = rv < lowlossRatio;

            % TR 38.901 Table 7.4.3-1. Penetration loss for wood is not
            % required because Table 7.4.3-2 is used below
            f = fc / 1e9;
            L_glass = 2 + 0.2*f;
            L_IIRglass = 23 + 0.3*f;
            L_concrete = 5 + 4*f;

            % TR 38.901 Table 7.4.3-2
            PL_in = 0.5 * d_2D_in;
            if (lowloss)
                PL_tw = 5 - 10*log10(0.3*10^(-L_glass/10) + 0.7*10^(-L_concrete/10));
                sigma_P = 4.4;
            else
                PL_tw = 5 - 10*log10(0.7*10^(-L_IIRglass/10) + 0.3*10^(-L_concrete/10));
                sigma_P = 6.5;
            end

        end

        if (spatialConsistency)
            if (sub6GHz)
                % For TR 38.901 Table 7.4.3-3 "single-frequency simulations
                % <6 GHz", sigma_P = 0 so 'rv' will have no effect. Also,
                % PL_tw has a fixed value so the use of different spatially
                % consistent random variables for the different building
                % types (low/high loss) is not relevant
                rv = 0;
            else
                % TR 38.901 Section Table 7.6.3.4-2, standard deviation for
                % O2I penetration loss for different building types is
                % uncorrelated, so two matrices exist, one for high loss
                % and one for low loss
                m = sc.O2I.AutoCorrMatrices(:,:,2 + lowloss);
                % TR 38.901 Section 7.6.3.3 says "penetration loss
                % deviation sigma_P ... can be modeled as a spatially
                % consistent random variable with correlation distance
                % 10m"; the interpretation here is that the term
                % N(0,sigma_P^2) in Equation 7.4-2 is formed by scaling a
                % spatially-consistent N(0,1) RV by the sigma_P given in TR
                % 38.901 Table 7.4.3-2
                rv = normalAutoCorrRVs(m,sc.O2I.FirstCoord,uepos(1:2));
            end
        else
            rv = rs.randn;
        end

        PL = PL_tw + PL_in + rv*sigma_P;

    end

end

% Environment height, TR 38.901 Table 7.4.1-1 Note 1
function h_E = environmentHeight(rs,s,d_2D,h_UT)

    h_E_UMi = 1;
    h_E_UMa = @()environmentHeight_UMa(rs,d_2D,h_UT);
    % h_E is not relevant for RMa, InH and InF
    h_E_RMa = 0;
    h_E_InH = 0;
    h_E_InF = 0;
    h_E = scenarioSwitch(s,h_E_UMi,h_E_UMa,h_E_RMa,h_E_InH,h_E_InF);

    function h_E = environmentHeight_UMa(rs,d_2D,h_UT)

        if (h_UT < 13)
            C = 0;
        else
            if (d_2D <= 18)
                % Note that the minimum 2-D distance for UMa in TR 38.901
                % Table 7.2-1 is 35m
                g = 0;
            else
                g = 5 / 4 * (d_2D/100)^3 * exp(-d_2D/150);
            end
            C = ((h_UT - 13)/10)^1.5 * g;
        end
        if (rs.rand < (1/(1+C)))
            h_E = 1;
        else
            v = 12:3:(h_UT-1.5);
            h_E = v(rs.randi(numel(v)));
        end

    end

end

function [LOS,pLOS,LOS_soft] = losCondition(chCfg,rs,d_2D_out,h_BS,h_UT,sc)

    persistent c;
    if (isempty(c))
        c = physconst('LightSpeed');
    end

    if (isempty(chCfg.LOSProbability))
        % TR 38.901 Table 7.4.2-1
        pLOS = losProbability(chCfg,d_2D_out,h_BS,h_UT);
    else
        % User-specified LOS probability
        pLOS = chCfg.LOSProbability;
    end

    spatialConsistency = ~isempty(sc);
    if (spatialConsistency)
        % TR 38.901 Section 7.6.3.3 "realization of a random variable
        % generated with [distance dependence] ... if the realization is
        % less than the LOS probability, the state is LOS; otherwise NLOS"
        if (~isfield(sc,'LOS'))
            validateSCRVSetup([]);
        end
        nSites = size(sc.LOS.AutoCorrMatrices,3);
        ai = mod(sc.NodeSubs(1) - 1,nSites) + 1;
        uepos = sc.UEPosition;
        % G is the "spatially consistent Gaussian random variable" for soft
        % LOS state
        G = normalAutoCorrRVs(sc.LOS.AutoCorrMatrices(:,:,ai),sc.LOS.FirstCoord,uepos(1:2));
        % Soft LOS state defined in TR 38.901 Section 7.6.3.3
        F = sqrt(2) * erfinv(2*pLOS - 1);
        lambda = c / chCfg.CenterFrequency;
        LOS_soft = 1/2 + 1/pi*atan(sqrt(20/lambda)*(G + F));
        % Transform G to the uniform distribution U(0,1)
        rv = probabilityIntegralTransform(G);
    else
        LOS_soft = [];
        rv = rs.rand;
    end

    LOS = rv < pLOS;

end

% LOS probability, TR 38.901 Section 7.4.2
function pLOS = losProbability(chCfg,d_2D_out,h_BS,h_UT)

    % Note that for InH scenarios, although Table 7.4.2-1 uses the variable
    % d_2D_in, its value is equal to d_2D_out here because this
    % implementation treats the link as outdoor i.e. it has no O2I
    % penetration loss. Similarly, d_2D_out is used in place of d_2D for
    % the InF scenario
    d_2D_in = d_2D_out;
    d_2D = d_2D_out;

    % For InH, augment the scenario name with the office scenario type
    s = chCfg.Scenario;
    if (s=="InH")
        if (chCfg.OfficeType=="Mixed")
            s = s + "-M";
        else % Open
            s = s + "-O";
        end
    end

    pLOS_UMi = @()losProbability_UMi(d_2D_out);
    pLOS_UMa = @()losProbability_UMa(d_2D_out,h_UT);
    pLOS_RMa = @()losProbability_RMa(d_2D_out);
    f = @losProbability_InH;
    pLOS_InH = {
        @()f(d_2D_in,[1.2 6.5],[4.7 32.6],0.32) % Mixed
        @()f(d_2D_in,[5 49],[70.8 211.7],0.54)  % Open
        };
    f = @(isL)losProbability_InF_SorD(chCfg,d_2D,isL,h_BS,h_UT);
    pLOS_InF = {
        @()f(true)  % SL
        @()f(true)  % DL
        @()f(false) % SH
        @()f(false) % DH
        @()1        % HH
        };
    pLOS = scenarioSwitch(s,pLOS_UMi,pLOS_UMa,pLOS_RMa,pLOS_InH,pLOS_InF);

    % If the LOS probability is not in [0,1], produce a warning
    if (pLOS<0 || pLOS>1)
        n = ["2-D outdoor distance"; "2-D indoor distance"; "2-D distance";
            "BS height"; "UE height";
            "clutter size"; "clutter density"; "clutter height"];
        v = [d_2D_out; d_2D_in; d_2D;
            h_BS; h_UT;
            chCfg.ClutterSize; chCfg.ClutterDensity; chCfg.ClutterHeight];
        nvstr = @(i)sprintf('%s (%0.3f)',n(i),v(i));
        nvfn = @(x)strjoin(arrayfun(nvstr,x,UniformOutput=false),', ');
        nUMi = @()1;
        nUMa = @()[1 5];
        nRMa = @()1;
        nInH = @()2;
        nInF = @()[3 4 5 6 7 8];
        warning('nr5g:h38901Channel:InvalidPLOS', ...
            'The value of the LOS probability (%0.4f) is outside the interval [0,1] for the specified scenario (%s), %s.', ...
            pLOS,s,nvfn(scenarioSwitch(s,nUMi,nUMa,nRMa,nInH,nInF)));
    end

    function pLOS = losProbability_UMi(d_2D_out)

        if (d_2D_out<=18)
            pLOS = 1;
        else
            pLOS = (18/d_2D_out + exp(-d_2D_out/36)*(1-18/d_2D_out));
        end

    end

    function pLOS = losProbability_UMa(d_2D_out,h_UT)

        if (d_2D_out<=18)
            % Note that the minimum 2-D distance for UMa in TR 38.901 Table
            % 7.2-1 is 35m
            pLOS = 1;
        else
            if (h_UT<=13)
                c = 0;
            else
                c = ((h_UT-13)/10)^(1.5);
            end
            pLOS = (18/d_2D_out + exp(-d_2D_out/63)*(1-18/d_2D_out)) * (1 + c*5/4*(d_2D_out/100)^3*exp(-d_2D_out/150));
        end

    end

    function pLOS = losProbability_RMa(d_2D_out)

        if (d_2D_out<=10)
            % Note that the minimum 2-D distance for RMa in TR 38.901 Table
            % 7.2-3 is 35m
            pLOS = 1;
        else
            pLOS = exp(-(d_2D_out-10)/1e3);
        end

    end

    function pLOS = losProbability_InH(d_2D_in,d_BP,den,m)

        if (d_2D_in < d_BP(1))
            pLOS = 1;
        elseif (d_2D_in < d_BP(2))
            pLOS = exp(-(d_2D_in - d_BP(1))/den(1));
        else
            pLOS = exp(-(d_2D_in - d_BP(2))/den(2)) * m;
        end

    end

    function pLOS = losProbability_InF_SorD(chCfg,d_2D,isL,h_BS,h_UT)

        d_clutter = chCfg.ClutterSize;
        r = chCfg.ClutterDensity;
        if (isL)
            k = -d_clutter / log(1 - r);
        else
            h_c = chCfg.ClutterHeight;
            k = -d_clutter / log(1 - r) * (h_BS - h_UT) / (h_c - h_UT);
        end
        pLOS = exp(-d_2D/k);

    end

end

% Configure array orientations
function [txorientation,rxorientation] = configureArrayOrientations(chCfg,uers,nodesubs,nodesiz)

    % Configure transmit array orientation in degrees
    if (~isempty(chCfg.TransmitArrayOrientation))
        txorientation = chCfg.TransmitArrayOrientation;
    else
        sectorIndex = nodesubs(2);
        numSectors = nodesiz(2);
        alphastep = 360/numSectors;
        alpha = (sectorIndex-1)*alphastep + 30;
        txorientation = [alpha; 0; 0];
    end

    % Configure receive array orientation in degrees
    if (~isempty(chCfg.ReceiveArrayOrientation))
        % A call is made to 'uers.rand' here to match the call in the
        % 'else' branch, even though the value generated here is not used.
        % The purpose of the call here is to produce the same RNG state
        % after this if/else, regardless of whether the 'if' or 'else'
        % branch is followed
        uers.rand;
        rxorientation = chCfg.ReceiveArrayOrientation;
    else
        rxorientation = [(uers.rand*360 - 180); 0; 0];
    end

end

% Configure transmit antenna array
function tx = configureTransmitAntennaArray(nodesiz,Nt,virtualization)

    % Set up structure for array, using nrCDLChannel to provide default
    % values
    persistent txarray;
    if (isempty(txarray))
        cdl = nrCDLChannel;
        txarray = rmfield(cdl.TransmitAntennaArray,'Orientation');
    end
    tx = txarray;

    % Establish the BS transmit array size
    numCSIRSPorts = Nt;
    [Ng,M,N] = panelDimensions(numCSIRSPorts);
    M = M * virtualization.K;
    N = N * virtualization.L;
    P = min(2,numCSIRSPorts);
    txsiz = [M N P 1 Ng];

    % If sectorization is configured, use the BS antenna element defined in
    % TR 38.901, otherwise use an isotropic antenna element for the BS
    numSectors = nodesiz(2);
    if (numSectors > 1)
        bsElement = '38.901';
    else
        bsElement = 'isotropic';
    end

    % Configure antenna array
    tx.Size = txsiz;
    tx.ElementSpacing = [0.5 0.5 0.5*M 0.5*N];
    if (P==1)
        tx.PolarizationAngles = 0;
    else
        tx.PolarizationAngles = [45 -45];
    end
    tx.Element = bsElement;

end

% Configure receive antenna array
function rx = configureReceiveAntennaArray(Nr)

    % Set up structure for array, using nrCDLChannel to provide default
    % values
    persistent rxarray;
    if (isempty(rxarray))
        cdl = nrCDLChannel;
        rxarray = rmfield(cdl.ReceiveAntennaArray,'Orientation');
    end
    rx = rxarray;

    % Configure array 
    numSRSPorts = Nr;
    rx.Size = [1 ceil(numSRSPorts/2) min(2,numSRSPorts) 1 1];

end

% Configure speed and direction
function [maximumDopplerShift,utDirectionOfTravel,movingScattererProportion,maximumScattererSpeed] = configureSpeedAndDirection(uers,s,centerFrequency,d_2D_in)

    persistent c;
    if (isempty(c))
        c = physconst('LightSpeed');
    end

    % TR 38.901 Section 7.6.10 (dual mobility) is not considered. Mobile
    % scatterers are not configured, and the BS has zero Doppler shift
    movingScattererProportion = 0;
    maximumScattererSpeed = 0;

    % TR 38.901 Table 7.2-1 for UMi and UMa. For RMa, assume that indoor
    % users are 3 km/h and users in cars are 120 km/h (inferred from Rep.
    % ITU-R M.2135-1). TR 38.901 Table 7.2-2 for InH. Assume 0 km/h for
    % InF
    indoor = (d_2D_in~=0);
    if (indoor)
        vRMa = 3.0;
    else
        vRMa = 120.0;
    end
    velocity = scenarioSwitch(s,3.0,3.0,vRMa,3.0,0.0);
    maximumDopplerShift = [(velocity*1000/3600)/c*centerFrequency 0];

    % TR 38.901 Table 7.2-1 for UMi and UMa "UT distribution
    % (horizontal) = Uniform" and assuming the same for RMa
    utDirectionOfTravel = [(uers.rand*360 - 180) 0; 90 90];

end

function def = channelLinkDefaults(chcfg)

    def.TXRUVirtualization = struct(K=1,Tilt=0,L=1,Pan=0);
    def.Seed = 0;
    def.d_2D_in = 0;
    def.CenterFrequency = 4e9;
    def.HasSmallScale = true;
    def.FastFading = true;
    def.n_fl = 1;
    def.Wrapping = false;
    def.SpatialConsistency = false;
    def.ScenarioExtents = [];
    if (~isCellular(chcfg.Scenario))
        def.InterSiteDistance = NaN;
        if (startsWith(chcfg.Scenario,"InF"))
            xy = chcfg.HallSize([1 2]);
            def.ScenarioExtents = [-xy/2 xy];
        end
    end
    def.TransmitArrayOrientation = [];
    def.ReceiveArrayOrientation = [];
    def.LOSProbability = [];

end

function chcfg = setDefaults(chcfg,def)

    for f = string(fieldnames(def).')

        if (~isfield(chcfg,f) || isempty(chcfg.(f)))
            chcfg.(f) = def.(f);
        end

    end

end

%% ========================================================================

% Prepare the TR 38.901 Section 7.5 fast fading model
function [channel,chinfo,chCfg,siteuers,SCRVs_out] = fastFadingChannelModel(chCfg)

    persistent c;
    if isempty(c)
        c = physconst('LightSpeed');
    end

    % Validate scenario
    persistent sv;
    if (isempty(sv))
        sv = nrPathLossConfig.Scenario_Values;
    end
    chCfg.Scenario = validatestring(chCfg.Scenario,sv,'h38901Channel','Scenario');

    % Prepare RNG for parameters that are specific to site and UE
    % (regardless of sector)
    nodesiz = chCfg.NodeSiz;
    nodesubs = chCfg.NodeSubs;
    siteueseed = sub2ind(nodesiz([1 3]),nodesubs(1),nodesubs(3)) - 1 + chCfg.Seed;
    siteuers = RandStream('mt19937ar','Seed',siteueseed);

    % Prepare RNG for parameters that are UE-specific (regardless of site
    % and sector)
    ueseed = nodesubs(3) - 1 + chCfg.Seed;
    uers = RandStream('mt19937ar','Seed',ueseed);

    % Prepare RNG for parameters that are site-specific (regardless of
    % sector and UE)
    siteseed = nodesubs(1) - 1 + chCfg.Seed;
    siters = RandStream('mt19937ar','Seed',siteseed);

    % Nonempty site positions indicate a new scenario configuration.
    % Reset any state related to the site positions
    persistent SCRVs;
    persistent clusterSCRVfn;
    if (~isempty(chCfg.SitePositions))
        if (spatialConsistency(chCfg)) % "Static", "ProcedureA" or "ProcedureB"
            % A system-wide RNG is used here:
            % * For O2I, a single autocorrelation matrix is created, 
            %   corresponding to "All-correlated" in TR 38.901
            %   Table 7.6.3.4-1
            % * For LOS / NLOS, an autocorrelation matrix per site is 
            %   created, corresponding to "Site-specific" in TR 38.901
            %   Table 7.6.3.4-1
            % * For cluster-specific RVs, an autocorrelation matrix per
            %   cluster is created and is re-used across sites under the 
            %   assumption that the intersite distance is larger than the
            %   distance at which the spatial correlation is considered to 
            %   be negligible
            % * For absolute time of arrival in InF scenarios, an
            %   autocorrelation matrix per site is created as TR 38.901
            %   Section 7.6.9 states that "delta tau is generated
            %   independently for links between the same UT and different
            %   BS sites"
            seed = prod(nodesiz) + chCfg.Seed;
            systemrs = RandStream('mt19937ar','Seed',seed);
            % Set up structure containing autocorrelation matrices for
            % spatially consistent RVs, and function for creating spatially
            % consistent cluster-specific RVs once number of clusters is
            % known
            [SCRVs,clusterSCRVfn] = setupSCRVs(chCfg,systemrs);
            % Store the system-wide RNG for creation of spatially
            % consistent cluster-specific RVs during mobility modelling
            if (isProcedureA(chCfg))
                setupSCRVsForProcedureA(systemrs);
            end
            if (isProcedureB(chCfg))
                setupSCRVsForProcedureB(systemrs);
            end
        else
            SCRVs = [];
            % Clear variables which store system-wide RNG and
            % cluster-specific RVs for use during mobility modelling
            softLOSState();
            setupSCRVsForProcedureA();
            setupSCRVsForProcedureB();
            generatePowersForProcedureA();
        end
    end

    % If spatial consistency is enabled, update node positions and
    % subscripts in structure containing variables related to spatial
    % consistency
    if (spatialConsistency(chCfg)) % "Static", "ProcedureA" or "ProcedureB"
        SCRVs.BSPosition = chCfg.BSPosition;
        SCRVs.UEPosition = chCfg.UEPosition;
        SCRVs.NodeSubs = chCfg.NodeSubs;
    end

    % Perform geographical distance based wrapping
    % (Rec. ITU-R M.2101-0 Attachment 2 to Annex 1)
    bspos = chCfg.BSPosition;
    uepos = chCfg.UEPosition;
    if (isCellular(chCfg.Scenario) && chCfg.Wrapping)
        if isfield(chCfg, 'NumCellSites')
            numCellSites = chCfg.NumCellSites;
        else
            numCellSites = chCfg.NodeSiz(1);
        end
        offsets = h38901Channel.wrappingOffsets(chCfg.InterSiteDistance,numCellSites,chCfg.NodeSiz(2));
    else
        offsets = [0 0];
    end
    d = vecnorm((bspos(1:2) + offsets) - uepos(1:2),2,2);
    [~,idx] = min(d);
    bspos = bspos + [offsets(idx,:) 0];
    h_BS = bspos(3);
    d_3D = vecnorm(bspos - uepos);

    % Set up intermediate variables    
    v_2D = uepos(1:2) - bspos(1:2);
    d_2D = vecnorm(v_2D);
    % Do not allow negative outdoor distance
    d_2D_out = max(d_2D - chCfg.d_2D_in,0);
    h_UT = uepos(3);
    indoor = (chCfg.d_2D_in~=0);
    freq = chCfg.CenterFrequency;
    fastFading = chCfg.HasSmallScale && chCfg.FastFading;

    % NOTE: The next steps are from TR 38.901 Section 7.5

    % ---------------------------------------------------------------------
    % Step 1 d): Give BS and UT antenna field patterns and array geometries

    if (chCfg.HasSmallScale)
        if (isfield(chCfg,'TransmitAntennaArray'))
            txarray = chCfg.TransmitAntennaArray;
        else
            txarray = configureTransmitAntennaArray(nodesiz,chCfg.NumTransmitAntennas,chCfg.TXRUVirtualization);
        end
        if (isfield(chCfg,'ReceiveAntennaArray'))
            rxarray = chCfg.ReceiveAntennaArray;
        else
            rxarray = configureReceiveAntennaArray(chCfg.NumReceiveAntennas);
        end
    end

    % ---------------------------------------------------------------------
    % Step 1 e): Give BS and UT array orientations

    if (chCfg.HasSmallScale)
        [txorientation,rxorientation] = configureArrayOrientations(chCfg,uers,nodesubs,nodesiz);
    end

    % ---------------------------------------------------------------------
    % Step 1 f): Speed and direction of motion of UT

    if (fastFading)
        [maximumDopplerShift,utDirectionOfTravel,movingScattererProportion,maximumScattererSpeed] = configureSpeedAndDirection(uers,chCfg.Scenario,freq,chCfg.d_2D_in);
    end

    % ---------------------------------------------------------------------
    % Step 1 g): Specify system center frequency f_c and bandwidth B
    % - Center frequency is given by chCfg.CenterFrequency
    % - Bandwidth is not required for TR 38.901 Section 7.5

    % ---------------------------------------------------------------------
    % Step 2: Assign propagation condition (NLOS/LOS)
    % Calculate propagation condition (LOS/NLOS) for the current BS-UE link

    [LOS,pLOS] = losCondition(chCfg,siteuers,d_2D_out,h_BS,h_UT,SCRVs);

    if (~isempty(chCfg.SitePositions))
        if (isSpatiallyConsistentMobility(chCfg)) % "ProcedureA" or "ProcedureB"
            % Store the spatially consistent RVs for LOS, for use in
            % calculating soft LOS state during mobility modelling
            softLOSState(SCRVs.LOS);
        end
    end

    % ---------------------------------------------------------------------
    % Step 3: Calculate pathloss
    % Calculate path loss and shadow fading for the current BS-UE link
    
    % Create the path loss function configured with environment height, LOS
    % condition, O2I penetration loss for indoor UEs and shadow fading.
    % Note that the pathloss is not actually calculated here, instead a
    % function handle is prepared
    h_E = environmentHeight(siteuers,chCfg.Scenario,d_2D,h_UT);
    O2I = outdoorToIndoor(uers,chCfg.Scenario,chCfg.d_2D_in,freq,SCRVs);
    plfn = @(SF)@(plc)@(txpos,rxpos,freq)(nrPathLoss(setfield(plc,'EnvironmentHeight',h_E),freq,LOS,txpos.',rxpos.') + O2I - SF);

    % ---------------------------------------------------------------------
    % Step 4: Generate large scale parameters (LSPs)
    % Calculate large scale parameters for the current BS-UE link.
    % Large scale parameters are:
    % * Shadow fading (SF)
    % * Ricean K factor (K)
    % * Delay spread (DS)
    % * Angular spreads (ASD, ASA, ZSD, ZSA)

    % NOTE: The sign of the shadow fading is defined so that positive SF
    % means more received power at UE than predicted by the path loss
    % model. (NOTE 2 from Table 7.5-6)

    % Get the parameters from TR 38.901 Tables 7.5-6 to 7.5-11
    params = getChannelModelParameters(chCfg,indoor,LOS,freq,bspos,uepos);
    N = params.NumClusters;

    % Generate large scale parameters
    [SF,K,DS,ASD,ASA,ZSD,ZSA] = generateLSP(chCfg,siters,uepos,params,indoor,LOS,fastFading);

    % Limit azimuth spread values to 104 degrees and zenith spread values
    % to 52 degrees
    ASD = min(ASD,104);
    ASA = min(ASA,104);
    ZSD = min(ZSD,52);
    ZSA = min(ZSA,52);

    % ---------------------------------------------------------------------
    % Step 5: Generate delays

    % Set up spatially consistent cluster-specific RVs for delays
    setupSpatialConsistency = spatialConsistency(chCfg) && fastFading && ~isfield(SCRVs,'Delays');
    if (setupSpatialConsistency) % "Static", "ProcedureA" or "ProcedureB"
        if (~isProcedureB(chCfg)) % "Static" or "ProcedureA"
            SCRVs.Delays = clusterSCRVfn(N);
        else % "ProcedureB"
            distance = 2 * c * 10^(params.mu_lgDS + params.sigma_lgDS);
            SCRVs.Delays = setupSCRVsForProcedureB(chCfg,distance,N);
        end
    end

    % Path delays are returned as a vector of NumClusters elements for the
    % current BS-UE link
    if (fastFading)
        % Note that the path delays scaling factor for LOS is only used in
        % nrCDLChannel but not in the generation of the cluster powers in
        % Step 6.
        if (~isProcedureB(chCfg)) % "None", "Static" or "ProcedureA"
            [pathDelays,pathDelaysScaling,tau_prime] = generateDelays(siteuers,params.r_tau,DS,N,indoor,LOS,K,SCRVs);
        else % "ProcedureB"
            [pathDelays,pathDelaysScaling,tau_prime] = generateDelaysForProcedureB(params.mu_lgDS,params.sigma_lgDS,N,indoor,LOS,SCRVs);
        end
    else % LOS only
        tau_prime = 0;
    end

    % ---------------------------------------------------------------------
    % Steps 6 and 7
    % The order of next two steps depends on whether or not TR 38.901
    % Section 7.6.3.2 Procedure B is being implemented:
    % * If Procedure B is not being implemented, the order of steps in 
    %   Section 7.5 is followed: Step 6 generates powers and Step 7 
    %   generates angles. This order is followed when spatial consistency
    %   is disabled, for TR 38.901 Section 7.6.3.1 "Spatial consistency 
    %   procedure", and for TR 38.901 Section 7.6.3.2 Procedure A
    % * If Procedure B is being implemented, Step 6 generates angles and
    %   Step 7 generates powers

    % Set up spatially consistent cluster-specific RVs for powers. Note
    % that unlike delays and angles, the SCRVs for powers have the same
    % correlation distance for all forms of spatial consistency, so can be
    % set up in common here
    if (setupSpatialConsistency) % "Static", "ProcedureA" or "ProcedureB"
        SCRVs.Powers = clusterSCRVfn(N);
    end

    % Store spatially consistent cluster-specific RVs for powers for use
    % during Procedure A
    if (setupSpatialConsistency)
        if (isProcedureA(chCfg))
            generatePowersForProcedureA(SCRVs.Powers);
        end
    end

    if (~isProcedureB(chCfg)) % "None", "Static" or "ProcedureA"

        % -----------------------------------------------------------------
        % Step 6: Generate cluster powers

        % Path powers are returned as a vector of N elements for the
        % current BS-UE link, where N=NumClusters-numel(pathPowers<(-25dB))
        if (fastFading)
            isMobilityUpdate = false;
            [pathPowers,pathPowersLOS,pathDelays,KFactorFirstCluster,tau_prime] = generatePowers(chCfg,siteuers,pathDelays,params.r_tau,DS,params.zeta,K,indoor,LOS,SCRVs,tau_prime,isMobilityUpdate);
        else
            pathPowersLOS = [];
        end

        % -----------------------------------------------------------------
        % Step 7: Generate arrival and departure angles
        % Generate arrival and departure angles for each cluster. These
        % represent the AnglesAOA, AnglesAOD, AnglesZOA, and AnglesZOD
        % properties of nrCDLChannel and are the variables phi_{n,AOA},
        % phi_{n,AOD}, theta_{n,ZOA}, and theta_{n,ZOD} in TR 38.901
        % Equation 7.5-13, 7.5-18, and 7.5-20, respectively

        % Set up spatially consistent cluster-specific RVs for angles
        if (setupSpatialConsistency) % "Static" or "ProcedureA"
            SCRVs.AnglesAoA.Offset = clusterSCRVfn(N);
            SCRVs.AnglesAoD.Offset = clusterSCRVfn(N);
            SCRVs.AnglesZoA.Offset = clusterSCRVfn(N);
            SCRVs.AnglesZoD.Offset = clusterSCRVfn(N);
            SCRVs.AnglesAoA.Sign = clusterSCRVfn(N);
            SCRVs.AnglesAoD.Sign = clusterSCRVfn(N);
            SCRVs.AnglesZoA.Sign = clusterSCRVfn(N);
            SCRVs.AnglesZoD.Sign = clusterSCRVfn(N);
        end

        if (chCfg.HasSmallScale)
            [anglesAoA, anglesAoD, anglesZoA, anglesZoD] = generateAngles(siteuers,bspos,uepos,ASA,ASD,ZSA,ZSD,pathPowersLOS,indoor,LOS,K,params,fastFading,SCRVs);
        end

    else % "ProcedureB"

        % -----------------------------------------------------------------
        % Step 6: Generate arrival and departure angles

        % Set up spatially consistent cluster-specific RVs for angles
        if (setupSpatialConsistency) % "ProcedureB"
            distance_xoD = 2 * c * 10^(params.mu_lgDS + params.sigma_lgDS);
            SCRVs.AnglesAoD = setupSCRVsForProcedureB(chCfg,distance_xoD,N);
            SCRVs.AnglesZoD = setupSCRVsForProcedureB(chCfg,distance_xoD,N);
            distance_xoA = 50;
            SCRVs.AnglesAoA = setupSCRVsForProcedureB(chCfg,distance_xoA,N);
            SCRVs.AnglesZoA = setupSCRVsForProcedureB(chCfg,distance_xoA,N);
        end

        if (chCfg.HasSmallScale)
            [phi_prime_AoA,phi_prime_AoD,theta_prime_ZoA,theta_prime_ZoD,phiLOS_AoA,phiLOS_AoD,thetaLOS_ZoA,thetaLOS_ZoD] = generateAnglesForProcedureB(bspos,uepos,indoor,LOS,params,fastFading,SCRVs);
        end

        % -----------------------------------------------------------------
        % Step 7: Generate cluster powers
        % Note that cluster powers are calculated using the "primed" angles
        % 'phi_prime_AoA', 'phi_prime_AoD', 'theta_prime_ZoA',
        % 'theta_prime_ZoD' produced by Step 6 above

        % Path powers are returned as a vector of NumClusters elements for
        % the current BS-UE link
        if (fastFading)
            [pathPowers,pathPowersLOS,KFactorFirstCluster] = generatePowersForProcedureB(tau_prime,DS,phi_prime_AoA,phi_prime_AoD,theta_prime_ZoA,theta_prime_ZoD,ASA,ASD,ZSA,ZSD,params.zeta,K,indoor,LOS,SCRVs);
        else
            pathPowersLOS = [];
        end

        % -----------------------------------------------------------------
        % Step 7b: Apply offset angles
        % The variables 'anglesAoA', 'anglesAoD', 'anglesZoA', 'anglesZoD'
        % represent the AnglesAOA, AnglesAOD, AnglesZOA, and AnglesZOD
        % properties of nrCDLChannel and are equal to the angles defined in
        % TR 38.901 Eqs. (7.6-17d), (7.6-17e), and (7.6-17f) not including
        % the addition of the ray offset angles (term C*alpha_m where 'C'
        % is the cluster-wise angle spread). Addition of the ray offset
        % angles is implemented in nrCDLChannel, which completes the
        % implementation of Eqs. (7.6-17d), (7.6-17e), and (7.6-17f).
        % Therefore "Step 7b: Apply offset angles" here prepares the
        % relevant properties of nrCDLChannel rather than actually applying
        % the ray offset angles

        if (chCfg.HasSmallScale)
            [anglesAoA, anglesAoD, anglesZoA, anglesZoD] = completeAnglesForProcedureB(phi_prime_AoA,phi_prime_AoD,theta_prime_ZoA,theta_prime_ZoD,phiLOS_AoA,phiLOS_AoD,thetaLOS_ZoA,thetaLOS_ZoD,params.mu_offsetZOD,indoor);
        end

    end

    % ---------------------------------------------------------------------
    % Step 8: Perform random coupling of rays
    % Generate the array to randomly couple the rays per each cluster. This
    % step is already implemented in nrCDLChannel but TR 38.901 Section 7.5
    % explicitly highlights the fact that all the parameters generated
    % until step 9 must be the same for links between co-located sectors
    % and a UE. For this reason, the coupling of rays is done outside
    % nrCDLChannel so that we can ensure the equality in the case of
    % sectorized sites.

    if (fastFading)
        N = size(pathPowers,2);
        M = params.NumRaysPerCluster;
        rayCoupling = generateRayCoupling(siteuers,N,M);
    end

    % ---------------------------------------------------------------------
    % Step 9: Generate cross polarization power ratios (XPRs)
    % Generate log-normal values of cross polarization power ratios

    if (fastFading)
        xpr = params.mu_XPR + params.sigma_XPR*siteuers.randn(N,M);
    end

    % ---------------------------------------------------------------------
    % Step 10: Draw random initial phases

    if (fastFading)
        initialPhases = generateInitialPhases(uers,N,M);
    end

    % ---------------------------------------------------------------------
    % Prepare for Step 11: Generate channel coefficients
    % NOTE: Step 11 is not performed here, it is performed when the small
    % scale channel returned by this function is executed

    % Create the small scale part of the channel (CDL) if required
    if (chCfg.HasSmallScale)

        % Create nrCDLChannel (for fastFading=true) or channel structure
        % (for fastFading=false)
        if (fastFading)

            ss = nrCDLChannel();
            ss.DelayProfile = 'Custom';
            ss.CarrierFrequency = freq;
            ss.SampleRate = chCfg.SampleRate;
            ss.NormalizeChannelOutputs = false;
            ss.ChannelFiltering = false;
            ss.NumTimeSamples = floor(chCfg.SampleRate * 1e-3); % one subframe
            ss.OutputDataType = 'single';
            ss.TransmitAntennaArray = txarray;
            ss.TransmitArrayOrientation = txorientation;
            ss.ReceiveAntennaArray = rxarray;
            ss.ReceiveArrayOrientation = rxorientation;
            ss.MaximumDopplerShift = maximumDopplerShift;
            ss.UTDirectionOfTravel = utDirectionOfTravel;
            ss.MovingScattererProportion = movingScattererProportion;
            ss.MaximumScattererSpeed = maximumScattererSpeed;
            ss.NormalizePathGains = false;

            ss.PathDelays = pathDelays/pathDelaysScaling;
            % Note that for DelayProfile='Custom', nrCDLChannel does not
            % explicitly implement Eq. (7.5-30), the terms involving the
            % K-factor K_R are not implemented inside nrCDLChannel and
            % instead must be accounted for in the AveragePathGains
            % property. This is achieved here by using 'pathPowersLOS', the
            % output of Eq. (7.5-8), to set AveragePathGains. TR 38.901
            % says "these power values are used only in equations (7.5-9)
            % and (7.5-14), but not in equation (7.5-22)" and Eq. (7.5-22)
            % is implemented inside nrCDLChannel. However, given the way
            % that nrCDLChannel implements Eq. (7.5-30), the net effect of
            % using the output of Eq. (7.5-8) as the input to nrCDLChannel
            % is equivalent to TR 38.901's requirement to use the output of
            % Eq. (7.5-6) as the input to Eq. (7.5-22)
            ss.AveragePathGains = 10*log10(pathPowersLOS);
            ss.AnglesAoA = anglesAoA;
            ss.AnglesAoD = anglesAoD;
            ss.AnglesZoA = anglesZoA;
            ss.AnglesZoD = anglesZoD;
            ss.HasLOSCluster = (LOS && ~indoor);
            if ss.HasLOSCluster
                ss.KFactorFirstCluster = KFactorFirstCluster;
            end
            ss.AngleSpreads = [params.C_ASD params.C_ASA params.C_ZSD params.C_ZSA];
            ss.RayCoupling = rayCoupling;
            ss.XPR = xpr;
            ss.InitialPhases = initialPhases;
            ss.NumStrongestClusters = min(N,2);
            if ~isnan(params.C_DS)
                ss.ClusterDelaySpread = params.C_DS * 1e-9;
            end

        else

            % TR 38.901 Section 7.8.1 "Large scale calibration" where "fast
            % fading channel is not modelled"
            ss = struct;
            ss.CarrierFrequency = freq;
            ss.TransmitAntennaArray = txarray;
            ss.TransmitArrayOrientation = txorientation;
            ss.ReceiveAntennaArray = rxarray;
            ss.ReceiveArrayOrientation = rxorientation;
            ss.AnglesAoD = anglesAoD;
            ss.AnglesAoA = anglesAoA;
            ss.AnglesZoD = anglesZoD;
            ss.AnglesZoA = anglesZoA;
            ss.SampleRate = chCfg.SampleRate;
            ss.PathDelays = 0;

        end

        % Configure path filters and path delays
        [pathFilters,pathDelays] = getPathFiltersAndDelays(ss);

    else

        % Configure empty small scale channel, path filters and path delays
        ss = [];
        pathFilters = [];
        pathDelays = [];

    end

    % ---------------------------------------------------------------------
    % Prepare for Step 12: Apply pathloss and shadowing
    % NOTE: Step 12 is not performed here, it is performed when the large
    % scale channel returned by this function is executed

    % Incorporate the shadow fading into the path loss
    plfn = plfn(SF);

    % Create the large scale part of the channel: path loss (including O2I
    % and SF) with geographical distance based wrapping
    if (isCellular(chCfg.Scenario) && chCfg.Wrapping)
        wrapping = true;
        ISD = chCfg.InterSiteDistance;
    else
        wrapping = false;
        ISD = 0;
    end
    if isfield(chCfg, 'NumCellSites')
        numCellSites = chCfg.NumCellSites;
    else
        numCellSites = chCfg.NodeSiz(1);
    end
    numSectors = chCfg.NodeSiz(2);
    lsfn = @(plc)largeScaleWrapper(wrapping,ISD,numCellSites,numSectors,plfn(plc));
    
    % ---------------------------------------------------------------------

    % Create channel info
    chinfo = newChannelInfo();
    chinfo.O2I = O2I;
    chinfo.pLOS = pLOS;
    chinfo.LOS = LOS;
    chinfo.h_E = h_E;
    chinfo.SF = SF;
    chinfo.K = K;
    chinfo.d_3D = d_3D;

    % Create the overall channel structure
    channel = newChannel();
    channel.CenterFrequency = freq;
    if (chCfg.HasSmallScale)
        channel.NumTransmitAntennas = getNumSignals(txarray);
        channel.NumReceiveAntennas = getNumSignals(rxarray);
    else
        channel.NumTransmitAntennas = chCfg.NumTransmitAntennas;
        channel.NumReceiveAntennas = chCfg.NumReceiveAntennas;
    end
    channel.LargeScale = lsfn;
    channel.SmallScale = ss;
    channel.TXRUVirtualization = chCfg.TXRUVirtualization;
    channel.PathFilters = pathFilters;
    channel.NodeSubs = chCfg.NodeSubs;
    channel.NodeSiz = chCfg.NodeSiz;
    channel.PathDelays = pathDelays;

    if (isSpatiallyConsistentMobility(chCfg)) % "ProcedureA" or "ProcedureB"
        % Record intermediate variables that are needed when applying
        % spatially consistent mobility
        channel.Mobility.tau_prime = tau_prime;
        if (chCfg.HasSmallScale)
            channel.Mobility.tau = ss.PathDelays;
        end
        channel.Mobility.DS = DS;
        channel.Mobility.r_tau = params.r_tau;
        channel.Mobility.zeta = params.zeta;
        channel.Mobility.K = K;
        channel.Mobility.indoor = indoor;
        if (fastFading)
            channel.Mobility.UTDirectionOfTravel = ss.UTDirectionOfTravel;
        else
            channel.Mobility.UTDirectionOfTravel = [0 0; 90 90];
        end
    end

    function numSignals = getNumSignals(array)

        if isstruct(array)
            numSignals = prod(array.Size);
        else % PhAST object
            if isa(array,'phased.internal.AbstractSubarray')
                numSignals = array.getNumSubarrays;
            else % phased.internal.AbstractArray
                numSignals = array.getNumElements;
            end
        end
    
    end

    % Return SCRVs for use in absolute time of arrival step
    SCRVs_out = SCRVs;

end

% Section 7.6.9 "Absolute time of arrival" for InF scenarios
function [channel,chinfo] = absoluteTOA(chCfg,channel,chinfo,siteuers,sc)

    persistent c0;
    if isempty(c0)
        c0 = physconst('LightSpeed');
    end

    if (chCfg.HasSmallScale)

        % TR 38.901 Table 7.6.9-1
        mu_lgDeltatau = -7.5;
        sigma_lgDeltatau = 0.4;

        % TR 38.901 Section 7.6.9 "Delta tau is generated from a lognormal
        % distribution"
        if (spatialConsistency(chCfg))
            nSites = size(sc.ATOA.AutoCorrMatrices,3);
            ai = mod(sc.NodeSubs(1) - 1,nSites) + 1;
            uepos = chCfg.UEPosition;
            rv = normalAutoCorrRVs(sc.ATOA.AutoCorrMatrices(:,:,ai),sc.ATOA.FirstCoord,uepos(1:2));
        else
            rv = siteuers.randn;
        end
        Deltatau = 10^(mu_lgDeltatau + sigma_lgDeltatau*rv);

        % TR 38.901 Section 7.6.9 "Delta tau should further be upper
        % bounded by 2L/c, where L is the largest dimension of the factory
        % hall"
        L = max(chCfg.HallSize);
        Deltatau = min(Deltatau,2*L/c0);

        % TR 38.901 Eq. (7.6-43) and Eq. (7.6-44)
        ss = channel.SmallScale;
        pathDelays = ss.PathDelays;
        pathDelays = pathDelays + (chinfo.d_3D / c0);
        if (chCfg.FastFading)
            LOS = ss.HasLOSCluster;
            pathDelays((1+LOS):end) = pathDelays((1+LOS):end) + Deltatau;
        end

        % Update small scale channel path delays
        channel.SmallScale.PathDelays = pathDelays;

        % Update path filters and path delays, note that the updated path 
        % delays are used for channel filtering in the full PHY
        [channel.PathFilters,channel.PathDelays] = getPathFiltersAndDelays(channel.SmallScale);

        % Record Delta tau in channel output (used to detect if absolute
        % ToA needs considered when applying TR 38.901 Secton 7.6.3.2
        % "Spatially-consistent UT/BS mobility modelling")
        channel.Mobility.Deltatau = Deltatau;

    end

end

function [pathFilters,pathDelays] = getPathFiltersAndDelays(ss)

    % If small scale channel is not a structure, it is an nrCDLChannel
    isCDLChannel = ~isstruct(ss);

    if (isCDLChannel)
        ssinfo = info(ss);
        pathDelays = ssinfo.PathDelays;
    else
        pathDelays = ss.PathDelays;
    end
    sampleDelays = ss.SampleRate * pathDelays;
    if (~isequal(pathDelays,0))
        pathFilters = cachedPathFilters(sampleDelays);
    else
        pathFilters = 1;
    end

end

function [SF,K,DS,ASD,ASA,ZSD,ZSA] = generateLSP(chCfg,rs,uepos,p,indoor,LOS,fastFading)
    % Generate the large scale parameters (LSP) for the given configuration
    % and propagation condition, as discussed in step 4 of TR 38.901
    % Section 7.5.

    % Cache auto-correlation matrices for each site to ensure that all UEs
    % connected to a site have the same random auto-correlation matrix and
    % to decrease the computation time of all LSPs
    persistent autoCorrMatrices firstCoord;
    persistent minpos maxpos;

    % Nonempty site positions indicate a new scenario configuration. Reset
    % any state related to the site positions
    allsitepos = chCfg.SitePositions;
    if ~isempty(allsitepos)
        % Get bounding box of scenario
        minpos = chCfg.ScenarioExtents(1:2);
        maxpos = minpos + chCfg.ScenarioExtents(3:4);
        % Initialize empty autocorrelation matrices, one row for each site
        % and columns for NLOS, LOS, and for each floor of O2I (floors
        % 1...8) if applicable. Note that if more than 8 floors are
        % encountered during simulation, the allocation here will be
        % extended automatically
        if (isCellular(chCfg.Scenario))
            nFloorsO2I = 8;
        else
            nFloorsO2I = 0;
        end
        autoCorrMatrices = cell(size(allsitepos,1),2+nFloorsO2I);
        firstCoord = cell(size(allsitepos,1),2+nFloorsO2I);
    end

    if ~isnan(chCfg.InterSiteDistance)
        % Generate auto-correlation matrices and extract the values related
        % to the UE position. This procedure is discussed in Section 3.3.1
        % of IST-4-027756 WINNER II Deliverable 1.1.2 v1.2, "WINNER II
        % Channel Models", IST-WINNER2, Tech. Rep., 2007
        nodeIdx = chCfg.NodeSubs(1);
        [ai,aj] = getLSPAutoCorrMatrixSubscripts(size(autoCorrMatrices,1),nodeIdx,indoor,LOS,chCfg.n_fl);
        if any(size(autoCorrMatrices)<[ai aj]) || isempty(autoCorrMatrices{ai,aj})
            [autoCorrMatrices{ai,aj},firstCoord{ai,aj}] = getLSPAutoCorrMatrices(rs,minpos,maxpos,p);
        end
        % Generate UE indices from grid
        [ui,uj] = pixelsubs(uepos(1:2),firstCoord{ai,aj});
        % Get the auto-correlated random values for each LSP, expressed as:
        % s_M = [SF K DS ASD ASA ZSD ZSA]'
        s_M = permute(autoCorrMatrices{ai,aj}(ui,uj,:),[3 1 2]);
    else
        % The user did not provide the value of the intersite distance so
        % the LSP are assumed without auto-correlation
        s_M = randn(rs,[7,1]);
    end

    % Generate cross-correlation matrix
    crossCorrMatrix = getCrossCorrelationMatrix(p);

    % Generate LSP
    correlatedLSP = crossCorrMatrix*s_M;
    SF = correlatedLSP(1).*p.sigma_SF; % Shadow fading in dB
    if fastFading
        K   = p.mu_K + correlatedLSP(2)*p.sigma_K; % K-factor in dB
        DS  = 10^(p.mu_lgDS + correlatedLSP(3)*p.sigma_lgDS);   % Delay spread in s
        ASD = 10^(p.mu_lgASD + correlatedLSP(4)*p.sigma_lgASD); % AOD spread in deg
        ASA = 10^(p.mu_lgASA + correlatedLSP(5)*p.sigma_lgASA); % AOA spread in deg
        ZSD = 10^(p.mu_lgZSD + correlatedLSP(6)*p.sigma_lgZSD); % ZOD spread in deg
        ZSA = 10^(p.mu_lgZSA + correlatedLSP(7)*p.sigma_lgZSA); % ZOA spread in deg
    else
        % In the case where fast fading is not required, only shadow fading
        % is needed
        [K,DS,ASD,ASA,ZSD,ZSA] = deal([]);
    end

    function [ai,aj] = getLSPAutoCorrMatrixSubscripts(nodeSiz,nodeIdx,indoor,LOS,n_fl)

        % Row subscript is node index
        ai = nodeIdx;

        if (~indoor)
            % Column subscript 1 = NLOS
            % Column subscript 2 = LOS
            aj = LOS + 1;
        else % O2I
            % If there are more nodes than floors
            if (nodeSiz > n_fl)
                % Under the assumption that the intersite distance is
                % larger than the distance at which the spatial correlation
                % is considered to be negligible, uncorrelated parameters
                % across the different floors can be obtained by simply
                % re-using the autocorrelation matrix for a different site
                aj = 3;
                ai = mod(ai - 1 + n_fl,nodeSiz) + 1;
            else
                % Column subscript 3 = floor 1 for O2I
                % ...
                % Column subscript N = floor N-2 for O2I
                aj = n_fl + 2;
            end
        end

    end

    function [autoCorrMatrices,firstCoord] = getLSPAutoCorrMatrices(rs,minpos,maxpos,p)
        % Generate the auto-correlation matrices for LSPs. Also return
        % 'firstCoord', the 2-D position of the first element of the pixel
        % grid (that is, the position corresponding to
        % autoCorrMatrices(1,1,:))

        % LSP correlation distances in horizontal plane, defined in Table 7.5-6
        % of TR 38.901
        distances = [p.corr_SF p.corr_K p.corr_DS p.corr_ASD p.corr_ASA p.corr_ZSD p.corr_ZSA];

        % Generate the auto-correlation matrices
        [autoCorrMatrices,firstCoord] = h38901Channel.createAutoCorrMatrices(rs,minpos,maxpos,distances);

        % Make sure to not propagate NaN for the NLOS/O2I cases
        autoCorrMatrices(isnan(autoCorrMatrices)) = 0;
    end

    function C0_sqrt = getCrossCorrelationMatrix(p)
        % Generate the square root of the cross-correlation matrix
        % C_{MxM}(0)

        C0 = zeros(7);
        % Construct the 7-by-7 cross-correlation matrix C0
        x = [.5, ...
             p.SFvsK,   .5, ...
             p.DSvsSF,  p.DSvsK,  .5, ...
             p.ASDvsSF, p.ASDvsK, p.ASDvsDS, .5, ...
             p.ASAvsSF, p.ASAvsK, p.ASAvsDS, p.ASDvsASA, .5, ...
             p.ZSDvsSF, p.ZSDvsK, p.ZSDvsDS, p.ZSDvsASD, p.ZSDvsASA, .5, ...
             p.ZSAvsSF, p.ZSAvsK, p.ZSAvsDS, p.ZSAvsASD, p.ZSAvsASA, p.ZSDvsZSA, .5];
        upperTriangleIndices = triu(true(7));
        C0(upperTriangleIndices) = x; % Fill in the upper piece
        C0 = C0 + C0.';

        % Generate sqrt(C0) such that C0_sqrt*C0_sqrt = C0
        C0_sqrt = sqrtm(C0);
    end
end

function [autoCorrMatrices,firstCoord] = createAutoCorrMatricesLocal(rs,minpos,maxpos,distances,M)
    % Generate auto-correlation matrices for specified correlation
    % distances. Also return 'firstCoord', the 2-D position of the first
    % element of the pixel grid (that is, the position corresponding to
    % autoCorrMatrices(1,1,:)). A decimation factor 'M' can optionally be
    % specified which creates the matrices at a lower resolution and then
    % interpolates back to the default resolution given by
    % getGridResolution(), this is a quicker approach for large correlation
    % distances

    % If decimation factor is not specified, use M = 1
    if (nargin==4)
        M = 1;
    end

    % Permute to put the vector of distances in the third dimension, to
    % make it easy to apply different distance scalings to the distance
    % matrix created below
    distances = permute(distances,[3 1 2]);

    % Get intermediate variables that depend on grid resolution, including
    % the effect of the decimation factor
    [firstCoord,imax,jmax,e,res] = getResolutionDependentVariables(minpos,maxpos,distances,M);

    % Generate grids with iid normal random variables ~N(0,1)
    D = numel(distances);
    normgrid = randn(rs,[imax jmax D]);

    % Calculate kernel distance matrix
    kernelmax = pos2pixel(2*e,0,res);
    kernelcenter = pos2pixel(2*e,e,res);
    offsets = repmat((1:kernelmax) - kernelcenter, kernelmax, 1);
    distance_matrix = sqrt((res * offsets).^2 + (res * offsets.').^2);

    % Calculate filter kernel impulse response
    H = exp(-distance_matrix./distances);

    % Compute auto-correlation by filtering the normal grids for each
    % random variable
    autoCorrMatrices = nan(size(normgrid));
    c = kernelcenter;
    for rv = 1:D
        filtered = filter2(H(:,:,rv),normgrid(:,:,rv));
        autoCorrMatrices(:,:,rv) = filtered / std(filtered(c:end-c+1,c:end-c+1),0,'all');
    end

    % If decimation is configured, perform interpolation back to the
    % default grid resolution
    if (M > 1)

        % Get intermediate variables corresponding to the default grid
        % resolution
        [firstCoordq,imaxq,jmaxq,~,res] = getResolutionDependentVariables(minpos,maxpos,distances);

        % Create subscript arrays which describe the locations of the
        % decimated auto-correlation matrices in terms of the default grid
        % resolution
        [x0,y0] = splitpixel(round((firstCoord - firstCoordq) / res));
        x = x0 + (1:M:imax*M);
        y = y0 + (1:M:jmax*M);
        z = 1:D;
        [X,Y,Z] = ndgrid(x,y,z);

        % Create interpolator
        gi = griddedInterpolant(X,Y,Z,autoCorrMatrices,'linear','linear');

        % Create subscript arrays of all elements in the auto-correlation
        % matrices for the default grid resolution, these are used as query
        % points for the interpolator
        xq = 1:imaxq;
        yq = 1:jmaxq;
        zq = z;
        [XQ,YQ,ZQ] = ndgrid(xq,yq,zq);

        % Perform interpolation
        autoCorrMatrices = gi(XQ,YQ,ZQ);

    end
end

function [firstCoord,imax,jmax,e,res] = getResolutionDependentVariables(minpos,maxpos,distances,M)
    % Calculate variables used for generation of auto-correlation matrices
    % that are dependent on the resolution implied by decimation factor 'M'

    % If decimation factor is not specified, use M = 1
    if (nargin==3)
        M = 1;
    end

    % Add extra distance 'e' to consider UEs on the edge of the map.
    % Consider 4.6 correlation distances (for the random variable with the
    % largest correlation distance). Beyond this distance the correlation
    % is < 1% (exp(-4.6)~=0.01) so considered to be negligible. This value
    % is also used for the filter kernel half size
    e = 4.6 * max(distances);

    % Ensure that extra distance is a multiple of the grid resolution
    res = getGridResolution() * M;
    e = e + mod(-e,res);

    % Calculate size of random matrix to be generated
    firstCoord = minpos + [-e -e];
    maxpos = maxpos + [e e];
    [imax,jmax] = splitpixel(pos2pixel(maxpos,firstCoord,res));
end

function res = getGridResolution()
    % Resolution of the auto-correlation pixel grid in meters

    res = 5;
end

function [i,j] = pixelsubs(xy,firstCoord)
    % [I,J] = pixelsubs(XY,FIRSTCOORD) returns row and column subscripts I
    % and J of the nearest autocorrelaton matrix element to given
    % coordinate pair, XY, given the coordinate pair of the first
    % autocorrelaton matrix element, FIRSTCOORD.

    [i,j] = splitpixel(pos2pixel(xy,firstCoord));
end

function px = pos2pixel(pos,firstCoord,res)
    % Convert current position in m to pixel in the auto-correlation grid

    if (nargin==2)
        res = getGridResolution();
    end
    px = floor((pos - firstCoord) / res) + 1;
end

function [i,j] = splitpixel(px)
    % Split columns of a 2-D pixel matrix into separate subscripts

    i = px(:,1);
    j = px(:,2);
end

function [SCRVs,clusterSCRVfn] = setupSCRVs(chCfg,systemrs)
    
    % O2I
    % TR 38.901 Section 7.6.3.3, building type, correlation
    % distance 50 m
    % TR 38.901 Section 7.6.3.3, penetration loss deviation
    % sigma_p, correlation distance 10 m
    % TR 38.901 Section Table 7.6.3.4-2, standard deviation for O2I
    % penetration loss for different building types is uncorrelated, so two
    % matrices are generated, one for high loss and one for low loss
    o2iDelta = [50 10 10];
    % The system-wide RNG corresponds to "All-correlated" in TR 38.901
    % Table 7.6.3.4-1
    minpos = chCfg.ScenarioExtents(1:2);
    maxpos = minpos + chCfg.ScenarioExtents(3:4);
    [SCRVs.O2I.AutoCorrMatrices,SCRVs.O2I.FirstCoord] = h38901Channel.createAutoCorrMatrices(systemrs,minpos,maxpos,o2iDelta);

    % LOS / NLOS state
    % TR 38.901 Table 7.6.3.1-2
    % Note that scenarioSwitch order is UMi, UMa, RMa, InH, InF
    s = chCfg.Scenario;
    losDelta = scenarioSwitch(s,50,50,60,10,@()chCfg.ClutterSize/2);
    % The system-wide RNG is used to create site-specific spatially
    % consistent parameters by creating one autocorrelation matrix per site
    nSites = size(chCfg.SitePositions,1);
    distances = repmat(losDelta,1,nSites);
    [SCRVs.LOS.AutoCorrMatrices,SCRVs.LOS.FirstCoord] = h38901Channel.createAutoCorrMatrices(systemrs,minpos,maxpos,distances);

    if (startsWith(chCfg.Scenario,"InF") && chCfg.AbsoluteTOA)
        % Absolute time of arrival TR 38.901 Table 7.6.9-1 "Correlation
        % distance in the horizontal plane". The table does not provide a
        % correlation distance for InF-HH, 11 m is assumed here
        atoaDelta = scenarioSwitch(s,NaN,NaN,NaN,NaN,{6 6 11 11 11});
        distances = repmat(atoaDelta,1,nSites);
        % The system-wide RNG is used to create site-specific spatially
        % consistent parameters by creating one autocorrelation matrix per
        % site
        [SCRVs.ATOA.AutoCorrMatrices,SCRVs.ATOA.FirstCoord] = h38901Channel.createAutoCorrMatrices(systemrs,minpos,maxpos,distances);
    end

    % Set up function for creating spatially consistent cluster-specific
    % RVs once number of clusters is known
    clusterSCRVfn = @(numClusters)createClusterAutoCorrMatricesStatic(chCfg,systemrs,numClusters);

end

function ssp = createClusterAutoCorrMatricesStatic(chCfg,rs,nClusters)
    % Generate the auto-correlation matrices for cluster-specific RVs in
    % Section 7.6.3.1 (i.e. chCfg.SpatialConsistency="Static")

    % TR 38.901 Table 7.6.3.1-2
    % Note that scenarioSwitch order is UMi, UMa, RMa, InH, InF
    s = chCfg.Scenario;
    deltaNLOS = scenarioSwitch(s,15,50,60,10,10);
    deltaLOS = scenarioSwitch(s,12,40,50,10,10);
    deltaO2I = scenarioSwitch(s,15,15,15,NaN,NaN);
    distances = [deltaNLOS deltaLOS deltaO2I];

    % Generate the auto-correlation matrices
    ssp = createClusterAutoCorrMatrices(chCfg,rs,nClusters,distances);
    
end

function ssp = createClusterAutoCorrMatrices(chCfg,rs,nClusters,distances,varargin)
    % Generate the auto-correlation matrices for cluster-specific RVs. Also
    % return 'firstCoord', the 2-D position of the first element of the
    % pixel grid (that is, the position corresponding to
    % autoCorrMatrices(1,1,:))

    % Get bounding box of the scenario
    minpos = chCfg.ScenarioExtents(1:2);
    maxpos = minpos + chCfg.ScenarioExtents(3:4);

    % Generate the auto-correlation matrices
    nSites = size(chCfg.SitePositions,1);
    D = size(distances,2);
    distances = repmat(distances,max(nSites,nClusters),1);
    [autoCorrMatrices,firstCoord] = createAutoCorrMatricesLocal(rs,minpos,maxpos,distances(:).',varargin{:});
    [X,Y,~] = size(autoCorrMatrices);
    autoCorrMatrices = reshape(autoCorrMatrices,X,Y,[],D);

    % Create output structure
    ssp.AutoCorrMatrices = autoCorrMatrices;
    ssp.FirstCoord = firstCoord;
end

function rvs = normalAutoCorrRVs(autoCorrMatrix,firstCoord,pos)
    % Generate normally distributed spatially correlated random variables

    [ui,uj] = pixelsubs(pos,firstCoord);
    ind = sub2ind(size(autoCorrMatrix),ui,uj);
    rvs = autoCorrMatrix(ind);
end

function rvs = uniformClusterSCRVs(thisSSP,sc,N,varargin)
    % Create uniformly distributed spatially correlated cluster-specific
    % random variables

    % Get normally distributed spatially cluster-specific RVs
    rvs = normalClusterSCRVs(thisSSP,sc,N,varargin{:});

    % Transform to uniform distribution U(0,1) using the probability
    % integral transform
    rvs = probabilityIntegralTransform(rvs);
end

function rvs = normalClusterSCRVs(thisSSP,sc,N,varargin)
    % Create normally distributed spatially correlated cluster-specific
    % random variables
    
    % Autocorrelation matrix subscript corresponding to UE position
    pos = sc.UEPosition(1:2);
    [ui,uj] = pixelsubs(pos,thisSSP.FirstCoord);
    ui = repmat(ui,1,N);
    uj = repmat(uj,1,N);

    % Autocorrelation matrix array subscript corresponding to cluster
    % Under the assumption that the intersite distance is larger than the
    % distance at which the spatial correlation is considered to be
    % negligible, uncorrelated parameters across the different clusters can
    % be obtained by simply re-using the autocorrelation matrix for a
    % different site
    autoCorrMatrices = thisSSP.AutoCorrMatrices;
    K = size(autoCorrMatrices,3);
    nodeIdx = sc.NodeSubs(1);
    uk = mod(nodeIdx-1 + (0:N-1),K) + 1;

    if (nargin==5)
        % Autocorrelation matrix array subscript corresponding to link type
        % (NLOS / LOS / O2I)
        indoor = varargin{1};
        LOS = varargin{2};
        if (indoor)
            ul = 3;
        else
            ul = LOS + 1;
        end
    else % nargin==4
        % Autocorrelation matrix array subscript is directly specfiied
        ul = varargin{1};
    end
    ul = repmat(ul,1,N);

    % Get normally distributed spatially cluster-specific RVs 
    ind = sub2ind(size(autoCorrMatrices),ui,uj,uk,ul);
    rvs = autoCorrMatrices(ind);
end

function rvs = probabilityIntegralTransform(rvs)
    % Transform to uniform distribution U(0,1) using the probability
    % integral transform

    rvs = 1/2 * (1 + erf(rvs / sqrt(2)));
end

function p = getChannelModelParameters(chCfg,indoor,LOS,freq,bspos,uepos)
    % Get the parameters listed in TR 38.901 Tables 7.5-6 to 7.5-11 for the
    % selected scenario, propagation condition (LOS/NLOS/O2I), carrier
    % frequency, BS and UE positions.

    persistent tables38901_fastFading;
    persistent pmap;

    if isempty(tables38901_fastFading)
        tables38901_fastFading = get38901FastFadingTables();
        pmap = cell(height(tables38901_fastFading),1);
    end

    % Carrier frequency in GHz
    fc = freq/1e9;

    % If scenario is InF-xx, remove the clutter and height specification 
    % part "-xx"
    scenario = chCfg.Scenario;
    if (startsWith(scenario,"InF"))
        xx = extractAfter(scenario,'-');
        scenario = extractBefore(scenario,"-");
    end

    % Get the parameters structure from the table row specified by
    % scenario, O2I state and LOS state.
    RowName = scenario + "_";
    if indoor
        RowName = RowName + "O2I";
    elseif LOS
        RowName = RowName + "LOS";
    else % NLOS
        RowName = RowName + "NLOS";
    end
    rowIdx = matches(tables38901_fastFading.Properties.RowNames,RowName);
    if (isempty(pmap{rowIdx}))
        p = table2struct(tables38901_fastFading(rowIdx,:));
        pmap{rowIdx} = p;
    else
        p = pmap{rowIdx};
    end

    if (scenario=="UMa" && fc<6)
        % Note 6 in TR 38.901 Table 7.5-6 Part-1
        % Note 4 in TR 38.901 Table 7.5-7
        fc = 6;
    elseif (scenario=="UMi" && fc<2)
        % Note 7 in TR 38.901 Table 7.5-6 Part-1
        fc = 2;
    elseif (scenario=="InH" && fc<6)
        % Note 6 in TR 38.901 Table 7.5-6 Part-2
        fc = 6;
    elseif (scenario=="InF" && fc<6)
        % Note 4 in TR 38.901 Table 7.5-10
        fc = 6;
    end

    % Get frequency-dependent values of Table 7.5-6
    if ~indoor && scenario=="UMa"
        p.mu_lgDS = p.mu_lgDS(fc);
        p.mu_lgASD = p.mu_lgASD(fc);
        p.C_DS = p.C_DS(fc);
        if ~LOS
            p.mu_lgASA = p.mu_lgASA(fc);
            p.mu_lgZSA = p.mu_lgZSA(fc);
        end
    elseif (~indoor && scenario=="UMi")
        p.mu_lgDS = p.mu_lgDS(fc);
        p.mu_lgASD = p.mu_lgASD(fc);
        p.mu_lgASA = p.mu_lgASA(fc);
        p.sigma_lgASA = p.sigma_lgASA(fc);
        p.mu_lgZSA = p.mu_lgZSA(fc);
        p.sigma_lgZSA = p.sigma_lgZSA(fc);
        if ~LOS
            p.sigma_lgDS = p.sigma_lgDS(fc);
            p.sigma_lgASD = p.sigma_lgASD(fc);
        end
    elseif (scenario=="InH")
        p.mu_lgDS = p.mu_lgDS(fc);
        p.mu_lgASA = p.mu_lgASA(fc);
        p.sigma_lgASA = p.sigma_lgASA(fc);
        p.mu_lgZSA = p.mu_lgZSA(fc);
        p.sigma_lgZSA = p.sigma_lgZSA(fc);
        if ~LOS
            p.sigma_lgDS = p.sigma_lgDS(fc);
        end
    elseif (scenario=="InF")
        if (LOS)
            p.mu_lgASA = p.mu_lgASA(fc);
            p.sigma_lgASA = p.sigma_lgASA(fc);
        end
        p.mu_lgZSA = p.mu_lgZSA(fc);
    end

    % Get shadow fading
    if ~indoor
        if (scenario=="InF")
            % Provide clutter and height specification part of scenario
            p.sigma_SF = p.sigma_SF(xx);
        end
        p.sigma_SF = p.sigma_SF(fc,bspos(:),uepos(:));
    end

    % Get volume and surface area-dependent values of Table 7.5-6 part 3
    if (scenario=="InF")
        [V,S] = getHallProperties(chCfg);
        p.mu_lgDS = p.mu_lgDS(V,S);
    end

    % Get dependent values of Tables 7.5-7 to 7.5-11
    h_UT = uepos(3);
    h_BS = bspos(3);
    v_2D = uepos(1:2) - bspos(1:2);
    d_2D = vecnorm(v_2D);
    if scenario=="UMa" % Table 7.5-7
        if ~indoor
            p.mu_lgZSD = p.mu_lgZSD(d_2D,h_UT);
            if ~LOS
                p.mu_offsetZOD = p.mu_offsetZOD(d_2D,h_UT,fc);
            end
        else
            p.mu_lgZSD = p.mu_lgZSD(LOS,d_2D,h_UT);
            p.sigma_lgZSD = p.sigma_lgZSD(LOS);
            p.mu_offsetZOD = p.mu_offsetZOD(LOS,d_2D,h_UT,fc);
        end
    elseif scenario=="UMi" % Table 7.5-8
        if ~indoor
            p.mu_lgZSD = p.mu_lgZSD(d_2D,h_UT,h_BS);
            if ~LOS
                p.mu_offsetZOD = p.mu_offsetZOD(d_2D);
            end
        else
            p.mu_lgZSD = p.mu_lgZSD(LOS,d_2D,h_UT,h_BS);
            p.mu_offsetZOD = p.mu_offsetZOD(LOS,d_2D);
        end
    elseif scenario=="RMa" % Table 7.5-9
        p.mu_lgZSD = p.mu_lgZSD(d_2D,h_UT);
        if indoor || ~LOS
            p.mu_offsetZOD = p.mu_offsetZOD(d_2D);
        end
    elseif scenario=="InH" % Table 7.5-10
        if (LOS)
            p.mu_lgZSD = p.mu_lgZSD(d_2D,h_UT,h_BS,fc);
            p.sigma_lgZSD = p.sigma_lgZSD(fc);
        end
    elseif scenario=="InF" % Table 7.5-11
        % Nothing required here, all values are constants
    end

    % Generate C_ZSD as per equation (7.5-20) of TR 38.901
    p.C_ZSD = (3/8)*10^p.mu_lgZSD;

    % Replace with 0s all the empty values which correspond to N/A in Table
    % 7.5-6. These values are related to K for the NLOS and O2I cases.
    f = fieldnames(p);
    f = f(structfun(@(x)isempty(x),p));
    for idx = 1:numel(f)
        p.(f{idx}) = 0;
    end

end

function [pathDelays,pathDelaysScaling,tau_prime] = generateDelays(rs,r_tau,ds,N,indoor,LOS,K,sc)
    % Generate cluster delays for the current BS-UE link, as discussed in Step 5 of TR 38.901 Section 7.5.
    % rs     - Random stream
    % r_tau  - Delay distribution proportionality factor
    % ds     - Delay spread
    % N      - Number of clusters
    % indoor - O2I flag
    % LOS    - LOS condition
    % K      - Ricean K-factor

    spatialConsistency = ~isempty(sc);
    if (spatialConsistency)
        rvs = uniformClusterSCRVs(sc.Delays,sc,N,indoor,LOS);
    else
        rvs = rs.rand(1,N);
    end
    tau_prime = -r_tau*ds*log(rvs); % Row vector, as per nrCDLChannel requirement - Eq. (7.5-1)
    pathDelays = sort(tau_prime-min(tau_prime)); % Eq. (7.5-2)
    pathDelaysScaling = 1; % By default, no scaling

    if (LOS && ~indoor)
        % Compute the path delays scaling factor the delays in case of LOS
        % condition
        pathDelaysScaling = 0.7705 - 0.0433*K + 0.0002*K^2 + 0.000017*K^3; % Scaling constant - Eq. (7.5-3)
    end
end

function [pathPowers,pathPowersLOS,pathDelays,KFactorFirstCluster,tau_prime] = generatePowers(chCfg,rs,pathDelays,r_tau,ds,zeta,K,indoor,LOS,sc,tau_prime,isMobilityUpdate)
    % Generate cluster powers for the current BS-UE link, as discussed in step 6 of TR 38.901 Section 7.5.
    % rs         - Random stream
    % pathDelays - Cluster delays (without LOS scaling)
    % r_tau      - Delay distribution proportionality factor
    % ds         - Delay spread
    % zeta       - Per cluster shadowing std in dB

    if (spatialConsistency(chCfg)) % "Static" or "ProcedureA"
        rvs = normalClusterSCRVs(sc.Powers,sc,size(pathDelays,2),indoor,LOS);
    else
        rvs = rs.randn(size(pathDelays));
    end
    shadowing = zeta*rvs;
    P_prime = exp(-pathDelays*(r_tau-1)/(r_tau*ds)).*10.^(-shadowing/10); % Eq. (7.5-5)

    % Eqs. (7.5-6), (7.5-7) and (7.5-8)
    [pathPowers,pathPowersLOS,KFactorFirstCluster] = normalizePowers(P_prime,K,indoor,LOS);

    if (~isProcedureA(chCfg) || ~isMobilityUpdate) % "Static" or initial link creation for "ProcedureA"
        % Remove clusters with less than -25 dB power compared to the
        % maximum cluster power
        pathPowers_dB = 10*log10(pathPowers);
        threshold = max(pathPowers_dB) - 25;
        toRemove = (pathPowers_dB < threshold);
        pathPowers(toRemove) = [];
        pathPowersLOS(toRemove) = [];
        pathDelays(toRemove) = [];
        tau_prime(toRemove) = [];
    end
end

function [pathPowers,pathPowersLOS,KFactorFirstCluster] = normalizePowers(P_prime,K,indoor,LOS)

    % Normalize the cluster powers
    % Eq. (7.5-6) or Eq. (7-6.17a) from Procedure B
    pathPowers = P_prime/sum(P_prime);

    if (LOS && ~indoor)
        Kr = 10^(K/10); % Ricean K-factor in linear scale
        % Power of LOS cluster
        % Eq. (7.5-7) or Eq. (7.6-17b) from Procedure B
        P1los = Kr/(Kr+1);
        % Eq. (7.5-8) or Eq. (7.6-17c) from Procedure B
        pathPowersLOS = 1/(Kr+1)*pathPowers + P1los*[1 zeros(1,size(pathPowers,2)-1)];
        % Calculate the K-factor for the LOS cluster P_1
        P1nlos = pathPowersLOS(1) - P1los;
        if (P1nlos==0)
            % Protect against infinite KFactorFirstCluster, the value here
            % is finite but will result in zero power in the NLOS part of
            % the first cluster in nrCDLChannel
            KFactorFirstCluster = fix(10*log10(realmax));
        else
            KFactorFirstCluster = 10*log10(P1los/P1nlos);
        end
    else
        pathPowersLOS = pathPowers;
        KFactorFirstCluster = -Inf;
    end

end

function [phi_AoA, phi_AoD, theta_ZoA, theta_ZoD] = generateAngles(rs,bspos,uepos,asa,asd,zsa,zsd,pathPowers,indoor,LOS,K,params,fastFading,sc)
    % Generate arrival and departure angles for each cluster for the
    % current BS-UE link. These angles represent the AnglesAOA, AnglesAOD,
    % AnglesZOA, and AnglesZOD properties of nrCDLChannel and are the
    % variables phi_{n,AOA}, phi_{n,AOD}, theta_{n,ZOA}, and theta_{n,ZOD}
    % in TR 38.901 Eqs. (7.5-13), (7.5-18), and (7.5-20), respectively.
    % That is, the output angles here do not include the addition of the
    % ray offset angles described in the equations (the ray offset angles
    % are the term C*alpha_m where 'C' is the cluster-wise angle spread.
    % Addition of the ray offset angles is implemented in nrCDLChannel,
    % which completes the implementation of Eqs. (7.5-13), (7.5-18), and
    % (7.5-20). If LOS, the input pathPowers has the values of the cluster
    % powers defined in Eq. (7.5-8), otherwise the values are those from
    % Eq. (7.5-6).

    % Determine LOS angles
    [phiLOS_AoA, phiLOS_AoD, thetaLOS_ZoA, thetaLOS_ZoD] = generateLOSAngles(bspos,uepos);

    if ~fastFading
        % Only LOS angles are needed
        phi_AoA = phiLOS_AoA;
        phi_AoD = phiLOS_AoD;
        theta_ZoA = thetaLOS_ZoA;
        theta_ZoD = thetaLOS_ZoD;
    else
        % Get number of clusters
        N = params.NumClusters;

        % Generate AoA
        phi_AoA = generateAzimuthAngle(rs,asa,phiLOS_AoA,pathPowers,indoor,LOS,K,N,sc,'AnglesAoA');

        % Generate AoD
        phi_AoD = generateAzimuthAngle(rs,asd,phiLOS_AoD,pathPowers,indoor,LOS,K,N,sc,'AnglesAoD');

        % Generate ZoA
        theta_ZoA = generateZenithAngle(rs,zsa,thetaLOS_ZoA,pathPowers,indoor,LOS,K,N,[],sc,'AnglesZoA');

        % Generate ZoD
        theta_ZoD = generateZenithAngle(rs,zsd,thetaLOS_ZoD,pathPowers,indoor,LOS,K,N,params.mu_offsetZOD,sc,'AnglesZoD');
    end
end

function [phiLOS_AoA, phiLOS_AoD, thetaLOS_ZoA, thetaLOS_ZoD] = generateLOSAngles(bspos,uepos)

    % Determine LOS angles, in degrees, between BS and UE
    d = uepos-bspos;
    if (norm(d)~=0)
        d = d/norm(d);
    end
    thetaLOS_ZoD = zenithAngle(d);
    phiLOS_AoD = azimuthAngle(d);
    thetaLOS_ZoA = zenithAngle(-d);
    phiLOS_AoA = azimuthAngle(-d);

end

function theta = zenithAngle(rho_hat)

    theta = real(acosd(rho_hat(3)));

end

function phi = azimuthAngle(rho_hat)

    theta = zenithAngle(rho_hat);
    phi = atan2d(rho_hat(2),rho_hat(1));
    phi(theta==0) = 0; % Set ambiguous phi to 0

end

function phi = generateAzimuthAngle(rs,angleSpread,phiLOS,pathPowers,indoor,LOS,K,numClusters,sc,sspName)
    % Generate azimuth angle phi in degrees

    % Scaling factors for the generation of AnglesAoA and AnglesAoD
    C_phi_NLOS_table = [4,     5,     8,     10,    11,    12,    14,    15,    16,    19,    20,    25; ...
                        0.779, 0.860, 1.018, 1.090, 1.123, 1.146, 1.190, 1.211, 1.226, 1.273, 1.289, 1.358]; % Table 7.5-2
    C_phi_NLOS = C_phi_NLOS_table(2, C_phi_NLOS_table(1,:)==numClusters);
    if (LOS && ~indoor)
        C_phi = C_phi_NLOS*(1.1035-0.028*K-0.002*K^2+0.0001*K^3); % Eq. (7.5-10)
    else % NLOS/O2I
        C_phi = C_phi_NLOS; % Eq. (7.5-10)
    end

    % Compute azimuth angle
    phi1 = (2*(angleSpread/1.4)*sqrt(-log(pathPowers/max(pathPowers))))/C_phi; % Eq. (7.5-9)
    [X,Y] = getRandomAngleComponents(rs,angleSpread,size(pathPowers,2),indoor,LOS,sc,sspName);
    if (LOS && ~indoor)
        phi = X.*phi1 + Y - (X(1)*phi1(1) + Y(1) - phiLOS); % Eq. (7.5-12)
    else % NLOS/O2I
        phi = X.*phi1 + Y + phiLOS; % Eq. (7.5-11)
    end
end

function theta = generateZenithAngle(rs,angleSpread,thetaLOS,pathPowers,indoor,LOS,K,numClusters,ZOD_offset,sc,sspName)
    % Generate zenith angle theta in degrees

    % Scaling factors for the generation of AnglesZoA and AnglesZoD
    C_theta_NLOS_table = [8,     10,    11,    12,    15,     19,    20,    25; ...
                          0.889, 0.957, 1.031, 1.104, 1.1088, 1.184, 1.178, 1.282]; % Table 7.5-4
    C_theta_NLOS = C_theta_NLOS_table(2, C_theta_NLOS_table(1,:)==numClusters);
    if (LOS && ~indoor)
        C_theta = C_theta_NLOS*(1.3086+0.0339*K-0.0077*K^2+0.0002*K^3); % Eq. (7.5-15)
    else % NLOS/O2I
        C_theta = C_theta_NLOS; % Eq. (7.5-15)
    end

    % Compute zenith angle
    theta1 = -angleSpread*log(pathPowers/max(pathPowers))/C_theta; % Eq. (7.5-14)
    [X,Y] = getRandomAngleComponents(rs,angleSpread,size(pathPowers,2),indoor,LOS,sc,sspName);
    thetaTmp = X.*theta1 + Y;
    if (LOS && ~indoor)
        theta = thetaTmp - (thetaTmp(1) - thetaLOS); % Eq. (7.5-17)
    else
        if isempty(ZOD_offset) % ZOA
            if indoor
                theta_bar = 90;
            else
                theta_bar = thetaLOS;
            end
            theta = thetaTmp + theta_bar; % Eq. (7.5-16)
        else % ZOD
            theta = thetaTmp + thetaLOS + ZOD_offset; % Eq. (7.5-19)
        end
    end
end

function [X,Y] = getRandomAngleComponents(rs,angleSpread,N,indoor,LOS,sc,sspName)
    % Get the random variables X and Y needed in the computation of the
    % azimuth and zenith angles in Step 7.

    % Random uniformly-distributed variable in the set {-1, 1}
    spatialConsistency = ~isempty(sc);
    if (spatialConsistency)
        thisSSP = sc.(sspName).Sign;
        rvs = uniformClusterSCRVs(thisSSP,sc,N,indoor,LOS);
    else
        rvs = rs.rand(1,N);
    end
    X = 2*(rvs>0.5)-1;

    % Random normally-distributed variable with variance (angleSpread/7)^2
    if (spatialConsistency)
        thisSSP = sc.(sspName).Offset;
        rvs = normalClusterSCRVs(thisSSP,sc,N,indoor,LOS);
    else
        rvs = rs.randn(1,N);
    end
    Y = angleSpread/7*rvs;
end

function coupling = generateRayCoupling(rs,L,M)
    % Generate random ray coupling matrix, as discussed in Step 8 of TR 38.901 Section 7.5.

    siz = [L M 3];
    rn = rs.rand(siz);
    [~,coupling] = sort(rn,2);
end

function initPhase = generateInitialPhases(rs,L,M)
    % Generate random initial phases matrix, as discussed in Step 10 of TR 38.901 Section 7.5.
    % Note that nrCDLChannel accepts initial phases in degrees.
    
    siz = [L M 4];
    initPhase = (rs.rand(siz)*360) - 180; % degrees
end

function tables38901_fastFading = get38901FastFadingTables()
    % Implementation of Tables 7.5-6 to 7.5-11 of TR 38.901. The output
    % here is a single table containing the union of all the tables from TR
    % 38.901. The output is transposed wrt the tables in TR 38.901 (that
    % is, scenario and LOS/NLOS/O2I state are the row headings, parameters
    % are the column headings) for faster access.

    %                       UMi LOS                            UMi NLOS             UMi O2I              UMa LOS                          UMa NLOS            UMa O2I          RMa LOS         RMa NLOS  RMa O2I                InH LOS                          InH NLOS                       InF LOS                             InF NLOS 
    % Table 7.5-6
    % Delay spread (DS)
    mu_lgDS     = {@(fc)(-0.240*log10(1+fc)-7.14);  @(fc)(-0.240*log10(1+fc)-6.83);  -6.62;  @(fc)(-6.955-0.0963*log10(fc));  @(fc)(-6.280-0.2040*log10(fc));  -6.62;           -7.49;            -7.43;  -7.47;  @(fc)(-0.01*log10(1+fc)-7.692);  @(fc)(-0.28*log10(1+fc)-7.173);  @(V,S)(log10(26*(V/S)+14)-9.35);  @(V,S)(log10(30*(V/S)+32)-9.44)}; 
    sigma_lgDS  = {                          0.38;  @(fc)( 0.160*log10(1+fc)+0.28);   0.32;                            0.66;                            0.39;   0.32;            0.55;             0.48;   0.24;                            0.18;   @(fc)(0.10*log10(1+fc)+0.055);                             0.15;                             0.19};
    % AOD spread (ASD)
    mu_lgASD    = {@(fc)(-0.050*log10(1+fc)+1.21);  @(fc)(-0.230*log10(1+fc)+1.53);   1.25;  @(fc)( 1.060+0.1114*log10(fc));  @(fc)( 1.500-0.1144*log10(fc));   1.25;            0.90;             0.95;   0.67;                            1.60;                            1.62;                             1.56;                             1.57};
    sigma_lgASD = {                          0.41;  @(fc)( 0.110*log10(1+fc)+0.33);   0.42;                            0.28;                            0.28;   0.42;            0.38;             0.45;   0.18;                            0.18;                            0.25;                             0.25;                              0.2};
    % AOA spread (ASA)
    mu_lgASA    = {@(fc)(-0.080*log10(1+fc)+1.73);  @(fc)(-0.080*log10(1+fc)+1.81);   1.76;                            1.81;  @(fc)( 2.080-0.2700*log10(fc));   1.76;            1.52;             1.52;   1.66;  @(fc)(-0.19*log10(1+fc)+1.781);  @(fc)(-0.11*log10(1+fc)+1.863);    @(fc)(-0.18*log10(1+fc)+1.78);                             1.72};
    sigma_lgASA = {@(fc)( 0.014*log10(1+fc)+0.28);  @(fc)( 0.050*log10(1+fc)+0.30);   0.16;                            0.20;                            0.11;   0.16;            0.24;             0.13;   0.21;   @(fc)(0.12*log10(1+fc)+0.119);   @(fc)(0.12*log10(1+fc)+0.059);      @(fc)(0.12*log10(1+fc)+0.2);                              0.3};
    % ZOA spread (ZSA)
    mu_lgZSA    = {@(fc)(-0.100*log10(1+fc)+0.73);  @(fc)(-0.040*log10(1+fc)+0.92);   1.01;                            0.95;  @(fc)( 1.512-0.3236*log10(fc));   1.01;            0.47;             0.58;   0.93;   @(fc)(-0.26*log10(1+fc)+1.44);  @(fc)(-0.15*log10(1+fc)+1.387);      @(fc)(-0.2*log10(1+fc)+1.5);    @(fc)(-0.13*log10(1+fc)+1.45)};
    sigma_lgZSA = {@(fc)(-0.040*log10(1+fc)+0.34);  @(fc)(-0.070*log10(1+fc)+0.41);   0.43;                            0.16;                            0.16;   0.43;            0.40;             0.37;   0.22;  @(fc)(-0.04*log10(1+fc)+0.264);  @(fc)(-0.09*log10(1+fc)+0.746);                             0.35;                             0.45};
    % Shadow fading (SF) [dB]
    SF = @(s,LOS)@(fc,bspos,uepos)getShadowFading(s,LOS,fc,bspos,uepos);
    sigma_SF    = {                SF("UMi",true);                 SF("UMi",false);      7;                  SF("UMa",true);                 SF("UMa",false);      7;  SF("RMa",true);  SF("RMa",false);      8;                  SF("InH",true);                 SF("InH",false);          @(xx)SF("InF-"+xx,true);         @(xx)SF("InF-"+xx,false)};
    % K-factor (K) [dB]
    mu_K        = {                             9;                              [];     [];                               9;                              [];     [];               7;               [];     [];                               7;                              [];                                7;                               []};
    sigma_K     = {                             5;                              [];     [];                             3.5;                              [];     [];               4;               [];     [];                               4;                              [];                                8;                               []};
    % Cross-Correlations
    ASDvsDS     = {                           0.5;                               0;    0.4;                             0.4;                             0.4;    0.4;               0;             -0.4;      0;                             0.6;                             0.4;                                0;                                0};
    ASAvsDS     = {                           0.8;                             0.4;    0.4;                             0.8;                             0.6;    0.4;               0;                0;      0;                             0.8;                               0;                                0;                                0};
    ASAvsSF     = {                          -0.4;                            -0.4;      0;                            -0.5;                               0;      0;               0;                0;      0;                            -0.5;                            -0.4;                                0;                                0};
    ASDvsSF     = {                          -0.5;                               0;    0.2;                            -0.5;                            -0.6;    0.2;               0;              0.6;      0;                            -0.4;                               0;                                0;                                0};
    DSvsSF      = {                          -0.4;                            -0.7;   -0.5;                            -0.4;                            -0.4;   -0.5;            -0.5;             -0.5;      0;                            -0.8;                            -0.5;                                0;                                0};
    ASDvsASA    = {                           0.4;                               0;      0;                               0;                             0.4;      0;               0;                0;   -0.7;                             0.4;                               0;                                0;                                0};
    ASDvsK      = {                          -0.2;                              [];     [];                               0;                              [];     [];               0;               [];     [];                               0;                              [];                             -0.5;                               []};
    ASAvsK      = {                          -0.3;                              [];     [];                            -0.2;                              [];     [];               0;               [];     [];                               0;                              [];                                0;                               []};
    DSvsK       = {                          -0.7;                              [];     [];                            -0.4;                              [];     [];               0;               [];     [];                            -0.5;                              [];                              0.7;                               []};
    SFvsK       = {                           0.5;                              [];     [];                               0;                              [];     [];               0;               [];     [];                             0.5;                              [];                                0;                               []};
    ZSDvsSF     = {                             0;                               0;      0;                               0;                               0;      0;            0.01;            -0.04;      0;                             0.2;                               0;                                0;                                0};
    ZSAvsSF     = {                             0;                               0;      0;                            -0.8;                            -0.4;      0;           -0.17;            -0.25;      0;                             0.3;                               0;                                0;                                0};
    ZSDvsK      = {                             0;                              [];     [];                               0;                              [];     [];               0;               [];     [];                               0;                              [];                                0;                               []};
    ZSAvsK      = {                             0;                              [];     [];                               0;                              [];     [];           -0.02;               [];     [];                             0.1;                              [];                                0;                               []};
    ZSDvsDS     = {                             0;                            -0.5;   -0.6;                            -0.2;                            -0.5;   -0.6;           -0.05;            -0.10;      0;                             0.1;                           -0.27;                                0;                                0};
    ZSAvsDS     = {                           0.2;                               0;   -0.2;                               0;                               0;   -0.2;            0.27;            -0.40;      0;                             0.2;                           -0.06;                                0;                                0};
    ZSDvsASD    = {                           0.5;                             0.5;   -0.2;                             0.5;                             0.5;   -0.2;            0.73;             0.42;   0.66;                             0.5;                            0.35;                                0;                                0};
    ZSAvsASD    = {                           0.3;                             0.5;      0;                               0;                            -0.1;      0;           -0.14;            -0.27;   0.47;                               0;                            0.23;                                0;                                0};
    ZSDvsASA    = {                             0;                               0;      0;                            -0.3;                               0;      0;           -0.20;            -0.18;  -0.55;                               0;                           -0.08;                                0;                                0};
    ZSAvsASA    = {                             0;                             0.2;    0.5;                             0.4;                               0;    0.5;            0.24;             0.26;  -0.22;                             0.5;                            0.43;                                0;                                0};
    ZSDvsZSA    = {                             0;                               0;    0.5;                               0;                               0;    0.5;           -0.07;            -0.27;      0;                               0;                            0.42;                                0;                                0};
    % Delay scaling parameter
    r_tau       = {                             3;                             2.1;    2.2;                             2.5;                             2.3;    2.2;            3.80;             1.70;   1.70;                             3.6;                               3;                              2.7;                                3};
    % XPR [dB]
    mu_XPR      = {                             9;                               8;      9;                               8;                               7;      9;              12;                7;      7;                              11;                              10;                               12;                               11};
    sigma_XPR   = {                             3;                               3;      5;                               4;                               3;      5;               4;                3;      3;                               4;                               4;                                6;                                6};
    % Number of clusters N
    NumClusters = {                            12;                              19;     12;                              12;                              20;     12;              11;               10;     10;                              15;                              19;                               25;                               25};
    % Number of rays per cluster M
    NumRaysPerCluster = {                      20;                              20;     20;                              20;                              20;     20;              20;               20;     20;                              20;                              20;                               20;                               20};
    % Cluster DS [ns]
    C_DS        = {                             5;                              11;     11;         @(fc)(getClusterDS(fc));         @(fc)(getClusterDS(fc));     11;              [];               [];     [];                              [];                              [];                               [];                               []};
    % Cluster ASD [deg]
    C_ASD       = {                             3;                              10;      5;                               5;                               2;      5;               2;                2;      2;                               5;                               5;                                5;                                5};
    % Cluster ASA [deg]
    C_ASA       = {                            17;                              22;      8;                              11;                              15;      8;               3;                3;      3;                               8;                              11;                                8;                                8};
    % Cluster ZSA [deg]
    C_ZSA       = {                             7;                               7;      3;                               7;                               7;      3;               3;                3;      3;                               9;                               9;                                9;                                9};
    % Pe cluster shadowing zeta [dB]
    zeta        = {                             3;                               3;      4;                               3;                               3;      4;               3;                3;      3;                               6;                               3;                                4;                                3};
    % Correlation distance in the horizontal plane [m]
    corr_DS     = {                             7;                              10;     10;                              30;                              40;     10;              50;               36;     36;                               8;                               5;                               10;                               10};
    corr_ASD    = {                             8;                              10;     11;                              18;                              50;     11;              25;               30;     30;                               7;                               3;                               10;                               10};
    corr_ASA    = {                             8;                               9;     17;                              15;                              50;     17;              35;               40;     40;                               5;                               3;                               10;                               10};
    corr_SF     = {                            10;                              13;      7;                              37;                              50;      7;              37;              120;    120;                              10;                               6;                               10;                               10};
    corr_K      = {                            15;                              [];     [];                              12;                              [];     [];              40;               [];     [];                               4;                              [];                               10;                               []};
    corr_ZSA    = {                            12;                              10;     25;                              15;                              50;     25;              15;               50;     50;                               4;                               4;                               10;                               10};
    corr_ZSD    = {                            12;                              10;     25;                              15;                              50;     25;              15;               50;     50;                               4;                               4;                               10;                               10};

    % Tables 7.5-7, 7.5-8, 7.5-9, 7.5-10, 7.5-11
    % ZOD spread (ZSD)
    mu_lgZSD = {
                                        % LOS                                                          NLOS                                                             O2I
              @(d_2D,h_UT,h_BS)getZODSpread_mu("UMi",0,1,d_2D,h_UT,h_BS);  @(d_2D,h_UT,h_BS)getZODSpread_mu("UMi",0,0,d_2D,h_UT,h_BS);  @(LOS,d_2D,h_UT,h_BS)getZODSpread_mu("UMi",0,LOS,d_2D,h_UT,h_BS);   % UMi
                        @(d_2D,h_UT)getZODSpread_mu("UMa",0,1,d_2D,h_UT);            @(d_2D,h_UT)getZODSpread_mu("UMa",0,0,d_2D,h_UT);            @(LOS,d_2D,h_UT)getZODSpread_mu("UMa",0,LOS,d_2D,h_UT);   % UMa
                        @(d_2D,h_UT)getZODSpread_mu("RMa",0,1,d_2D,h_UT);            @(d_2D,h_UT)getZODSpread_mu("RMa",0,0,d_2D,h_UT);                  @(d_2D,h_UT)getZODSpread_mu("RMa",1,1,d_2D,h_UT);   % RMa
        @(d_2D,h_UT,h_BS,fc)getZODSpread_mu("InH",0,1,d_2D,h_UT,h_BS,fc);                                                        1.08;  %                                -                                  % InH
                                                                    1.35;                                                         1.2   %                                -                                  % InF
        };
    sigma_lgZSD = {
                           % LOS                        NLOS                                    O2I
                                       0.35;                         0.35;                                 0.35; % UMi
                getZODSpread_sigma("UMa",1);  getZODSpread_sigma("UMa",0);  @(LOS)getZODSpread_sigma("UMa",LOS); % UMa
                                       0.34;                         0.30;                                 0.30; % RMa
        @(fc)getZODSpread_sigma("InH",1,fc);                         0.36;  %                    -               % InH
                                       0.35;                         0.55   %                    -               % InF
        };
    % ZOD offset
    mu_offsetZOD = {
                % LOS                                       % NLOS                                                      % O2I
        getZODOffset("UMi",0,1);                  @(d_2D)getZODOffset("UMi",0,0,d_2D);                  @(LOS,d_2D)getZODOffset("UMi",0,LOS,d_2D); % UMi
        getZODOffset("UMa",0,1);  @(d_2D,h_UT,fc)getZODOffset("UMa",0,0,d_2D,h_UT,fc);  @(LOS,d_2D,h_UT,fc)getZODOffset("UMa",0,LOS,d_2D,h_UT,fc); % UMa
        getZODOffset("RMa",0,1);                  @(d_2D)getZODOffset("RMa",0,0,d_2D);                        @(d_2D)getZODOffset("RMa",1,0,d_2D); % RMa
                              0;                                                    0;  %                                  -                       % InH
                              0;                                                    0   %                                  -                       % InF
        };

    % Transpose the table so that each row represents a scenario and
    % LOS/NLOS/O2I state
    tables38901_fastFading = table( ...
        mu_lgDS,sigma_lgDS,mu_lgASD,sigma_lgASD,mu_lgASA,sigma_lgASA, ...
        mu_lgZSA,sigma_lgZSA,sigma_SF,mu_K,sigma_K,ASDvsDS,ASAvsDS, ...
        ASAvsSF,ASDvsSF,DSvsSF,ASDvsASA,ASDvsK,ASAvsK,DSvsK,SFvsK, ...
        ZSDvsSF,ZSAvsSF,ZSDvsK,ZSAvsK,ZSDvsDS,ZSAvsDS,ZSDvsASD, ...
        ZSAvsASD,ZSDvsASA,ZSAvsASA,ZSDvsZSA,r_tau,mu_XPR,sigma_XPR, ...
        NumClusters,NumRaysPerCluster,C_DS,C_ASD,C_ASA,C_ZSA,zeta, ... 
        corr_DS,corr_ASD,corr_ASA,corr_SF,corr_K,corr_ZSA,corr_ZSD, ...
        mu_lgZSD,sigma_lgZSD,mu_offsetZOD ...
        );
    tables38901_fastFading.Properties.RowNames = [
        "UMi_LOS", "UMi_NLOS", "UMi_O2I", ...
        "UMa_LOS", "UMa_NLOS", "UMa_O2I", ...
        "RMa_LOS", "RMa_NLOS", "RMa_O2I", ...
        "InH_LOS", "InH_NLOS", ...
        "InF_LOS", "InF_NLOS"
      ];

    function sigma_SF = getShadowFading(scenario,LOS,fc,bspos,uepos)

        persistent plc;

        if (isempty(plc))
            plc = nrPathLossConfig;
        end

        freq = fc*1e9; % nrPathLoss wants frequency in Hz
        [~,sigma_SF] = nrPathLoss(setfield(plc,'Scenario',scenario),freq,LOS,bspos,uepos); %#ok<SFLD>

    end

    function C_DS = getClusterDS(fc)
        % Get the value of the cluster delay spread for UMa with LOS/NLOS
        % propagation
        C_DS = max(0.25, 6.5622 - 3.4084*log10(fc));
    end

    function mu_lgZSD = getZODSpread_mu(scenario,indoor,LOS,d_2D,h_UT,h_BS,fc)
        % Get the value of mu_lgZSD from Tables 7.5-7, 7.5-8, 7.5-9, 7.5-10
        d_2D_km = d_2D/1000; % Convert to km
        switch scenario
            case "UMi"
                if LOS
                    mu_lgZSD = max(-0.21, -14.8*d_2D_km+0.01*abs(h_UT-h_BS) + 0.83);
                else
                    mu_lgZSD = max(-0.5, -3.1*d_2D_km+0.01*max(h_UT-h_BS,0) + 0.2);
                end
            case "UMa"
                if LOS
                    mu_lgZSD = max(-0.5, -2.1*d_2D_km-0.01*(h_UT-1.5) + 0.75);
                else
                    mu_lgZSD = max(-0.5, -2.1*d_2D_km-0.01*(h_UT-1.5) + 0.9);
                end
            case "RMa"
                if indoor
                    mu_lgZSD = max(-1, -0.19*d_2D_km-0.01*(h_UT-1.5) + 0.28);
                elseif LOS
                    mu_lgZSD = max(-1, -0.17*d_2D_km-0.01*(h_UT-1.5) + 0.22);
                else
                    mu_lgZSD = max(-1, -0.19*d_2D_km-0.01*(h_UT-1.5) + 0.28);
                end
            case "InH"
                if LOS
                    mu_lgZSD = -1.43*log10(1+fc) + 2.228;
                end
        end
    end

    function sigma_lgZSD = getZODSpread_sigma(scenario,LOS,fc)
        % Get the value of sigma_lgZSD from Tables 7.5-7 and 7.5-10
        switch scenario
            case "UMa"
                if LOS
                    sigma_lgZSD = 0.40;
                else
                    sigma_lgZSD = 0.49;
                end
            case "InH"
                if LOS
                    sigma_lgZSD = 0.13*log10(1+fc) + 0.30;
                end
        end
    end

    function mu_offsetZOD = getZODOffset(scenario,indoor,LOS,d_2D,h_UT,fc)
        % Get the value of mu_offsetZOD from Tables 7.5-7, 7.5-8, 7.5-9
        % fc in GHz, d_2D in m
        switch scenario
            case "UMi"
                if LOS
                    mu_offsetZOD = 0;
                else
                    mu_offsetZOD = 10^(-1.5*log10(max(10,d_2D)) + 3.3);
                end
            case "UMa"
                if LOS
                    mu_offsetZOD = 0;
                else
                    a = @(fc)(0.208*log10(fc) - 0.782);
                    b = 25;
                    c = @(fc)(-0.13*log10(fc) + 2.03);
                    e = @(fc)(7.66*log10(fc) - 5.96);
                    mu_offsetZOD = e(fc) - 10^(a(fc)*log10(max(b,d_2D)) + c(fc) - 0.07*(h_UT-1.5));
                end
            case "RMa"
                if indoor || ~LOS
                    mu_offsetZOD = atan((35-3.5)/d_2D) - atan((35-1.5)/d_2D);
                else % LOS
                    mu_offsetZOD = 0;
                end
        end
    end
end

function [pathGains,sampleTimes] = spanSlot(obj,bsID,ueID,pathGains,sampleTimes)

    if (~isscalar(sampleTimes))

        % Get OFDM information and calculate sample time Ts
        if (~isnan(bsID))
            node = obj.theBSMap(bsID);
        else
            node = obj.theUEMap(ueID);
        end
        ofdmInfo = node.OFDMInfo;
        Ts = 1 / ofdmInfo.SampleRate;

        % Calculate the whole number of subframes elapsed up to the start
        % of the packet, and the corresponding number of samples
        wholeSubframes = floor((sampleTimes(1)+Ts) / 1e-3);
        samplesPerSubframe = 1e-3 / Ts;
        wholeSubframeSamples = wholeSubframes * samplesPerSubframe;

        % Calculate the start and end sample index of each slot of the
        % subframe
        samplesPerSlot = sum(reshape(ofdmInfo.SymbolLengths,ofdmInfo.SymbolsPerSlot,[]),1);
        slotStartSamples = wholeSubframeSamples + cumsum([0 samplesPerSlot(1:end-1)]);
        slotEndSamples = slotStartSamples + samplesPerSlot;

        % Get the start and end sample index of the packet
        packetStartSample = sampleTimes(1) / Ts;
        packetEndSample = sampleTimes(end) / Ts;

        % Find the start and end sample index of the slot which contains
        % the packet
        deltaStart = abs(packetStartSample - slotStartSamples);
        deltaEnd = abs(packetEndSample - slotEndSamples);
        [deltaStart_min,slotIdxStart] = min(deltaStart);
        [deltaEnd_min,slotIdxEnd] = min(deltaEnd);
        if (deltaStart_min < deltaEnd_min)
            slotIdx = slotIdxStart;
        else
            slotIdx = slotIdxEnd;
        end
        slotStartSample = slotStartSamples(slotIdx);
        slotEndSample = slotEndSamples(slotIdx);

        % Extend the path gains and sample times if they do not encompass
        % the start or end of the slot
        if (packetStartSample > slotStartSample)
            slotStartTime = slotStartSample * Ts;
            if (slotStartTime~=sampleTimes(1))
                sampleTimes = [slotStartTime; sampleTimes];
                pathGains = [pathGains(1,:,:,:); pathGains];
            end
        end
        if (packetEndSample < slotEndSample)
            slotEndTime = slotEndSample * Ts;
            if (slotEndTime~=sampleTimes(end))
                sampleTimes = [sampleTimes; slotEndTime];
                pathGains = [pathGains; pathGains(end,:,:,:)];
            end
        end

    end

end

function y = isCellular(s)

    y = any(s==["UMi","UMa","RMa"]);

end

function [V,S] = getHallProperties(chCfg)

    siz = chCfg.HallSize;
    V = prod(siz);
    S = 2 * sum(prod(siz([1 2; 1 3; 2 3]),2));

end

function h = cachedPathFilters(sampleDelays)

    % Set up and cache path filters for a set of different fractional
    % delays
    persistent fh fd;
    if (isempty(fh))
        cdl = nrCDLChannel;
        cdl.SampleRate = 1;
        cdl.DelayProfile = 'Custom';
        fd = 0:0.02:1.00;
        cdl.PathDelays = fd;
        cdl.AveragePathGains = zeros(size(fd));
        cdl.AnglesAoA = zeros(size(fd));
        cdl.AnglesAoD = zeros(size(fd));
        cdl.AnglesZoA = zeros(size(fd));
        cdl.AnglesZoD = zeros(size(fd));
        fh = getPathFilters(cdl).';
    end

    % Separate the path sample delays into fractional and integer parts
    sampleDelays = sampleDelays.';
    fracDelays = mod(sampleDelays,1);
    intDelays = sampleDelays - fracDelays;

    % For each path, establish which cached fractional delay filter has the
    % nearest delay to the path delay, and create a matrix of these filters
    % for all the paths
    [~,fdidx] = min(abs(fracDelays - fd),[],2);
    hfd = fh(fdidx,:);

    % For each path, establish the position of the last fractional delay
    % filter tap when also considering the required integer delay, and use
    % the maximum of these positions across the paths to set the size of
    % the output path filter matrix
    np = size(hfd,1);
    ntaps = max(double(hfd~=0) .* (1:size(hfd,2)),[],2);
    lasttap = intDelays + ntaps;
    nh = max(lasttap);

    % For each path, position the fractional delay filter in the output
    % path filter matrix according to the required integer delay
    h = zeros(np,nh);
    for i = 1:np
        h(i,intDelays(i)+(1:ntaps(i))) = hfd(i,1:ntaps(i));
    end

end

% Calculate path gains corresponding to LOS channel ray
function pathgains = noFastFadingPathGainsLocal(ch)

    % Create persistent variables
    persistent models;
    persistent dbar_tx;
    persistent dbar_rx;
    persistent prevch;

    % If the function is called with no input arguments, or the current
    % input arguments imply different antenna arrays from the previous
    % call, clear the persistent variables
    if (nargin==0 || ~sameAntennaArray(ch,prevch))
        models = [];
        dbar_tx = [];
        dbar_rx = [];
    end
    if (nargin==0)
        pathgains = [];
        return;
    end
    prevch = ch;

    % Get the small scale part of the channel (will be a structure rather
    % than an nrCDLChannel, because it was created with fastFading=false)
    ss = ch.SmallScale;

    % Calculate antenna element positions and orientations
    % NOTE: assumes that the antenna arrays are the same for all links 
    % except that the transmit antenna array orientation differs for each 
    % sector in a site, and that all sites have the same sectorization
    if (isempty(models))
        models = cell(1,ch.NodeSiz(2));
    end
    sector = ch.NodeSubs(2);
    if (isempty(models{sector}))
        models{sector} = nr5g.internal.nrCDLChannel.makeCDLChannelAntennaArrayStructure(ss.TransmitAntennaArray,ss.ReceiveAntennaArray,ss.TransmitArrayOrientation,ss.ReceiveArrayOrientation,ss.CarrierFrequency);
    end
    model = models{sector};

    % Calculate location vectors for antenna elements
    % NOTE: assumes that element positions are the same for all links
    if (isempty(dbar_tx))
        dbar_tx = locationVectors(model.TransmitAntennaArray);
        dbar_rx = locationVectors(model.ReceiveAntennaArray);
    end

    % Calculate field terms and spherical unit vectors corresponding to LOS
    % azimuth / elevation angles, for the transmit and receive antenna
    % arrays
    S = model.NumInputSignals;
    [F_tx,rhat_tx] = angleDependentTerms(model.TransmitAntennaArray,ss.AnglesAoD,ss.AnglesZoD,S);
    U = model.NumOutputSignals;
    [F_rx,rhat_rx] = angleDependentTerms(model.ReceiveAntennaArray,ss.AnglesAoA,ss.AnglesZoA,U);

    % Calculate the LOS path gains for all transmit and receive antennas
    pathgains = zeros(S,U);
    for s = 1:S
        for u = 1:U
            fieldTerm = ((F_rx(1,u)) .* F_tx(1,s) - (F_rx(2,u)) .* F_tx(2,s)).';
            pathgains(s,u) = fieldTerm .* locationTerm(rhat_rx,dbar_rx,u) .* locationTerm(rhat_tx,dbar_tx,s);
        end
    end
    pathgains = permute(pathgains,[3 4 1 2]);

end

% Determine if transmit and receive antenna arrays in 'x' and 'y' will
% result in the same antenna location vectors and antenna orientations
function same = sameAntennaArray(x,y)

    if (isempty(y))
        same = false;
        return;
    end

    xss = x.SmallScale;
    yss = y.SmallScale;
    same = ...
        (x.NodeSiz(2)==y.NodeSiz(2)) && ...
        isequal(xss.TransmitAntennaArray,yss.TransmitAntennaArray) && ...
        isequal(xss.ReceiveAntennaArray,yss.ReceiveAntennaArray) && ...
        (xss.CarrierFrequency==yss.CarrierFrequency);
    % Transmit antenna array orientation is not relevant because
    % 'noFastFadingPathGainsLocal' assumes that the orientation only
    % depends on the sector, and variables are cached for each sector
    raa = xss.ReceiveAntennaArray;
    if (same && (~isstruct(raa) || ~strcmpi(raa.Element,'isotropic')))
        % Receive antenna array orientation is only relevant if the receive
        % antenna element pattern is a PhAST object or is not isotropic
        same = isequal(xss.ReceiveArrayOrientation,yss.ReceiveArrayOrientation);
    end

end

% Vector of locations of all elements in an antenna array
function dbar = locationVectors(a)

    if ~isempty(a.SubarrayPositions) % Partitioned/Replicated subarrays
        radiatorPositions = a.SubarrayPositions;
    else
        radiatorPositions = a.ElementPositions;
    end

    % NOTE: assumes a.Position = [0; 0; 0]
    dbar = reshape(radiatorPositions,3,[]);

end

% Field terms and spherical unit vectors corresponding to LOS azimuth /
% elevation angles, for given antenna array, angles, and number of antennas
function [F,rhat] = angleDependentTerms(aa,AoX,ZoX,N)

    F = arrayfun(@(s)wireless.internal.channelmodels.getFieldTerm(aa,ZoX,AoX,s),1:N,'UniformOutput',false);
    F = cat(2,F{:});
    rhat = wireless.internal.channelmodels.getRhoHat(AoX,ZoX);

end

% Location-related channel term, a function of the spherical unit vectors
% corresponding to LOS azimuth / elevation angles and element locations
function l = locationTerm(rhat,dbar,ant)

    l = exp(1i*2*pi*rhat.'*dbar(:,ant));

end

% Apply phase term in Eq. 7.5-29 related to d_3D
function pathGains = applyPhaseLOS_d_3D(ch,pathGains,los)

    persistent c;
    if (isempty(c))
        c = physconst('LightSpeed');
    end

    lambda_0 = c / ch.CenterFrequency;
    d_3D = getDistance3D(ch.ChannelConfiguration);
    pathGains(:,los,:,:) = pathGains(:,los,:,:) * exp(-1i*2*pi*d_3D/lambda_0);

end

function d_3D = getDistance3D(chCfg)

    bspos = chCfg.BSPosition;
    uepos = chCfg.UEPosition;
    d_3D = vecnorm(bspos - uepos);

end

function checkUEPosition(extents,uepos)

    minpos = extents(1:2);
    maxpos = minpos + extents(3:4);
    if (any(uepos(1:2) < minpos) || any(uepos(1:2) > maxpos))
        error('nr5g:h38901Channel:UEOutsideSystem','UE position x=%0.3f, y=%0.3f is outside the system boundary x=%0.3f...%0.3f, y=%0.3f...%0.3f.',uepos(1),uepos(2),minpos(1),maxpos(1),minpos(2),maxpos(2));
    end

end

function y = spatialConsistency(s)

    if (isnumeric(s.SpatialConsistency) || islogical(s.SpatialConsistency))
        y = s.SpatialConsistency;
    else
        y = ~strcmpi(s.SpatialConsistency,"None");
    end

end

function y = isProcedureA(s)

    y = spatialConsistency(s) && strcmpi(s.SpatialConsistency,"ProcedureA");

end

function y = isProcedureB(s)

    y = spatialConsistency(s) && strcmpi(s.SpatialConsistency,"ProcedureB");

end

function y = isSpatiallyConsistentMobility(s)

    y = spatialConsistency(s) && any(strcmpi(s.SpatialConsistency,["ProcedureA" "ProcedureB"]));

end

function obj = mobilityManager(cfgin)

    cfg = cfgin;

    % These variables apply to update t_k:
    k = [];       % index of update
    d = [];       % distance accumulator
    state_k = []; % state

    obj.update = @update;
    obj.getState = @getState;

    function varargout = update(ch,t,tx,rx)

        % Check that spatial consistency of channel link configuration
        % matches spatial consistency procedure configured here
        chSpatialConsistency = ch.ChannelConfiguration.SpatialConsistency;
        if (~strcmpi(chSpatialConsistency,cfg.SpatialConsistency))
            error('nr5g:h38901Channel:ChangedSpatialConsistency','Spatial consistency configured for channel link (%s) must match spatial consistency procedure configured when calling h38901Channel.spatiallyConsistentMobility (%s).',chSpatialConsistency,cfg.SpatialConsistency);
        end

        % Get update distance. Note that TR 38.901 only describes an
        % "update distance" for Procedure A. However, the implementation
        % here also applies the "update distance" to Procedure B, in order
        % that the channel is not updated excessively frequently in the
        % case of low speed mobility
        updateDistance = cfg.UpdateDistance;
        absoluteTOA = ~isempty(ch.Mobility.Deltatau);
        if (absoluteTOA)
            % "In case the absolute time of arrival modeling of clause 
            % 7.6.9 is simultaneously used, the update distance should be 
            % the minimum of 1 meter and d_3D/10"
            d_3D = getDistance3D(ch.ChannelConfiguration);
            updateDistance = min(updateDistance,d_3D/10);
        end

        if (isempty(ch.SmallScale))

            % If the small scale channel is not enabled, return without
            % updating the channel
            updated = false;

        elseif (isempty(k)) % first input, t_0

            % NOTE: TR 38.901 defines that t_0 (i.e. t_k with k = 0) occurs
            % at t = 0, but the current input start time may be non-zero if
            % this first input occurs after t = 0. For consistency the code
            % here lets t_0 occur at that time

            % Initialize 'k' and the distance accumulator 'd'
            k = 0;
            d = 0;

            % Store the current input state as the state for update t_0
            validateInfo(tx,'TX');
            validateInfo(rx,'RX');
            state_k.in = makeInputState(t,tx,rx);

            % Perform any channel setup required at t_0
            t_update = 0;
            [ch,state_k] = updateSmallScaleChannel(cfg,ch,t,k,state_k,t_update);
            updated = true;

        else % inputs subsequent to t_0

            % NOTE: TR 38.901 talks about a current update t_k and previous
            % update t_{k-1}, but the code here establishes if the distance
            % travelled is sufficient to update 'k'. Therefore the process
            % is described in terms of the most recent update t_k and the
            % next update t_{k+1}

            % 't' is the current input time, check that it is greater than
            % or equal to the time of the most recent update
            if (t < state_k.in.time)
                error('nr5g:h38901Channel:NegativeTimeDifference','The current simulation time (%0.9gs) must be greater than or equal to the previous simulation time (%0.9gs).',t,state_k.in.time);
            end

            % 'state_k' is:
            % * the state at exactly t_k if the accumulated distance 'd' is
            %   zero
            % -or-
            % * the state at some point between t_k and t_{k+1} if 'd' is
            %   non-zero

            % Calculate Δt between the input and the most recent update
            delta_t = t - state_k.in.time;

            % Δd is the "distance of UT/BS" covered since the most recent
            % update (that update occurred at or after t_k and before
            % t_{k+1})
            s = relative_speed(state_k.in);
            delta_d = s * delta_t;

            % If (d + Δd) equals or exceeds the update distance, an update
            % is required i.e. 'k' needs to be increased
            if (updateDistance~=0)
                updated = (d + delta_d) >= updateDistance;
            else
                updated = (delta_d~=0);
            end
            if (updated)

                % Get the current input state
                inputstate = makeInputState(t,tx,rx);

                % Add Δd to the accumulated distance
                d = d + delta_d;

                % Calculate Δk, the integer number of update distances
                % travelled - the input arrival rate might be low enough
                % that more than one update distance has been travelled
                if (updateDistance~=0)
                    delta_k = floor(d / updateDistance);
                else
                    delta_k = 1;
                end

                % Get the distance and time since the last update i.e
                % between t_k and t_{k+Δk}. Note that the time is
                % calculated using the relative speed from the most recent
                % update, which is assumed to be the relative speed at t_k
                % (although as noted above, that update and therefore its
                % velocities may have occurred after t_k but before
                % t_{k+1})
                if (updateDistance~=0)
                    d_update = delta_k * updateDistance;
                else
                    d_update = delta_d;
                end
                t_update = d_update / s;

                % update index 'k'
                k = k + delta_k;

                % Update small scale channel
                [ch,state_k] = updateSmallScaleChannel(cfg,ch,t,k,state_k,t_update);

                % Retain any "unprocessed distance" in the accumulated
                % distance; if 'd' is non-zero then the current time is
                % beyond t_{k+Δk}, but has not reached t_{k+Δk+1}. If 'd'
                % is zero, the current time is exactly t_{k+Δk}
                d = d - d_update;

                % Store the current input state as the state for update
                % t_k, with non-zero 'd' accounting for the case that the
                % current input time is beyond the new t_k
                state_k.in = inputstate;

            end

        end

        if (updated)
    
            % Get the path filters and path delays
            [ch.PathFilters,ch.PathDelays] = getPathFiltersAndDelays(ch.SmallScale);

        end

        varargout{1} = ch;
        varargout{2} = updated;
        if (nargout==3)
            % Calculate soft LOS state, TR 38.901 Section 7.6.3.3 Eq.
            % (7.6-18). Note that Eq. (7.6-19) is not implemented here. The
            % variable 'LOS_soft' is an input to that equation, along with
            % the channel coefficients matrices H^{LOS} and H^{NLOS}
            % produced by evaluating suitably-configured small scale
            % channels. Note that the current node positions are used to
            % calculate 'LOS_soft', not the positions at the last update.
            % This allows for calculation of Eq. (7.6-19) for node
            % positions between updates
            chCfg = ch.ChannelConfiguration;
            chCfg.LOSProbability = [];
            bspos = tx.Position;
            uepos = rx.Position;
            LOS_soft = softLOSState(chCfg,bspos,uepos);
            varargout{3} = LOS_soft;
        end

    end

    function [out_k,out_d,out_state_k] = getState()

        out_k = k;
        out_d = d;
        out_state_k = state_k;

    end

end

function s = makeInputState(t,tx,rx)

    s.time = t;

    s.txpos = tx.Position;
    s.v_bar_tx = tx.Velocity.';
    s.omega_tx = tx.RotationVelocity;

    s.rxpos = rx.Position;
    s.v_bar_rx = rx.Velocity.';
    s.omega_rx = rx.RotationVelocity;

end

function validateInfo(x,name)

    validateattributes(x.Position,{'numeric'},{'row','ncols',3},[name '.Position']);
    validateattributes(x.Velocity,{'numeric'},{'row','ncols',3},[name '.Velocity']);
    validateattributes(x.RotationVelocity,{'numeric'},{'column','nrows',3},[name '.RotationVelocity']);

end

function vn = relative_speed(in)

    vn = norm(in.v_bar_rx - in.v_bar_tx);

end

function [ch,state_k] = updateSmallScaleChannel(cfg,ch,t,k,state_k,delta_t)

    % Apply cluster-wise spatial consistency procedure
    if (isProcedureA(cfg))
        [ch.SmallScale,state_k] = procedureA(ch,k,state_k,delta_t);
    else % isProcedureB
        [ch.SmallScale,state_k] = procedureB(ch,k,state_k,delta_t);
    end

    % Update properties that aren't cluster-wise
    ch = updateNonClusterProperties(ch,t,state_k.in,delta_t);

end

% TR 38.901 Section 7.6.3.2 Procedure A
function [ss,s_k] = procedureA(ch,k,s_k,delta_t)

    persistent c;
    if (isempty(c))
        c = physconst('LightSpeed');
    end

    ss = ch.SmallScale;

    % If small scale channel is not a structure, it is an nrCDLChannel
    isCDLChannel = ~isstruct(ss);

    if (isCDLChannel)
        LOS = ss.HasLOSCluster;
    else
        LOS = true;
    end

    persistent ssp_X;
    if (isempty(ssp_X))
        % Set up cluster-specific SCRVs for random variable X_n now that
        % the number of clusters is known
        chCfg = ch.ChannelConfiguration;
        N = numel(ss.PathDelays);
        ssp_X = setupSCRVsForProcedureA(chCfg,N);
    end

    % Set up elements of state that can be re-used for each 'k'
    if (k==0)

        % Random variable X_n for Eq. (7.6-10b) and Eq. (7.6-10c)
        % "is not changed during the UT mobility within a drop"
        N = numel(ss.PathDelays);
        sc.UEPosition = s_k.in.rxpos;
        sc.NodeSubs = ch.NodeSubs;
        ul = 1;
        rvs = uniformClusterSCRVs(ssp_X,sc,N-LOS,ul);
        X = 2*(rvs>0.5)-1;
        if (LOS)
            X = [1 X];
        end
        s_k.X = X;

    end

    absoluteTOA = ~isempty(ch.Mobility.Deltatau);

    if (k==0) % "for k = 0" case in Eq. (7.6-9)

        % Calculate cluster delays
        tau_Delta = zeros(size(ss.PathDelays));
        if (absoluteTOA)
            tau_Delta((1+LOS):end) = ch.Mobility.Deltatau;
        else
            tau_prime = ch.Mobility.tau_prime; % the delays from Eq. (7.5-1)
            tau_Delta((1+LOS):end) = min(tau_prime);
        end
        tau = ch.Mobility.tau; % cluster delays as in Step 11 in Clause 7.5
        d_3D = getDistance3D(ch.ChannelConfiguration);
        tau_tilde = tau + tau_Delta + d_3D/c;

    else % "for k > 0" case in Eq. (7.6-9)

        % Calculate the change in cluster delays and angles based on the
        % state for the most recent update 's_k' (note that the cluster
        % delays depend on the angles, and vice versa)
        delta_tau_tilde = calculate_delta_tau_tilde(ss,s_k.in,delta_t);
        [delta_ZoD,delta_AoD,delta_ZoA,delta_AoA] = calculate_delta_angles(ss,s_k,delta_t);

        % Update cluster delays
        tau_tilde = s_k.tau_tilde + delta_tau_tilde;

        % Update cluster angles
        ss.AnglesZoD = ss.AnglesZoD + delta_ZoD;
        ss.AnglesAoD = ss.AnglesAoD + delta_AoD;
        ss.AnglesZoA = ss.AnglesZoA + delta_ZoA;
        ss.AnglesAoA = ss.AnglesAoA + delta_AoA;

    end

    % Store cluster delays
    s_k.tau_tilde = tau_tilde;

    % Normalize delays according to Eq (7.6-10a), replacement for
    % Eq. (7.5-2)
    tau = tau_tilde - min(tau_tilde);

    % Update cluster powers according to Section 7.5 step 6 based on
    % updated delays from Eq (7.6-10a)
    if (isCDLChannel)
        uepos = s_k.in.rxpos;
        [~,pathPowersLOS] = generatePowersForProcedureA(ch,tau,uepos);
    end

    % Update small scale channel properties related to cluster delays and
    % powers
    if (isCDLChannel)
        ss.AveragePathGains = 10*log10(pathPowersLOS);
    end
    ss.PathDelays = tau;

    % If absolute time of arrival is enabled, re-instate the propagation
    % time delay
    if (absoluteTOA)
        ss.PathDelays = ss.PathDelays + min(tau_tilde);
    end

end

% TR 38.901 Section 7.6.3.2 Procedure B
function [ss,s_k] = procedureB(ch,~,s_k,~)

    % Update channel link configuration. Note that 'n_fl' and 'd_2D_in' are
    % not updated because TR 38.901 Section 7.6.3.3 states "It is noted
    % that soft indoor/outdoor states are not modelled in this TR"
    chCfg = ch.ChannelConfiguration;
    chCfg.BSPosition = s_k.in.txpos;
    chCfg.UEPosition = s_k.in.rxpos;

    % Remove site positions from channel link configuration, it is assumed
    % they were present in the call to h38901Channel.createChannelLink that
    % originally created the channel link and therefore spatially
    % consistent random variables (SCRVs) have already been set up
    chCfg.SitePositions = [];

    % Create new channel link
    originalch = ch;
    ch = h38901Channel.createChannelLink(chCfg);
    ss = ch.SmallScale;

    % If small scale channel is not a structure, it is an nrCDLChannel
    isCDLChannel = ~isstruct(ss);

    % Restore properties of original small scale channel if needed
    if (isCDLChannel)
        ss = restoreProperty(ss,originalch,'ChannelFiltering');
    end
    ss = restoreProperty(ss,originalch,'TransmitArrayOrientation');
    ss = restoreProperty(ss,originalch,'ReceiveArrayOrientation');

end

function ss = restoreProperty(ss,originalch,name)

    ss.(name) = originalch.SmallScale.(name);

end

function ch = updateNonClusterProperties(ch,t,in,delta_t)

    persistent c;
    if (isempty(c))
        c = physconst('LightSpeed');
    end

    ss = ch.SmallScale;

    % If small scale channel is not a structure, it is an nrCDLChannel
    isCDLChannel = ~isstruct(ss);

    % TR 38.901 Section 7.6.10 dual tx and rx mobility is supported, the
    % maximum Doppler shift and direction of travel of both link ends are
    % uppdated, but scatterer mobility is not configured

    if (isCDLChannel)
        % Update maximum Doppler shift
        lambda = c / ss.CarrierFrequency;
        v_rx = norm(in.v_bar_rx);
        v_tx = norm(in.v_bar_tx);
        ss.MaximumDopplerShift = [v_rx v_tx] / lambda;
    end

    % Update direction of travel
    direction = [azimuthElevation(in.v_bar_rx) azimuthElevation(in.v_bar_tx)];
    ch.Mobility.UTDirectionOfTravel = direction;
    if (isCDLChannel)
        ss.UTDirectionOfTravel = direction;
    end

    % Update array orientations
    ss.TransmitArrayOrientation = rotateAntennaArray(ss.TransmitArrayOrientation,in.omega_tx,delta_t);
    ss.ReceiveArrayOrientation = rotateAntennaArray(ss.ReceiveArrayOrientation,in.omega_rx,delta_t);

    if (isCDLChannel)
        % Update current time
        ss.InitialTime = t;
    end

    ch.SmallScale = ss;

end

function out = softLOSState(varargin)
    % Generate soft LOS state 'LOS_soft' defined in TR 38.901 Section
    % 7.6.3.3 Eq. (7.6-18)

    persistent fn;

    if (nargin==0) % called to reset state

        fn = [];
        out = fn;

    elseif (nargin==1) % called to store SCRVs for LOS

        lsp = varargin{1};

        % Create a function that binds the spatially consistent RVs for LOS
        % to the function for generating the soft LOS state, and store it
        % as a persistent variable. This is required as these SCRVs are not
        % otherwise available during mobility modelling
        rs = []; % not used when spatial consistency is enabled
        % Do not allow negative outdoor distance
        d_2D_out = @(chCfg,bspos,uepos)max(vecnorm(uepos(1:2)-bspos(1:2)) - chCfg.d_2D_in,0);
        h_BS = @(bspos)bspos(3);
        h_UT = @(uepos)uepos(3);
        scfn = @(chCfg,bspos,uepos)struct(NodeSubs=chCfg.NodeSubs,LOS=lsp,BSPosition=bspos,UEPosition=uepos);
        fn = @(chCfg,bspos,uepos)losCondition(chCfg,rs,d_2D_out(chCfg,bspos,uepos),h_BS(bspos),h_UT(uepos),scfn(chCfg,bspos,uepos));
        out = fn;

    else % nargin==3, called to generate soft LOS state 

        validateSCRVSetup(fn);
        chCfg = varargin{1};
        bspos = varargin{2};
        uepos = varargin{3};

        % Create the soft LOS state by calling the 'losCondition' function
        [~,~,LOS_soft] = fn(chCfg,bspos,uepos);
        out = LOS_soft;

    end

end

function out = setupSCRVsForProcedureA(varargin)
    % Generate the auto-correlation matrices for cluster-specific RV 'X'
    % used in rotation matrices in spatial consistency Procedure A (i.e.
    % chCfg.SpatialConsistency="ProcedureA")

    persistent fn;

    if (nargin==0) % called to reset state

        fn = [];
        out = fn;

    elseif (nargin==1) % called to store system-wide RNG

        systemrs = varargin{1};

        % Create a function that binds the system-wide RNG to the function
        % that creates autocorrelation matrices, and store it as a
        % persistent variable. This is required as this RNG is not
        % otherwise available during mobility modelling
        fn = @(chCfg,distance,numClusters)createClusterAutoCorrMatrices(chCfg,systemrs,numClusters,distance);
        out = fn;

    else % nargin==2, called to create SCRVs

        validateSCRVSetup(fn);        
        chCfg = varargin{1};
        numClusters = varargin{2};

        % TR 38.901 Section 7.6.3.2 Procedure A, after Eq. (7.6-10c): "The
        % cluster specific decorrelation distances are 60m, 15m, 50m and
        % 10m for RMa, UMi, UMa and Indoor scenarios, respectively"
        % Note that scenarioSwitch order is UMi, UMa, RMa, InH, InF
        s = chCfg.Scenario;
        distance = scenarioSwitch(s,15,50,60,10,10);

        % Create the autocorrelation matrices by calling the
        % 'createClusterAutoCorrMatrices' function
        out = fn(chCfg,distance,numClusters);

    end

end

function varargout = generatePowersForProcedureA(varargin)
    % Generate cluster powers as discussed in step 6 of TR 38.901 Section
    % 7.5, when updating cluster powers in spatial consistency Procedure A
    % (i.e. chCfg.SpatialConsistency="ProcedureA")

    persistent fn;

    if (nargin==0) % called to reset state

        fn = [];
        varargout{1} = fn;

    elseif (nargin==1) % called to store SCRVs for powers

        ssp = varargin{1};

        % Create a function that binds the spatially consistent
        % cluster-specific RVs for powers to the function for generating
        % powers, and store it as a persistent variable. This is required
        % as these SCRVs are not otherwise available during mobility
        % modelling
        chcfg = struct(SpatialConsistency="ProcedureA");
        rs = []; % not used when spatial consistency is enabled
        scfn = @(uepos,ns)struct(Powers=ssp,UEPosition=uepos,NodeSubs=ns);
        isMobilityUpdate = true;
        fn = @(tau,m,LOS,uepos,ns)generatePowers(chcfg,rs,tau,m.r_tau,m.DS,m.zeta,m.K,m.indoor,LOS,scfn(uepos,ns),m.tau_prime,isMobilityUpdate);
        varargout{1} = fn;

    else % nargin==3, called to generate powers

        validateSCRVSetup(fn);
        ch = varargin{1};
        tau = varargin{2};
        uepos = varargin{3};

        % Generate the powers by calling the 'generatePowers' function
        m = ch.Mobility;
        LOS = ch.SmallScale.HasLOSCluster;
        ns = ch.NodeSubs;
        [varargout{1:nargout}] = fn(tau,m,LOS,uepos,ns);

    end

end

function out = setupSCRVsForProcedureB(varargin)
    % Generate the auto-correlation matrices for cluster-specific RVs used
    % in spatial consistency Procedure B
    % (i.e. chCfg.SpatialConsistency="ProcedureB")

    persistent fn;

    if (nargin==0) % called to reset state

        fn = [];
        out = fn;

    elseif (nargin==1) % called to store system-wide RNG

        systemrs = varargin{1};

        % Create a function that binds the system-wide RNG to the function
        % that creates autocorrelation matrices, and store it as a
        % persistent variable. This is required as this RNG is not
        % otherwise available during mobility modelling. For correlation
        % distances greater than 50 m, the matrices are created at a lower
        % resolution and then interpolated back to the default resolution
        % given by getGridResolution(), this is a quicker approach for
        % large correlation distances
        getDecimationFactor = @(d)min(ceil(d / 50));
        fn = @(chCfg,distance,numClusters)createClusterAutoCorrMatrices(chCfg,systemrs,numClusters,distance,getDecimationFactor(distance));
        out = fn;

    else % nargin==3, called to create SCRVs

        validateSCRVSetup(fn);
        chCfg = varargin{1};
        distance = varargin{2};
        numClusters = varargin{3};

        % Separate SCRVs for NLOS, LOS and O2I
        distances = [distance distance distance];

        % Create the autocorrelation matrices by calling the
        % 'createClusterAutoCorrMatrices' function
        out = fn(chCfg,distances,numClusters);

    end

end

function validateSCRVSetup(x)

    if (isempty(x))
        error('nr5g:h38901Channel:NoSitePositions','When spatial consistency is configured, the first call to h38901Channel/createChannelLink must include non-empty SitePositions in order to set up spatially consistent random variables (SCRVs).')
    end

end

function [pathDelays,pathDelaysScaling,tau_prime] = generateDelaysForProcedureB(mu_lgDS,sigma_lgDS,N,indoor,LOS,sc)
    % Generate cluster delays according to Step 5 of TR 38.901 Section 7.6.3.2 Procedure B

    % Delays are drawn randomly from a uniform distribution
    rvs = uniformClusterSCRVs(sc.Delays,sc,N,indoor,LOS);
    tau_prime = rvs * 2 * 10^(mu_lgDS + sigma_lgDS);

    % Normalise the delays by subtracting the minimum delay
    pathDelays = tau_prime - min(tau_prime);

    % In the case of LOS, set the delay of the first cluster to 0
    if (LOS && ~indoor)
        pathDelays(1) = 0;
    end

    % No path delay scaling is specified by Step 5 for Procedure B
    pathDelaysScaling = 1;

end

function [pathPowers,pathPowersLOS,KFactorFirstCluster] = generatePowersForProcedureB(tau_prime,ds,phi_prime_AoA,phi_prime_AoD,theta_prime_ZoA,theta_prime_ZoD,asa,asd,zsa,zsd,zeta,K,indoor,LOS,sc)

    % Per cluster shadowing term
    rvs = normalClusterSCRVs(sc.Powers,sc,size(tau_prime,2),indoor,LOS);
    shadowing = zeta*rvs;

    % Ricean K-factor in linear scale
    Kr = 10^(K/10);

    % In case of LOS condition, delay spreads and angle spreads are
    % adjusted to preserve the spread when Eq. (7.6-17c) is applied in
    % 'normalizePowers' below
    if (LOS && ~indoor)
        ds = sqrt(1 + Kr/2) * ds;
        asa = sqrt(1 + Kr) * asa;
        asd = sqrt(1 + Kr) * asd;
        zsa = sqrt(1 + Kr) * zsa;
        zsd = sqrt(1 + Kr) * zsd;
    end

    % Eq. (7.6-17)
    expangle = @(ang,spread)exp(-sqrt(2)*abs(ang) / spread);
    P_prime = exp(-tau_prime/ds);
    P_prime = P_prime .* expangle(phi_prime_AoA,asa) .* expangle(phi_prime_AoD,asd);
    P_prime = P_prime .* expangle(theta_prime_ZoA,zsa) .* expangle(theta_prime_ZoD,zsd);
    P_prime = P_prime .* 10.^(-shadowing/10);

    % Eqs. (7.6-17a), (7.6-17b) and (7.6-17c)
    [pathPowers,pathPowersLOS,KFactorFirstCluster] = normalizePowers(P_prime,K,indoor,LOS);

end

function [phi_prime_AoA,phi_prime_AoD,theta_prime_ZoA,theta_prime_ZoD,phiLOS_AoA,phiLOS_AoD,thetaLOS_ZoA,thetaLOS_ZoD] = generateAnglesForProcedureB(bspos,uepos,indoor,LOS,params,fastFading,sc)
    % Generate arrival and departure angles for each cluster for the
    % current BS-UE link. 

    % Determine LOS angles
    [phiLOS_AoA,phiLOS_AoD,thetaLOS_ZoA,thetaLOS_ZoD] = generateLOSAngles(bspos,uepos);
     
    if (~fastFading)
        % Only LOS is needed when fastFading=false. Procedure B Step 6 says
        % "In case of LOS, set the angles of the first cluster to 0", so
        % the output angles will be zero here. The final cluster angles
        % will be the LOS angles produced by 'generateLOSAngles' above once
        % 'completeAnglesForProcedureB' has been executed in Step 7b
        phi_prime_AoA = 0;
        phi_prime_AoD = 0;
        theta_prime_ZoA = 0;
        theta_prime_ZoD = 0;
    else
        N = params.NumClusters;

        % Generate AoA
        phi_prime_AoA = generateAngleForProcedureB(params.mu_lgASA,params.sigma_lgASA,N,indoor,LOS,sc,'AnglesAoA');

        % Generate AoD
        phi_prime_AoD = generateAngleForProcedureB(params.mu_lgASD,params.sigma_lgASD,N,indoor,LOS,sc,'AnglesAoD');

        % Generate ZoA
        theta_prime_ZoA = generateAngleForProcedureB(params.mu_lgZSA,params.sigma_lgZSA,N,indoor,LOS,sc,'AnglesZoA');

        % Generate ZoD
        theta_prime_ZoD = generateAngleForProcedureB(params.mu_lgZSD,params.sigma_lgZSD,N,indoor,LOS,sc,'AnglesZoD');
    end

end

function a_prime = generateAngleForProcedureB(mu_lg,sigma_lg,N,indoor,LOS,sc,sspName)

    % Eq. (7.6-16)
    rvs = uniformClusterSCRVs(sc.(sspName),sc,N,indoor,LOS); % unif(0,1)
    rvs = (rvs - 0.5) * 2; % unif(-1,1)
    a_prime = rvs * (2 * 10^(mu_lg + sigma_lg));

    % In case of LOS, set the angles of the first cluster to 0
    if (LOS && ~indoor)
        a_prime(1) = 0;
    end

end

function [anglesAoA, anglesAoD, anglesZoA, anglesZoD] = completeAnglesForProcedureB(phi_prime_AoA,phi_prime_AoD,theta_prime_ZoA,theta_prime_ZoD,phiLOS_AoA,phiLOS_AoD,thetaLOS_ZoA,thetaLOS_ZoD,mu_offsetZOD,indoor)

    % Note: Eqs. (7.6-17e) and (7.6-17f) in TR 38.901 refer to phi'_ZOA and
    % phi'_ZOD respectively. However, zenith angles are typically denoted
    % by the variable theta. Step 6 states that Eq. (7.6-16), which defines
    % phi'_AoA, "is repeated independently for AOD, AOA, ZOD, and ZOA with
    % corresponding maximum angles". The most literal interpretation would
    % result in all four variables being differently subscripted phi', but
    % the implementation here uses theta for the zenith angles

    % Eq. (7.6-17d)
    anglesAoA = phi_prime_AoA + phiLOS_AoA;
    anglesAoD = phi_prime_AoD + phiLOS_AoD;

    % Eq. (7.6-17e)    
    if indoor
        theta_bar = 90;
    else
        theta_bar = thetaLOS_ZoA;
    end
    anglesZoA = theta_prime_ZoA + theta_bar;

    % Eq. (7.6-17f)
    anglesZoD = theta_prime_ZoD + thetaLOS_ZoD + mu_offsetZOD;

end

function ao = rotateAntennaArray(ao,omega,delta_t)

    % 'delta_t' is in seconds
    delta_t_minutes = delta_t / 60;

    % 'omega' is in revolutions per minute (RPM)
    revs = omega * delta_t_minutes;
    degrees = revs * 360;
    ao = ao + degrees;

end

function azel = azimuthElevation(v_bar)

    v = norm(v_bar);
    if (v~=0)
        v_bar = v_bar / v;
    end
    phi = azimuthAngle(v_bar);
    theta = zenithAngle(v_bar);

    azel = [phi; theta];

end


function delta_tau_tilde = calculate_delta_tau_tilde(ss,in,delta_t)

    persistent c;
    if (isempty(c))
        c = physconst('LightSpeed');
    end

    % Eq. (7.6-10)
    r_hat_rx = permute(wireless.internal.channelmodels.getRhoHat(ss.AnglesAoA,ss.AnglesZoA),[1 3 2]);

    % Eq. (7.6-10aa)
    r_hat_tx = permute(wireless.internal.channelmodels.getRhoHat(ss.AnglesAoD,ss.AnglesZoD),[1 3 2]);

    % Change in delay in Eq. (7.6-9)
    delta_tau_tilde = -((r_hat_rx.'*in.v_bar_rx + r_hat_tx.'*in.v_bar_tx) / c) * delta_t;
    delta_tau_tilde = delta_tau_tilde.';

end

function [delta_ZoD,delta_AoD,delta_ZoA,delta_AoA] = calculate_delta_angles(ss,s_k,delta_t)

    persistent c;
    if (isempty(c))
        c = physconst('LightSpeed');
    end

    % If small scale channel is not a structure, it is an nrCDLChannel
    isCDLChannel = ~isstruct(ss);

    % Setup
    in = s_k.in;
    N = numel(s_k.X);
    R = repmat(eye(3),[1 1 N]);
    R(2,2,:) = s_k.X;
    phi_AoD = ss.AnglesAoD;
    phi_AoA = ss.AnglesAoA;
    theta_ZoD = ss.AnglesZoD;
    theta_ZoA = ss.AnglesZoA;
    if (isCDLChannel)
        LOS = ss.HasLOSCluster;
    else
        LOS = true;
    end
    if (LOS)
        % Set up values that ensure R_Z and R_Y return identity matrices
        % for the LOS path, in order to follow the "for LOS" case of
        % Eq. (7.6-10b) and Eq. (7.6-10c)
        phi_AoD(1) = -180;
        phi_AoA(1) = 0;
        theta_ZoD(1) = 90;
        theta_ZoA(1) = 90;
    end
    fR_Z = @(alpha)wireless.internal.channelmodels.rotationMatrix([alpha; 0; 0]);
    fR_Y = @(beta)wireless.internal.channelmodels.rotationMatrix([0; beta; 0]);
    R_Z = @(alpha)d3fn(fR_Z,alpha);
    R_Y = @(beta)d3fn(fR_Y,beta);
    R = pagemtimes(R_Y(90 - theta_ZoD),R);
    R = pagemtimes(R,R_Y(90 - theta_ZoA));
    dot = @(v,fn,theta,phi)permute(pagemtimes(v,'transpose',permute(fn(theta,phi),[1 3 2]),'none'),[1 3 2]);

    % Eq. (7.6-10b)
    R_rx = pagemtimes(R_Z(phi_AoD + 180),R);
    R_rx = pagemtimes(R_rx,R_Z(-phi_AoA));
    v_bar_prime_rx = pagemtimes(R_rx,in.v_bar_rx) - in.v_bar_tx;

    % Eq. (7.6-10c)
    R_tx = pagemtimes(R_Z(-phi_AoD),R);
    R_tx = pagemtimes(R_tx,R_Z(phi_AoA + 180));
    v_bar_prime_tx = pagemtimes(R_tx,in.v_bar_tx) - in.v_bar_rx;

    % Change in angle in Eq. (7.6-11)
    delta_AoD = rad2deg(dot(v_bar_prime_rx,@phi_hat,ss.AnglesZoD,ss.AnglesAoD) ./ (c*s_k.tau_tilde.*sind(ss.AnglesZoD)) * delta_t);
    delta_AoD(~isfinite(delta_AoD)) = 0; % Set ambiguous azimuth to zero

    % Change in angle in Eq. (7.6-12)
    delta_ZoD = rad2deg(dot(v_bar_prime_rx,@theta_hat,ss.AnglesZoD,ss.AnglesAoD) ./ (c*s_k.tau_tilde) * delta_t);

    % Change in angle in Eq. (7.6-13)
    delta_AoA = rad2deg(dot(v_bar_prime_tx,@phi_hat,ss.AnglesZoA,ss.AnglesAoA) ./ (c*s_k.tau_tilde.*sind(ss.AnglesZoA)) * delta_t);
    delta_AoA(~isfinite(delta_AoA)) = 0; % Set ambiguous azimuth to zero

    % Change in angle in Eq. (7.6-14)
    delta_ZoA = rad2deg(dot(v_bar_prime_tx,@theta_hat,ss.AnglesZoA,ss.AnglesAoA) ./ (c*s_k.tau_tilde) * delta_t);

end

function out = d3fn(fn,x)

    out = arrayfun(fn,x,UniformOutput=false);
    out = cat(3,out{:});

end

% Eq. (7.1-13)
function out = theta_hat(theta,phi)

    out = [cosd(theta).*cosd(phi); cosd(theta).*sind(phi); -sind(theta)];

end

% Eq. (7.1-14)
function out = phi_hat(~,phi)

    out = [-sind(phi); cosd(phi); zeros(size(phi))];

end
