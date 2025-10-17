classdef hPlotPacketTransitions < handle
    %hPlotPacketTransitions Plots the live packet communication across
    %time and frequency for Bluetooth, WLAN, and coexistence nodes
    %
    %   hPlotPacketTransitions(NODES,SIMULATIONTIME,Name=Value) creates
    %   the live packet communication visualization across time and frequency
    %   for Bluetooth, WLAN, interfering WLAN and coexistence nodes. The first
    %   subplot is time vs nodes. The second subplot is time vs frequency. You
    %   can specify additional name-value arguments in any order as
    %   (Name1=Value1, ..., NameN=ValueN).
    %
    %   NODES is the configured wireless nodes to visualize. Specify this as an
    %   array of objects or cell array of scalar objects of type bluetoothNode,
    %   bluetoothLENode, wlanNode, helperCoexNode, or
    %   helperInterferingWLANNode.
    %
    %   SIMULATIONTIME is a finite positive scalar indicating the simulation
    %   time in seconds.
    %
    %   hPlotPacketTransitions properties (Configurable as NV pair in
    %   constructor):
    %
    %   FrequencyPlotFlag - Display the time vs frequency plot VisualizeAtEnd
    %   - Display visualization only at the end of simulation
    %
    %   Note:
    %     * The figure shows packet communication across time and frequency for
    %       all the configured nodes.
    %     * The top plot shows state transition across time for all the
    %       nodes.
    %     * The bottom plot shows the packet communication in frequency
    %       across time. This plot can be disabled using the FrequencyPlotFlag
    %       NV pair.
    %     * In frequency domain if two or more packets operate in same
    %       frequency, the packet might be with higher contrast indicating
    %       overlap.
    %     * The legend indicating the states and node type can be seen on the
    %       right for both the plots.
    %     * The time slider and edit box between both the plot can be used to
    %       go to a specific time period.
    %     * The full view check box can be used to view all the nodes and the
    %       complete simulation duration.
    %     * The export button can be used to save a snip of the visualization.
    %     * Interactions in the UI figure are applicable only at the end of
    %       simulation.

    %   Copyright 2023-2024 The MathWorks, Inc.

    %% Configurable properties as NV pair in constructor
    properties (SetAccess=private)
        %FrequencyPlotFlag Enable the time vs frequency plot
        %   FrequencyPlotFlag is a flag indicating whether to enable or disable the
        %   frequency plot. Specify this flag as a logical scalar. The default
        %   value is true.
        FrequencyPlotFlag (1,1) {logical} = true

        %VisualizeAtEnd Display visualization only at the end of simulation
        %   VisualizeAtEnd is a flag indicating whether to display the
        %   visualization only at the end of simulation. If this flag is disabled,
        %   the packet communication can be visualized during the run time of the
        %   simulation. If this flag is enabled, the visualization will be
        %   displayed only after the simulation is over. Irrespective of the state
        %   of this flag, packet communication will be available for the complete
        %   simulation time. For long running simulation enable this flag for
        %   better performance. Specify this flag as a logical scalar. The default
        %   value is false.
        VisualizeAtEnd (1,1) {logical} = false
    end

    properties (SetAccess=private,Hidden) % Figure Based properties
        %PacketCommUIFigure UI Figure handle for the packet communication
        %visualization
        PacketCommUIFigure

        %StateAxes UI Axes handle for time vs nodes displaying the state
        %transitions
        StateAxes

        %FreqAxes UI Axes handle for time vs frequency displaying the packet
        %communication over frequency
        FreqAxes

        %pMinTimeSliderMin Handle for horizontal slider to vary the minimum time
        %value of the plots
        pMinTimeSlider

        %pMaxTimeSlider Handle for horizontal slider to vary the maximum time value
        %of the plots
        pMaxTimeSlider

        %pEditMinHandle Handle for the edit box to vary the minimum time value
        %of the plots
        pEditMinHandle

        %pEditMaxHandle Handle for the edit box to vary the maximum time value
        %of the plots
        pEditMaxHandle

        %pTimeLine1 Handle of the line to display in the state transition figure
        pTimeLine1

        %pTimeLine2 Handle of the line to display in the time vs frequency figure
        pTimeLine2
    end

    properties (Access=private) % Properties related to the figure display
        %pMaxNodesToDisplay Maximum number of nodes to display in Y axis of state
        %transitions plot. If the frequency plot is enabled the maximum number of
        %nodes to display is set as 6. If the frequency plot is disabled maximum
        %number of nodes to display is set as 20.
        pMaxNodesToDisplay

        %pMinMaxFreqLimits Minimum and maximum Y limits in the frequency domain
        pMinMaxFreqLimits

        %pWLANFreqLimits Minimum and maximum frequency limits for WLAN nodes
        pWLANFreqLimits = []

        %pYTickOrder List of Y tick information of the state transition plot. First
        %column indicates ID of the node, second column indicates device ID incase
        %of WLAN, and the third column indicates the type of technology. 1
        %indicates WLAN, 3 indicates Bluetooth LE, and 4 indicates Bluetooth
        %BR/EDR. This is used to plot the packets in their respective positions in
        %the state transition plot.
        pYTickOrder

        %pXLimits X axis limit values stored for usage when reverting from full
        %view
        pXLimits

        %pYLimits Y axis limit values of state transition plot stored for usage
        %when reverting from full view
        pYLimits

        %pFreqLimits Y axis limit values of frequency plot stored for usage when
        %reverting from full view
        pFreqLimits

        %pIsInitialized Flag to check whether the visualization is initialized
        pIsInitialized = false

        %pObjIdxsOnDisp Object indices on display
        pObjIdxsOnDisp

        %pLastEndTimePlotted End time of the last packet plotted during live
        %visualization
        pLastEndTimePlotted = 0

        %pLastCalledTime Time for which the last update was called to draw packets
        pLastCalledTime = [0 10]

        %pLastCalledNodes Nodes for which the last update was called to draw
        %packets
        pLastCalledNodes = [1 2]
    end

    properties (Access=private) % Network simulation related information
        %CommunicationInfo Complete information from all the events
        CommunicationInfo

        %InfoCount Number of information received from events or callbacks
        InfoCount = 0

        %StartTimeList Start time of all the information stored in
        %CommunicationInfo
        StartTimeList

        %EndTimeList End time of all the information stored in CommunicationInfo
        EndTimeList

        %YTicksList Y tick of all the information stored in CommunicationInfo
        YTicksList

        %pMaxEndTime Maximum end time among the list of information received
        pMaxEndTime = 0

        %pStatePatch List of all the patch handles in the state axes. The order is
        %same as the pStateColors.
        pStatePatch

        %pFreqPatch List of all the patch handles in the frequency axes. The order
        %is same as the pFreqColors
        pFreqPatch

        %pSimulationTime Run time of the network simulation
        pSimulationTime = 0

        %pIsSimInProgress Flag indicating whether the network simulation is running
        pIsSimInProgress = true

        %CommunicationInfoStruct Default structure for the information to be stored
        CommunicationInfoStruct = struct("Type",0, ...
            "StartTime",-1, ...
            "Duration",-1, ...
            "FillColor","", ...
            "YTickVal",-1, ...
            "FreqPlot",false, ...
            "FreqStart",-1, ...
            "ChannelBandwidth",-1, ...
            "FreqColor","")
    end

    properties (Access=private)
        % List of state colors for the pre existing states
        pTxColor = "--mw-color-success"
        pIdleEIFSColor = "--mw-backgroundColor-tertiary"
        pContendColor = "--mw-color-warning"
        pRxColorUs = "--mw-graphics-colorOrder-9-quaternary"
        pRxColorOthers = "--mw-graphics-colorOrder-2-tertiary"
        pRxColorSuccess = "--mw-graphics-colorOrder-11-secondary"
        pRxColorFailure = "--mw-color-error"

        %pLineColor Color of the vertical line shown in both state and frequency
        %axes to map the time
        pLineColor = "--mw-graphics-colorOrder-5-secondary"

        %pStatePacketHeight Height of the packet plotted in state axes
        pBarHeight = 1

        %pStateColors List of colors for all the states displayed
        pStateColors

        %pFreqColors List of colors for all the wireless technologies
        pFreqColors

        %pNumPacketsToShow Number of packets that can be visible in the screen
        pNumPacketsToShow = 10^4
    end

    %% Standard node relevant properties
    properties (Access=private)
        %pCommunicationInfoWLAN Packet type and packet color of WLAN packets
        pCommunicationInfoWLAN = {1,"--mw-graphics-colorOrder-8-primary"}

        %pCommunicationInfoLE Packet type and packet color of Bluetooth LE packets
        pCommunicationInfoLE = {3,"--mw-graphics-colorOrder-4-quaternary"}

        %pCommunicationInfoBREDR Packet type and packet color of Bluetooth BR/EDR
        %packets
        pCommunicationInfoBREDR = {4,"--mw-graphics-colorOrder-11-quaternary"}

        %pCommunicationInfoLE6GHz Packet type and packet color of Bluetooth LE 6 GHz packets
        pCommunicationInfoLE6GHz = {5,"--mw-graphics-colorOrder-2-quaternary"}

        %pLECenterFrequencies Bluetooth LE channel center frequencies in Hz
        pLECenterFrequencies = [2404:2:2424 2428:2:2478 2402 2426 2480]*1e6
    end

    properties (SetAccess=private)
        %WLANNodes List of WLAN nodes to visualize. The value is an array of
        %objects of type wlanNode.
        WLANNodes

        %BluetoothNodes List of Bluetooth nodes to visualize. The value is an array
        %of objects of type bluetoothNode or bluetoothLENode object.
        BluetoothNodes

        %CoexNodes List of coexistence nodes to visualize. The value is an array of
        %objects of type helperCoexNode object.
        CoexNodes

        %InterferingNodes List of interfering WLAN or LE nodes to
        %visualize. The value is an array of objects of type
        %helperInterferingWLANNode or helperInterferingBluetoothNode object.
        InterferingNodes
    end

    properties (Access=private)
        %pPPDUDuration PPDU duration of the WLAN packet received
        pPPDUDuration = zeros(1,3)

        %pNumBREDRNodes Number of Bluetooth BR/EDR nodes configured
        pNumBREDRNodes = 0

        %pNumLENodes Number of Bluetooth LE nodes configured
        pNumLENodes = 0

        %pNumWLANNodes Number of WLAN nodes configured
        pNumWLANNodes = 0

        %pNumInterferingNodes Number of interfering nodes configured
        pNumInterferingNodes = 0

        %pNumLE6GHzNodes Number of Bluetooth LE 6 GHz nodes configured
        pNumLE6GHzNodes
    end

    %% Constructor
    methods
        function obj = hPlotPacketTransitions(nodes,simulationTime,varargin)

            % Validate and store the simulation time
            validateattributes(simulationTime,{'numeric'},{'nonempty','scalar','positive','finite'},mfilename,"simulationDuration");
            obj.pSimulationTime = max(round(simulationTime,9),1e-9);

            % Name-value pair check
            if mod(nargin-2,2)==1
                error("Invalid property/value pair arguments");
            end

            % Store the NV pair information
            for idx = 1:2:nargin-2
                obj.(char(varargin{idx})) = varargin{idx+1};
            end

            % Validate the nodes
            if iscell(nodes)
                for idx = 1:numel(nodes)
                    validateattributes(nodes{idx},["bluetoothLENode","bluetoothNode","wlanNode","helperCoexNode","helperInterferingWLANNode","helperInterferingBluetoothNode"], ...
                        {'scalar'},mfilename,"nodes");
                end
            else
                validateattributes(nodes(1),["bluetoothLENode","bluetoothNode","wlanNode","helperCoexNode","helperInterferingWLANNode","helperInterferingBluetoothNode"], ...
                    {'scalar'},mfilename,"nodes");
                nodes = num2cell(nodes);
            end

            % Use weak-references for cross-linking handle objects
            objWeakRef = matlab.lang.WeakReference(obj);

            % Segregate the nodes based on their type and attach necessary listeners
            % for plotting the visualization
            [countCoex,countBluetooth,countWLAN,countInt] = deal(1);
            [numBREDRNodes,numLENodes,numWLANNodes,numIntWLANNodes,numHelperLENodes] = deal(0);
            for idx = 1:numel(nodes)
                if isa(nodes{idx},"helperCoexNode")
                    obj.CoexNodes{countCoex} = nodes{idx};
                    % Add listener to Bluetooth device
                    addlistener(nodes{idx},"PacketTransmissionStarted",@(srcNode,eventData) objWeakRef.Handle.bluetoothPlotCallback(srcNode,eventData));
                    addlistener(nodes{idx},"PacketReceptionEnded",@(srcNode,eventData) objWeakRef.Handle.bluetoothPlotCallback(srcNode,eventData));

                    % Add listener to WLAN device
                    addlistener(nodes{idx},"StateChanged",@(srcNode,eventData) objWeakRef.Handle.wlanPlotCallback(srcNode,eventData,obj.CoexNodes));
                    addlistener(nodes{idx},"MPDUDecoded",@(srcNode,eventData) objWeakRef.Handle.wlanPlotCallback(srcNode,eventData,obj.CoexNodes));

                    % Increment device count based on the added devices in coexistence node
                    if ~isempty(nodes{idx}.BluetoothDevice.DeviceName)
                        numBREDRNodes = numBREDRNodes+1;
                    end
                    if ~isempty(nodes{idx}.BluetoothLEDevice.DeviceName)
                        numLENodes = numLENodes+1;
                    end
                    if ~isempty(nodes{idx}.WLANDevice.DeviceName)
                        numWLANNodes = numWLANNodes+1;
                    end
                    countCoex = countCoex+1;
                elseif isa(nodes{idx},"bluetoothNode") || isa(nodes{idx},"bluetoothLENode")
                    obj.BluetoothNodes{countBluetooth} = nodes{idx};
                    % Add listener to Bluetooth device
                    addlistener(obj.BluetoothNodes{countBluetooth},"PacketTransmissionStarted",@(srcNode,eventData) objWeakRef.Handle.bluetoothPlotCallback(srcNode,eventData));
                    addlistener(obj.BluetoothNodes{countBluetooth},"PacketReceptionEnded",@(srcNode,eventData) objWeakRef.Handle.bluetoothPlotCallback(srcNode,eventData));
                    if isa(nodes{idx},"helperBluetoothLE6GHzNode")
                        addlistener(obj.BluetoothNodes{countBluetooth},"ChannelAccessEnded",@(srcNode,eventData) objWeakRef.Handle.bluetoothPlotCallback(srcNode,eventData));
                    end

                    % Increment Bluetooth node count
                    if isa(nodes{idx},"helperBluetoothLE6GHzNode")
                        numHelperLENodes = numHelperLENodes+1;
                    elseif isa(nodes{idx},"bluetoothNode")
                        numBREDRNodes = numBREDRNodes+1;
                    elseif isa(nodes{idx},"bluetoothLENode")
                        numLENodes = numLENodes+1;
                    end
                    countBluetooth = countBluetooth+1;
                elseif isa(nodes{idx},"wlanNode")
                    obj.WLANNodes{countWLAN} = nodes{idx};
                    % Add listener to WLAN device
                    addlistener(obj.WLANNodes{countWLAN},"StateChanged",@(srcNode,eventData) objWeakRef.Handle.wlanPlotCallback(srcNode,eventData,obj.WLANNodes));
                    addlistener(obj.WLANNodes{countWLAN},"MPDUDecoded",@(srcNode,eventData) objWeakRef.Handle.wlanPlotCallback(srcNode,eventData,obj.WLANNodes));

                    % Increment WLAN node count
                    numWLANNodes = numWLANNodes+1;
                    countWLAN = countWLAN+1;
                elseif isa(nodes{idx},"helperInterferingWLANNode") || isa(nodes{idx},"helperInterferingBluetoothNode")
                    obj.InterferingNodes{countInt} = nodes{idx};
                    addlistener(obj.InterferingNodes{countInt},"PacketTransmissionStarted",@(srcNode,eventData) objWeakRef.Handle.interferingPlotCallback(srcNode,eventData));
                    numIntWLANNodes = numIntWLANNodes+1;
                    countInt = countInt+1;
                end
            end
            obj.pNumBREDRNodes = numBREDRNodes;
            obj.pNumLENodes = numLENodes;
            obj.pNumWLANNodes = numWLANNodes;
            obj.pNumInterferingNodes = numIntWLANNodes;
            obj.pNumLE6GHzNodes = numHelperLENodes;

            % Schedule actions in the network simulator
            networkSimulator = wirelessNetworkSimulator.getInstance;
            if ~obj.VisualizeAtEnd
                scheduleAction(networkSimulator,@(varargin) objWeakRef.Handle.drawCommunicationInfo([],[],[],[]),[],0,0.005);
            end
        end
    end

    %% Figure Based Methods
    methods (Access=private)
        function initializeVisualization(obj)
            %initializeVisualization Initialize the visualization

            import matlab.graphics.internal.themes.specifyThemePropertyMappings

            % Use weak-references for cross-linking handle objects
            objWeakRef = matlab.lang.WeakReference(obj);

            % Schedule action in the network simulator at the simulation end time to
            % update the visualization
            networkSimulator = wirelessNetworkSimulator.getInstance;
            simulationTime = obj.pSimulationTime;
            scheduleAction(networkSimulator,@(varargin) objWeakRef.Handle.updateFigureAtSimulationEnd(),[],simulationTime);

            % Get screen resolution and calculate the figure dimension
            set(0,"units","pixels");
            resolution = get(0,"screensize");
            screenWidth = resolution(3);
            screenHeight = resolution(4);
            figureWidth = screenWidth*0.85;
            figureHeight = screenHeight*0.8;

            % Create the UI figure
            if obj.FrequencyPlotFlag
                figUIName = "Packet Communication Over Time And Frequency";
            else
                figUIName = "Packet Communication Over Time";
            end
            obj.PacketCommUIFigure = uifigure("Name",figUIName,"Tag","Packet Communication UI");
            if obj.VisualizeAtEnd
                obj.PacketCommUIFigure.Visible = "off";
            end
            obj.PacketCommUIFigure.Position = [60 60 figureWidth figureHeight];
            obj.PacketCommUIFigure.WindowButtonDownFcn =  @(~,eventData) objWeakRef.Handle.updateLine(eventData,zeros(1,0));
            figureGrid = uigridlayout(obj.PacketCommUIFigure,[11 4], ... % 11 rows and 4 columns
                "ColumnWidth",{'fit','fit','fit','fit','5x','5x','fit','fit','fit'}, ...
                "RowHeight",{'fit','0.2x','1.8x','1x','1x','fit','1x','0.2x','1.8x','1x','1x'}, ...
                "Tag","Figure Grid");

            % Initialize the colors based on the theme
            currTheme = obj.PacketCommUIFigure.Theme;
            if isempty(currTheme)
                currTheme = matlab.graphics.internal.themes.lightTheme;
            end

            getThemeColor = @(currTheme,color) matlab.graphics.internal.themes.getAttributeValue(currTheme, color);

            obj.pContendColor = getThemeColor(currTheme, obj.pContendColor);
            obj.pTxColor = getThemeColor(currTheme, obj.pTxColor);
            obj.pIdleEIFSColor = getThemeColor(currTheme, obj.pIdleEIFSColor);
            obj.pRxColorUs = getThemeColor(currTheme, obj.pRxColorUs);
            obj.pRxColorOthers = getThemeColor(currTheme, obj.pRxColorOthers);
            obj.pRxColorSuccess = getThemeColor(currTheme, obj.pRxColorSuccess);
            obj.pRxColorFailure = getThemeColor(currTheme, obj.pRxColorFailure);
            obj.pLineColor = getThemeColor(currTheme, obj.pLineColor);
            obj.pCommunicationInfoWLAN{2} = getThemeColor(currTheme, obj.pCommunicationInfoWLAN{2});
            obj.pCommunicationInfoLE{2} = getThemeColor(currTheme, obj.pCommunicationInfoLE{2});
            obj.pCommunicationInfoBREDR{2} = getThemeColor(currTheme, obj.pCommunicationInfoBREDR{2});
            obj.pCommunicationInfoLE6GHz{2} = getThemeColor(currTheme, obj.pCommunicationInfoLE6GHz{2});

            % Add title
            title = createLabel(obj,figureGrid,figUIName,"Packet Communication Figure Title",18,"Bold",1,[2 8]);
            title.HorizontalAlignment = "center";

            % If the frequency plot is disabled, display the state transition for the
            % complete screen and place the slider below
            if ~obj.FrequencyPlotFlag
                figureGrid.RowHeight{10} = "1.5x";
                stateFigRow = [2 9];
                sliderRow = 11;
                obj.pMaxNodesToDisplay = 20;
            else
                stateFigRow = [2 5];
                sliderRow = 7;
                obj.pMaxNodesToDisplay = 6;
            end

            % Axes for time vs nodes: State transition plot
            panelState = createPanel(obj,figureGrid,"State Panel",stateFigRow,[2 8]);
            panelState.Scrollable = "off";
            panelState.BorderColor = panelState.BackgroundColor;
            panelState.AutoResizeChildren = "off";
            obj.StateAxes = uiaxes("Parent",panelState,"Tag","State Transition Axes");
            obj.StateAxes.XLabel.String = "Time (seconds)";
            obj.StateAxes.YLabel.String = "Node Names";
            obj.StateAxes.FontSize = 14;
            obj.StateAxes.Title.String = "State Transitions of Nodes";
            obj.StateAxes.Title.FontSize = 16;
            obj.StateAxes.Position(4) = panelState.InnerPosition(4)-panelState.BorderWidth-10;
            obj.StateAxes.Position(3) =  panelState.InnerPosition(3)-10;
            xlim(obj.StateAxes,[0 simulationTime]);
            panelState.SizeChangedFcn = @(~,~) objWeakRef.Handle.resizeAxes(panelState,obj.StateAxes,true);
            hold(obj.StateAxes,"on")

            % Axes for time vs frequency: Channel occupancy plot
            if obj.FrequencyPlotFlag
                panelFreq = createPanel(obj,figureGrid,"Freq Panel",[8 11],[2 8]);
                panelFreq.Scrollable = "off";
                panelFreq.BorderColor = panelFreq.BackgroundColor;
                panelFreq.AutoResizeChildren = "off";
                obj.FreqAxes = uiaxes("Parent",panelFreq,"Tag","Frequency Axes");
                obj.FreqAxes.XLabel.String = "Time (seconds)";
                obj.FreqAxes.YLabel.String = "Frequency (MHz)";
                obj.FreqAxes.FontSize = 14;
                obj.FreqAxes.Title.String = "Packet Communication over Frequency";
                obj.FreqAxes.Position(4) = panelFreq.InnerPosition(4)-panelFreq.BorderWidth-10;
                obj.FreqAxes.Position(3) =  obj.StateAxes.Position(3);
                panelFreq.SizeChangedFcn = @(~,~) objWeakRef.Handle.resizeAxes(panelFreq,obj.FreqAxes,false);
                obj.FreqAxes.Title.FontSize = 16;
                hold(obj.FreqAxes,"on")
            end

            % Initialize the limit mode for axes
            props = ["X" "Y" "Z" "A" "C"]+"LimMode";
            for prop=props
                obj.StateAxes.(prop)="manual";
                if obj.FrequencyPlotFlag
                    obj.FreqAxes.(prop)="manual";
                end
            end

            % Initialize the axes interactivity options
            props = ["RotateSupported" "DatatipsSupported" "BrushSupported"];
            for prop=props
                obj.StateAxes.InteractionOptions.(prop)="off";
                if obj.FrequencyPlotFlag
                    obj.FreqAxes.InteractionOptions.(prop)="off";
                end
            end

            % Calculate the tick list and corresponding labels
            tickBase = 1;
            tickIdx = 1;
            yTicksList = zeros(0,1);
            yTickLabels = cell(0,1);

            % Update the Y Ticks of coexistence nodes with node and device for Y axis
            % labels
            wlanFreqInfoCoex = [];
            for idx = 1:numel(obj.CoexNodes)
                if ~isempty(obj.CoexNodes{idx}.WLANDevice.DeviceName)
                    [tickIdx,tickBase,yTicksList,yTickLabels,wlanFreqInfoCoex] = ...
                        yTicksWLAN(obj,tickIdx,tickBase,yTicksList,yTickLabels,obj.CoexNodes(idx));
                    tickBase = tickBase - obj.pBarHeight*0.75;
                end

                for btIdx = [obj.pCommunicationInfoLE{1} obj.pCommunicationInfoBREDR{1}]
                    if (btIdx==4 && ~isempty(obj.CoexNodes{idx}.BluetoothDevice.DeviceName))
                        [tickIdx,tickBase,yTicksList,yTickLabels] = ...
                            yTicksBluetooth(obj,btIdx,["BR" "EDR"],tickIdx,tickBase, ...
                            yTicksList,yTickLabels,obj.CoexNodes{idx},obj.CoexNodes{idx}.BluetoothDevice.DeviceID);
                        tickBase = tickBase - obj.pBarHeight*0.75;
                    elseif(btIdx==3 && ~isempty(obj.CoexNodes{idx}.BluetoothLEDevice.DeviceName))
                        [tickIdx,tickBase,yTicksList,yTickLabels] = ...
                            yTicksBluetooth(obj,btIdx,"LE",tickIdx,tickBase,yTicksList, ...
                            yTickLabels,obj.CoexNodes{idx},obj.CoexNodes{idx}.BluetoothLEDevice.DeviceID);
                        tickBase = tickBase - obj.pBarHeight*0.75;
                    end
                end
                tickBase = tickBase + obj.pBarHeight*0.75;
                tickBase = tickBase + obj.pBarHeight*1.25;
            end
            tickBase = tickBase - obj.pBarHeight*1.25;
            tickIdx = numel(yTickLabels)+1;

            % Update the Y Ticks of WLAN nodes
            [tickIdx,tickBase,yTicksList,yTickLabels,wlanFreqInfoWLAN] = ...
                yTicksWLAN(obj,tickIdx,tickBase,yTicksList,yTickLabels,obj.WLANNodes);

            % Update the Y Ticks of Bluetooth nodes
            for idx = 1:numel(obj.BluetoothNodes)
                [tickIdx,tickBase,yTicksList,yTickLabels] = ...
                    yTicksBluetooth(obj,0,"",tickIdx,tickBase,yTicksList, ...
                    yTickLabels,obj.BluetoothNodes{idx},1);
            end

            % Update the Y Ticks of interfering nodes
            wlanFreqInfoIntWLAN = [];
            for idx = 1:numel(obj.InterferingNodes)
                if isa(obj.InterferingNodes{idx},"helperInterferingWLANNode")
                    wlanFreqInfoIntWLAN = zeros(numel(obj.InterferingNodes),3);
                    [tickIdx,tickBase,yTicksList,yTickLabels,wlanFreqInfoIntWLAN(idx,:)] = ...
                        yTicksIntWLAN(obj,0,tickIdx,tickBase,yTicksList, ...
                        yTickLabels,obj.InterferingNodes{idx});
                else
                    [tickIdx,tickBase,yTicksList,yTickLabels] = ...
                        yTicksBluetooth(obj,0,"",tickIdx,tickBase,yTicksList, ...
                        yTickLabels,obj.InterferingNodes{idx},1);
                end
            end

            % Update the Y Ticks of WLAN nodes
            if ~isempty(wlanFreqInfoCoex)
                wlanFreqInfo = [wlanFreqInfoCoex;wlanFreqInfoIntWLAN;wlanFreqInfoWLAN];
            else
                wlanFreqInfo = [wlanFreqInfoIntWLAN;wlanFreqInfoWLAN];
            end

            % If the number of nodes in Y axis is less than maximum nodes to display
            % update Y Axis
            if numel(yTickLabels)<obj.pMaxNodesToDisplay
                obj.pMaxNodesToDisplay = numel(yTickLabels);
            end
            ymax = max(yTicksList(1:obj.pMaxNodesToDisplay))+obj.pBarHeight;
            ylim(obj.StateAxes,[0 ymax]);

            % Update the Y ticks
            yticks(obj.StateAxes,yTicksList);
            obj.StateAxes.YTickLabel = yTickLabels;
            set(obj.StateAxes,"TickLength",[0 0])

            % Legend for state transition plot
            stateTransitionHeading = "Operating States";
            stateLegendPanel = createPanel(obj,figureGrid,"State Transition Legend Panel",[1 4],9);
            stateLegendPanel.Clipping = "off";
            stateLegendGrid = uigridlayout(stateLegendPanel,[100 4]);
            for idx = 1:numel(stateLegendGrid.ColumnWidth)
                stateLegendGrid.ColumnWidth{idx} = 55;
            end

            % Find the operating states and node types based on the nodes added
            statesList = string.empty;
            typeList = string.empty;
            if obj.pNumWLANNodes+obj.pNumInterferingNodes>0 && obj.pNumBREDRNodes+obj.pNumLENodes+obj.pNumLE6GHzNodes==0 % Only WLAN nodes
                statesList = ["Transmission"; "Idle/EIFS/SIFS"; "Contention"; ...
                    "Reception (Destined to node)";"Reception"; "Reception Failure"];
                obj.pStateColors = [obj.pTxColor;obj.pIdleEIFSColor;obj.pContendColor; ...
                    obj.pRxColorUs;obj.pRxColorOthers;obj.pRxColorFailure];
                typeList = "WLAN";
                obj.pFreqColors = obj.pCommunicationInfoWLAN{2};
            elseif obj.pNumWLANNodes+obj.pNumInterferingNodes==0 && obj.pNumBREDRNodes+obj.pNumLENodes+obj.pNumLE6GHzNodes>0 % Only Bluetooth nodes
                statesList = ["Transmission"; "Idle/EIFS/SIFS"; "Reception Success"; "Reception Failure"];
                obj.pStateColors = [obj.pTxColor;obj.pIdleEIFSColor;obj.pRxColorSuccess;obj.pRxColorFailure];
                if obj.pNumLENodes==0 && obj.pNumBREDRNodes>0
                    typeList = "Bluetooth BR/EDR";
                    obj.pFreqColors = obj.pCommunicationInfoBREDR{2};
                elseif obj.pNumBREDRNodes==0 && obj.pNumLENodes>0
                    typeList = "Bluetooth LE";
                    obj.pFreqColors = obj.pCommunicationInfoLE{2};
                elseif obj.pNumBREDRNodes>0 && obj.pNumLENodes>0
                    typeList = ["Bluetooth LE"; "Bluetooth BR/EDR"];
                    obj.pFreqColors = [obj.pCommunicationInfoLE{2};obj.pCommunicationInfoBREDR{2}];
                end
            elseif obj.pNumWLANNodes+obj.pNumInterferingNodes>0 && obj.pNumBREDRNodes+obj.pNumLENodes+obj.pNumLE6GHzNodes>0
                statesList = ["Transmission"; "Idle/EIFS/SIFS"; "Contention"; "Reception (Destined to node)"; ...
                    "Reception"; "Reception Success"; "Reception Failure"];
                obj.pStateColors = [obj.pTxColor;obj.pIdleEIFSColor;obj.pContendColor; ...
                    obj.pRxColorUs;obj.pRxColorOthers;obj.pRxColorSuccess;obj.pRxColorFailure];
                if obj.pNumLENodes==0 && obj.pNumBREDRNodes>0
                    typeList = ["WLAN"; "Bluetooth BR/EDR"];
                    obj.pFreqColors = [obj.pCommunicationInfoWLAN{2};obj.pCommunicationInfoBREDR{2}];
                elseif obj.pNumBREDRNodes==0 && obj.pNumLENodes>0
                    typeList = ["WLAN"; "Bluetooth LE"];
                    obj.pFreqColors = [obj.pCommunicationInfoWLAN{2};obj.pCommunicationInfoLE{2}];
                elseif obj.pNumBREDRNodes>0 && obj.pNumLENodes>0
                    typeList = ["WLAN"; "Bluetooth LE"; "Bluetooth BR/EDR"];
                    obj.pFreqColors = [obj.pCommunicationInfoWLAN{2};obj.pCommunicationInfoLE{2}; ...
                        obj.pCommunicationInfoBREDR{2}];
                elseif obj.pNumBREDRNodes==0 && obj.pNumLENodes==0
                    typeList = "WLAN";
                    obj.pFreqColors = obj.pCommunicationInfoWLAN{2};
                end
            end
            if obj.pNumLE6GHzNodes>0
                if ~any(statesList=="Contention")
                    statesList(end+1) = "Contention";
                    obj.pStateColors(end+1,:) = obj.pContendColor;
                end
                typeList(end+1) = "Bluetooth LE 6 GHz";
                obj.pFreqColors(end+1,:) = obj.pCommunicationInfoLE6GHz{2};
            end
            if obj.pNumInterferingNodes>0 && isa(obj.InterferingNodes{1},"helperInterferingBluetoothNode")
                typeList(end+1) = "Interfering Bluetooth LE";
                obj.pFreqColors(end+1,:) = obj.pCommunicationInfoLE{2};
            end
            obj.pStatePatch = cell(0,size(obj.pStateColors,1));
            obj.pFreqPatch = cell(0,size(obj.pFreqColors,1));

            % Write the labels for state
            rowIdxStart = 0;
            stateLabel = createLabel(obj,stateLegendGrid,stateTransitionHeading, ...
                strjoin([stateTransitionHeading "Label"]),14,"Bold",[rowIdxStart+1 rowIdxStart+3],[1 4]);
            rowIdxStart = stateLabel.Layout.Row(2);
            rowIdxStartTemp = rowIdxStart;
            for idx = 1:numel(statesList)
                rowIdxEnd = rowIdxStart+3;
                stateLabel = createLabel(obj,stateLegendGrid,statesList(idx),strjoin([statesList(idx) "Label"]), ...
                    14,"normal",[rowIdxStart+1 rowIdxEnd],[2 4]);
                rowIdxStart = rowIdxEnd;
            end

            % Draw box for legend to show colors
            rowIdxStart = rowIdxStartTemp;
            for idx = 1:numel(statesList)
                rowIdxEnd = rowIdxStart+3;
                stateBoxPanel = createPanel(obj,stateLegendGrid,"State Legend Panel",[rowIdxStart+1 rowIdxEnd],1);
                stateBoxPanel.BackgroundColor = obj.pStateColors(idx,:);
                rowIdxStart = rowIdxEnd;
            end

            if obj.FrequencyPlotFlag
                % If frequency plot is enabled create legend for frequency axes
                freqLegendPanel = createPanel(obj,figureGrid,"Freq Comm Legend Panel",[9 11],9);
                freqLegendPanel.Clipping = "off";
                freqLegendGrid = uigridlayout(freqLegendPanel,[100 4]);
                for idx = 1:numel(freqLegendGrid.ColumnWidth)
                    freqLegendGrid.ColumnWidth{idx} = 55;
                end

                rowIdxStart = 0;
                for idx = 1:numel(typeList)
                    rowIdxEnd = rowIdxStart+3;
                    freqColorPanel = createPanel(obj,freqLegendGrid,"Freq Legend Color Panel",[1+rowIdxStart rowIdxEnd],1);
                    freqColorPanel.BackgroundColor = obj.pFreqColors(idx,:);
                    stateLabel = createLabel(obj,freqLegendGrid,typeList(idx),strjoin([typeList(idx) "Label"]), ...
                        14,"normal",[rowIdxStart+1 rowIdxEnd],[2 4]);
                    rowIdxStart = rowIdxEnd;
                end
            end

            % Implements the slider for time minimum
            editLabelHandle = createLabel(obj,figureGrid,"Min Time","Time Min Label",14,"bold",sliderRow-1,3);
            obj.pEditMinHandle = uieditfield(figureGrid,"numeric","Tag","Time Min Edit Field","Limits",[0 simulationTime]);
            assignRowColumnToLayout(obj,obj.pEditMinHandle,sliderRow-1,4);
            obj.pEditMinHandle.ValueChangedFcn = @(~,eventData) objWeakRef.Handle.updateXAxisLimits(eventData.Value,true);
            obj.pMinTimeSlider = uislider(figureGrid,"Tag","Time Min Slider");
            assignRowColumnToLayout(obj,obj.pMinTimeSlider,sliderRow,[3 5]);
            obj.pMinTimeSlider.Limits = [0 simulationTime];
            obj.pMinTimeSlider.ValueChangedFcn = @(~,eventData) objWeakRef.Handle.updateXAxisLimits(eventData.Value,true);

            % Implements the slider for time maximum
            editLabelHandle = createLabel(obj,figureGrid,"Max Time","Time Max Label",14,"bold",sliderRow-1,7);
            obj.pEditMaxHandle = uieditfield(figureGrid,"numeric","Tag","Time Max Edit Field","Limits",[0 simulationTime]);
            assignRowColumnToLayout(obj,obj.pEditMaxHandle,sliderRow-1,8);
            obj.pEditMaxHandle.ValueChangedFcn = @(edtHandle,eventData) objWeakRef.Handle.updateXAxisLimits(eventData.Value,false);
            obj.pMaxTimeSlider = uislider(figureGrid,"Tag","Time Max Slider");
            assignRowColumnToLayout(obj,obj.pMaxTimeSlider,sliderRow,[6 8]);
            obj.pMaxTimeSlider.Limits = [0 simulationTime];
            obj.pMaxTimeSlider.ValueChangedFcn = @(~,eventData) objWeakRef.Handle.updateXAxisLimits(eventData.Value,false);

            % Implements the drop down for frequency
            maxLimitsYFreq = [];
            wlanFreqLimits = [];
            if obj.FrequencyPlotFlag
                limitsYFreq = [2400 2500];
                maxLimitsYFreq = limitsYFreq; % For full view
                ddItems = [];
                if ~isempty(wlanFreqInfo)
                    wlanFreqLimits(1,:) = (wlanFreqInfo(:,2)-wlanFreqInfo(:,3))/1e6;
                    wlanFreqLimits(2,:) = (wlanFreqInfo(:,2)+wlanFreqInfo(:,3))/1e6;
                end
                % Based on the presence of nodes, find the possible frequency bands
                if obj.pNumBREDRNodes+obj.pNumLENodes==0
                    if all(wlanFreqInfo(:,2)/1e6<=2500)
                    elseif all(wlanFreqInfo(:,2)/1e6>=5000) && all(wlanFreqInfo(:,2)/1e6<=5900)
                        limitsYFreq = [5000 6000];
                        maxLimitsYFreq = limitsYFreq;
                    elseif all(wlanFreqInfo(:,2)/1e6>=5900)
                        limitsYFreq = [5900 7200];
                        maxLimitsYFreq = limitsYFreq;
                    elseif any(wlanFreqInfo(:,2)/1e6>=2400) && any(wlanFreqInfo(:,2)/1e6>=5000) && ...
                            all(wlanFreqInfo(:,2)/1e6<=5900)
                        ddItems = ["2.4 GHz band";"5 GHz band"];
                        maxLimitsYFreq(2) = 6000;
                    elseif any(wlanFreqInfo(:,2)/1e6>=2400) && any(wlanFreqInfo(:,2)/1e6>=5900) && ...
                            ~(any(wlanFreqInfo(:,2)/1e6>=5000) && any(wlanFreqInfo(:,2)/1e6<=5900))
                        ddItems = ["2.4 GHz band";"6 GHz band"];
                        maxLimitsYFreq(2) = 7200;
                    elseif all(wlanFreqInfo(:,2)/1e6>=5000)
                        ddItems = ["5 GHz band";"6 GHz band"];
                        maxLimitsYFreq = [5000 7200];
                    else
                        ddItems = ["2.4 GHz band";"5 GHz band";"6 GHz band"];
                        maxLimitsYFreq(2) = 7200;
                    end
                else
                    if obj.pNumWLANNodes+obj.pNumInterferingNodes>0
                        if (any(wlanFreqInfo(:,2)/1e6>=5000) && any(wlanFreqInfo(:,2)/1e6<=5900)) && ...
                                any(wlanFreqInfo(:,2)/1e6>=5900)
                            ddItems = ["2.4 GHz band";"5 GHz band";"6 GHz band"];
                            maxLimitsYFreq(2) = 7200;
                        elseif any(wlanFreqInfo(:,2)/1e6>=5000) && all(wlanFreqInfo(:,2)/1e6<=5900)
                            ddItems = ["2.4 GHz band";"5 GHz band"];
                            maxLimitsYFreq(2) = 6000;
                        elseif any(wlanFreqInfo(:,2)/1e6>=5900)
                            ddItems = ["2.4 GHz band";"6 GHz band"];
                            maxLimitsYFreq(2) = 7200;
                        end
                    end
                end
                if obj.pNumLE6GHzNodes>0
                    ddItems = ["6 GHz band";"2.4 GHz band";"5 GHz band"];
                    limitsYFreq = [5900 7200];
                    maxLimitsYFreq = limitsYFreq;
                end
                obj.pMinMaxFreqLimits = maxLimitsYFreq;
                obj.pWLANFreqLimits = wlanFreqLimits;
                if numel(ddItems)>1
                    % More than 1 frequency band exists for frequency axes, add a drop down
                    stateLabel = createLabel(obj,freqLegendGrid,"Freq Band","Freq Drop Down Label", ...
                        14,"normal",[rowIdxStart+1 rowIdxStart+5],[1 3]);
                    freqWindowHandle = uidropdown(freqLegendGrid,"Items",ddItems,"Value",ddItems(1),"Tag","Freq Drop Down");
                    assignRowColumnToLayout(obj,freqWindowHandle,[rowIdxStart+1 rowIdxStart+5],[4 6]);
                    freqWindowHandle.ValueChangedFcn = @(freqVarHandle,eventData) objWeakRef.Handle.freqDDCallback(eventData);
                end
                calcYTicks(obj,limitsYFreq);
            end

            % Panel for axes configurations
            cfgPanel = createPanel(obj,figureGrid,"Config Panel",[5 8],9);
            cfgGrid = uigridlayout(cfgPanel,[2 2],"ColumnWidth",{'fit','fit'});

            % Check box for full view
            checkBoxHandle = uicheckbox(cfgGrid,"Text","Full view","FontSize",14,"Tag","Full View");
            assignRowColumnToLayout(obj,checkBoxHandle,1,1);
            checkBoxHandle.ValueChangedFcn = @(checkHandle,eventData) objWeakRef.Handle.fullViewUpdate();

            % Save Snip button
            saveImageHandle = uibutton(cfgGrid,"Text","Export","FontSize",14,"Tag","Export");
            assignRowColumnToLayout(obj,saveImageHandle,2,1);
            saveImageHandle.ButtonPushedFcn = @(btnHandle,eventData) objWeakRef.Handle.saveSnipCallback(obj.PacketCommUIFigure);

            % Implements scrollbar for more number of nodes
            if numel(yTickLabels)>obj.pMaxNodesToDisplay
                verticalSliderHandle = uislider(figureGrid,"Orientation","vertical","Limits",[1 numel(yTickLabels)], ...
                    "MajorTicks",1:numel(yTickLabels),"MinorTicks",[],"Tag","Node Slider");
                if obj.FrequencyPlotFlag
                    assignRowColumnToLayout(obj,verticalSliderHandle,[3 5],1);
                else
                    assignRowColumnToLayout(obj,verticalSliderHandle,[3 8],1);
                end
                verticalSliderHandle.ValueChangedFcn = @(sliderUI,eventData) objWeakRef.Handle.nodeSliderCallback(eventData);
            end

            % Adds the line at t=0
            updateLine(obj,[],0)

            % Add callbacks to X axis of the state and frequency axes
            obj.StateAxes.XAxis.LimitsChangedFcn = @(src,evt) objWeakRef.Handle.xLimChangedFcn(src,evt);
            if obj.FrequencyPlotFlag
                obj.FreqAxes.XAxis.LimitsChangedFcn = @(src,evt) objWeakRef.Handle.xLimChangedFcn(src,evt);
            end

            % Disable interactions during the run time of simulation
            updateInteractivityTime(obj,false)
            updateInteractivityAxes(obj,false)

            % Initialize the properties
            preAllocValue = 1e6;
            obj.StartTimeList = zeros(1,preAllocValue);
            obj.EndTimeList = zeros(1,preAllocValue);
            obj.YTicksList = zeros(1,preAllocValue);
            obj.CommunicationInfo = repmat(obj.CommunicationInfoStruct,1,preAllocValue);
            obj.pXLimits = [0 simulationTime];

            % Pauses for all the actions to reflect for live visualization
            if ~obj.VisualizeAtEnd
                drawnow
                waitfor(obj.PacketCommUIFigure, "FigureViewReady", "on")
            end
        end

        function resizeAxes(obj,panelHandle,axess,flag)
            %resizeAxes Updates the size of the axes to make it have starting position
            %in UI figure

            axess.Position(4) = panelHandle.InnerPosition(4)-panelHandle.BorderWidth-10;
            axess.Position(3) = panelHandle.Position(3)-10;
            if ~flag
                if obj.FrequencyPlotFlag % Frequency axes enabled
                    obj.FreqAxes.InnerPosition(1) = obj.StateAxes.InnerPosition(1);
                    obj.FreqAxes.InnerPosition(3) = obj.StateAxes.InnerPosition(3);
                    obj.FreqAxes.Position(4) = panelHandle.InnerPosition(4)-panelHandle.BorderWidth-10;
                end
            else
                obj.StateAxes.Position(4) = panelHandle.InnerPosition(4)-panelHandle.BorderWidth-10;
                obj.StateAxes.Position(3) = panelHandle.Position(3)-10;
            end
        end

        function panelHandle = createPanel(obj,parentHandle,tagName,row,column)
            %createPanel Creates UI panel

            panelHandle = uipanel(parentHandle,"Scrollable","on","Tag",tagName);
            assignRowColumnToLayout(obj,panelHandle,row,column);
        end

        function labelHandle = createLabel(obj,parentHandle,labelName,tagName,fontSize,fontWeight,row,column)
            %createLabel Creates UI label

            labelHandle = uilabel(parentHandle,"Text",labelName,"WordWrap","on","FontWeight",fontWeight',"Tag",tagName);
            assignRowColumnToLayout(obj,labelHandle,row,column)
            if ~isempty(fontSize)
                labelHandle.FontSize = fontSize;
            end
        end

        function assignRowColumnToLayout(~,uiHandle,row,column)
            %assignRowColumnToLayout Assigns row and column to the specified UI handle

            if ~isempty(row)
                uiHandle.Layout.Row = row;
            end
            if ~isempty(column)
                uiHandle.Layout.Column = column;
            end
        end

        function calcYTicks(obj,limitsYFreq)
            %calcYTicks Based on the Y Frequency limits returns the Y Ticks

            % Get the Y Frequency limits
            maxFreqLimits = obj.pMinMaxFreqLimits;
            wlanFreqLimits = obj.pWLANFreqLimits;
            if ~isempty(wlanFreqLimits)
                minWLANLim = min(maxFreqLimits(1),min(wlanFreqLimits(1,:)));
                maxWLANLim = max(maxFreqLimits(2),max(wlanFreqLimits(2,:)));
            end

            % Update frequency axes limits
            obj.FreqAxes.YLim = limitsYFreq;

            % Calculate the nominal Y frequency limits
            limitsYFreq = [ceil(limitsYFreq(1)/100)*100 floor(limitsYFreq(2)/100)*100];

            % Calculate the Y Tick labels and check if the minimum and maximum Y limits
            % are consistent with WLAN frequency limits and update if necessary
            switch string(limitsYFreq)
                case ["2400" "2500"]
                    yTickLabels = {'2400';'2410';'2420';'2430';'2440';'2450'; ...
                        '2460';'2470';'2480';'2490';'2500'};
                    if ~isempty(wlanFreqLimits)
                        if any(minWLANLim<2400)
                            yTickLabels = [{num2str(minWLANLim(minWLANLim<2400))};yTickLabels];
                            obj.FreqAxes.YLim(1) = minWLANLim(minWLANLim<2400);
                        end
                        if any(maxWLANLim>2500 & maxWLANLim<6000)
                            yTickLabels{end+1} = num2str(maxWLANLim(maxWLANLim>2500 & maxWLANLim<6000));
                            obj.FreqAxes.YLim(2) = maxWLANLim(maxWLANLim>2500 & maxWLANLim<6000);
                        end
                    end
                case ["2400" "6000"]
                    yTickLabels = {'2400';'2500';'3000';'3500';'4000'; ...
                        '4500';'5000';'5500';'6000'};
                case ["2400" "7200"]
                    yTickLabels = {'2400';'2500';'3000';'4000';'5000'; ...
                        '5500';'6000';'6500';'7000';'7200'};
                case ["5000" "6000"]
                    yTickLabels = {'5000';'5100';'5200';'5300';'5400'; ...
                        '5500';'5600';'5700';'5800';'5900';'6000'};
                    if ~isempty(wlanFreqLimits)
                        if any(minWLANLim<5000 & minWLANLim>2500)
                            yTickLabels = [{num2str(minWLANLim(minWLANLim<5000))};yTickLabels];
                            obj.FreqAxes.YLim(1) = minWLANLim(minWLANLim<5000);
                        end
                        if any(maxWLANLim>6000 & maxWLANLim<7200)
                            yTickLabels{end+1} = num2str(maxWLANLim(maxWLANLim>2500 & maxWLANLim<6000));
                            obj.FreqAxes.YLim(2) = maxWLANLim(maxWLANLim>2500 & maxWLANLim<6000);
                        end
                    end
                case ["5900" "7200"]
                    yTickLabels = {'6000';'6100';'6200';'6300';'6400';'6500';'6600';'6700';'6800';'6900';'7000';'7100';'7200'};
                    if ~isempty(wlanFreqLimits)
                        if any(minWLANLim<6000 & minWLANLim>5800)
                            yTickLabels = [{num2str(minWLANLim(minWLANLim<6000))};yTickLabels];
                            obj.FreqAxes.YLim(1) = minWLANLim(minWLANLim<6000);
                        end
                        if any(maxWLANLim>7200)
                            yTickLabels{end+1} = num2str(maxWLANLim(maxWLANLim>7200));
                            obj.FreqAxes.YLim(2) = maxWLANLim(maxWLANLim>7200);
                        end
                    end
            end

            % Update Y Ticks and labels
            obj.FreqAxes.YTick = str2double(string(cell2mat(yTickLabels)));
            obj.FreqAxes.YTickLabel = yTickLabels;
        end
    end

    %% Standard Based:Figure relevant Methods
    methods (Access=private)
        function [tickIdx,tickBase,yTicksList,yTickLabels] = yTicksBluetooth(obj,btIdx,techName,tickIdx,tickBase,yTicksList,yTickLabels,node,deviceID)
            %yTicksBluetooth Update and return Y ticks of Bluetooth nodes

            obj.pYTickOrder(tickIdx,:) = [node.ID deviceID btIdx];
            yTickLabels{tickIdx} = strrep(node.Name,"_","\_");
            if isa(node,"helperCoexNode")
                tickBase = tickBase + obj.pBarHeight*0.5;
                if ~contains(node.Name,techName)
                    nodeName = strrep(node.Name,"_","\_");
                    yTickLabels{tickIdx} = strjoin([nodeName techName]);
                end
            else
                tickBase = tickBase + obj.pBarHeight*1.25;
            end
            tickBase = tickBase + obj.pBarHeight*1.25;
            yTicksList(tickIdx) = tickBase;
            tickIdx = tickIdx+1;
            tickBase  = tickBase+obj.pBarHeight;
        end

        function [tickIdx,tickBase,yTicksList,yTickLabels,wlanFreqInfo] = yTicksIntWLAN(obj,intWLANIdx,tickIdx,tickBase,yTicksList,yTickLabels,node)
            %yTicksIntWLAN Update and return Y ticks of interfering WLAN nodes

            obj.pYTickOrder(tickIdx,:) = [node.ID 0 intWLANIdx];
            nodeName = strrep(node.Name,"_","\_");
            yTickLabels{tickIdx} = strjoin([nodeName num2str(node.CenterFrequency/1e6) "MHz"]);
            tickBase = tickBase + obj.pBarHeight*2.5;
            yTicksList(tickIdx) = tickBase;
            tickIdx = tickIdx+1;
            tickBase  = tickBase+obj.pBarHeight;
            wlanFreqInfo = [node.ID node.CenterFrequency node.Bandwidth];
        end

        function [tickIdx,tickBase,yTicksList,yTickLabels,wlanFreqInfo] = yTicksWLAN(obj,tickIdx,tickBase,yTicksList,yTickLabels,nodes)
            %yTicksWLAN Update and return the WLAN frequency information and tick
            %information

            % WLAN information for Y Ticks
            numWLAN = numel(nodes);
            wlanFreqInfo = zeros(numWLAN*3,2);
            if numWLAN>0
                isCoex = isa(nodes{1},"helperCoexNode");
                if isCoex
                    tickBase = tickBase + obj.pBarHeight*0.5;
                    node = nodes{1}.WLANDevice;
                    node.ID = nodes{1}.ID;
                    node.Name = nodes{1}.Name;
                    nodeName = strrep(node.Name,"_","\_");
                else
                    obj.pPPDUDuration = zeros(numWLAN,3);
                end
                wlanFreqInfo = zeros(numWLAN*3,3);
                count = 1;

                for idx = 1:numWLAN
                    if ~isCoex
                        node = nodes{idx};
                        freq = node.ReceiveFrequency;
                        devCfg = getDeviceConfig(obj,node);
                        bw = [devCfg.ChannelBandwidth];
                        nodeName = strrep(node.Name,"_","\_");
                        deviceIDs = 1:numel(freq);
                    else
                        freq = node.ReceiveFrequency;
                        devCfg = getDeviceConfig(obj,node);
                        bw = [devCfg.ChannelBandwidth];
                        deviceIDs = [node.DeviceID];
                    end
                    numFreq =  numel(freq);
                    wlanFreqInfo(count:count+numFreq-1,:) = [node.ID*ones(numFreq,1) freq' bw']; % Node ID Frequency Bandwidth
                    obj.pYTickOrder(size(obj.pYTickOrder,1)+1:size(obj.pYTickOrder,1)+numFreq,:) =  ...
                        [node.ID*ones(numFreq,1) deviceIDs' ones(numFreq,1)];
                    count = count+numFreq;
                    tickBase  = tickBase + obj.pBarHeight*1.25;
                    for inIdx = 1:numFreq
                        yTicksList(tickIdx) = tickBase;
                        if isCoex
                            yTickLabels{tickIdx} = strjoin([nodeName "WLAN" num2str(freq(inIdx)/1e6) "MHz"]);
                        else
                            yTickLabels{tickIdx} = strjoin([nodeName num2str(freq(inIdx)/1e6) "MHz"]);
                        end
                        tickIdx = tickIdx+1;
                        tickBase  = tickBase + obj.pBarHeight*1.5;
                    end
                    tickBase  = tickBase - obj.pBarHeight*0.5;
                end
                wlanFreqInfo = wlanFreqInfo(1:count-1,:);
            end
        end
    end

    %% Live Update Functions
    methods (Access=private)
        function drawCommunicationInfo(obj,xmin,xmax,startNodeIdx,endNodeIdx)
            %drawCommunicationInfo Draw the communication information stored on both
            %axes

            % If previously and currently called time and nodes are same, return
            if isequal(obj.pLastCalledTime,[xmin xmax]) && isequal(obj.pLastCalledNodes,[startNodeIdx endNodeIdx])
                return;
            end

            % Draw the information
            packetCount = obj.InfoCount;
            if packetCount>0
                % During run time, calculate time and nodes to be plotted
                if isempty(xmin)
                    % Calculate the start and end time
                    xmin = obj.pLastEndTimePlotted;
                    xmax = obj.CommunicationInfo(packetCount).StartTime+obj.CommunicationInfo(packetCount).Duration;
                    if xmin>=xmax
                        xmin = obj.pLastCalledTime(1);
                    end
                end
                if isempty(startNodeIdx)
                    % Find the Y axis tick index
                    yticksList = find(obj.StateAxes.YTick>=obj.StateAxes.YLim(1));
                    startNodeIdx = yticksList(1);
                    yticksList = find(obj.StateAxes.YTick<=obj.StateAxes.YLim(2));
                    endNodeIdx = yticksList(end);
                end

                % Find the indices of all packets that fall in the time and node interval
                startTimeList = obj.StartTimeList(1:packetCount);
                endTimeList = obj.EndTimeList(1:packetCount);
                yTicksList = obj.YTicksList(1:packetCount);
                objIdxsToDraw = find(yTicksList>=startNodeIdx & yTicksList<=endNodeIdx & ...
                    endTimeList>=xmin & startTimeList<=xmax);

                if numel(objIdxsToDraw)>0
                    if ~isequal(obj.pObjIdxsOnDisp,objIdxsToDraw)
                        % Draw the patch objects
                        drawStateTransitionPatch(obj,objIdxsToDraw);
                        if obj.FrequencyPlotFlag
                            drawChannelOccupancyPatch(obj,objIdxsToDraw);
                        end
                        obj.pLastCalledTime = [xmin xmax];
                        obj.pLastCalledNodes = [startNodeIdx endNodeIdx];
                    end
                    % Update the axes limits
                    set(obj.StateAxes,"xlim",[xmin,xmax]);
                    if obj.FrequencyPlotFlag
                        set(obj.FreqAxes,"xlim",[xmin,xmax]);
                    end
                    updateXLim(obj,xmin,xmax)
                    updateYLim(obj,startNodeIdx,endNodeIdx)
                    drawnow
                else
                    set(obj.StateAxes,"xlim",[xmin,xmax]);
                    set(obj.FreqAxes,"xlim",[xmin,xmax]);
                    updateXLim(obj,xmin,xmax)
                    updateYLim(obj,startNodeIdx,endNodeIdx)
                end

                % Initialize the last axes limits plotted and displayed
                obj.pLastCalledTime = [xmin xmax];
                obj.pLastCalledNodes = [startNodeIdx endNodeIdx];
                ns = wirelessNetworkSimulator.getInstance;
                obj.pLastEndTimePlotted = ns.CurrentTime;
                obj.pObjIdxsOnDisp = objIdxsToDraw;
            end
        end

        function drawStateTransitionPatch(obj,objIdxsToDraw)
            %drawStateTransitionPatch Draw the patch object for state transition axes

            % Segregate the packet indices based on the state colors
            toDrawState = cell(1,size(obj.pStateColors,1));
            for packetIdx = objIdxsToDraw
                commInfo  = obj.CommunicationInfo(packetIdx);
                packetColorIdx = find(commInfo.FillColor(2)==obj.pStateColors(:,2));

                % Assign the packet index to the particular state color
                if isempty(toDrawState{packetColorIdx}) || size(toDrawState,2)<packetColorIdx
                    toDrawState{packetColorIdx} = packetIdx;
                else
                    toDrawState{packetColorIdx} = [toDrawState{packetColorIdx} packetIdx];
                end
            end

            % Calculate the vertices and faces of the rectangles for the patches
            for colorIdx = size(toDrawState,2):-1:1
                packetIndices = toDrawState{colorIdx};
                if numel(packetIndices)>0
                    verticesCount = 0;
                    vertices = [0 0];
                    nodePositions = ones(1,numel(packetIndices));
                    for idx = 1:numel(packetIndices)
                        packetIdx = packetIndices(idx);
                        commInfo = obj.CommunicationInfo(packetIdx);
                        startTime = commInfo.StartTime;
                        duration = commInfo.Duration;
                        nodeYTickPos = obj.StateAxes.YTick(commInfo.YTickVal) - obj.pBarHeight/2;
                        nodePositions(idx) = nodeYTickPos;
                        height = obj.pBarHeight;
                        vertices(verticesCount+1,:) = [startTime nodeYTickPos];
                        vertices(verticesCount+2,:) = [startTime nodeYTickPos+height];
                        vertices(verticesCount+3,:) = [startTime+duration nodeYTickPos+height];
                        vertices(verticesCount+4,:) = [startTime+duration nodeYTickPos];
                        verticesCount = verticesCount+4;
                    end
                    faces = reshape(1:verticesCount,4,[])';

                    % Draw the patch objects
                    if isempty(obj.pStatePatch) || size(obj.pStatePatch,2)<colorIdx || ...
                            isempty(obj.pStatePatch{colorIdx})
                        obj.pStatePatch{colorIdx} = patch(obj.StateAxes,Vertices=vertices,Faces=faces, ...
                            FaceColor=obj.pStateColors(colorIdx,:),EdgeColor="none");
                        patchObj = obj.pStatePatch{colorIdx};
                    else
                        patchObj = obj.pStatePatch{colorIdx};
                        patchObj.Vertices = vertices;
                        patchObj.Faces = faces;
                    end
                    alpha = ones(size(faces,1),1);
                    patchObj.FaceAlpha = "flat";
                    patchObj.EdgeAlpha = "flat";
                    patchObj.FaceVertexAlphaData = alpha;
                    patchObj.AlphaDataMapping = "none";
                end
            end
        end

        function drawChannelOccupancyPatch(obj,objIdxsToDraw)
            %drawChannelOccupancyPatch Draw the patch object for channel occupancy axes

            % Segregate the packet indices based on the frequency colors
            toDrawState = cell(1,size(obj.pFreqColors,1));
            for packetIdx = objIdxsToDraw
                commInfo = obj.CommunicationInfo(packetIdx);
                packetColorIdx = find(commInfo.FreqColor(3)==obj.pFreqColors(:,3));

                % Assign the packet index to the particular frequency color
                if isempty(toDrawState{packetColorIdx}) || size(toDrawState,2)<packetColorIdx
                    toDrawState{packetColorIdx} = packetIdx;
                else
                    toDrawState{packetColorIdx} = [toDrawState{packetColorIdx} packetIdx];
                end
            end

            % Calculate the vertices and faces of the rectangles for the patches
            for colorIdx = size(toDrawState,2):-1:1
                packetIndices = toDrawState{colorIdx};
                if numel(packetIndices)>0
                    verticesCount = 0;
                    vertices = [0 0];
                    freqFaceVal = {};
                    freqFaceIdx = 1;
                    for idx = 1:numel(packetIndices)
                        packetIdx = packetIndices(idx);
                        commInfo = obj.CommunicationInfo(packetIdx);
                        if commInfo.FreqPlot
                            startTime = commInfo.StartTime;
                            duration = commInfo.Duration;
                            channelFrequencyStart = commInfo.FreqStart;
                            channelBandwidth = commInfo.ChannelBandwidth;
                            numFreq = numel(channelFrequencyStart);
                            nodeYTickPos = obj.StateAxes.YTick(commInfo.YTickVal) - obj.pBarHeight/2;
                            for freqIdx = 1:numFreq
                                vertices(verticesCount+1,:) = [startTime channelFrequencyStart(freqIdx)];
                                vertices(verticesCount+2,:) = [startTime ...
                                    channelFrequencyStart(freqIdx)+channelBandwidth(freqIdx)];
                                vertices(verticesCount+3,:) = [startTime+duration ...
                                    channelFrequencyStart(freqIdx)+channelBandwidth(freqIdx)];
                                vertices(verticesCount+4,:) = [startTime+duration ...
                                    channelFrequencyStart(freqIdx)];
                                verticesCount = verticesCount+4;
                            end
                            freqFaceVal{idx} = freqFaceIdx:freqFaceIdx+numFreq-1;
                            freqFaceIdx = freqFaceIdx+numFreq;
                        end
                    end
                    faces = reshape(1:verticesCount,4,[])';

                    % Draw the patch objects
                    if isempty(obj.pFreqPatch) || size(obj.pFreqPatch,2)<colorIdx || ...
                            isempty(obj.pFreqPatch{colorIdx})
                        obj.pFreqPatch{colorIdx} = patch(obj.FreqAxes,Vertices=vertices,Faces=faces, ...
                            FaceColor=obj.pFreqColors(colorIdx,:),EdgeColor="none");
                        patchObj = obj.pFreqPatch{colorIdx};
                    else
                        patchObj = obj.pFreqPatch{colorIdx};
                        patchObj.Vertices = vertices;
                        patchObj.Faces = faces;
                    end
                    alpha = ones(size(faces,1),1)*0.5;
                    patchObj.FaceAlpha = "flat";
                    patchObj.EdgeAlpha = "flat";
                    patchObj.FaceVertexAlphaData = alpha;
                    patchObj.AlphaDataMapping = "none";
                end
            end
        end

        function updateXLim(obj,xLower,xHigher)
            %updateXLim Update X axis limits of the figure and the interactivity
            %options

            % Calculate the X limits
            xLower = max(xLower,0);
            xHigher = min(xHigher,obj.pMaxEndTime);

            % Update the edit box value and slider value
            obj.pEditMinHandle.Value = xLower;
            obj.pEditMaxHandle.Value = xHigher;
            obj.pMinTimeSlider.Value = xLower;
            obj.pMaxTimeSlider.Value = xHigher;

            % Set the axes limits
            set(obj.StateAxes,"xlim",[xLower,xHigher]);
            if obj.FrequencyPlotFlag
                set(obj.FreqAxes,"xlim",[xLower,xHigher]);
            end
        end

        function updateYLim(obj,startNode,endNode)
            %updateYLim Update Y axis limits of the figure and the interactivity
            %options

            sliderHandle = findobj(obj.PacketCommUIFigure,Tag="Node Slider");
            if ~isempty(sliderHandle)
                % Calculate the Y limits
                yLower = obj.StateAxes.YTick(startNode)-1.25;
                yHigher = obj.StateAxes.YTick(endNode)+1.25;

                % Update the slider limits
                sliderHandle.Value = endNode;

                % Set the axes limits
                set(obj.StateAxes,"ylim",[yLower,yHigher]);
            end
        end

        function updateFigureAtSimulationEnd(obj,varargin)
            %updateFigureAtSimulationEnd Update the visualization at the end of
            %simulation

            % Simulation is over
            obj.pIsSimInProgress = false;

            % Calculate maximum end time to be displayed
            obj.pMaxEndTime = max(obj.pMaxEndTime,obj.pSimulationTime);

            % Enable interactions of the axes after the simulation end
            updateInteractivityAxes(obj,true)

            % Calculate the minimum and maximum limits
            numPacketsToShow = obj.pNumPacketsToShow;
            endTime = obj.pMaxEndTime;
            if obj.InfoCount+obj.InfoCount>numPacketsToShow
                startTime = obj.CommunicationInfo(end-numPacketsToShow/2+1).StartTime;

                % Draw the information for the limits and update the axes
                drawCommunicationInfo(obj,startTime,endTime,[],[]);
            else
                checkBoxHandle = findobj(obj.PacketCommUIFigure,Tag="Full View");
                checkBoxHandle.Value = true;
                fullViewUpdate(obj);
            end

            % Add line to the state figure at simulation end time
            obj.pTimeLine1 = xline(obj.StateAxes,endTime,"LineWidth",1,"Tag","Time line 1","Color",obj.pLineColor);
            if obj.FrequencyPlotFlag
                obj.pTimeLine2 = xline(obj.FreqAxes,endTime,"LineWidth",1,"Tag","Time line 2","Color",obj.pLineColor);
            end

            % Update the slider and edit box limits
            obj.pMinTimeSlider.Limits(2) = endTime;
            obj.pMaxTimeSlider.Limits(2) = endTime;
            obj.pEditMaxHandle.Limits = [0 endTime];

            % Display the figure
            obj.PacketCommUIFigure.Visible = "on";
        end

        %% Callbacks to the UI Components
        function xLimChangedFcn(obj,~,eventData)
            %xLimChangedFcn Callback function when X limits has been updated in the
            %axes

            % Get minimum and maximum X limits
            xLower = eventData.NewLimits(1);
            xHigher = eventData.NewLimits(2);

            % Calculate minimum and maximum Y limits
            yLower = find(obj.StateAxes.YTick > obj.StateAxes.YLim(1)-1.25);
            yHigher = find(obj.StateAxes.YTick < obj.StateAxes.YLim(end)+1.25);
            yLower = yLower(1);
            yHigher = yHigher(end);

            % Calculate minimum and maximum X limits when X limits are not within the
            % limits
            numPacketsToShow = obj.pNumPacketsToShow;
            if xLower<0
                xLower = 0;
                if xHigher<0
                    if obj.InfoCount*2>=numPacketsToShow
                        xHigher = obj.CommunicationInfo(numPacketsToShow).StartTime;
                    else
                        xHigher = obj.pMaxEndTime;
                    end
                end
            end
            if xHigher>obj.pMaxEndTime
                xHigher = obj.pMaxEndTime;
                if xLower>obj.pMaxEndTime
                    if obj.InfoCount*2>=numPacketsToShow
                        xLower = obj.CommunicationInfo(end-numPacketsToShow+1).StartTime;
                    else
                        xLower = 0;
                    end
                end
            end

            % Draw the information for the limits and update the axes
            drawCommunicationInfo(obj,xLower,xHigher,yLower,yHigher);
        end

        function updateXAxisLimits(obj,newTimeVal,isMin)
            %updateXAxisLimits Update the time limits of the figures

            % Calculate the limits for the time axis
            xlimit_lower = obj.StateAxes.XLim(1);
            xlimit_higher = obj.StateAxes.XLim(2);
            if isMin
                validateattributes(newTimeVal,'double',{'scalar','<=',xlimit_higher},"Minimum Time Value");
                xlimit_lower = newTimeVal;
            else
                validateattributes(newTimeVal,'double',{'scalar','>=',xlimit_lower},"Maximum Time Value");
                xlimit_higher = newTimeVal;
            end

            % Draw the information for the limits and update the axes
            drawCommunicationInfo(obj,xlimit_lower,xlimit_higher,[],[]);
        end

        function fullViewUpdate(obj)
            %fullViewUpdate Callback for the full view of the X Axis

            checkBoxHandle = findobj(obj.PacketCommUIFigure,Tag="Full View");

            % Get the X and Y axis limits
            if checkBoxHandle.Value
                timeLimit = [0 obj.pMaxEndTime];
                obj.pXLimits = obj.StateAxes.XLim;
                yLower = find(obj.StateAxes.YTick > obj.StateAxes.YLim(1)-1.25);
                yHigher = find(obj.StateAxes.YTick < obj.StateAxes.YLim(end)+1.25);
                obj.pYLimits(1,:) = [yLower(1) yHigher(end)];
                nodeLimit = [1 numel(obj.StateAxes.YTick)];
                if obj.FrequencyPlotFlag
                    obj.pFreqLimits = obj.FreqAxes.YLim;
                end
                freqLimits = obj.pMinMaxFreqLimits;
            else
                timeLimit = obj.pXLimits;
                nodeLimit = [obj.pYLimits(1,1) obj.pYLimits(1,2)];
                freqLimits = obj.pFreqLimits;
            end

            % Draw the information for the limits and update the axes
            drawCommunicationInfo(obj,timeLimit(1),timeLimit(2),nodeLimit(1),nodeLimit(2));

            % Update the Y limits
            if obj.FrequencyPlotFlag
                obj.FreqAxes.YLim = freqLimits;
                calcYTicks(obj,freqLimits);
            end

            % Enable or disable other figure handles accordingly
            updateInteractivityTime(obj,~checkBoxHandle.Value)
            handlesToUpdate = ["Freq Drop Down";"Node Slider"];
            for idx = 1:numel(handlesToUpdate)
                objHandle = findobj(obj.PacketCommUIFigure,Tag=handlesToUpdate(idx));
                if ~isempty(objHandle)
                    objHandle.Enable = ~checkBoxHandle.Value;
                end
            end
        end

        function nodeSliderCallback(obj,eventData)
            %nodeSliderCallback Callback for the vertical slider if there are more
            %nodes

            % Calculate the nominal Y axis limits
            maxYLimitBase = max(obj.StateAxes.YTick(1:obj.pMaxNodesToDisplay)) + obj.pBarHeight;
            value = round(eventData.Value);

            % Check nominal values with specified value and update
            if value<numel(obj.StateAxes.YTick)
                ylimit_higher = obj.StateAxes.YTick(value)+obj.pBarHeight;
            else
                ylimit_higher = obj.StateAxes.YTick(end)+obj.pBarHeight;
            end

            % Reset if max Y limit if necessary
            if ylimit_higher<maxYLimitBase
                ylimit_higher = maxYLimitBase;
            end

            % Calculate minimum Y limit
            ylimit_lower = ylimit_higher-maxYLimitBase-obj.pBarHeight;
            if ylimit_lower<1
                ylimit_lower = 1;
            end

            % Calculate the start and end node index
            yticksList = find(obj.StateAxes.YTick>=ylimit_lower);
            if isempty(yticksList)
                return;
            end
            startNodeIdx = yticksList(1);
            yticksList = find(obj.StateAxes.YTick<=ylimit_higher);
            if isempty(yticksList)
                return;
            end
            endNodeIdx = yticksList(end);

            % Draw the information for the limits and update the axes
            drawCommunicationInfo(obj,obj.StateAxes.XLim(1),obj.StateAxes.XLim(2),startNodeIdx,endNodeIdx);

            % Update the Y Axis limits of state transition axes
            obj.StateAxes.YLim = [ylimit_lower,ylimit_higher];
        end

        function freqDDCallback(obj,eventData)
            %freqDDCallback Callback for the frequency drop down

            % Calculate Y axis limits and ticks based on the frequency limits
            objHandle = findobj(obj.PacketCommUIFigure,Tag="Freq Drop Down");
            if ~isempty(objHandle)
                eventVal = char(eventData.Value);

                % Calculate Y axis limits and ticks based on the frequency limits
                switch eventVal
                    case "2.4 GHz band"
                        freqLimits = [2400 2500];
                    case "5 GHz band"
                        freqLimits = [5000 6000];
                    case "6 GHz band"
                        freqLimits = [5900 7200];
                end

                % Update the frequency figure Y limits, Ticks, and Labels
                calcYTicks(obj,freqLimits)
            end
        end

        function saveSnipCallback(~,figHandle)
            %saveSnipCallback Save the UI figure as a image

            % Open dialog to get the file name and file format
            filter = {'*.jpg';'*.jpeg';'*.png';'*.tif';'*.tiff';'*.pdf'};
            [filename,filepath] = uiputfile(filter);

            % Once file name is entered, save the UI figure with the specified format
            if ischar(filename)
                exportapp(figHandle,[filepath filename]);
            end
        end

        function updateLine(obj,eventData,value)
            %updateLine Add a line handle to the figure

            % No action when simulation is running
            if obj.pIsSimInProgress
                return;
            end

            % Return of UI figure was closed
            if isempty(obj.PacketCommUIFigure) || ~isvalid(obj.PacketCommUIFigure)
                return;
            end

            % Get the value when clicked on any axes
            if isempty(value) && isa(eventData.Source.CurrentObject,"matlab.ui.control.UIAxes")
                % Delete the old line and draw at the new position clicked
                delete(obj.pTimeLine1)
                value = obj.StateAxes.CurrentPoint(1);
                obj.pTimeLine1 = xline(obj.StateAxes,value,"LineWidth",1,"Tag","Time line 1","Color",obj.pLineColor);
                if obj.FrequencyPlotFlag
                    delete(obj.pTimeLine2)
                    obj.pTimeLine2 = xline(obj.FreqAxes,value,"LineWidth",1,"Tag","Time line 2","Color",obj.pLineColor);
                end
            end
        end

        function updateInteractivityTime(obj,flag)
            %updateInteractivityTime Update interactivity of the controls for time axes

            % Enable or disable the UI controls
            handlesToUpdate = ["Time Min Label";"Time Max Label";...
                "Time Min Edit Field";"Time Max Edit Field"; ...
                "Time Min Slider";"Time Max Slider"];
            for idx = 1:numel(handlesToUpdate)
                objHandle = findobj(obj.PacketCommUIFigure,Tag=handlesToUpdate(idx));
                if ~isempty(objHandle)
                    objHandle.Enable = flag;
                end
            end
        end

        function updateInteractivityAxes(obj,flag)
            %updateInteractivityAxes Update interactivity of the controls for the
            %entire axes

            % Enable or disable the full view UI controls
            handlesToUpdate = ["Full View";"Export";"Freq Drop Down";"Node Slider"];
            for idx = 1:numel(handlesToUpdate)
                objHandle = findobj(obj.PacketCommUIFigure,Tag=handlesToUpdate(idx));
                if ~isempty(objHandle)
                    objHandle.Enable = flag;
                end
            end

            % Enable or disable the axes interactions
            if ~flag
                disableDefaultInteractivity(obj.StateAxes);
                if obj.FrequencyPlotFlag
                    disableDefaultInteractivity(obj.FreqAxes);
                end
            else
                enableDefaultInteractivity(obj.StateAxes)
                if obj.FrequencyPlotFlag
                    enableDefaultInteractivity(obj.FreqAxes);
                end
            end
            obj.StateAxes.Toolbar.Visible = flag;
            if obj.FrequencyPlotFlag
                obj.FreqAxes.Toolbar.Visible = flag;
            end
        end
    end
    %% Standard Based Callbacks for Events
    methods (Access=private)
        function wlanPlotCallback(obj,srcNode,eventData,nodes)
            %wlanPlotCallback Calculates the information necessary for live state
            %transition of WLAN node and plots the transition

            if ~obj.pIsInitialized
                % If visualization is not initialized, create the visualization
                obj.pIsInitialized = true;
                initializeVisualization(obj);
            end

            % Return of UI figure was closed
            if isempty(obj.PacketCommUIFigure) || ~isvalid(obj.PacketCommUIFigure)
                return;
            end

            % Get the WLAN node information
            if isa(srcNode,"helperCoexNode")
                sourceNode = srcNode.WLANDevice;
                sourceNode.ID = srcNode.ID;
                sourceNode.Name = srcNode.Name;
            else
                sourceNode = srcNode;
            end

            notificationData = eventData.Data;

            % Get the position in the WLAN nodes and in the Y tick list
            wlanNodeYTickIdx = find(sourceNode.ID==obj.pYTickOrder(:,1));
            nodeIDs = zeros(1,numel(nodes));
            for idx = 1:numel(nodes)
                nodeIDs(idx) = nodes{idx}.ID;
            end
            nodeIdx = find(nodeIDs==sourceNode.ID);
            deviceIdx = notificationData.DeviceID;
            if isa(nodes{nodeIdx},"helperCoexNode")
                nodes{nodeIdx} = nodes{nodeIdx}.WLANDevice;
                deviceIdx = notificationData.DeviceIndex;
            end

            if ~isscalar(wlanNodeYTickIdx)
                % If more than one Y tick is found, check based on device ID
                wlanDeviceYTickIdx = notificationData.DeviceID==obj.pYTickOrder(wlanNodeYTickIdx,2);
                wlanNodeYTickIdx = wlanNodeYTickIdx(wlanDeviceYTickIdx);
            end

            channelBandwidth = [];
            channelFrequencyStart = [];
            fillColor = [];
            freqPlot = false;
            % Update information based on the state of received information
            if isfield(notificationData,"State")
                switch notificationData.State
                    case "Contention"
                        fillColor = obj.pContendColor;
                        startTime = notificationData.CurrentTime-notificationData.Duration;
                        duration = notificationData.Duration;
                    case "Transmission"
                        fillColor = obj.pTxColor;
                        startTime = notificationData.CurrentTime;
                        duration = notificationData.Duration;
                        freqPlot = true;
                        centerFrequency = notificationData.Frequency/1e6;
                        channelBandwidth = notificationData.Bandwidth/1e6;
                        channelFrequencyStart = centerFrequency-(channelBandwidth/2);
                end

            else % MPDUDecoded event. This event is triggered after decoding of the received packet
                if all(notificationData.FCSFail) || notificationData.PHYDecodeFail
                    % Received packet failed decoding at PHY or during FCS check at MAC
                    fillColor = obj.pRxColorFailure;
                else % Successful decoding at MAC

                    % Find the address of the received packet
                    if ~isstruct(notificationData.MPDU)
                        % Find first MPDU with FCS pass
                        passIndex = find(~notificationData.FCSFail,1,"first");
                        hex_packet = dec2hex(notificationData.MPDU{passIndex});

                        % Get the MAC address of the received packet and address where the packet
                        % is destined to.
                        destinationMACAddress = reshape(hex_packet(5:10,:).',12,1).';
                        receivedMACAddress = nodes{nodeIdx}.MAC(deviceIdx).MACAddress;
                    else
                        % Get the MAC address of the received packet and address where the packet
                        % is destined to.
                        destinationMACAddress = notificationData.MPDU.Address1;
                        receivedMACAddress = nodes{nodeIdx}.MAC(deviceIdx).MACAddress;
                    end

                    % Set the state of the reception depending on the received address of the
                    % packet and the address the packet is destines to.
                    if strcmp(destinationMACAddress,receivedMACAddress)
                        % Expected MAC address is same, hence received packet is successful.
                        fillColor = obj.pRxColorUs;
                    elseif strcmp(destinationMACAddress,"FFFFFFFFFFFF") % Broadcast
                        % Initialize a flag indicating if received packet is destined to received
                        % node
                        isDestinedToUs = true;

                        % Trigger frames and Multi-STA BA are also broadcast frames. But they can
                        % be classified as destined to us or others. These frames are supported
                        % only in abstracted MAC, where the MPDU is a structure.
                        if isstruct(notificationData.MPDU)
                            senderMACAddress = notificationData.MPDU.Address2;
                            if strcmp(notificationData.MPDU.FrameType,"Trigger")
                                % Trigger frame is destined to a STA if it is sent by its associated AP and
                                % has the STA AID in any of the AID12 subfields of the user info fields.
                                isDestinedToUs = strcmp(senderMACAddress, nodes{nodeIdx}.MAC(deviceIdx).BSSID) && ...
                                    any(nodes{nodeIdx}.MAC(deviceIdx).AID == notificationData.MPDU.TriggerConfig.AID12);
                            elseif strcmp(notificationData.MPDU.FrameType,"Multi-STA-BA")
                                % Multi-STA BA frame is destined to a STA if it is sent by its associated
                                % AP and has the STA AID in any of the AID11 subfields of the per AID-TID
                                % info fields.
                                isDestinedToUs = strcmp(senderMACAddress, nodes{nodeIdx}.MAC(deviceIdx).BSSID) && ...
                                    any(nodes{nodeIdx}.MAC(deviceIdx).AID == notificationData.MPDU.AID);
                            end
                        end

                        % Check if the packet is destined to the received node. If true, then
                        % packet was successfully received. Else packet was successfully decoded
                        % but not destined to the node, hence Reception Others state.
                        if isDestinedToUs
                            fillColor = obj.pRxColorUs;
                        else
                            fillColor = obj.pRxColorOthers;
                        end
                    else
                        % Expected MAC address is different, hence received packet is successfully
                        % decoded but was destined for some other node.
                        fillColor = obj.pRxColorOthers;
                    end
                end

                % Packet reception event is triggered only at the end of packet reception
                % and decoding. Hence the duration of the packet received is current time -
                % PPDU start time.
                duration = notificationData.CurrentTime-notificationData.PPDUStartTime;

                % Get the start time of the PPDU
                startTime = notificationData.PPDUStartTime;

                % For channel occupancy plot
                freqPlot = true;
                centerFrequency = notificationData.Frequency/1e6;
                channelBandwidth = notificationData.Bandwidth/1e6;
                channelFrequencyStart = centerFrequency-(channelBandwidth/2);
            end

            if ~isempty(fillColor)
                % Calculate the Y tick position
                if ~isscalar(wlanNodeYTickIdx)
                    for idx = 1:numel(wlanNodeYTickIdx)
                        if obj.pYTickOrder(wlanNodeYTickIdx(idx),3)==obj.pCommunicationInfoWLAN{1}
                            wlanNodeYTickIdx = wlanNodeYTickIdx(idx);
                            break;
                        end
                    end
                end

                % Store the packet information
                storeCommunicationInfo(obj,obj.pCommunicationInfoWLAN,startTime,duration,fillColor, ...
                    wlanNodeYTickIdx,freqPlot,channelFrequencyStart,channelBandwidth)
            end
        end

        function devCfg = getDeviceConfig(~, node)
            %getDeviceConfig Returns the object holding MAC/PHY configuration

            if isa(node.DeviceConfig, 'wlanMultilinkDeviceConfig')
                devCfg = node.DeviceConfig.LinkConfig;
            else
                devCfg = node.DeviceConfig;
            end
        end

        function bluetoothPlotCallback(obj,srcNode,eventData)
            %bluetoothPlotCallback Calculates the information necessary for live state
            %transition of Bluetooth nodes and plots the transition

            if ~obj.pIsInitialized
                % If visualization is not initialized, create the visualization
                obj.pIsInitialized = true;
                initializeVisualization(obj);
            end

            % Return of UI figure was closed
            if isempty(obj.PacketCommUIFigure) || ~isvalid(obj.PacketCommUIFigure)
                return;
            end

            notificationData = eventData.Data;
            % Update information based on the state of received information
            if strcmp(eventData.EventName,"PacketTransmissionStarted")
                fillColor = obj.pTxColor;
                startTime = notificationData.CurrentTime;
                duration = notificationData.PacketDuration;
            elseif strcmp(eventData.EventName,"PacketReceptionEnded")
                startTime = notificationData.CurrentTime-notificationData.PacketDuration;
                if notificationData.SuccessStatus
                    fillColor = obj.pRxColorSuccess;
                else
                    fillColor = obj.pRxColorFailure;
                end
                duration = notificationData.PacketDuration;
            elseif strcmp(eventData.EventName,"ChannelAccessEnded")
                startTime = notificationData.CurrentTime-notificationData.Duration;
                fillColor = obj.pContendColor;
                duration = notificationData.Duration;
            end

            % Get the frequency for the channel number received
            if isa(srcNode,"bluetoothNode") || (isa(srcNode,"helperCoexNode") && any(strcmp(notificationData.PHYMode,["BR","EDR2M","EDR3M"])))
                commInfo = obj.pCommunicationInfoBREDR;
                channelBandwidth = 1;
                channelNumber = notificationData.ChannelIndex;
                channelFrequencyStart = 2402+channelNumber-0.5;
            else
                channelBandwidth = 2;
                channelFrequencyStart = [];
                if isa(srcNode,"helperBluetoothLE6GHzNode")
                    commInfo = obj.pCommunicationInfoLE6GHz;
                    if ~strcmp(eventData.EventName,"ChannelAccessEnded")
                        centerFrequency6000 = (5945:2:6423)*1e6;
                        channelNumber = notificationData.ChannelIndex;
                        channelFrequencyStart = ((centerFrequency6000(channelNumber+1))/1e6)-channelBandwidth/2;
                    else
                        channelBandwidth = [];
                    end
                else
                    commInfo = obj.pCommunicationInfoLE;
                    channelNumber = notificationData.ChannelIndex;
                    channelFrequencyStart = ((obj.pLECenterFrequencies(channelNumber+1))/1e6)-channelBandwidth/2;
                end
            end

            % Calculate the Y tick position
            bluetoothNodeYTickIdx = find(srcNode.ID==obj.pYTickOrder(:,1));
            if ~isscalar(bluetoothNodeYTickIdx)
                % If more than one Y tick is found, check based on device ID
                btDeviceYTickIdx = notificationData.DeviceID==obj.pYTickOrder(bluetoothNodeYTickIdx,2);
                bluetoothNodeYTickIdx = bluetoothNodeYTickIdx(btDeviceYTickIdx);
            end

            % Store the packet information
            storeCommunicationInfo(obj,commInfo,startTime,duration,fillColor, ...
                bluetoothNodeYTickIdx,~isempty(channelFrequencyStart),channelFrequencyStart,channelBandwidth)
        end

        function interferingPlotCallback(obj,srcNode,eventData)
            %interferingPlotCallback Calculates the information necessary
            %for live state transition of interfering WLAN or Bluetooth
            %nodes and plots the transition

            if ~obj.pIsInitialized
                % If visualization is not initialized, create the visualization
                obj.pIsInitialized = true;
                initializeVisualization(obj);
            end

            % Return of UI figure was closed
            if isempty(obj.PacketCommUIFigure) || ~isvalid(obj.PacketCommUIFigure)
                return;
            end

            notificationData = eventData.Data;

            % Update information based on the state of received information
            fillColor = obj.pTxColor;
            startTime = notificationData.CurrentTime;

            % Get the frequency for the channel number received
            if isa(srcNode,"helperInterferingWLANNode")
                commInfo = obj.pCommunicationInfoWLAN;
                channelBandwidth = srcNode.Bandwidth/1e6;
                channelFrequencyStart = (srcNode.CenterFrequency/1e6)-channelBandwidth/2;
                duration = notificationData.PacketDuration;
            elseif isa(srcNode,"helperInterferingBluetoothNode")
                commInfo = obj.pCommunicationInfoLE;
                channelBandwidth = srcNode.Bandwidth/1e6;
                channelFrequencyStart = (notificationData.CenterFrequency-notificationData.Bandwidth/2)/1e6;
                duration = notificationData.Duration;
            end

            % Calculate the Y tick position
            intNodeYTickIdx = find(srcNode.ID==obj.pYTickOrder(:,1));
            if ~isscalar(intNodeYTickIdx)
                % If more than one Y tick is found, check based on node type
                intDeviceYTickIdx = commInfo{1}==obj.pYTickOrder(intNodeYTickIdx,3);
                intNodeYTickIdx = intNodeYTickIdx(intDeviceYTickIdx);
            end

            % Store the packet information
            storeCommunicationInfo(obj,commInfo,startTime,duration,fillColor, ...
                intNodeYTickIdx,true,channelFrequencyStart,channelBandwidth)
        end

        function storeCommunicationInfo(obj,CommunicationInfoTech,startTime,duration,fillColor, ...
                yTickIdx,freqPlot,channelFrequencyStart,channelBandwidth)
            %storeCommunicationInfo Store packet information

            % Update the maximum end time and X limits
            if startTime+duration>obj.pMaxEndTime
                obj.pMaxEndTime = startTime+duration;
                if obj.pMaxEndTime>obj.pEditMaxHandle.Limits(2)
                    obj.pEditMaxHandle.Limits(2) = obj.pMaxEndTime;
                    obj.pMaxTimeSlider.Limits(2) = obj.pMaxEndTime;
                    obj.pXLimits(2) = obj.pMaxEndTime;
                end
            end

            % Update the packet information in the structure
            commInfo = obj.CommunicationInfoStruct;
            commInfo.Type = CommunicationInfoTech{1};
            commInfo.StartTime = startTime;
            commInfo.Duration = duration;
            commInfo.FillColor = fillColor;
            commInfo.YTickVal = yTickIdx;
            commInfo.FreqPlot = freqPlot;
            commInfo.FreqStart = channelFrequencyStart;
            commInfo.ChannelBandwidth = channelBandwidth;
            commInfo.FreqColor = CommunicationInfoTech{2};

            % Add the packet information to properties for plotting
            obj.InfoCount = obj.InfoCount+1;
            packetCount = obj.InfoCount;
            obj.StartTimeList(1,packetCount) = startTime;
            obj.EndTimeList(1,packetCount) = startTime+duration;
            obj.YTicksList(1,packetCount) = yTickIdx;
            obj.CommunicationInfo(1,packetCount) = commInfo;
        end
    end
end