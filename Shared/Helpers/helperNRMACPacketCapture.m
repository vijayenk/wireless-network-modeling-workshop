classdef helperNRMACPacketCapture < handle
    %helperNRMACPacketCapture Captures the NR MAC packets
    %
    %   helperNRMACPacketCapture methods:
    %
    %   close - Close the packet writer

    %   Copyright 2024-2025 The MathWorks, Inc.

    properties(SetAccess=private)
        %PCAPObj Holds the nrPCAPWriter object
        PCAPObj
    end

    methods
        function obj = helperNRMACPacketCapture(nodes,fileName)
            %helperNRMACPacketCapture Captures the NR MAC packets
            %
            % nodes - List of nodes of type nrGNB or nrUE
            %
            % fileName - Name of the packet capture (PCAP) file

            % Create an nrPCAPWriter object to capture packets into a PCAP file
            obj.PCAPObj = nrPCAPWriter(FileName=fileName, FileExtension='pcap');
            for idx=1:numel(nodes)
                addlistener(nodes(idx), 'PacketTransmissionStarted', @(src, eventData) logPackets(obj,eventData.Data));
                addlistener(nodes(idx), 'PacketReceptionEnded', @(src, eventData) logPackets(obj,eventData.Data));
            end
        end

        function close(obj)
            %close Close the packet writer

            delete(obj.PCAPObj);
        end
    end

    methods(Hidden)
        function logPackets(obj, data)
            %logPackets Process the event data and write it into a PCAP file
            %
            % logPackets(obj, data) Process the event data and write it into a PCAP file
            %
            % OBJ is a helperNRMACPacketCapture object
            %
            % DATA is a structure and contains information about the packet
            % transmission or reception event

            if data.SignalType ~= "CSIRS"  && data.SignalType ~= "SRS"
                % Determine the radio type (duplex mode)
                if data.DuplexMode =="TDD"
                    radioType = obj.PCAPObj.RadioTDD;
                else
                    radioType = obj.PCAPObj.RadioFDD;
                end

                % Determine the link type
                if data.LinkType == "UL"
                    linkType = obj.PCAPObj.Uplink;
                else
                    linkType = obj.PCAPObj.Downlink;
                end

                % Prepare the Meta data
                timingInfo = data.TimingInfo;
                packetMetaData = struct('RadioType', radioType, 'RNTIType', obj.PCAPObj.CellRNTI, 'RNTI', data.RNTI, ...
                    'HARQID', data.HARQID, 'SystemFrameNumber', mod(timingInfo(1), 1024), 'SlotNumber', timingInfo(2), 'LinkDir', linkType);

                % Log the packet based on the event time in seconds and
                % convert it to the microseconds
                write(obj.PCAPObj, data.PDU, round(data.CurrentTime*1e6), PacketInfo=packetMetaData);
            end
        end
    end
end