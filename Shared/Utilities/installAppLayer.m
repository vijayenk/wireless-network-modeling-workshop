function installAppLayer(networkSimulator, appDataRate)

gNB = networkSimulator.Nodes{1};
UE = networkSimulator.Nodes{2};

% Install the DL application traffic on gNB for the UE node
dlApp = networkTrafficOnOff(GeneratePacket=true,OnTime=Inf,OffTime=0,DataRate=appDataRate);
addTrafficSource(gNB,dlApp,DestinationNode=UE);

% Install the UL application traffic on the UE node for the gNB
ulApp = networkTrafficOnOff(GeneratePacket=true,OnTime=Inf,OffTime=0,DataRate=appDataRate);
addTrafficSource(UE,ulApp)
end