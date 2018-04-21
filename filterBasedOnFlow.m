function filtered = filterBasedOnFlow(labeled, flowMagnitude, maxFlow)
% filtered = FILTERBASEDONFLOW(labeled, flowMagnitude, maxFlow) filters the
% provided labeled image, one that can be retrieved using BlobAnalysis, and
% filters it using FLOWMAGNITUDE and maxFlow. LABELED and FLOWMAGNITUDE
% must have the same size and MAXFLOW must be a scalar.
%
% RETURNS a labeled image, FILTERED, that contains only the blobs that (1)
% have a flow, i.e. flow>0, and a flow<maxFlow. FILTERED is of the same
% size as LABELED.
    temp = logical(labeled).*flowMagnitude;
    temp2 = temp>0;
    temp3 = temp<maxFlow;
    filtered = labeled.*uint8(temp3).*uint8(temp2);
end