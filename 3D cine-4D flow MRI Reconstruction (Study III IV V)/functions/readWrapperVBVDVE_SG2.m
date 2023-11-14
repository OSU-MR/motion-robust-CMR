function [data, samp, param, cMaps, rawData, ref] = readWrapperVBVDVE_SG2(p, version)

switch version
    case 'VD'
        [data, samp, param, cMaps, rawData, ref] = readWrapper3(p);
        
    case 'VE'
        [data, samp, param, cMaps, rawData] = readWrapper4(p);
        rawData = rawData{2};
        param.TR = param.TR * 1000;
        ref = [];

end

