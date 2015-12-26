function dNdS = simdNdS(subS, subN, nTrials)
%input: subS {substitution:array of frequencies}, subN {substitution:array of freq}, 
%output: null distribution (dNdS) of dN/dS
    dNdS = [];
    
    for i=1:nTrials
        dS = 0;
        dN = 0;
        for c = intersect(subS.keys,subN.keys)
            s = char(c);
            freqpool = [subS(s) subN(s)];
            if (~strcmp(s,'CA') && ~strcmp(s,'GT') && ~strcmp(s,'TC') && ~strcmp(s,'GA'))
                dS = dS+filteredMean(randsample(freqpool,length(subS(s)),true));
                dN = dN+filteredMean(randsample(freqpool,length(subN(s)),true));
            end
        end
        dNdS = [dNdS dN/dS];
    end
    
end

