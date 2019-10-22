function [feat_ilev] = CONVERT2ILEV(feat_lev, slicePts, altPts)
% converts TIEGCM matrix on the midpoint levels to matrix on interface
% levels

% feat_lev = matrix containing features on levs with size of [slicePts, altPts]
% slicePts, can be either numbers lat points or lon points. whatever you want
% altPts =  number of altitude points
    
feat_ilev = zeros( size(feat_lev) );

% find average between elements to get on the interfaces
    if slicePts ~= 0
        for alt = 1:altPts
            if alt ~= 1           
                for slc = 1:slicePts
                    if slc ~= 1           
                        feat_ilev(slc,alt) = 0.5*(feat_lev(slc,alt) + feat_lev(slc-1, alt-1));
                    else       
                        feat_ilev(slc,alt) = feat_lev(slc,alt);        
                    end
                end

            else
                for slc = 1:slicePts
                    if slc ~= 1   
                        feat_ilev(slc,alt) = 0.5*(feat_lev(slc,alt) + feat_lev(slc-1, alt));
                    else       
                        feat_ilev(slc,alt) = feat_lev(slc,alt);          
                    end
                end
            end
        end

    else    % only converting an array
        for alt = 1:altPts
            if alt ~= 1           
                feat_ilev(alt) = 0.5*(feat_lev(alt) + feat_lev(alt-1));
            else
                feat_ilev(alt) = feat_lev(alt);
            end
        end
    end

    
end



