function mac=MAC2(phi)
% This function is to compute and plot Modal Assurance Criterion (MAC) matrix between identified mode shapes
% rectangle around the peaks.
% Programmer: Mohammad Farshchin, Ph.D candidate at The UofM
% Email: Mohammad.Farshchin@gmail.com
% Last modified: 9/23/2015
% phi: matrix of the identified mode shapes
% mac: MAC matrix
% Example: load ModeShapes.mat;MAC_Matrix=MAC(phi)
    for I=1:size(phi,2)
        for J=1:size(phi,2)
            mac(I,J)=Mac(phi(:,I),phi(:,J));
        end
    end
end

