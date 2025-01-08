function [rgb maxval] = DKL2RGB(dkl,monxyY,shrink)
% ---- Inputs ----
% dkl:              Nx3 matrix of DKL values
%                   [Lum LM S]
%                   along axes, range within monitor gamut is -.5 to .5
% monxyY:           3x3 matrix of xyY values of R, G, B channels at maximum
%                   [[ Rx  Ry  RY  ]
%                    [ Gx  Gy  GY  ]
%                    [ Bx  By  BY  ]]
% shrink:           0 or 1
%                   1 means to shrink the RGB values to fit the 0-1 gamut
%                       (outputs corrected RGB to console)
%                   0 means resulting RGB values are kept as is
% -- Outputs --
% rgb:              Nx3 matrix of RGB values
%                   [R G B]
%                   values range from 0 to 1
% 
% 
% Derived from deprecated function fromDKL:
% -  no longer uses global variables
% -  uses updated initmon file (initmon_ver2)
% 
% 10/2023: changes made by LH


if nargin <3
    shrink=0;
end

% dkl(:,3)=-dkl(:,3);% change S axis - MT 01/12/2015
[RGBCOL COLRGB] = initmon_ver2(monxyY);
rgb = COLRGB * dkl';


maxval = max(max(abs(rgb)));
if maxval > 0.5
    if shrink==1
        if maxval < .5 + .000001
            fprintf('very small epsilon\n')
        end
        rgb = rgb/(2*maxval);
        fprintf('\tCORRECTED rgb: %g %g %g\n', rgb(1)+.5, rgb(2)+.5, rgb(3)+.5);
        fprintf(2,'\tALL inputted RGBs have been scaled, not just out-of-gamut RGBs.\n')
        % warning: note that this corrects ALL the rgbs, not just the out-of-gamut one
    else
        warning('One or more values are outside the gamut.')
    end
end


rgb = rgb'+0.5;

