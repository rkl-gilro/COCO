function [dkl] = RGB2DKL(rgb,monxyY)
% ---- Inputs ----
% rgb:              Nx3 matrix of RGB values
%                   [R G B]
%                   values range from 0 to 1
% monxyY:           3x3 matrix of xyY values of R, G, B channels at maximum
%                   [[ Rx  Ry  RY  ]
%                    [ Gx  Gy  GY  ]
%                    [ Bx  By  BY  ]]
% 
% ---- Outputs ----
% dkl:              Nx3 matrix of DKL values
%                   [Lum LM S]
%                   along axes, range within monitor gamut is -.5 to .5
% 
% 
% Derived from deprecated function toDKL:
% -  no longer uses global variables
% -  uses updated initmon file (initmon_ver2)
% 
% 10/2023: changes made by LH


[RGBCOL COLRGB] = initmon_ver2(monxyY);

dkl = (rgb - 0.5) * RGBCOL';
% % change S axis - MT 01/12/2015
% dkl(:,3)=-dkl(:,3);

