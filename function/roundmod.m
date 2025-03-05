% SYNTAX:
%   [rx] = roundmod(x,y);
%
% INPUT:
%   x = values to be rounded
%   y = resolution
%
% OUTPUT:
%   rx = rounded values
%
% DESCRIPTION:
%   Rounds the input values to the nearest float determined by the
%   resolution in input.
%% ------------------------------------------------------------------------
function [rx] = roundmod(x,y)
rx = round(x./y).*y;

