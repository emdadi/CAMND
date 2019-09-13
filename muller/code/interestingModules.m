function [ vmods, fmods ] = interestingModules( model, flux, modules, nutrient, notshow )
%PARTITIONMODULES find the interesting modules of a given network
%   This method filters out the interesting modules of reactions with
%   variable flux rates, i.e. modules without any net-production are dropped 
%   (they correspond to thermodynamically infeasible cycles or are
%   uncoupled to the whole system).
%   Additionally, the set of fixed reactions is grouped into modules. After
%   not-show metabolites and all reactions with variable flux rate are
%   removed from the system, the remaining network may decompose into
%   several connected components. These connected components are then
%   identified as the fixed modules.
%
%   Input:
%    model    - the model (COBRA format)
%    flux     - a feasible flux through the network
%    modules  - a decomposition of the network (variable reactions) into modules
%    nutrient - a list of reaction IDs that should be explicitely depcited
%    notshow  - a list of metabolite IDs of metabolites that should be
%                   hidden (e.g. h2o or h)
%   Output:
%    vmods    - selected modules with reactions of variable flux rates
%    fmods    - selected modules with reactions of fixed flux rates


% Copyright (C) 2013  Arne Müller <arne.mueller@fu-berlin.de>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.




epsilon = 1e-5;

keep_module = false(1,size(modules,2));
for i=1:size(modules,2)
    keep_module(i) = norm(model.S(:,modules(:,i)) * flux(modules(:,i)),'inf') > epsilon;
end

% compute set of interesting reactions R that we want to visualize in modules
% etc.
R = logical(modules(:,keep_module) * ones(sum(keep_module),1));

bio = find(model.c);
assert(length(bio) == 1);

R(bio) = true;
if ~isempty(nutrient)
    R(nutrient) = true;
end

% compute set M of metabolites that should be shown
% these are metabolites that interface with the interesting reactions
M = logical(model.S(:,R) ~= 0) * ones(sum(R),1);
M = logical(M);
M(notshow) = false;

% identify blocked reactions (will be removed)
F = flux > epsilon | flux < -epsilon;

% find connected components of non-interesting reactions
C = findConnectedComponents(model.S(~M,F & ~R)');

% generate a module for each component of non-interesting reactions
% (these connect the interesting modules)
fmodules = false(size(model.S,2), size(C,2));
fmodules(F & ~R,:) = C;

% generate modules for nutrients
nutrient_mods = false(length(model.rxns), length(nutrient));
for i=1:length(nutrient)
    nutrient_mods(nutrient(i),i) = true;
end

% generate modules for biomass (only 1 modules)
biomass_mods = false(length(model.rxns), 1);
biomass_mods(bio) = true;

% generate output matrices
vmods = modules(:,keep_module);
fmods = [nutrient_mods, fmodules, biomass_mods];

end

