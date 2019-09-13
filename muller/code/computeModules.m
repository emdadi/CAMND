function [modules, var] = computeModules( model )
%COMPUTEMODULES Computes the modules of this model
% Input:
%		model 			model to compute modules for (COBRA format)
%	
%	Output:
%		modules 		logical matrix with columns corresponding to modules and rows corresponding to reactions
%								each column has true for each reaction contained in the module and false otherwise
%								the matrix only contains the modules of reactions with variable flux rate
%		var					matrix with 2 columns. The first column gives minimal flux rates, the second gives maximal flux rates for each reaction


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




% This is ugly (but we have to choose some error tolerances)!
% If the algorithm does not work for you, try to play with these values. Numerical instabilities can cause weird results.
epsilon = 1e-7;
model.precision = 1e-8; % only important for method 1,2
model.dual_precision = 1e-8; % only important for method 1,2

% run FVA to find fixed reactions
% Choose the method that you want to use. 

% method 1,2: 	use the fva feature of the metaopt toolbox. 
%								This allows us to solve FVA with higher precision and faster than the cobra toolbox.
%								The metaopt toolbox is however rather a pain to install.
% method 2:			uses a variant of the metaopt toolbox that sometimes can circumvent some numerical troubles
%
% method 3:			use the cobra toolbox for FVA. This method is slower and the numerical precision cannot be tuned.

%By default we use method 3, because it only needs the Cobra Toolbox
method = 3;

if method == 1
	% method 1
	var = metaopt('fva', model);
elseif method == 2
	% method 2
	model.int_rxns = {};
	var = metaopt('tfva', model);
elseif method == 3
	% method 3
	[minvar, maxvar] = fluxVariability(model, 0, 'max', model.rxns, true);
	var = [minvar, maxvar];
else
	error('not a valid method');
end

% find fixed reactions (a small amount of false positives is ok, since these usually turn out to be modules containing a single reaction)
V = var(:,2) - var(:,1) > epsilon;

% restrict only to reactions with variable flux rate
S = model.S(:,V);

% compute fundamental circuits
K = findCircuitBase(S);
% compute connected components of the fundamental circuits
m = logical(findConnectedComponents(K));

% S did not contain all the reactions, so make sure we get the indices right.
modules = false(size(model.S,2), size(m,2));
modules(V,:) = m;

end

