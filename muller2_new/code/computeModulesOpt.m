function [modules, var, flux] = computeModulesOpt( model )
%COMPUTEMODULESOPT Computes the modules of the optimal flux space of this model
% Input:
%		model 			model to compute modules for (COBRA format)
%	
%	Output:
%		modules 		logical matrix with columns corresponding to modules and rows corresponding to reactions
%								each column has true for each reaction contained in the module and false otherwise
%								the matrix only contains the modules of reactions with variable flux rate
%		var					matrix with 2 columns. The first column gives minimal flux rates, the second gives maximal flux rates for each reaction.
%		flux				one optimal flux vector. This is useful for computing interface fluxes of modules.


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



% we use FBA to compute the optimal flux rate. Since we use the result of this computation as a constraint later on for FVA. We must solve it with high numerical precision. If we do not do this, it can happen that some of the LPs in the FVA run are cosidered infeasible due to numerical instabilities. Note however, that matlab is only working with doubles and hence even if we solve with very high precision, the output cannot store this precision. Hence, using a precision higher that 1e-20 will likely not yield any improvement. If your model is so numerically unstable that you need higher precision, you'll have to reimplement the module detection code in rational arithmetic (likely outside matlab).

model.precision = 1e-20; % only for method 1
model.dual_precision = 1e-20; % only for method 1

% we run FBA using the metaopt toolbox or using the Cobra toolbox.
% method 1: use the metaopt toolbox. This allows us to solve FBA with arbitrarily high precision. 
% method 2: use the cobra toolbox. Here we just have to hope that the default precision is good enough.

% Default method uses cobra toolbox, so that we don't have the dependency to the metaopt toolbox
method = 2;
if method == 1
	[opt, flux] = metaopt('fba', model);
elseif method == 2
	sol = optimizeCbModel(model);
	opt = sol.f;
	flux = sol.x;
else
	error('no such method');	
end

% check if the model has an objective function. In the current implementation it is only allowed to have one objective reaction.
% Adaptations to arbitrary objective functions are however easy to realize, but require the addition of artificial reactions.
obj = find(model.c);

assert(length(obj) == 1);
assert(model.c(obj) == 1);

% fix the objective reaction to its optimal value.
model.ub(obj) = opt;
model.lb(obj) = opt;

% we have now restriced the flux space to only contain optimal fluxes. Hence we can now compute the modules.
[modules, var] = computeModules(model);

end

