function [ ] = drawModule( model, module )
%DRAWMODULE draws a module using graphviz
% This function writes graphviz code to a file and then runs graphviz. Finally, the figure is shown in a pdf viewer.
%
% Input:
%   model  - an cobra model structure
%   module - logical or index array specifying reactions in the module


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




% don't show arcs for stoichiometric coefficients smaller than epsilon.
epsilon = 1e-4;

% we are only interested in the reactions in the module.
S = model.S(:,module);
if islogical(module)
    rxns = find(module);
else
    rxns = module;
end

% identify metabolites used in the module (and ignore the rest)
mets = find(sum((S ~= 0),2) ~= 0);
S = S(mets, :);

% write graphviz file
fid = fopen('module.gv','w'); % this is the file to which we will write the output
    if fid == -1
        error('unable to open file for writing graphviz output');
    end
    
    fprintf(fid, 'digraph G {\n');

    for i=1:length(rxns)
        for j=1:length(mets)
            
            if S(j,i) > epsilon
                if S(j,i) == 1
                    fprintf(fid,'R%d -> M%d\n', i, j);
                else
                    fprintf(fid,'R%d -> M%d [label="%2.2g"]\n', i, j, full(S(j,i)));
                end
            elseif S(j,i) < -epsilon
                if S(j,i) == -1
                    fprintf(fid,'M%d -> R%d\n', j, i);
                else
                    fprintf(fid,'M%d -> R%d [label="%2.2g"]\n', j, i, -full(S(j,i)));
                end
            end
        end
    end
    
    for i=1:length(rxns)
        fprintf(fid,'R%d [shape = box, label = "%s"];\n', i, model.rxns{rxns(i)});
    end
    for j=1:length(mets)
        fprintf(fid,'M%d [label = "%s"];\n', j, model.mets{mets(j)});
    end
    
    fprintf(fid, '}\n');
    fclose(fid);
    
    % filename of the output figure
    resultfile = sprintf('module.pdf');
    
    % run graphviz
    system(['dot -Tpdf module.gv -o' resultfile]);
    
    % show drawn figure
    system('evince module.pdf &');
end

