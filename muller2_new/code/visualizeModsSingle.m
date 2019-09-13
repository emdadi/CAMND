function [metG, numMods, numMets] = visualizeModsSingle( model, notshow, fmods, vmods, flux, nutrient )
% creates graphviz output for vizualizing the interaction of modules.
% The graphviz output is written in a file.
%
% Input: 
% model     - metabolic model to visualize
% notshow   - metabolite indicies of metabolites that should not be shown (e.g. h2o)
% fmods     - matrix whose columns contain modules with only fixed
%               reactions
% vmods     - matrix whose columns contain modules with variable reactions
% flux      - a feasible flux
% nutrient  - indices of nutrient input reactions to be shown separately
%
% Output:
%		metG		- metabolites that were grouped together
%		numMods - number of modules
%		numMets	- number of metabolites


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




% we don't draw arcs for fluxes with a flow value less than 1e-4
epsilon = 1e-4;

M = true(1,length(model.mets));
M(notshow) = false;

num_mods = 0;
num_mods = num_mods + size(vmods,2);

  vrxns = false(1,length(model.rxns));
  for i=1:size(vmods,2);
    vrxns(vmods(:,i)) = true;
  end


% build stoichiometric matrix of modules
%Spar = sparse(sum(M),0);
% also build matrix that encodes all the metabolites that are produced /
% consumed in a module but are in steady-state inside the module
Si = sparse(sum(M), size(fmods,2) + num_mods);
j = 0;
for i=1:size(fmods,2)
    j = j+1;
    if any(fmods(vrxns,i))
        % modules is not used, because there is a variable flux module
        % so just don't initialize
    else
        Si(:,j) = model.S(M, fmods(:,i)) * flux(fmods(:,i));
    end
end
for k=1:size(vmods,2)
    j = j+1;
        Si(:,j) = model.S(M, vmods(:,k)) * flux(vmods(:,k));
end


% also build matrix that encodes all the metabolites that are produced /
% consumed in a module but are in steady-state inside the module
Shidden = sparse(sum(M),size(fmods,2) + num_mods);
for i=1:size(fmods,2)
    Shidden(:,i) = any(model.S(M, fmods(:,i)) ~= 0,2);
end
j = size(fmods,2);
for k=1:size(vmods,2)
    j = j+1;
    Shidden(:,j) = any(model.S(M, vmods(:,k)) ~= 0,2);
end

Shidden(abs(Si) > epsilon) = false;

Mfind = find(M);
model.mets(Mfind(sum(Shidden,2) >= 2))

for i=1:size(fmods,2)
    fprintf('hidden %d: \n',i);
    disp(model.mets(Mfind(Shidden(:,i) & any(Si,2))))
end

% misuse this old method to find "parallel" metabolites
[fcv, lambda] = findParallelReactions(Si', epsilon);

metgroups = false(0,sum(M));
metG = zeros(length(model.mets), 0);

j = 0;
for i=1:length(fcv)
    if fcv(i) == i
        % by the construction of fcv, we know that the other parallel
        % metabolites only have indices greater than i
        parallel = (fcv == i);
        assert(all(find(parallel) >= i));
        % in contrast to the coupled method, we can just take the
        % stoichiometry of one of the metabolites.
        % we want to take the metabolite with biggest stoichiometry
        % since its parallel over Spar, we can just look at one reaction
        % that contains i
        
        parI = find(parallel);
        
        test_rxn = find(ones(1,sum(parallel))*(abs(Si(parallel,:)) > epsilon*10) >= 1);
        
        if isempty(test_rxn)
            % we may have missed it because our epsilon was too big
            test_rxn = find(ones(1,sum(parallel))*(abs(Si(parallel,:)) > epsilon) >= 1);
        end
        
        if ~isempty(test_rxn)
            % only create metabolite nodes for metabolites that are really
            % used

            [valselect, iselect] = max(abs(Si(parallel,test_rxn(1))));

            assert(~isempty(iselect));

            j = j+1;
            metreps(j) = parI(iselect);
            metgroups(j,parallel) = true;
            
            lambdaj2 = zeros(sum(M),1);
            lambdaj2(parallel) = lambda(parI(iselect)) ./ lambda(parallel);
            metG(M,j) = lambdaj2;
        end
    end
end

    S = Si(metreps,:);
            
    Spos = false(length(metreps), size(S,2));
    Sneg = false(length(metreps), size(S,2));
    
    Spos = S > epsilon;
    Sneg = S < -epsilon;
    
    used_metrep = any(Spos | Sneg,2);
    used_mods = any(Spos | Sneg, 1);
    
        
    fid = fopen('network.gv','w'); % this is the file to which we will write the output
    if fid == -1
        error('unable to open file for writing graphviz output');
    end
    
    fprintf(fid, 'digraph G {\n');
    
    for i=1:size(fmods,2)
        for j=1:size(S,1)
            metname = sprintf('M%d',j);
            rxnname = sprintf('F%d',i);
            
		        if Sneg(j,i)
	            	fprintf(fid,'%s -> %s [label="%2.2f"]\n',metname, rxnname, -full(S(j,i)));
		        end
		        if Spos(j,i)
	            	fprintf(fid,'%s -> %s [label="%2.2f"]\n',rxnname, metname, full(S(j,i)));
		        end
        end
    end
    
     for k=1:size(vmods,2)
        i = k + size(fmods,2);
        for j=1:size(S,1)
            metname = sprintf('M%d',j);
            rxnname = sprintf('V%d',k);
            
		        if Sneg(j,i)
	            	fprintf(fid,'%s -> %s [label="%2.2f"]\n',metname, rxnname, -full(S(j,i)));
		        end
		        if Spos(j,i)
	            	fprintf(fid,'%s -> %s [label="%2.2f"]\n',rxnname, metname, full(S(j,i)));
		        end
        end
    end

		

    
    
    % also build matrix that encodes all the metabolites that are produced /
    % consumed in a module but are in steady-state inside the module
    Shidden = sparse(sum(M),size(fmods,2) + num_mods);
    for i=1:size(fmods,2)
        Shidden(:,i) = any(abs(model.S(M, fmods(:,i))) > epsilon,2);
    end
    j = size(fmods,2);

    for k=1:size(vmods,2)
        j = j+1;fprintf(fid,'%s [shape = box];\n', rxnname);
        Shidden(:,j) = any(abs(model.S(M, vmods(:,k))) > epsilon,2);
    end
 
    Shidden(abs(Si) > epsilon) = false;
    Shidden(abs(Si) > epsilon) = false;
    
    Mfind = find(M);
    hiddenmets = find(sum(Shidden,2) >= 2);
    
    model.mets(Mfind(hiddenmets))

    for i=1:size(fmods,2)
        fprintf('hidden F%d: \n',i);        
        disp(model.mets(Mfind(Shidden(:,i) & any(abs(Si) > epsilon,2))))
        disp(model.mets(Mfind(hiddenmets(Shidden(hiddenmets,i) ~= 0))))
    end

    for i=1:size(vmods,2)
        j = size(fmods,2)+i;
        fprintf('hidden V%d: \n',i);
        disp(model.mets(Mfind(Shidden(:,j) & any(abs(Si) > epsilon,2))))
        disp(model.mets(Mfind(hiddenmets(Shidden(hiddenmets,j) ~= 0))))
    end
    
    % TODO: draw hidden interactions
    
    mets_show = model.mets(M);
    for j=1:size(S,1)
        metname = sprintf('M%d',j);
        if sum(metgroups(j,:)) == 1
        	label = mets_show{metreps(j)};
        	fprintf(fid,'%s [label = "%s"];\n', metname, label);
        end
    end
    for i=1:size(fmods,2)
        rxnname = sprintf('F%d',i);
        if sum(fmods(:,i)) == 1
        	label = model.rxns{fmods(:,i)};
        	fprintf(fid,'%s [shape = box, label = "%s"];\n', rxnname, label);
        else
            fprintf(fid,'%s [shape = box];\n', rxnname);
        end
    end
    
    for i=1:size(vmods,2)
        j = size(fmods,2)+i;
        rxnname = sprintf('V%d',i);
        if used_mods(j)
            fprintf(fid,'%s [shape = box];\n', rxnname);
        end
    end

    bio = find(model.c);
		if length(bio) > 0
			for i=1:size(fmods,2)
				if fmods(bio,i)
					fprintf(fid,'{rank = sink; F%d}\n',i);
				end
				if any(fmods(nutrient, i))
					fprintf(fid,'{rank = source; F%d}\n',i);
				end
			end
		end
    
    fprintf(fid, '}\n');
    fclose(fid);
    
    resultfile = sprintf('network.pdf');

		% run graphviz to create figure (pdf)    
    system(['dot -Tpdf network.gv -o' resultfile]);
    %system('evince network_poster.pdf &');
    
    numMets = sum(used_metrep);
    numMods = sum(used_mods);
end



% detect columns with proportional stoichiometry. We use the flux coupling code from F2C2 (Larhlimi et al., F2C2: a fast tool for the computation of flux coupling in genome-scale metabolic networks, BMC Bioinformatics 2012, 13:57  doi:10.1186/1471-2105-13-57)
function [fcv, lambda] = findParallelReactions(S, tol)
    S(abs(S)<=tol) = 0;
    Sl= logical(S);
    [n m] = size(S);
    fcv = zeros(m, 1);
    lambda = zeros(m, 1);
    for i=1:m
        if fcv(i)~=0
            continue;
        end
        if any(Sl(:,i))
            % reaction is not blocked
            % reaction is trivially coupled to itself
            fcv(i) = i;
            lambda(i) = 1;
            % detect other reactions its coupled to
            if i < m
                for j=(i+1):m
                    if nnz(Sl(:,i)-Sl(:,j))==0
                        v = S(Sl(:,i),i)./S(Sl(:,j),j);
                        if isempty(v)
                            continue
                        end
                        if all(abs(v-v(1))<=tol)
                            fcv(j) = i;
                            lambda(j) = v(1);
                        end
                    end
                end
            end
        else
            % reaction is blocked
            fcv(i) = 0;
            lambda(i) = 0;
        end
    end
end
