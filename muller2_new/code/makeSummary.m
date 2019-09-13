function [ ] = makeSummary( model, notshow, nutrient, name )
%MAKESUMMARY create a folder with the big network, the modules (except big
%fixed ones) and a latex file giving extra information
%
% This method computes the modules, generates the figures and creates a summary file.
%
%	Input:
%		model			model (cobra format) for which to run the analysis
%		notshow		indices of metabolites that should not be shown in the interaction figure (e.g. h20)
%		nutrient	indices of nutrients (reactions) that should be shown explicitely
%		name			name of the network (specifies the name of the folder where the results are stored)


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




[modules, var, flux] = computeModulesOpt( model );

[modv, modf] = interestingModules( model, flux, modules, nutrient, notshow );

[metG, numMods] = visualizeModsSingle( model, notshow, modf, modv, flux, nutrient );

% the directory where the results will be stored. YOU MAY WANT TO MODIFY THIS
outdir = ['../../results/decomposition/matroids/' name '/viz/'];
mkdir(outdir);

% create latex document that gives a summary and list metabolite groups

    fid = fopen('network_summary.tex','w'); % this is the file to which we will write the output
    if fid == -1
        error('unable to open file for writing latex document');
    end
    
    fprintf(fid, '\\documentclass[a4paper,12pt]{article}\n\n');
    fprintf(fid, '\\usepackage{geometry}\n');
    fprintf(fid, '\\usepackage{fancyhdr}\n');
    fprintf(fid, '\\geometry{left=2cm, top=3cm, right=2cm, bottom=2.5cm}\n');
    fprintf(fid, '\\pagestyle{fancy}\n');
    
    fprintf(fid, '\\title{Summary Information on %s}\n', strrep(name,'_',' '));
    fprintf(fid, '\\author{Arne C. M\\"uller, Frank J. Bruggeman, Brett G. Olivier, Leen Stougie}\n');
    
    fprintf(fid, '\\begin{document}\n');
    fprintf(fid, '\\maketitle\n');
    fprintf(fid, '\\section*{Statistics}\n');
    fprintf(fid, '\\begin{description}\n');
    fprintf(fid, '\\item[No. Metabolites:] %d \n', size(metG,2));
    fprintf(fid, '\\item[No. Modules:] %d \n', full(numMods));
    fprintf(fid, '\\end{description}\n');
    
    fprintf(fid, '\\section*{Metabolite Groups}\n');
    for i=1:size(metG,2)
        if sum(metG(:,i) ~= 0) > 1
            fprintf(fid, '\\paragraph{Metabolite Group M%d:}\n', i);
            fprintf(fid, '\\begin{itemize}\n');
            for j=1:size(metG,1)
                if metG(j,i) ~= 0
                    fprintf(fid, '\\item $%2.4f * \\mathtt{%s}$ (%s) \n',metG(j,i), model.mets{j}, model.metNames{j});
                end
            end
            fprintf(fid, '\\end{itemize}\n');
        end
    end

    fprintf(fid, '\\end{document}\n');
    fclose(fid);
    
    system('pdflatex network_summary.tex');
    system(['mv network_summary.pdf ' outdir]);
    system(['mv network.pdf ' outdir]);
    
    % also create figures for modules and move them to target directory
    
    for i=1:size(modv,2)
        if sum(modv(:,i)) < 100
            drawModule( model, modv(:,i) );
            tName = sprintf('V%d.pdf',i);
            system(['mv module.pdf ' outdir tName]);
        end
    end
    
    for i=1:size(modf,2)
        if sum(modf(:,i)) < 100 && sum(modf(:,i)) > 1
            drawModule( model, modf(:,i) );
            tName = sprintf('F%d.pdf',i);
            system(['mv module.pdf ' outdir tName]);
        end
    end
end

