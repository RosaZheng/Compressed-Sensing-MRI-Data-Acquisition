% Generates a coefficient mask of a given size using the given parameters
% The mask will be origin mirrored on the left half if symmetric = 1
% If incoherent, does Monte Carlo algorithm with masktrials for minimum
% peak interference

function [mask npdf] = cs_generatemask(imgsize, ratio, uniformity, coherence, samplingtype, symmetry, density, masktrials)
% imgsize, 
% ratio, 
% uniformity, 'Uniform'
% coherence, string: 'Incoherent', 'Coherent'
% samplingtype, 'Line' or 'Point' (not setup yet)
% symmetry, 'Symmetric' or 'Asymmetric'
% density, 
% masktrials
%

    row = imgsize(1);
    col = imgsize(2);
    
    
    if(mod(col,2)==0) 
        startcol = (col / 2) + 1;
    else
        startcol = (col + 1) / 2;
    end
    if(mod(row,2)==0) 
        dcrow = (row / 2) + 1;
    else
        dcrow = (row + 1) / 2;
    end
    
    % can be replaced by 
    startcol=ceil((col+1)/2);
    dcrow = ceil((row+1)/2);
    
    npdf = ones(row,col);
    mininterference = 1e99;
    
    if(~strcmp(coherence,'Incoherent'))
       masktrials = 1; %force to 1 if 'Incoherent'
    end
    % generate a new random mask each time
    for q = 1:masktrials
    
        if(strcmp(uniformity,'Uniform'))

            if(strcmp(coherence,'Incoherent'))

                if(strcmp(samplingtype,'Line'))

                    if(strcmp(symmetry,'Symmetric')||strcmp(symmetry,'symmetry') 

                        % Uniform, incoherent, symmetric, line-based sampling
                        n = ceil(col/2);
                        m = ceil(col*ratio/2);
                        perm = randperm(n)-1+startcol; % uniform permutation
                        perm = perm(1:m);
                        if(sum(find(perm==startcol))==0) % keep DC line
                            perm(1) = startcol;
                        end
                        mask = zeros(row,col);
                        mask(:,perm) = 1;
                        mask = cs_restoresymmetry(mask);

                    else

                        % Uniform, incoherent, asymmetric, line-based sampling
                        n = ceil(col);
                        m = ceil(col*ratio);
                        perm = randperm(n);
                        perm = perm(1:m);
                        if(sum(find(perm==startcol))==0) % keep DC line
                            perm(1) = startcol;
                        end
                        mask = zeros(row,col);
                        mask(:,perm) = 1;

                    end

                else

                    if(strcmp(symmetry,'Symmetric'))

                        % Uniform, incoherent, symmetric point-based sampling
                        m = ceil(row*col*ratio/2);
                        perm = randperm(((col-startcol+1)*row));
                        perm = perm(1:m);
                        permcol = floor((perm-1)/row);
                        permrow = perm - permcol*row;
                        permcol = permcol + startcol;
                        mask = zeros(row,col);
                        for i=1:length(perm)
                            mask(permrow(i),permcol(i)) = 1;
                        end
                        mask(dcrow,startcol) = 1; % keep DC point
                        mask = cs_restoresymmetry(mask);

                    else

                        % Uniform, incoherent, asymmetric point-based sampling
                        m = ceil(row*col*ratio);
                        perm = randperm(col*row);
                        perm = perm(1:m);
                        permcol = floor((perm-1)/row)+1;
                        permrow = perm - permcol*row;
                        mask = zeros(row,col);
                        for i=1:length(perm)
                            mask(permrow(i),permcol(i)) = 1;
                        end
                        mask(dcrow,startcol) = 1; % keep DC point

                    end

                end

            else

                if(strcmp(samplingtype,'Line'))

                    if(strcmp(symmetry,'Symmetric'))

                        % Uniform, coherent, symmetric line-based sampling
                        n = ceil(col/2);
                        m = ceil(col*ratio/2);
                        colpos = startcol:(n/m):col;
                        mask = zeros(row,col);
                        for i=1:length(colpos)
                            mask(:,round(colpos(i))) = 1;
                        end
                        mask = cs_restoresymmetry(mask);

                    else

                        % Uniform, coherent, asymmetric line-based sampling
                        n = ceil(col);
                        m = ceil(col*ratio);
                        colpos = 1:(n/m):col;
                        mask = zeros(row,col);
                        for i=1:length(colpos)
                            mask(:,round(colpos(i))) = 1;
                        end

                    end

                else

                    if(strcmp(symmetry,'Symmetric'))

                        % Uniform, coherent, symmetric point-based sampling
                        n = row*ceil(col/2);
                        m = ceil(row*col*ratio/2);
                        spacing = sqrt(n/m);
                        rowpos = 1:spacing:row;
                        colpos = startcol:spacing:col;
                        mask = zeros(row,col);
                        for a=1:length(rowpos)
                            for b=1:length(colpos)
                                mask(round(rowpos(a)),round(colpos(b))) = 1;
                            end
                        end
                        mask(dcrow,startcol) = 1; % keep DC point
                        mask = cs_restoresymmetry(mask);

                    else

                        % Uniform, coherent, asymmetric point-based sampling
                        n = row*ceil(col);
                        m = ceil(row*col*ratio);
                        spacing = sqrt(n/m);
                        rowpos = 1:spacing:row;
                        colpos = 1:spacing:col;
                        mask = zeros(row,col);
                        for a=1:length(rowpos)
                            for b=1:length(colpos)
                                mask(round(rowpos(a)),round(colpos(b))) = 1;
                            end
                        end
                        mask(dcrow,startcol) = 1; % keep DC point

                    end

                end

            end

        else

            if(strcmp(coherence,'Incoherent'))

                if(strcmp(samplingtype,'Line'))

                    if(strcmp(symmetry,'Symmetric'))

                        % Nonuniform, incoherent, symmetric line-based sampling
                        m = ceil(col*ratio/2);
                        stddevedge = 3*density/sqrt(ratio);
                        mask = zeros(row,col);
                        mask(:,startcol) = 1; % keep DC line
                        for i=2:m
                            b = startcol;
                            while(mask(1,b)==1)
                                b = startcol + round(abs(randn*((col/2)/stddevedge)));
                                if(b>col)
                                    b = startcol;
                                end
                            end
                            mask(:,b) = 1;
                        end
                        mask = cs_restoresymmetry(mask);

                    else

                        % Nonuniform, incoherent, asymmetric line-based sampling
                        m = ceil(col*ratio);
                        stddevedge = 3*density/sqrt(ratio*2);
                        mask = zeros(row,col);
                        mask(:,startcol) = 1; % keep DC line
                        for i=2:m
                            b = startcol;
                            while(mask(1,b)==1)
                                b = startcol + round(randn*((col/2)/stddevedge));
                                if(b>col||b<1)
                                    b = startcol;
                                end
                            end
                            mask(:,b) = 1;
                        end

                    end

                    npdf = pdf('norm',1:col,col/2,(col/2)/stddevedge);
                    npdf = ones(row,1)*(npdf/max(npdf));

                else

                    if(strcmp(symmetry,'Symmetric'))

                        % Nonuniform, incoherent, symmetric point-based sampling
                        m = ceil(row*col*ratio/2);
                        stddevedge = 2*density/sqrt(ratio);
                        mask = zeros(row,col);
                        mask(dcrow,startcol) = 1; % keep DC point
                        for i=2:m
                            a = dcrow;
                            b = startcol;
                            while(mask(a,b)==1)
                                a = dcrow + round(randn*((row/2)/stddevedge));
                                b = startcol + round(abs(randn*((col/2)/stddevedge)));
                                if(a<1||a>row)
                                    a = dcrow;
                                end
                                if(b>col)
                                    b = startcol;
                                end
                            end
                            mask(a,b) = 1;
                        end
                        mask = cs_restoresymmetry(mask);

                    else

                        % Nonuniform, incoherent, asymmetric point-based sampling
                        m = ceil(row*col*ratio);
                        stddevedge = 2*density/sqrt(ratio*2);
                        mask = zeros(row,col);
                        mask(dcrow,startcol) = 1; % keep DC point
                        for i=2:m
                            a = dcrow;
                            b = startcol;
                            while(mask(a,b)==1)
                                a = dcrow + round(randn*((row/2)/stddevedge));
                                b = startcol + round(randn*((col/2)/stddevedge));
                                if(a<1||a>row)
                                    a = dcrow;
                                end
                                if(b>col)
                                    b = startcol;
                                end
                            end
                            mask(a,b) = 1;
                        end

                    end

                    npdf = pdf('norm',1:col,col/2,(col/2)/stddevedge);
                    npdf = (npdf'*npdf);

                end

            else

                % Nonuniform, coherent sampling

            end

        end
    
        % compute total interference and see how it compares
        spectrum = fft2(mask);
        newinterference = max(abs(spectrum(2:end)));
        if(newinterference<mininterference)
            finalmask = mask;
            mininterference = newinterference;
        end
        
    end
    
    mask = finalmask;
    
end

