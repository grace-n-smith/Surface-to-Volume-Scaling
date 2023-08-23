%FIGURES?
set(0,'DefaultFigureVisible','off')
%set(0,'DefaultFigureVisible','on')

%%% PROVIDE PATH TO FOLDER WITH IMAGES 
%path_to_ims = "E:\Light microscopy\normal_vs_polyploidy\onlyfullimg\osmotic_shock_rep1\15x13\";
%path_to_ims = "E:\Light microscopy\normal_vs_polyploidy\onlyfullimg\osmotic_shock_rep1\15x13\";
path_to_ims = "E:\Light microscopy\normal_vs_polyploidy\08-08-2022_L1210_FUCCI_Membrane_Barasertib_osmotic-shock_rep2\";

%%% PROVIDE PATH TO FOLDER FOR ANALYSES
% path_to_analyses = "E:\Dropbox (MIT)\Surface to Volume Scaling Project Data Sharing\UROP_code_analyses\Analyses\"
% path_to_analyses = "C:\Users\Lenovo\Dropbox (MIT)\Surface to Volume Scaling Project Data Sharing\UROP_code_analyses\Analyses\";
path_to_analyses = "E:\Analyses\";

path_to_drive = "E:\Analyses\10-May-2023\";


% scan directory for *.dv
image_files = dir(strcat(path_to_ims,'\*R3D_D3D.dv'));
[num_of_files, ~] = size(image_files); % count files (for iterable)

%DANGER ZONE
dz = 60;
z_size = 15;
t_size = 13;


%%% import data from csv %%%

%what analysis to use?
date = datetime(2023,5,10);

% column names (for personal reference, should be specified in github README)
varNames = ["imageID", "im_num", "cellID", "radius", "z_position", "t_position", "x_pos", "y_pos", "curated", "omit", "edited", "comments", "true radius", "inner baseline", "outer baseline", "avg peak height", "FWHM", "Integration"];

% SUMMON THE DATA
datatable = readcell("E:\Analyses\10-May-2023\datatable_20-Mar-2023_curation rep2 4.xlsx");

%%
% find indices of each image to enable iteration
[vert, horz] = size(datatable);
indices = zeros(num_of_files, 1);
%% 

% this array finds where the beginning of the data is for each image
prev = 0;
for row  = 1:vert
    if datatable{row, 2} ~= prev
        indices(datatable{row, 2}, 1) = row;
    end
    prev = datatable{row, 2};
    %disp(indices)
end
indices = indices(1:30,:);
indices(30, 1) = 195;

%disp(indices)
% scan directory for *.mat (storing eclipse images produced by cookiecutter.m)
eclipse_files = dir(strcat(path_to_drive,'\*.mat'));
[num_of_eclipses, ~] = size(eclipse_files); % count files (for iterable)

%% 

for im_num = 1:num_of_files

    for cellnum = indices(im_num):indices(im_num+1)-1 %iterate through every cell  in every image
        
        %each cell has it's own eclipse matrix to pull from

        if ismissing(datatable{cellnum, 10}) %only analyze cells that were included in the manual curation
            [cellslice, new_filepath] = find_cellslice(image_files, im_num, cellnum, path_to_drive, datatable);
            mkdir(new_filepath);
            disp(new_filepath); %displays to show progress

            %take cell data from table
            xcenter = datatable{cellnum, 7};
            ycenter = datatable{cellnum, 8}; 
            radius  = datatable{cellnum, 4};
            zlayer  = datatable{cellnum, 5};
            tlayer  = datatable{cellnum, 6};

            [X, Y] = make_endpoints(xcenter, ycenter, radius*1.5); % creates X, Y, arrays containing the 
                                                                   % endpoints of all of the line
                                                                   % profiles for plotting
            
            % Display line profiles on image
            fig = figure();
            imshow(transpose(cellslice), [])
            hold on
            for i = 1:360
                plot(X(i,:),Y(i, :),'Color','r','LineWidth', 0.0002)
                hold on
            end
            text(0, 14, extractAfter(image_files(im_num).name, strlength(image_files(im_num).name)-15) + "   z = " + string(zlayer) + ", t = " + string(tlayer) + ", cell #" + string(datatable{cellnum, 3}), 'Color', 'red', 'FontSize', 10, 'Interpreter', 'none');
            saveas(fig, new_filepath + "/profilevisualization.jpg")
            hold off

            % each line profile is its own column
            c = double.empty(0, 360);
            for i = 1:360
                ci = improfile(transpose(cellslice), X(i,:),Y(i, :),radius*1.5);
                c = [c ci];
            end

            fig2 = figure();
            plot(c)
            saveas(fig2, new_filepath + "/allprofsnotaligned.jpg")


            %%% currently the line profiles are not aligned, and before
            %%% alignment we need to make sure to toss the line profiles
            %%% that have true zeros withing the found circle (intersecting
            %%% circles).

            % make all true zeros into NaNs
            indices2 = find(abs(c)==0.0);
            c(indices2) = NaN;
            
            profcontainsnan = zeros(360, 1);
            for i = 1:360
               if sum(isnan(c(1:radius, i)), "all") ~= 0
                    c(:, i) = NaN(size(c(:, i)));
                    profcontainsnan(i) = 1;
                end
            end
            
            
            cleaned_profiles = c;
            % maxima within a certain distance of the identified radius get pinged
            % find index of maximum on each row???
            zero_sec = NaN(round(radius-radius*0.30), 360);
            cleaned_profiles(1:round(radius-radius*0.30), 1:360) = zero_sec;
             
%                 for i = 1:360
%                     [peaks, pos] = findpeaks(cleaned_profiles(:, i));
%                     display(cleaned_profiles(:, i));
%                     [M,I] = max(peaks);
%                     loc(i, 1) = pos(I);
%                 end
            [peaks, pos] = max(cleaned_profiles);
            loc = pos';
            xrelativecir = zeros(360, 1);
            yrelativecir = zeros(360, 1);
            for i = 1:360
                xrelativecir(i, 1) = cos(i*2*pi/360)*loc(i) +xcenter;%sin(2*pi*i/360)*loc(i);
                yrelativecir(i, 1) = sin(i*2*pi/360)*loc(i) +ycenter;%cos(2*pi*i/360)*loc(i);
            end
            fig4 = figure();
            imshow(transpose(cellslice), [])
            hold on
            
            [prof_i, prof_j] = size(cleaned_profiles);
            aligned_profiles = NaN(prof_i, prof_j);
            amount_moved = NaN(360, 1);

            for j = 1:prof_j

                if profcontainsnan(j) == 0
                    col = cleaned_profiles(:, j);
                    amount_moved(j, 1) = -(radius - int32(loc(j)));
                    colprime = circshift(col, radius - int32(loc(j)));
                    if loc(j) - radius < 0
                        %colprime(abs(int32(loc(j))-radius):end,1) = NaN(abs(int32(loc(j))-radius), 1);
                        colprime = vertcat(NaN(abs(int32(loc(j))-radius), 1), colprime(abs(int32(loc(j))-radius-1):size(colprime, 1)));
                    elseif loc(j) - radius > 0
                        colprime((length(colprime) - abs(int32(loc(j))-radius) + 1):length(colprime)) = NaN(abs(int32(loc(j))-radius), 1);
                    end
                    aligned_profiles(:, j) = colprime;
                end
                
            end

%             %EXCLUDE SHIT THAT IS IN THE ANGULAR BLACKOUT ZONES
%             if datatable{cellnum, 12???} ~= empty???
%                 
%             end

            %AVG THE DISPLACEMRNT TO GET TRUE RADIUS
            radshift = mean(amount_moved, "omitnan");
            truerad = radius + radshift;
            datatable{cellnum, 13} = truerad;
            colors = zeros(size(amount_moved, 1), 3);
            colors(:, 1) = 1; 
            colors(amount_moved>=0, 1) = 0;
            colors(amount_moved>=0, 2) = 1;

            scatter(xrelativecir, yrelativecir, [], colors) %red if decrease in radius, green if increase
            saveas(fig4, new_filepath + "\maxalongprofile.jpg")

            
            
        
            

            %PLOT ALIGNED PROFILES
            fig5 = figure();
            plot(aligned_profiles)
            saveas(fig5, new_filepath + "/allprofsaligned.jpg")
            
            %MEAN OF ALIGNED PROFILES
            fig6 = figure();
            meanprofile = mean(aligned_profiles, 2, 'omitnan');
            plot(meanprofile)
            saveas(fig6, new_filepath + "/meanaligned.jpg")

            %SAVE BASELINES
            innerbaselinepos = 0.90;
            outerbaselinepos = 1.20;
            datatable{cellnum, 14} = meanprofile(int32(radius*innerbaselinepos));
            datatable{cellnum, 15} = meanprofile(int32(radius*outerbaselinepos));
            
            %SAVE MAX OF MEAN PROFILE
            datatable{cellnum, 16} = max(meanprofile);

            %FWHM
            [datatable{cellnum, 17}, mover2, leftind, rightind] = findFWHM(meanprofile, datatable{cellnum, 14}, radius);
            fig8 = figure();
            plot(meanprofile);
            hold on;
            line = line([leftind rightind], [mover2 mover2]);
            line.LineStyle = '-';
            line.Color = 'red';
            hold off;
            saveas(fig8, new_filepath + "/meanFWHM.jpg");


            %INTEGRAL
            ind = find(meanprofile(round(radius*0.8):end)>=datatable{cellnum, 14});	% Find indicies where profile>innerbaseline
            leftindex = min(ind) + round(radius*0.8);
            rightindex = max(ind) + round(radius*0.8);
            intwithbase = trapz(meanprofile(leftindex:rightindex, 1));
            base = (rightindex - leftindex)*datatable{cellnum, 14};
            int = intwithbase - base;
            datatable{cellnum, 18} = int;

            fig9 = figure();
            patchx = [leftindex:(rightindex-1); leftindex:(rightindex-1); (leftindex+1):rightindex; (leftindex+1):rightindex];
            patchy = [datatable{cellnum, 14}*ones(1, rightindex-leftindex); (meanprofile(leftindex:(rightindex-1)))'; (meanprofile((leftindex+1):rightindex))'; datatable{cellnum, 14}*ones(1, rightindex-leftindex)];
%             size(patchx)
%             size(patchy)

            integralpatch = patch(patchx, patchy, 'r');
            hold on;
            meanplot = plot(meanprofile);
            
            set(integralpatch, 'facecolor', [1 0.8 0.8], 'edgecolor', 'none');
            hold off;
            saveas(fig9, new_filepath + "/meanintegral.jpg");



            fig7 = figure();
            %plot(nanmean(adjusted_line_profs))
            e = transpose(std(aligned_profiles, 0, 2,"omitnan"));
            [jmax, imax] = size(aligned_profiles);
            avg = transpose(mean(aligned_profiles, 2, "omitnan"));
            x4 = [1:jmax-1; 1:jmax-1; 2:jmax; 2:jmax];
            y4 = [avg(1:end-1) + e(1:end-1); avg(1:end-1) - e(1:end-1); avg(2:end) - e(2:end); avg(2:end) + e(2:end)];
            
            hp = patch(x4, y4, 'r');
            hold on;
            meanplot2 = plot(meanprofile);
            
            set(hp, 'facecolor', [1 0.8 0.8], 'edgecolor', 'none');
            %set(hl, 'color', 'r');
            hold off;
            saveas(fig7, new_filepath + "/meanalignedwithsd.jpg")

            clearvars cellslice amount_moved X Y x4 y4 patchx patchy e colprime colors col aligned_profiles cleaned_profiles indices2 zero_sec c avg meanprofile fig fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9 line
        end
        
    end
    
end
%% 
mask = cellfun(@(x) any(isa(x,'missing')), datatable); % using isa instead of ismissing allows white space through
datatable(mask) = {[]}
writecell(datatable, path_to_analyses +"10-May-2023\datatable_profiledata_rep2_timepoint1_ver2" + ".xlsx" );
%%
%FUNCTIONS
function [cellslice, new_filepath] = find_cellslice(image_files, im_num, cellnum, path_to_drive, datatable)
    tag = image_files(im_num).name;
    tag2 = tag(strlength(tag)-15:strlength(tag)-11);
    %need to summon this cell in particular
    new_filepath = path_to_drive + "im" + tag2 + "_cell" + string(datatable{cellnum, 3});
    filepath = string(path_to_drive + "im" + tag2 + "_cell" + string(datatable{cellnum, 3}) + ".mat");
    cellslice = load(filepath).cellslice;
    %display(cellslice)
end

function [X, Y] = make_endpoints(xcenter, ycenter, radius)
    
    X1 = xcenter * ones(360,1);
    Y1 = ycenter * ones(360,1);
    
    temp = transpose(0:(2*pi/360):(2*pi - 2*pi/360));
    X2 = radius*cos(temp) + X1;
    Y2 = radius*sin(temp) + Y1;

    X = [X1 X2];
    Y = [Y1 Y2];

end


function [FWHM, mover2, leftind, rightind] = findFWHM(fx, baseline, radius)
    %findFWHM: Finds the full width half max (FWHM) of a function
    [m, n] = max(fx);		%	Find maximum value and index
    mover2 = (m-baseline)/2 + baseline;
    ind = find(fx(round(radius*0.8):end)>=((m-baseline)/2 + baseline));	%	Find indicies where I>=max(I)/2
    leftind = min(ind) + round(radius*0.8);			%	Leftmost index
    rightind = max(ind) + round(radius*0.8);		%	Rightmost index
    %	Get FWHM
    FWHM = abs(rightind-leftind);
end

function c = conditional(condition , a , b)
    if condition
        c = a;
    else
        c = b;
    end
end
