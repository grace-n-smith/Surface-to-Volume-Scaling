
%FIGURES?
set(0,'DefaultFigureVisible','off')

%%% PROVIDE PATH TO FOLDER WITH IMAGES 
%path_to_ims = "E:\Light microscopy\normal_vs_polyploidy\onlyfullimg\osmotic_shock_rep1\15x13\";
path_to_ims = "E:\Light microscopy\normal_vs_polyploidy\08-08-2022_L1210_FUCCI_Membrane_Barasertib_osmotic-shock_rep2\";

%%% PROVIDE PATH TO FOLDER FOR ANALYSES
% path_to_analyses = "E:\Dropbox (MIT)\Surface to Volume Scaling Project Data Sharing\UROP_code_analyses\Analyses\"
path_to_analyses = "E:\Analyses\";

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

% column names
varNames = ["imageID", "im_num", "cellID", "radius", "z_position", "t_position", "x_pos", "y_pos", "curated", "omit", "edited", "comment"];

% SUMMON THE DATA
datatable = readcell("E:\Analyses\10-May-2023\datatable_20-Mar-2023_curation rep2.xlsx");

%%


% need to do something to find indices of each image
[vert, horz] = size(datatable);
indices = zeros(num_of_files, 1);

% this array says where the beginning of the data is for each image
prev = 0;
for row  = 2:vert %starts at 2 because of title row
    if datatable{row, 2} ~= prev
        indices(datatable{row, 2}, 1) = row;
    end
    prev = datatable{row, 2};
end

%%

mkdir(path_to_analyses, string(date));
path_to_all_cutouts = path_to_analyses + string(date);

count = 2;

for im_num = 1:num_of_files
    norm_im = normalizeim(path_to_ims, image_files, im_num, z_size, t_size);

    while datatable{count, 2} == im_num
        %create copy of image to make cookie cutter thing from. this should
        %be the layer that the cell equator is on
        if ismissing(datatable{count, 10})
            cutout = savecut(datatable, image_files, im_num, path_to_all_cutouts, count, norm_im);
            
            
        end
        count = count + 1;
    end
    clearvars norm_im
end


% %%
% for cellnum = 1:size(datatable, 1)
%     %open image
%     path = strcat(path_to_ims, datatable{cellnum, 1});
%     data = bfopen(char(path));
%     
%     image = data{1, 1}; %this cell contains images + path
%     A = image(:, 1); %just images
%     multiDim = reshape(A, z_size, t_size); %changing from linear to matrix for comfortable indexing
%     
%     %cell array of matrices
%     repackaged_im = repackage(z_size, t_size, multiDim);
% 
%     %will normalize to 4000? %%% ONLY NORMALIZES t=1 %%%
%     norm_im = cell(z_size, t_size);
%     for i = 1:z_size
%         for j = 1:t_size
%             im = repackaged_im{i, j};
%             normalized = 4000/mean2(im) * im;
%             norm_im{i, j} = normalized;
%         end
%     end
%     %disp(norm_im)
%     
%     %cut
%     for i = 1:size(datatable, 1)
%         if ((datatable{i, 2} == datatable{cellnum, 2})) && (cellnum ~= i)
%             x = datatable{i, 5};
%             y = datatable{i, 6};
%             r = datatable{i, 4};
%             %x^2 + y^2 = r^2
%             
%             %iterate through image i guess
%             save = cell(z_size, 1)
%             for j = 1:z_size
%                 cut = norm_im{j, 1};
%                 %for k = 1:size(norm_im(1))
%                     for l = 1:size(norm_im(2))
%                         if (distance([x, y], [j, l]) < r + dz)
%                             cut(j,l) = 0;
%                             %display("NaN")
%                         end
%                     end
%                 %end
%                 save{j, 1} = cut;
%                 %imshow(b, [])
%             end
%             
%         end
%     end
% end
% %%
% imshow(save{8, 1}, [])
% %%
% 
% for im_num = 3:3
%     
%     for cellnum = 1:size(datatable, 1)
% 
%     
%     end
%     %cut = cookiecutter(datatable, image, sorted)
%     
%     %FOR EACH CELL
%     profiles = draw_line_profiles(cellnum, sorted, image, dz); %chooses max circle
%     %display(profiles) %each profile is a horizontal line?
%     aligned_profiles = alignPeaks(cellnum, sorted, profiles);
%     avg = averageProfile(aligned_profiles);
%     %display(aligned_profiles)
%     %negavg = -avg;
%     %[peakValues, indexesOfPeaks, w, p] = findpeaks(negavg)
%     %plot(profiles)
% 
%     inner_base = avg(r-50);
%     outer_base = avg(r+50);
% 
%     
%     lineprof = line([0, 360],[inner_base, inner_base]);
%     saveas(lineprof, path_to_analyses + string(datetime("today"))+"\"+image_files(im_num).name + '\line'+ cellnum + '.jpg')
% 
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function norm_im = normalizeim(path_to_ims, image_files, im_num, z_size, t_size)
    %open image
    path = strcat(path_to_ims, image_files(im_num).name);
    data = bfopen(char(path));
    
    image = data{1, 1}; %this cell contains images + path
    A = image(:, 1); %just images
    multiDim = reshape(A, z_size, t_size); %changing from linear to matrix for comfortable indexing
    
    %cell array of matrices
    repackaged_im = repackage(z_size, t_size, multiDim);

    %will normalize to 4000? %%% ONLY NORMALIZES t=1 %%%
    norm_im = cell(z_size, t_size);
    for i = 1:z_size
        for j = 1:t_size
            im = repackaged_im{i, j};
            normalized = 4000/mean2(im) * im;
            norm_im{i, j} = normalized;
        end
    end
    %disp(norm_im)
end

function [group, cutout] = savecut(datatable, image_files, im_num, path, count, norm_im)
    %display(selected_circle);
    z = datatable{count, 5};
    t = datatable{count, 6};
    frame = norm_im{z, t};

    x = datatable{count, 7};
    y = datatable{count, 8};
    r = datatable{count, 4};
    %COMMENT
    dz = 150;
    
    cropped = imcrop(frame, [x-(dz+r) y-(dz+r) 2*(dz+r) 2*(dz+r)]);
    cutout = figure();
    imshow(cropped, []);
    group = viscircles(size(cropped)./2, r,'EdgeColor','r');
    text(0, 14, extractAfter(image_files(im_num).name, strlength(image_files(im_num).name)-15) + "   z = " + string(z) + ", t = " + string(t) + ", cell #" + string(datatable{count, 3}), 'Color', 'red', 'FontSize', 10, 'Interpreter', 'none');
    tag = image_files(im_num).name;
    tag2 = tag(strlength(tag)-15:strlength(tag)-11);
    saveas(group , string(path + '\image' + tag2 +  '_cell' +string(datatable{count, 3}) +'.jpg'));
end

function repackaged_im = repackage(z_size, t_size, multiDim)
    % re-package data to be cell array of matrices
    repackaged_im = cell(z_size, t_size);
    for i = 1:z_size
        for j = 1:t_size
            sample = multiDim(i, j);
            repackaged_im{i, j} = sample{1, 1};
        end
    end
    %display(repackaged_im);
end