%%% SHOULD THIS CODE DISPLAY LIVE FIGURES?
set(0,'DefaultFigureVisible','off') % only save figures to filesystem
%set(0,'DefaultFigureVisible','on') % display figures as they are created

%%% PROVIDE PATH TO FOLDER WITH IMAGES 
pathToIms = "E:\Light microscopy\normal_vs_polyploidy\08-08-2022_L1210_FUCCI_Membrane_Barasertib_osmotic-shock_rep1\";
%%% PROVIDE PATH TO FOLDER FOR ANALYSES
pathToAnalyses = "E:\Analyses\";
%path_to_analyses = "C:\Users\Lenovo\Documents\Local UROP data\Analyses\";

%%% PROVIDE IMAGE DETAILS
Z_SIZE            = 15;                 % Z LAYERS
T_SIZE            = 13;                 % TIMEPOINTS
DZ                = 60;                 % DANGER ZONE PX (microscope edge artifacts)

% scan directory for *.dv
imageFiles = dir(fullfile(pathToIms,'\*R3D_D3D.dv')); %file ending that all relevant images have
[nFiles, ~] = size(imageFiles); % count files (for iterable)

% empty array for image paths
files = strings(nFiles, 1);

% fill `files` with image paths
for i = 1:nFiles
    a = strcat(pathToIms , "\" ,imageFiles(i).name);
    files(i, 1) = a;
end

% column names
varNames = ["imageID", "im_num", "cellID", "radius", "z_position", "t_position", "x_pos", "y_pos", "edited", "omit"];
datatable = cell(0, 10);

%loci.common.DebugTools.setRootLevel('DEBUG');
mkdir(pathToAnalyses, string(datetime("today")));


%%
for imageNo = 1:nFiles 

    %make directory for this image
    imDirPath = string(pathToAnalyses + string(datetime("today")) + "\");
    mkdir(imDirPath, imageFiles(imageNo).name);
    
    %%% TEMP DISPLAY %%%
    % display(image_files(im_num).name)
    %curate editing ignored
    %open image
    repack_im = repackageim(pathToIms, imageFiles, imageNo, Z_SIZE, T_SIZE);

    [~, ~, ~, ~, sorted] = identifycells(repack_im, 1); %run circle finder and sorter, specify t slice
    % sorted: [x_center, y_center, radius, z_slice, t_slice]

    representatives = pickLargest(sorted); % data for the single circles that represent the equator of a cell

    [key, group] = makeKeyImage(repack_im, representatives, imageFiles(imageNo).name);
    %keyim = imshow(key)
    %saveas(keyim,'C:\Users\Lenovo\Documents\Local UROP Data\Analyses\'+string(datetime("today"))+"_image"+string(im_num)+'\key'+string(im_num)+'.jpg');
    keyim = figure();
    %imshow(repackaged_im{7, 1})
    %plot(group, [0, 0, 1, 1])
    %axes('pos',[.1 .6 .5 .3])
    saveas(group , string(pathToAnalyses + string(datetime("today"))+"\"+imageFiles(imageNo).name + '\key.jpg'));


    %%%%% FOR EACH CELL %%%%%
    for cellnum = 1:size(sorted, 2)
        largest = pickLargest(sorted);
        %create empty table for each cell to append to the big 
        celldata = cell(1, 9);
        celldata{1, 1} = imageFiles(imageNo).name;
        celldata{1, 2} = imageNo;
        celldata{1, 3} = cellnum;
        celldata{1, 4} = largest(cellnum, 3);
        celldata{1, 5} = largest(cellnum, 4);
        celldata{1, 6} = largest(cellnum, 5);
        celldata{1, 7} = largest(cellnum, 1);
        celldata{1, 8} = largest(cellnum, 2);

        datatable = vertcat(datatable, celldata); % try [datatable, celldata] 'writemode', 'append'

        disp("Displaying cell #" + cellnum);
        [superimposed, croppedim, r, all_layers] = displayLargest(sorted, repack_im, cellnum, imageNo, imageFiles, pathToAnalyses);
        saveas(croppedim, pathToAnalyses + string(datetime("today"))+"\"+ imageFiles(imageNo).name + "\" + "cropped"+ cellnum + ".jpg");
        clearvars celldata
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clearvars keyim representatives path norm_im multiDim sorted
    
end
%% 
writecell(datatable, pathToAnalyses + string(datetime("today"))+"\datatable_" + string(datetime("today")) + ".csv" );
%%
%datatable = importdata("datatable2.csv")

%%

function repackagedIm = repackageim(pathToIms, imageFiles, imageNo, z_size, t_size)
    %open image
    path = strcat(pathToIms, imageFiles(imageNo).name);
    data = bfopen(convertStringsToChars(path));
    
    image = data{1, 1}; %this cell contains images + path
    A = image(:, 1); %just images
    multiDim = reshape(A, z_size, t_size); %changing from linear to matrix for comfortable indexing
    
    %cell array of matrices
    repackagedIm = repackage(z_size, t_size, multiDim);

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
end


function dis = distance(circle1, circle2)
    %dis = sqrt((center1(1) - center2(1))^2 + ((center1(2) - center2(2))^2));
    dis = sqrt((circle2(1, 1) - circle1(1, 1))^2  + (circle2(1, 2) - circle1(1, 2))^2 );
end

function circles = getCircles(im)
    imshow(imfill(im, "holes"), [])
    [centers, radii] = imfindcircles(imfill(im, "holes"),[90 300],'Sensitivity',0.97,'Method','TwoStage', "ObjectPolarity","bright");
    
    circles = [centers, radii];
end

function largest = pickLargest(sorted)
    %retrieve number of cells
    [~, len] = size(sorted);

    largest = zeros(len, 5); %there are five values for each cell
    %sorted
    for i = 1:len
        %for each identified cell
        allCirc = sorted{1, i};
        [a, ~] = size(allCirc);
        if a > 1
            [~, index] = max(allCirc);
            largest(i, :) = allCirc(index(3), :);
            %i
            %all_circ
        else
            largest(i, :) = allCirc;
        end
    end
end

function [superimposed, croppedim, r, allLayers] = displayLargest(sorted_circles, ztimage, id, imageNo, imageFiles, pathToAnalyses)
    selected_circle = sorted_circles{1, id};

    [a, ~] = size(selected_circle);
    if a > 1
        [~, index] = max(selected_circle);
        maxRadRow = index(3);
    else
        maxRadRow = 1;
    end


    %display(selected_circle);
    z = selected_circle(maxRadRow, 4);
    t = selected_circle(maxRadRow, 5);
    frame = ztimage(z, t);

    superimposed = figure();
    %imshow(frame{1,1}, [])

    x = selected_circle(maxRadRow, 1);
    y = selected_circle(maxRadRow, 2);
    r = selected_circle(maxRadRow, 3);
    %COMMENT
    text(0, 0, string("Slice z = " + string(z) + ", t = " + string(t)))
    viscircles([x,y], r,'EdgeColor','r');
    dz = 150;
    
    cropped = imcrop(frame{1,1}, [x-(dz+r) y-(dz+r) 2*(dz+r) 2*(dz+r)]);
    croppedim = figure();
    imshow(cropped, [])
    viscircles(size(cropped)./2, r,'EdgeColor','r');
    text(0, 14, extractAfter(imageFiles(imageNo).name, strlength(imageFiles(imageNo).name)-15), 'Color', 'red', 'FontSize', 10, 'Interpreter', 'none');

    r = selected_circle(maxRadRow, 3);
    
    allLayers = [];
    %group g
    for i = 1:size(ztimage, 1)
        slice = ztimage(i, t);
        cropped2 = imcrop(slice{1,1}, [x-(dz+r) y-(dz+r) 2*(dz+r) 2*(dz+r)]);
        if i == z
            
            disp("with marker")
            [vert, horz] = size(allLayers);
            allLayers = horzcat(allLayers, max(allLayers,[],'all')*ones(vert,20));
            allLayers = horzcat(allLayers, cropped2);%g = viscircles(size(cropped)./2, r,'EdgeColor','r');
            allLayers = horzcat(allLayers, max(allLayers,[],'all')*ones(vert,20));
            %imshow(all_layers, [])
        else
            allLayers = horzcat(allLayers, cropped2);%g = viscircles(size(cropped)./2, r,'EdgeColor','r');
            %imshow(all_layers, [])
        end
        
        %imshow(cropped2, [])
        
        %add slice # in corner
        %cropped3 = insertText(cropped2, [2*r, 2*r], i, FontSize=10, BoxColor="white", BoxOpacity=0, TextColor="red");
        %imshow(cropped3, [])
        %xy = [(i-1)*2*r*dz, r];
        %cropped3 = rgb2gray(insertText(cropped2, xy, [string(i)], FontSize=20, TextColor="white"));%BoxColor="white", BoxOpacity=0.0, );
        
        %all_layers = rgb2gray(insertText(all_layers, xy, [string(i)], FontSize=20, TextColor="white"));
    end
    %min(all_layers,[],'all')
    %max(all_layers,[],'all')
    grayim = mat2gray(allLayers);
    imwrite(grayim, pathToAnalyses + string(datetime("today"))+"\"+imageFiles(imageNo).name + '\cell_' +string(id) + '_all_layers.jpg')
    %imshow(all_layers, [])
end

function frame = getFrame(zt_image, z, t)
    arguments
        zt_image (:,:)
        z {mustBeInteger}
        t {mustBeInteger}
    end

    frame = zt_image{z, t};

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

function [keyfig, hgroup] = makeKeyImage(im, representatives, name)
    
    key = zeros(size(im{1,1}, 1), size(im{1,1}, 2), 3); % converting to "RGB"

    key(:,:,1) = im{1,1};
    key(:,:,2) = im{1,1};
    key(:,:,3) = im{1,1};
    [maxval, ~] = max(key);
    [maxval2, ~] = max(maxval);

    key = (key)/(maxval2(1, 1, 1))*1;
    
    key = insertText(key, representatives(:, 1:2), (1:size(representatives))', FontSize=60, BoxColor="white", BoxOpacity=0.7, TextColor="black");
    
    keyfig = figure();
    imshow(key)
    for i = 1:size(representatives)
        hgroup = viscircles([representatives(i, 1), representatives(i, 2)], representatives(i,3),'EdgeColor','r');
        hold on
    end
    hgroup2 = text(0, 25, extractAfter(name, strlength(name)-15), 'Color', 'red', 'FontSize', 10, 'Interpreter', 'none');
    hold off
   
end

function [all, friendly, lonely, not_danger, sorted] = identifycells(zt_image, t)
    arguments
        zt_image (:,:)
        t {mustBeInteger}
    end
    
    dz = 60;

    %%% retrieve size of image to use in iteration
    [zmax, tmax] = size(zt_image);

    %%% find circle data for all z layers
    disp("Finding Circle Data...")
    cellLocation = []; %each circle is defined by [x, y, r]
    for z = 1:zmax
        frame = getFrame(zt_image, z, t);
        circLoc = getCircles(frame);
        pos = ones(size(circLoc));
        pos(:,3) = [];
        pos(:, 1) = pos(:, 1)*z;
        pos(:, 2) = pos(:, 2)*t;
        cellLocation = vertcat(cellLocation, [circLoc, pos]); %[x_center, y_center, radius, z_slice, t_slice]
    end
    %display(cell_location)

    %%% weed out lonely circles
    disp("Removing lonely circles")
    friendlyCircles = [];
    lonelyCircles = [];
    for i = 1:size(cellLocation)[1];
        is_close = 0;
        for j = 1:size(cellLocation)[1];

            if distance(cellLocation(i, 1:2), cellLocation(j, 1:2)) >= dz/3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %lol
            else
                is_close = is_close + 1;
            end

        end
        
        if is_close > 2
            friendlyCircles = vertcat(friendlyCircles, cellLocation(i, :));
        else
            lonelyCircles = vertcat(lonelyCircles, cellLocation(i, :));
        end
    end
    %display(lonely_circles)
    %display(friendly_circles)
    
    all = {cellLocation};
    friendly = {friendlyCircles};
    lonely = {lonelyCircles};

    %toss border cells
    disp("Removing Danger Zone circles")
    notDangerZone = [];

    for i = 1:size(friendlyCircles, 1)
        
        if not((friendlyCircles(i, 1) <= dz + friendlyCircles(i, 3)) || (friendlyCircles(i, 1) >= (2040 - (dz + (friendlyCircles(i, 3))))) || (friendlyCircles(i, 2) <= (dz + friendlyCircles(i, 3))) || (friendlyCircles(i, 2) >= (2040 - (dz + friendlyCircles(i, 3)))))
            notDangerZone = vertcat(friendlyCircles(i, :), notDangerZone);
        end

    end

    not_danger = {notDangerZone};

    %make matrix of circles identified for each cell
    disp("Sorting circles...")
    zerosForCat = zeros(size(notDangerZone));
    forIdentification = [notDangerZone, zerosForCat]; %zeros represent whether the circle has been accounted for in classification
    
    %sorted circles (yes im a clown dont @ me)
    sorted = {};
    for j = 1:size(forIdentification, 1)

        %identify base circles
        for i = 1:size(forIdentification, 1)
            if forIdentification(i, 6) == 0
                %sorkels is defined to be a cell array containing an array
                %in each cell. It's length should be equal to the number of
                %circles in the image that are being classified (cells fit
                %for data collection)
                sorted{j} = [];
                sorted{j} = vertcat(forIdentification(i, 1:5), sorted{j});
                forIdentification(i, 6) = 1;
                break
            end
        end
    
        %identify close circles
        for i = 1:size(forIdentification, 1)
    
            if forIdentification(i, 6) == 0
    
                if distance(sorted{j}(1, 1:2), forIdentification(i, 1:2)) <= dz
                    sorted{j} = vertcat(forIdentification(i, 1:5), sorted{j});
                    forIdentification(i, 6) = 1;
                end
          
            end
            
        end
    
        if sum(forIdentification(:, 6)) == size(forIdentification(1))
            break
        end

    end

end