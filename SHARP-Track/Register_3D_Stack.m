% ------------------------------------------------------------------------
%          3D registration of serial slices
% ------------------------------------------------------------------------


%% ENTER FILE LOCATION
% note: use grayscale images of adjust code
% note: transform all images in the processed images folder or adjust code

% directory of histology
processed_images_folder = 'C:\Users\SWC\Downloads\test\test\processed'; 

% neuron data file to LOAD
neuron_file = 'C:\Drive\Histology\neuron_data.mat';

% translate the neuron coordinates to the cropped and downsampled histology
% reference frame with scaling and translation in X Y Z
rescaling_factors = [11.5439330543933	11.5567375886525	40];
translation_factors = [430 282 1];

% registered neuron data file to SAVE
registered_neuron_file = 'C:\Drive\Histology\neuron_data_registered.mat';

% affine transformation to SAVE
affine_transformation_file = 'C:\Drive\Histology\affine_transformation.mat';
transformation_matrix_file = 'C:\Drive\Histology\transformation_matrix.mat';

% plane to view ('coronal', 'sagittal', 'transverse') - only tested on coronal place so far!
plane = 'coronal';


%% GENERATE 3D STACK OF SLICE IMAGES

% find image names
processed_images = dir([processed_images_folder filesep '*.tif']);
processed_image_names = natsortfiles({processed_images.name});

% get the dimensions of the 3D stack
total_num_files = size(processed_images,1); disp(['found ' num2str(total_num_files) ' processed slice images']);
if strcmp(plane,'coronal')
    slice_size = [800 1140];
elseif strcmp(plane,'sagittal')
    slice_size = [800 1320];
elseif strcmp(plane,'transverse')
    slice_size = [1140 1320];
end
stack = zeros(slice_size(1), slice_size(2), total_num_files, 'uint8');

% fill the stack with the images
for i = 1:total_num_files
    processed_image_name = processed_image_names{i};
    current_slice_image = imread(fullfile(processed_images_folder, processed_image_name));
    if ~sum(size(current_slice_image)==slice_size)
        current_slice_image = current_slice_image(1:slice_size(1), 1:slice_size(2));
    end
    stack(:,:,i) = current_slice_image;
end


%% EXTRACT COORDINATES OF MATCHING POINTS

atlas_coordinates = zeros(0,3);
stack_coordinates = zeros(0,3);

for i = 1:total_num_files
    % load clicked transform points
    slice_name = processed_image_names{i}(1:end-4);
    folder_transformations = fullfile(processed_images_folder, ['transformations' filesep]);
    load([folder_transformations slice_name '_transform_data.mat']);
    
    % get atlas coordinates - get z using allen location from save_transform
    xy_atlas = save_transform.transform_points{1};
    z_atlas = save_transform.allen_location{1};
    
    % get offset due to angling
    slice_angle = save_transform.allen_location{2};
    offset_map = generate_offset_map(slice_size, slice_angle);
    indices_from_clicked_points = sub2ind(slice_size, xy_atlas(:, 2), xy_atlas(:, 1));
    offset = offset_map(indices_from_clicked_points);
    z_atlas = z_atlas + offset;
    
    % put into array
    current_atlas_coordinates = zeros(length(z_atlas), 3);
    current_atlas_coordinates(:, 1:2) = xy_atlas;
    current_atlas_coordinates(:, 3) = z_atlas;
    atlas_coordinates = cat(1, atlas_coordinates, current_atlas_coordinates);

    % get image stack coordinates
    xy_histology = save_transform.transform_points{2};
    z_histology = i;
    
    % put into array
    current_stack_coordinates = zeros(length(z_atlas), 3);
    current_stack_coordinates(:, 1:2) = xy_histology;
    current_stack_coordinates(:, 3) = z_histology;
    stack_coordinates = cat(1, stack_coordinates, current_stack_coordinates);
    
end


%% GENERATE 3D AFFINE TRANSFORMATION (least square solution)

% initialize arrays
M = size(atlas_coordinates,1);
X = [atlas_coordinates ones(M,1)];
U = stack_coordinates;
% solve for the first three columns of T (We know that X * T = U)
Tinv = X \ U;
% add fourth column
Tinv(:,4) = [0 0 0 1]';
T = inv(Tinv);
T(:,4) = [0 0 0 1]';
% generate affine transform object
affine_transform = affine3d(T);
% save the output
save(affine_transformation_file, "affine_transform")
save(transformation_matrix_file, "T")

%% APPLY AND VERIFY 3D AFFINE TRANSFORMATION

% Apply transform
R = imref3d([800 1140 1320]);
registered_stack = imwarp(stack, affine_transform, 'OutputView',R);
if strcmp(plane,'coronal')
    registered_stack = permute(registered_stack, [3 1 2]);
end

% test on points clicked during the transformation
stack_points_test = ones(length(U), 4);
stack_points_test(:, 1:3) = stack_coordinates;
transform_prediction = stack_points_test*T;
% RMSE = sqrt( mean((transform_prediction(:,1:3) - atlas_coordinates).^2, 'all'));
mean_absolute_error_in_microns = mean(abs(transform_prediction(:,1:3) - atlas_coordinates), 'all')*10


%% APPLY 3D AFFINE TRANSFORMATION TO A LABELLED NEURON

% create variable called neuron_data -- now rescaled_neuron_data which is
% in the reference frame of the downsampled and translated histology images
% ID, TYPE, X, Y, Z, R, P, SYN
load(neuron_file)

% extract position from neuron
xyz_neuron = neuron_data(:, 3:5);
xyz_neuron_rescaled = xyz_neuron ./ rescaling_factors + translation_factors;

% apply affine transformation to neuron coordinates
xyz_neuron_for_transform = ones(length(xyz_neuron_rescaled), 4);
xyz_neuron_for_transform(:, 1:3) = xyz_neuron_rescaled;
xyz_neuron_transformed = xyz_neuron_for_transform * T;
xyz_neuron_transformed = round(xyz_neuron_transformed(:,1:3));
save(registered_neuron_file, "xyz_neuron_transformed")

%% FUNCTIONS

function offset_map = generate_offset_map(slice_size, slice_angle)
  
  if slice_angle(1)==0; offset_DV = 0;
  else; offset_DV = -slice_angle(1):sign(slice_angle(1)):slice_angle(1);
  end; start_index_DV = 1; 
 
  % loop through AP offsets
  num_DV_iters_add_ind = floor( (slice_size(1) - floor( slice_size(1) / length(offset_DV))*length(offset_DV)) / 2);
  for curr_DV_iter = 1:length(offset_DV)
      cur_offset_DV = offset_DV(curr_DV_iter);
      if cur_offset_DV == slice_angle(1)
          end_index_DV = slice_size(1);
      elseif curr_DV_iter <= num_DV_iters_add_ind  || length(offset_DV - curr_DV_iter) < num_DV_iters_add_ind
          end_index_DV = start_index_DV + floor( slice_size(1) / length(offset_DV));
      else
          end_index_DV = start_index_DV + floor( slice_size(1) / length(offset_DV)) - 1;
      end
      
       if slice_angle(2)==0;  offset_ML = 0;
       else; offset_ML = -slice_angle(2):sign(slice_angle(2)):slice_angle(2);
       end; start_index_ML = 1;
      % nested: loop through ML offsets
      num_ML_iters_add_ind = floor( (slice_size(2) - floor( slice_size(2) / length(offset_ML))*length(offset_ML)) / 2);
      for curr_ML_iter = 1:length(offset_ML)
          cur_offset_ML = offset_ML(curr_ML_iter);
          if cur_offset_ML == slice_angle(2)
              end_index_ML = slice_size(2);
          elseif curr_ML_iter <= num_ML_iters_add_ind  || length(offset_ML - curr_ML_iter) < num_ML_iters_add_ind
              end_index_ML = start_index_ML + floor( slice_size(2) / length(offset_ML));
          else
              end_index_ML = start_index_ML + floor( slice_size(2) / length(offset_ML)) - 1;
          end

          offset_map(start_index_DV:end_index_DV, start_index_ML:end_index_ML) = cur_offset_DV + cur_offset_ML;

          start_index_ML = end_index_ML + 1;
       end
      start_index_DV = end_index_DV + 1;
  end  
  
end

image_folder = 
function NeuronDataToDownsampledStackSpace(image_folder, image_file_names, microns_per_pixel, microns_per_pixel_after_downsampling, atlas_reference_size)
    microns_between_slices = 40;
    z_translation_microns = 0;
    neuron_data_path = 'C:\Drive\Histology\neuron_data.mat';
    zero_indexing_in_neuron_data = true;
    % load first image
    image = imread(fullfile(image_folder,image_file_names{1}));
    original_image_size = size(image);
    downsampled_image_size = round(original_image_size*microns_per_pixel/microns_per_pixel_after_downsampling);
    XY_downsampling_factor = original_image_size ./ downsampled_image_size;
    Z_rescaling_factor = microns_between_slices;
    rescaling_factors = [XY_downsampling_factor Z_rescaling_factor];
    translation_factors = [floor((atlas_reference_size(2) - downsampled_image_size(2)) / 2) + mod(downsampled_image_size(2),2) ...
                            floor((atlas_reference_size(1) - downsampled_image_size(1)) / 2) + mod(downsampled_image_size(1),2) ...
                                z_translation_microns/Z_rescaling_factor];
    if zero_indexing_in_neuron_data
        translation_factors = translation_factors + [1 1 1];
    end
    neuron_data = load(neuron_data_path);
    neuron_data = neuron_data.neuron_data;
    colXYZ = [3, 4, 5]; 
    rescaled_neuron_data = neuron_data;
end