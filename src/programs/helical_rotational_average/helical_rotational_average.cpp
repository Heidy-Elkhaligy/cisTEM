#include "../../core/core_headers.h"
#include <memory>

class
        helical_rotational_average : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

void AverageAlongZ(Image* current_volume);
void AverageRotationallyPerSlice(Image* volume);

IMPLEMENT_APP(helical_rotational_average)

// override the DoInteractiveUserInput

void helical_rotational_average::DoInteractiveUserInput( ) {

    wxString input_filename;
    wxString output_filename;

    UserInput* my_input = new UserInput("helical_rotational_average", 1.0);
    input_filename      = my_input->GetFilenameFromUser("Input helical density map file", "Name of the helical structure file", "my_helical_structure.mrc", true);
    output_filename     = my_input->GetFilenameFromUser("Output azimuthal average helical map", "Name of the output azimuthal average helical structure file", "my_azimuthal_average.mrc", false);

    delete my_input;

    my_current_job.Reset(3); // + 1 the number shown in the string below
    my_current_job.ManualSetArguments("tt", input_filename.ToUTF8( ).data( ), output_filename.ToUTF8( ).data( ));
}

// override the do calculation method which will be what is actually run..

bool helical_rotational_average::DoCalculation( ) {
    // get the arguments for this job..
    wxString input_filename  = my_current_job.arguments[0].ReturnStringArgument( );
    wxString output_filename = my_current_job.arguments[1].ReturnStringArgument( );

    if ( ! DoesFileExist(input_filename.ToStdString( )) ) {
        SendError(wxString::Format("Error: Mask %s not found\n", input_filename.ToStdString( )));
        exit(-1);
    }
    // initiate I/O variables
    MRCFile my_input_filename(input_filename.ToStdString( ), false);
    MRCFile my_output_filename(output_filename.ToStdString( ), true);

    long  number_of_input_images = my_input_filename.ReturnNumberOfSlices( );
    float pixel_size             = my_input_filename.ReturnPixelSize( );
    Image my_image;
    Image my_volume;

    my_image.Allocate(my_input_filename.ReturnXSize( ), my_input_filename.ReturnYSize( ), true);
    my_volume.Allocate(my_input_filename.ReturnXSize( ), my_input_filename.ReturnYSize( ), my_input_filename.ReturnZSize( ), true, true);
    my_image.SetToConstant(0.0);
    my_volume.SetToConstant(0.0);

    // average along Z the input volume then do rotational average per slice
    my_volume.ReadSlices(&my_input_filename, 1, number_of_input_images);
    AverageAlongZ(&my_volume);
    AverageRotationallyPerSlice(&my_volume);
    // my_volume.ApplyRampFilter( );
    //wxPrintf("volume z dimension is %i, image z dimension is %i\n", my_volume.physical_address_of_box_center_z, my_image.physical_address_of_box_center_z);
    my_volume.WriteSlices(&my_output_filename, 1, number_of_input_images);

    Image               model_volume;
    Image               projection_volume_3d;
    ReconstructedVolume input_3d;
    Image               my_slice;
    Image               projection_volume_image;
    Image               padded_projection_volume_image;
    AnglesAndShifts     my_parameters;

    projection_volume_3d.Allocate(my_input_filename.ReturnXSize( ) * 2, my_input_filename.ReturnYSize( ) * 2, my_input_filename.ReturnZSize( ) * 2, true, true);
    my_slice.Allocate(my_input_filename.ReturnXSize( ), my_input_filename.ReturnYSize( ), true);
    model_volume.Allocate(2 * my_slice.logical_x_dimension, 2 * my_slice.logical_y_dimension, 2 * my_slice.logical_x_dimension, true);
    model_volume.SetToConstant(0.0);
    projection_volume_3d.SetToConstant(0.0);
    my_slice.SetToConstant(0.0);
    // copy one slice from the azimuthal average volume
    // will this work?
    my_slice.CopyFrom(&my_volume);

    float edge_value = my_slice.ReturnAverageOfRealValuesOnEdges( );
    my_slice.Resize(model_volume.logical_x_dimension, model_volume.logical_y_dimension, 1, edge_value);

    //my_slice.QuickAndDirtyWriteSlice("my_slice.mrc", 1, my_input_filename.ReturnZSize( ));

    // fill in the model volume with the azimuthal average slice
    long volume_counter = 0;
    for ( int z = 0; z < model_volume.logical_z_dimension; z++ ) {
        for ( int y = 0; y < model_volume.logical_y_dimension; y++ ) {
            for ( int x = 0; x < model_volume.logical_x_dimension; x++ ) {
                long pixel_coord_xy                      = my_slice.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
                model_volume.real_values[volume_counter] = my_slice.real_values[pixel_coord_xy];
                volume_counter++;
            }
            volume_counter += my_slice.padding_jump_value;
        }
    }

    input_3d.InitWithDimensions(my_input_filename.ReturnXSize( ) * 2, my_input_filename.ReturnYSize( ) * 2, my_input_filename.ReturnZSize( ) * 2, pixel_size);
    input_3d.density_map->CopyFrom(&my_volume);
    float mask_radius    = FLT_MAX; //100 - FLT_MAX
    input_3d.mask_radius = mask_radius;
    input_3d.PrepareForProjections(0.0, 2.0 * pixel_size); // 0.0, 2.0 * pixel_size float low resolution limit and high resolution limit, bool approximate bining = F and apply_bining = T

    projection_volume_3d.CopyFrom(input_3d.density_map);
    // deallocate the reconstruction volume
    input_3d.Deallocate( );

    projection_volume_image.Allocate(my_input_filename.ReturnXSize( ), my_input_filename.ReturnYSize( ), true);
    padded_projection_volume_image.Allocate(my_input_filename.ReturnXSize( ) * 2, my_input_filename.ReturnYSize( ) * 2, false); // as my volume now is already padded so no need to add extra padding

    my_parameters.Init(90.0, 90.0, 90.0, 0.0, 0.0);
    projection_volume_3d.ExtractSlice(padded_projection_volume_image, my_parameters);
    padded_projection_volume_image.SwapRealSpaceQuadrants( ); // must do this step as image is not centered in the box
    padded_projection_volume_image.BackwardFFT( );
    padded_projection_volume_image.object_is_centred_in_box = true;
    padded_projection_volume_image.ClipInto(&projection_volume_image);

    wxString projection_filename;
    size_t   ext_pos;

    ext_pos             = output_filename.find('.');
    output_filename     = output_filename.substr(0, ext_pos);
    projection_filename = output_filename + "_projection.mrc";
    MRCFile output_projection(projection_filename.ToStdString( ), true);
    projection_volume_image.WriteSlice(&output_projection, 1);

    //my_slice.QuickAndDirtyWriteSlice("my_slice.mrc", 1, my_input_filename.ReturnZSize( ));

    // Image Zaverage;
    // Zaverage.Allocate(my_input_filename.ReturnXSize( ), my_input_filename.ReturnYSize( ), true);
    // Zaverage.SetToConstant(0.0);

    // for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
    //     my_image.ReadSlice(&my_input_filename, image_counter + 1);
    //     my_image.AverageRotationally( );
    //     Zaverage.AddImage(&my_image);
    // }
    // Zaverage.DivideByConstant(number_of_input_images);

    // for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
    //     Zaverage.QuickAndDirtyWriteSlice("my_rotational_average_then_along_z_volume.mrc", image_counter + 1);
    // }

    wxPrintf("\nHelical rotational average program ended successfully!\n");

    return true;
}

void AverageAlongZ(Image* current_volume) {

    // Ensure image is in real space
    bool input_was_fourier = false;
    if ( ! current_volume->is_in_real_space ) {
        current_volume->BackwardFFT( );
        input_was_fourier = true;
    }

    long nx = current_volume->logical_x_dimension;
    long ny = current_volume->logical_y_dimension;
    long nz = current_volume->logical_z_dimension;

    long slice_size = nx * ny;

    // allocate 2D buffer
    std::vector<float> avg(slice_size, 0.0f);

    //----------------------------------------------------------------------
    // 1. Compute sum across all z slices
    //----------------------------------------------------------------------

    long index = 0;
    for ( long z = 0; z < nz; z++ ) {
        long row_index = 0;

        for ( long y = 0; y < ny; y++ ) {
            for ( long x = 0; x < nx; x++ ) {
                avg[row_index] += current_volume->real_values[index];
                row_index++;
                index++;
            }
            index += current_volume->padding_jump_value;
        }
    }

    // Convert sum → mean
    for ( long i = 0; i < slice_size; i++ )
        avg[i] /= float(nz);

    //----------------------------------------------------------------------
    // 2. Copy averaged slice back into every z plane
    //----------------------------------------------------------------------

    index = 0;
    for ( long z = 0; z < nz; z++ ) {
        long row_index = 0;

        for ( long y = 0; y < ny; y++ ) {
            for ( long x = 0; x < nx; x++ ) {
                current_volume->real_values[index] = avg[row_index];
                row_index++;
                index++;
            }
            index += current_volume->padding_jump_value;
        }
    }

    if ( input_was_fourier )
        current_volume->ForwardFFT( );
}

void AverageRotationallyPerSlice(Image* volume) {
    // image must be in real space
    bool input_image_in_fourier_space = false;
    if ( ! volume->is_in_real_space ) {
        volume->BackwardFFT( );
        input_image_in_fourier_space = true;
    }
    // max radius in real space is sqrt(2)*0.5*logical_dimension
    long number_of_rings = volume->logical_x_dimension;
    //float edge_value = current_image->ReturnAverageOfRealValues(std::min(current_image->physical_address_of_box_center_x - 2, current_image->physical_address_of_box_center_y - 2), true);
    float edge_value;
    auto  ring_axis   = std::make_unique<float[]>(number_of_rings); // size is the number of elements
    auto  ring_values = std::make_unique<float[]>(number_of_rings);
    auto  ring_weight = std::make_unique<float[]>(number_of_rings);
    ZeroArray(ring_values.get( ), number_of_rings);
    ZeroArray(ring_weight.get( ), number_of_rings);

    long central_x_pixel = volume->physical_address_of_box_center_x;
    long central_y_pixel = volume->physical_address_of_box_center_y;
    // Add third dimension
    long central_z_pixel = 0;

    double radius;
    double difference;
    long   index_of_bin;
    long   counter;

    // intialize values and weights (number of bins run from 0 to N-1)
    // TODO: remove this comment
    // number_of_rings == dimensions of images in original stack
    // Divide by number of rings - 1 because of 0 indexing
    for ( counter = 0; counter < number_of_rings; counter++ ) {
        ring_axis[counter] = 0.0 + counter * ((sqrt(pow(volume->physical_address_of_box_center_x, 2) + pow(volume->physical_address_of_box_center_y, 2))) - 0.0) / float(number_of_rings - 1); // Diagonal because we're using the physical address (0 is the upper left corner)
        // = counter * (sqrt(pow(physical_address_of_box_center_x, 2) + pow(physical_address_of_box_center_y, 2) + pow(physical_address_of_box_center_z, 2)));
    }

    // edge radius in real space is 0.5*logical_dimension
    long edge_bin = long((0.5 * volume->logical_x_dimension - ring_axis[0]) / (ring_axis[1] - ring_axis[0]));

    // now go through and work out the average;

    counter = 0;

    // z-dim first; that's each image in the volume (looking "through" the image rather than at the side)
    // Nested loops mean move across the image starting at the first image (z), the first row of pixels, going column by column (left to right across the image, going from top to bottom)
    for ( long z = 0; z < 1; z++ ) {
        for ( long y = 0; y < volume->logical_y_dimension; y++ ) {
            for ( long x = 0; x < volume->logical_x_dimension; x++ ) {
                radius       = sqrtf(powf(double(central_x_pixel - x), 2) + powf(double(central_y_pixel - y), 2)); // Here we're moving ring by ring closer to the corner of the image
                index_of_bin = long((radius - ring_axis[0]) / (ring_axis[1] - ring_axis[0]));

                // Only if the index of the current bin index is greater than that of the edge bin, add the current ring_value to the voxel (real_values array) and give it weight
                if ( index_of_bin >= edge_bin ) {
                    ring_values[index_of_bin] += volume->real_values[counter];
                    //ring_values[index_of_bin] += edge_value;
                    ring_weight[index_of_bin] += 1;
                }
                // Otherwise, calculate the difference between the radius and ring_axis (I guess this is finding the "true" radius of the ring?)
                else {
                    // Determines the weight of each ring that's less than the maximum of 1?
                    difference = (radius - ring_axis[index_of_bin]) / (ring_axis[index_of_bin + 1] - ring_axis[index_of_bin]);

                    ring_values[index_of_bin] += volume->real_values[counter] * (1 - difference);
                    ring_values[index_of_bin + 1] += volume->real_values[counter] * difference;
                    ring_weight[index_of_bin] += (1 - difference);
                    ring_weight[index_of_bin + 1] += difference;
                }

                counter++;
            }
            counter += volume->padding_jump_value;
        }
    }
    // divide by number of members...

    for ( counter = 0; counter < number_of_rings; counter++ ) {
        if ( ring_weight[counter] != 0.0 )
            ring_values[counter] /= ring_weight[counter];
    }

    // put the data back into the image

    counter = 0;

    for ( long z = 0; z < volume->logical_z_dimension; z++ ) {
        float z_radius_squared = powf(float(central_z_pixel - z), 2);
        for ( long y = 0; y < volume->logical_y_dimension; y++ ) {
            float y_radius_squared = powf(float(central_y_pixel - y), 2); // same
            for ( long x = 0; x < volume->logical_x_dimension; x++ ) {
                float x_radius_squared = powf(float(central_x_pixel - x), 2);
                radius                 = sqrtf(y_radius_squared + x_radius_squared);
                index_of_bin           = long((radius - ring_axis[0]) / (ring_axis[1] - ring_axis[0]));

                if ( index_of_bin >= edge_bin ) {
                    // set corner values to average at edge
                    volume->real_values[counter] = ring_values[edge_bin - 1];
                    //volume->real_values[counter] = edge_value;
                }
                else {
                    difference                   = (radius - ring_axis[index_of_bin]) / (ring_axis[index_of_bin + 1] - ring_axis[index_of_bin]);
                    volume->real_values[counter] = (ring_values[index_of_bin] * (1 - difference)) + (ring_values[index_of_bin + 1] * difference);
                }

                counter++;
            }
            counter += volume->padding_jump_value;
        }
    }

    // All below here moved to DoCalculation in azimuthal_average
    // padding radial average with edge values
    //edge_value = volume->ReturnAverageOfRealValuesOnEdges( );
    //volume->Resize(volume->logical_x_dimension, volume->logical_y_dimension, volume->logical_z_dimension, edge_value);

    // put the data back into the image and padding with helix in the z-direction
    /*
    long pixel_coord_xy  = 0;
    long pixel_coord_xyz = 0;
    counter              = 0;

    // volume only exists for putting info into the volume
    for ( z = 0; z < volume->logical_z_dimension; z++ ) {
        for ( y = 0; y < volume->logical_y_dimension; y++ ) {
            for ( x = 0; x < volume->logical_x_dimension; x++ ) {
                pixel_coord_xy = volume->ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
                //pixel_coord_xyz = current_volume->ReturnReal1DAddressFromPhysicalCoord(x, y, z);
                //current_volume->real_values[pixel_coord_xyz] = current_image->real_values[pixel_coord_xy];
                volume->real_values[counter] = volume->real_values[pixel_coord_xy];
                counter++;
            }
            counter += volume->padding_jump_value;
        }
    }*/

    if ( input_image_in_fourier_space )
        volume->ForwardFFT( );
}
