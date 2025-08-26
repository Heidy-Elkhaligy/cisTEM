#include "../../core/core_headers.h"

class
        create_sphere_mask : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

void create_white_sphere_mask(Image* mask_file, int x_sphere_center, int y_sphere_center, int z_spehere_center, float radius);
void create_black_sphere_mask(Image* mask_file, int x_sphere_center, int y_sphere_center, int z_spehere_center, float radius);
IMPLEMENT_APP(create_sphere_mask)

// override the DoInteractiveUserInput

void create_sphere_mask::DoInteractiveUserInput( ) {

    wxString output_filename;
    int      box_size;
    float    pixel_size;
    int      x_mask_center       = 1;
    int      y_mask_center       = 1;
    int      z_mask_center       = 1;
    int      sphere_mask_radius  = 1;
    bool     create_white_sphere = true;

    UserInput* my_input = new UserInput("create_sphere_mask", 1.0);
    output_filename     = my_input->GetFilenameFromUser("Output mask file name", "Name of the output mask file", "my_mask.mrc", false);
    box_size            = my_input->GetIntFromUser("Box size of the mask (pixels)", "The box size of the sphere mask, 1 = default cubic box size is 256 pixels", "1", 1);
    pixel_size          = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
    x_mask_center       = my_input->GetIntFromUser("X center of the mask (pixels)", "The X center of the provided mask or the created mask by the program, 1 = default is 0.75 * box size", "1", 1);
    y_mask_center       = my_input->GetIntFromUser("Y center pf the mask (pixels)", "The Y center of the provided mask or the created mask by the program, 1 = default is 0.5 * box size", "1", 1);
    z_mask_center       = my_input->GetIntFromUser("Z center of the mask (pixels)", "The Z center of the provided mask or the created mask by the program, 1 = default is 0.5 * box size", "1", 1);
    sphere_mask_radius  = my_input->GetIntFromUser("Sphere mask radius (pixels)", "The radius of the provided mask or the created mask by the program, 1 = default is 0.3 * box size", "1", 1);
    create_white_sphere = my_input->GetYesNoFromUser("Create a white sphere mask?", "If, Yes a white sphere mask inside the box will be create. Otherwise, a black sphere mask on white background will be created", "Yes");

    delete my_input;

    my_current_job.Reset(10); // + 1 the number shown in the string below
    my_current_job.ManualSetArguments("tifiiiib", output_filename.ToUTF8( ).data( ),
                                      box_size,
                                      pixel_size,
                                      x_mask_center,
                                      y_mask_center,
                                      z_mask_center,
                                      sphere_mask_radius,
                                      create_white_sphere);
}

// override the do calculation method which will be what is actually run..

bool create_sphere_mask::DoCalculation( ) {
    // get the arguments for this job..
    wxString output_filename     = my_current_job.arguments[0].ReturnStringArgument( );
    int      box_size            = my_current_job.arguments[1].ReturnIntegerArgument( );
    float    pixel_size          = my_current_job.arguments[2].ReturnFloatArgument( );
    int      x_mask_center       = my_current_job.arguments[3].ReturnIntegerArgument( );
    int      y_mask_center       = my_current_job.arguments[4].ReturnIntegerArgument( );
    int      z_mask_center       = my_current_job.arguments[5].ReturnIntegerArgument( );
    int      sphere_mask_radius  = my_current_job.arguments[6].ReturnIntegerArgument( );
    bool     create_white_sphere = my_current_job.arguments[7].ReturnBoolArgument( );
    // initiate I/O variables
    MRCFile my_output_filename(output_filename.ToStdString( ), true);

    if ( box_size == 1 ) {
        box_size = 256;
    }
    if ( x_mask_center == 1 ) {
        x_mask_center = 0.75 * box_size;
    }
    if ( y_mask_center == 1 ) {
        y_mask_center = 0.5 * box_size;
    }
    if ( z_mask_center == 1 ) {
        z_mask_center = 0.5 * box_size;
    }
    if ( sphere_mask_radius == 1 ) {
        sphere_mask_radius = 0.25 * box_size;
    }

    Image my_mask_smooth;
    my_mask_smooth.Allocate(box_size, box_size, box_size, true);
    my_mask_smooth.SetToConstant(0.0);
    if ( create_white_sphere ) {
        create_white_sphere_mask(&my_mask_smooth, x_mask_center, y_mask_center, z_mask_center, sphere_mask_radius);
    }
    else {
        create_black_sphere_mask(&my_mask_smooth, x_mask_center, y_mask_center, z_mask_center, sphere_mask_radius);
    }
    my_mask_smooth.ForwardFFT( );
    my_mask_smooth.GaussianLowPassFilter((pixel_size * 2) / 150);
    my_mask_smooth.BackwardFFT( );
    my_output_filename.my_header.SetDimensionsVolume(box_size, box_size, box_size);
    my_output_filename.my_header.SetPixelSize(pixel_size);
    my_mask_smooth.WriteSlices(&my_output_filename, 1, box_size);

    return true;
}

void create_white_sphere_mask(Image* mask_file, int x_sphere_center, int y_sphere_center, int z_spehere_center, float radius) {

    int   boxsize = mask_file->logical_x_dimension;
    int   i, j, k;
    int   dx, dy, dz;
    float d;
    // initialize the mask to 0.0
    mask_file->SetToConstant(0.0);

    long pixel_counter = 0;
    for ( k = 0; k < mask_file->logical_z_dimension; k++ ) {
        for ( j = 0; j < mask_file->logical_y_dimension; j++ ) {
            for ( i = 0; i < mask_file->logical_x_dimension; i++ ) {
                dx = i - x_sphere_center;
                dy = j - y_sphere_center;
                dz = k - z_spehere_center;
                d  = sqrtf(dx * dx + dy * dy + dz * dz);
                if ( d < radius ) {
                    mask_file->real_values[pixel_counter] = 1.0;
                }
                else {
                    mask_file->real_values[pixel_counter] = 0.0;
                }
                pixel_counter++;
            }
            pixel_counter += mask_file->padding_jump_value;
        }
    }
    //mask_file->QuickAndDirtyWriteSlices("make_white_sphere_mask_inside_function.mrc", 1, mask_file->logical_z_dimension);
}

void create_black_sphere_mask(Image* mask_file, int x_sphere_center, int y_sphere_center, int z_spehere_center, float radius) {

    int   boxsize = mask_file->logical_x_dimension;
    int   i, j, k;
    int   dx, dy, dz;
    float d;
    // initialize the mask file to be 1.0
    mask_file->SetToConstant(1.0);

    long pixel_counter = 0;

    for ( k = 0; k < mask_file->logical_z_dimension; k++ ) {
        for ( j = 0; j < mask_file->logical_y_dimension; j++ ) {
            for ( i = 0; i < mask_file->logical_x_dimension; i++ ) {
                dx = i - x_sphere_center;
                dy = j - y_sphere_center;
                dz = k - z_spehere_center;
                d  = sqrtf(dx * dx + dy * dy + dz * dz);
                if ( d < radius ) {
                    mask_file->real_values[pixel_counter] = 0.0;
                }
                else {
                    mask_file->real_values[pixel_counter] = 1.0;
                }
                pixel_counter++;
            }
            pixel_counter += mask_file->padding_jump_value;
        }
    }
}
