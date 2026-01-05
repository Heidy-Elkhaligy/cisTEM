#include "../../core/core_headers.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip> // for std::fixed and std::setprecision
#include <memory>

class
        align_classaverage_tubes : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

typedef struct ctf_parameters {
    float acceleration_voltage; // keV
    float spherical_aberration; // mm
    float amplitude_contrast;
    float defocus_1; // A
    float defocus_2; // A
    float astigmatism_angle; // degrees
    float lowest_frequency_for_fitting; // 1/A
    float highest_frequency_for_fitting; // 1/A
    float astigmatism_tolerance; // A
    float pixel_size; // A
    float additional_phase_shift; // rad
} ctf_parameters;

std::vector<float>  sum_image_columns(Image* current_image);
float               max_abs_column_sum(Image* current_image);
std::pair<int, int> FindOuterTubeEdges(const std::vector<float>& cols, float min_tube_diameter, float max_tube_diameter);
void                sum_image_direction(Image* current_image, int dim);
//Functions for the average images bins
void InitializeCTFSumOfSquares(Image& current_image, std::vector<float>& ctf_sum_of_squares);
void ApplyCTFAndReturnCTFSumOfSquares(Image& image, CTF ctf_to_apply, bool absolute, bool apply_beam_tilt, bool apply_envelope, std::vector<float>& ctf_sum_of_squares);
void divide_by_ctf_sum_of_squares(Image& current_image, std::vector<float>& ctf_sum_of_squares);

IMPLEMENT_APP(align_classaverage_tubes)

// override the DoInteractiveUserInput

void align_classaverage_tubes::DoInteractiveUserInput( ) {
    wxString input_images;
    wxString output_images;
    bool     class_average = true;
    wxString input_star_filename;
    float    acceleration_voltage;
    float    spherical_aberration;
    float    amplitude_contrast;
    float    pixel_size;
    float    min_tube_diameter   = 0.0;
    float    max_tube_diameter   = 0.0;
    float    outer_mask_radius   = 0;
    float    low_pass_resolution = 50.0;
    bool     use_auto_corr;

    // bool     save_new_angles = false;
    // wxString angles_filename;
    // bool     create_average_power_spectrum;
    // wxString average_power_spectrum_filename;
    bool use_memory = false;
    int  max_threads;

    UserInput* my_input = new UserInput("align_classaverage_tubes", 1.00);
    input_images        = my_input->GetFilenameFromUser("Input images file name", "Filename of helical tube stack", "helical_stack.mrc", true);
    output_images       = my_input->GetFilenameFromUser("Output images", "The top down 2D tubes images", "aligned_helical_stack.mrc", false);
    class_average       = my_input->GetYesNoFromUser("Input stack contains class average images?", "Default input is class average images. However, if this question is set to No, it will assume input stack is normal images and will add an extra cross-correlation step to determine correct psi angle", "Yes");
    if ( ! class_average ) {
        input_star_filename  = my_input->GetFilenameFromUser("Input star file", "The input star file", "my_parameters.star", true);
        acceleration_voltage = my_input->GetFloatFromUser("Acceleration voltage (keV)", "Acceleration voltage, in keV", "300.0", 0.0, 500.0);
        spherical_aberration = my_input->GetFloatFromUser("Spherical aberration (mm)", "Objective lens spherical aberration", "2.7", 0.0);
        amplitude_contrast   = my_input->GetFloatFromUser("Amplitude contrast", "Fraction of total contrast attributed to amplitude contrast", "0.07", 0.0);
    }
    pixel_size          = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
    min_tube_diameter   = my_input->GetFloatFromUser("Minimum tube diameter", "The minimum tube diameter for searching and bining tubes", "30.0", 0.0);
    max_tube_diameter   = my_input->GetFloatFromUser("Maximum tube diameter", "The maximum tube diameter for searching and bining tubes", "60.0", 0.0);
    outer_mask_radius   = my_input->GetFloatFromUser("Outer mask radius for masking the images during tube alignment (pixels)", "Outer mask radius to use when searching and aligning tubes in pixels, zero mean no masking should be applied", "0", 0);
    low_pass_resolution = my_input->GetFloatFromUser("Resolution limit for low pass filtering", "Resolution limit for low pass filter. Only this resolution or worse information will be retained", "50.0", 0.0);
    use_auto_corr       = my_input->GetYesNoFromUser("Use auto correlation to find helix axis psi rotation?", "If yes, will use auto correlation to find helix axis psi rotation. If no, FT will be used.", "NO");

    // save_new_angles = my_input->GetYesNoFromUser("Do you want to save the in-plane rotation angles used to align images?", "If yes, a text file will be generated with the psi angles", "NO");
    // if ( save_new_angles ) {
    //     angles_filename = my_input->GetStringFromUser("Filename for the saved angles file", "Filename for the in-plane rotation of the input images and rotation applied to align them vertically", "psi_angles.txt");
    // }
    // create_average_power_spectrum = my_input->GetYesNoFromUser("Do you want to create an average power spectrum mrc file?", "If yes, average power spectrum will be saved to an mrc file", "NO");
    // if ( create_average_power_spectrum == true ) {
    //     average_power_spectrum_filename = my_input->GetFilenameFromUser("Output average power spectrum file name ", "The output average power spectrum file", "average_power_spectrum_image.mrc", false);
    // }
    use_memory = my_input->GetYesNoFromUser("Allocate images to memory?", "Choice between memory allocation or using functions; no is recommended for systems with limited memory.", "NO");

#ifdef _OPENMP
    max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "when threading, what is the max threads to run", "1", 1);
#else
    max_threads = 1;
#endif

    delete my_input;

    my_current_job.Reset(16);
    my_current_job.ManualSetArguments("ttbtffffffffbbi", input_images.ToUTF8( ).data( ),
                                      output_images.ToUTF8( ).data( ),
                                      class_average,
                                      input_star_filename.ToUTF8( ).data( ),
                                      acceleration_voltage,
                                      spherical_aberration,
                                      amplitude_contrast,
                                      pixel_size,
                                      min_tube_diameter,
                                      max_tube_diameter,
                                      outer_mask_radius,
                                      low_pass_resolution,
                                      use_auto_corr,
                                      use_memory,
                                      max_threads); //update_star_file, input_star_filename.ToUTF8( ).data( ),
}

// override the do calculation method which will be what is actually run..

bool align_classaverage_tubes::DoCalculation( ) {
    wxString input_images         = my_current_job.arguments[0].ReturnStringArgument( );
    wxString output_images        = my_current_job.arguments[1].ReturnStringArgument( );
    bool     class_average        = my_current_job.arguments[2].ReturnBoolArgument( );
    wxString input_star_filename  = my_current_job.arguments[3].ReturnStringArgument( );
    float    acceleration_voltage = my_current_job.arguments[4].ReturnFloatArgument( );
    float    spherical_aberration = my_current_job.arguments[5].ReturnFloatArgument( );
    float    amplitude_contrast   = my_current_job.arguments[6].ReturnFloatArgument( );
    float    pixel_size           = my_current_job.arguments[7].ReturnFloatArgument( );
    float    min_tube_diameter    = my_current_job.arguments[8].ReturnFloatArgument( );
    float    max_tube_diameter    = my_current_job.arguments[9].ReturnFloatArgument( );
    float    outer_mask_radius    = my_current_job.arguments[10].ReturnFloatArgument( );
    float    low_pass_resolution  = my_current_job.arguments[11].ReturnFloatArgument( );
    bool     use_auto_corr        = my_current_job.arguments[12].ReturnBoolArgument( );
    bool     use_memory           = my_current_job.arguments[13].ReturnBoolArgument( );
    int      max_threads          = my_current_job.arguments[14].ReturnIntegerArgument( );

    MRCFile my_input_images(input_images.ToStdString( ), false);
    MRCFile my_output_images(output_images.ToStdString( ), true);
    // MRCFile* my_output_power_spectrum;
    // if ( create_average_power_spectrum ) {
    //     my_output_power_spectrum = new MRCFile(average_power_spectrum_filename.ToStdString( ), true, true);
    // }
    long number_of_input_images = my_input_images.ReturnNumberOfSlices( );

    ctf_parameters* ctf_parameters_stack;
    if ( ! class_average ) {
        ctf_parameters_stack = new ctf_parameters[number_of_input_images];
        //cisTEM star
        cisTEMParameters input_star_file;
        if ( (is_running_locally && ! DoesFileExist(input_star_filename.ToStdString( ))) ) {
            SendErrorAndCrash(wxString::Format("Error: Input star file %s not found\n", input_star_filename));
        }
        //CisTEM star
        input_star_file.ReadFromcisTEMStarFile(input_star_filename.ToStdString( ));
        for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
            ctf_parameters_stack[image_counter].acceleration_voltage          = acceleration_voltage;
            ctf_parameters_stack[image_counter].spherical_aberration          = spherical_aberration;
            ctf_parameters_stack[image_counter].amplitude_contrast            = amplitude_contrast;
            ctf_parameters_stack[image_counter].defocus_1                     = input_star_file.ReturnDefocus1(image_counter);
            ctf_parameters_stack[image_counter].defocus_2                     = input_star_file.ReturnDefocus2(image_counter);
            ctf_parameters_stack[image_counter].astigmatism_angle             = input_star_file.ReturnDefocusAngle(image_counter);
            ctf_parameters_stack[image_counter].lowest_frequency_for_fitting  = 0.0;
            ctf_parameters_stack[image_counter].highest_frequency_for_fitting = 0.5;
            ctf_parameters_stack[image_counter].astigmatism_tolerance         = 0.0;
            ctf_parameters_stack[image_counter].pixel_size                    = pixel_size;
            ctf_parameters_stack[image_counter].additional_phase_shift        = input_star_file.ReturnPhaseShift(image_counter);
            //wxPrintf("The current image is %li and its defocus is %f\n", image_counter+1, input_star_file.ReturnDefocus1(image_counter) );
        }
    }

    std::vector<float> tube_rotation(number_of_input_images, 0.0f);
    std::vector<float> best_sum_column(number_of_input_images, 0.0f);
    std::vector<float> x_shift_column(number_of_input_images, 0.0f);
    Image              current_image;
    Image              temp_image;
    Image              final_image;

    int x_dim;
    int y_dim;
    current_image.ReadSlice(&my_input_images, 1);
    float center_peak_index = current_image.logical_y_dimension / 2;
    x_dim                   = current_image.logical_x_dimension;
    y_dim                   = current_image.logical_y_dimension;

    // initiate default parameters for the ApplyCTFAndReturnCTFSumOfSquares function
    // (May be change that later to be expert options inputs???)
    bool absolute        = false;
    bool apply_beam_tilt = false;
    bool apply_envelope  = false;

    // CTF object
    CTF current_ctf;

    auto CTFSumOfSquares = std::make_shared<std::vector<float>>( );

    auto CTFSumOfSquaresFinal = std::make_shared<std::vector<float>>( );

    InitializeCTFSumOfSquares(current_image, *CTFSumOfSquares);
    InitializeCTFSumOfSquares(current_image, *CTFSumOfSquaresFinal);

    std::vector<std::vector<float>> all_columns_sum(number_of_input_images, std::vector<float>(x_dim, 0.0)); // why I am saving those values?????

    // parameters for masking the FT image
    float cosine_edge       = 10.0;
    float outside_weight    = 0.0;
    float filter_radius     = 0.0;
    float outside_value     = 0.0;
    bool  use_outside_value = false;
    long  image_counter;

    // rotation (degrees)
    float psi_min  = 0.0;
    float psi_max  = 360.0;
    float psi_step = 4.0;

    // input stack low pass filtered and masked
    Image* image_stack_filtered_masked;
    if ( use_memory )
        image_stack_filtered_masked = new Image[number_of_input_images];
    else
        image_stack_filtered_masked = nullptr;

    if ( use_memory ) {
        wxPrintf("\nLoading images to memory...\n\n");
        ProgressBar* loading_progress = new ProgressBar(number_of_input_images);
#pragma omp parallel num_threads(max_threads) shared(loading_progress, my_input_images, image_stack_filtered_masked, x_dim, y_dim, outer_mask_radius, low_pass_resolution, pixel_size)
        {
#pragma omp for schedule(static) ordered
            for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
// Read from disk
#pragma omp critical
                image_stack_filtered_masked[image_counter].ReadSlice(&my_input_images, image_counter + 1);
                // Normalize the image using cisTEM Normalize
                image_stack_filtered_masked[image_counter].Normalize( );
                // Here the masking is important as we want to only find the rotation of the tubes around the center or near the center
                if ( outer_mask_radius != 0 ) {
                    image_stack_filtered_masked[image_counter].CircleMask(outer_mask_radius);
                }
                // FT the image
                image_stack_filtered_masked[image_counter].ForwardFFT( );
                // convert the central pixel to zero (Is that done in real or Fouriier space??)
                image_stack_filtered_masked[image_counter].ZeroCentralPixel( );

                // will applying a low pass filter here improve finding the correct rotation in FT
                image_stack_filtered_masked[image_counter].GaussianLowPassFilter((pixel_size * 2) / low_pass_resolution);

// if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
#pragma omp ordered
                loading_progress->Update(image_counter + 1);
            }
        }
        delete loading_progress;
    }
    //////////// DEBUGGING
    // MRCFile lazy_output("rotated_images_based_on_initial_psi_after_gaussian_filter.mrc", true);
    // if ( ! lazy_output.IsOpen( ) ) {
    //     lazy_output.OpenFile("rotated_images_based_on_initial_psi_after_gaussian_filter.mrc", true);
    //     if ( ! lazy_output.IsOpen( ) ) {
    //         wxPrintf("ERROR: Could not open '%s' for writing\n", "rotated_images_based_on_initial_psi_after_gaussian_filter.mrc");
    //         DEBUG_ABORT;
    //     }
    // }
    // lazy_output.my_header.SetNumberOfImages(number_of_input_images);
    // lazy_output.my_header.SetDimensionsImage(x_dim, y_dim);
    // lazy_output.SetPixelSize(pixel_size);
    // lazy_output.WriteHeader( );

    wxPrintf("\nFinding psi rotation...\n\n");
    ProgressBar* my_progress = new ProgressBar(number_of_input_images);

#pragma omp parallel for schedule(dynamic, 1) num_threads(max_threads) default(none) shared(my_input_images, image_stack_filtered_masked, use_memory, best_sum_column, tube_rotation, number_of_input_images, low_pass_resolution, max_threads, cosine_edge, outside_weight, filter_radius, outside_value, use_outside_value, x_dim, y_dim, \
                                                                                            psi_step, pixel_size, psi_min, psi_max, my_progress, outer_mask_radius, min_tube_diameter, center_peak_index, x_shift_column, tube_rotation, use_auto_corr, all_columns_sum, max_tube_diameter) private(current_image, image_counter, final_image)

    for ( image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
        //wxPrintf("The image counter for slice reading is %li \n", image_counter + 1);
        // read the current image in the stack
        if ( use_memory ) {
            current_image.CopyFrom(&image_stack_filtered_masked[image_counter]);
        }
        else {
#pragma omp critical
            current_image.ReadSlice(&my_input_images, image_counter + 1);

            // Normalize the image using cisTEM Normalize
            current_image.Normalize( );
            //circle mask the image in real space using the outer mask radius to avoid the noise on the edges
            if ( outer_mask_radius != 0 ) {
                current_image.CircleMask(outer_mask_radius);
            }
            // FT the image
            current_image.ForwardFFT( );
            // convert the central pixel to zero (Is that done in real or Fouriier space??)
            current_image.ZeroCentralPixel( );

            // will applying a low pass filter here improve finding the correct rotation in FT
            current_image.GaussianLowPassFilter((pixel_size * 2) / low_pass_resolution);
        }

        if ( use_auto_corr ) {
            for ( long pixel_counter = 0; pixel_counter < current_image.real_memory_allocated / 2; pixel_counter++ ) {

                // calculating the amplitude is not needed by let's see and print its value
                float amplitude                             = abs(current_image.complex_values[pixel_counter]);
                current_image.complex_values[pixel_counter] = amplitude * amplitude + I * 0.0f;
            }
            // return the image to real space again and save them to see the correlation
            current_image.BackwardFFT( );
            current_image.SwapRealSpaceQuadrants( );
            // set the image is centered inside the box as true
            current_image.object_is_centred_in_box = true;
        }
        Image temp_image;
        Image power_image;
        // Allocate memory for the image that will be rotated
        if ( use_auto_corr ) {
            temp_image.Allocate(x_dim, y_dim, true);
            temp_image.SetToConstant(0.0);
        }
        else {
            temp_image.Allocate(x_dim, y_dim, false);
            temp_image.SetToConstant(0.0);
        }

        float              local_best_sum = -FLT_MAX;
        float              local_best_psi = 0.0f;
        std::vector<float> local_columns_sum_vector;

        for ( float psi = psi_min; psi <= psi_max; psi += psi_step ) {
            if ( use_auto_corr ) {
                temp_image.CopyFrom(&current_image);
                temp_image.Rotate2DInPlace(psi, FLT_MAX);
            }
            else {
                AnglesAndShifts rotation_angle;
                rotation_angle.Init(0.0, 0.0, psi, 0.0, 0.0);
                temp_image.CopyFrom(&current_image);
                temp_image.SwapRealSpaceQuadrants( );
                Image rotated_image;
                rotated_image.Allocate(x_dim, y_dim, false); // the rotated images real values will be the rotated FT
                rotated_image.SetToConstant(0.0);
                temp_image.RotateFourier2D(rotated_image, rotation_angle);
                // allocate memory for power image
                power_image.Allocate(x_dim, y_dim, true);
                power_image.SetToConstant(0.0);
                rotated_image.ComputeAmplitudeSpectrumFull2D(&power_image);
                //temp_image.Rotate2DInPlace(psi, FLT_MAX);
                //temp_image.QuickAndDirtyWriteSlice("temp_power_image.mrc", 1);
                rotated_image.Deallocate( );
            }

            float column_sum;

            // Use a threshold to find the line corresponding to tube axis angle of rotation

            if ( use_auto_corr ) {
                // Filter the autocorrelation image so that any pixel values above average are kept to ensure correct rotation calculation
                // should we make that in the advanced options something to be optimized by user
                Image binary_mask;
                float image_average;
                image_average = temp_image.ReturnAverageOfRealValues( );
                binary_mask.Allocate(x_dim, y_dim, true);
                binary_mask.SetToConstant(0.0);
                binary_mask.CopyFrom(&temp_image);
                binary_mask.Binarise(image_average);
                //binary_mask.QuickAndDirtyWriteSlice("binary_mask.mrc", 1);

                float filter_edge = 40.0;
                temp_image.ApplyMask(binary_mask, cosine_edge / pixel_size, outside_weight, pixel_size / filter_radius, pixel_size / filter_edge, outside_value, use_outside_value);
                binary_mask.Deallocate( );
                // Use a threshold to find the line corresponding to tube axis
                column_sum = max_abs_column_sum(&temp_image);
            }
            else {
                // Use a threshold to find the line corresponding to tube axis angle of rotation
                float image_average;
                float image_sd;
                float image_threshold;
                Image binary_mask;
                image_average = power_image.ReturnAverageOfRealValues( );
                image_sd      = sqrt(power_image.ReturnVarianceOfRealValues( ));
                // Threshold value is 2 sd away from mean to eliminate any outliers
                image_threshold = image_average + (3 * image_sd);
                binary_mask.CopyFrom(&power_image);
                binary_mask.Binarise(image_threshold);
                //binary_mask.QuickAndDirtyWriteSlice("binary_mask.mrc", 1);

                float filter_edge = 40.0;
                power_image.ApplyMask(binary_mask, cosine_edge / pixel_size, outside_weight, pixel_size / filter_radius, pixel_size / filter_edge, outside_value, use_outside_value);
                //power_image.QuickAndDirtyWriteSlice("binarized_power_image.mrc", 1);

                column_sum = max_abs_column_sum(&power_image);
            }

            // save the 1D vector to examine in case there is problems
            std::vector<float> column_sum_vector;
            if ( use_auto_corr ) {
                column_sum_vector = sum_image_columns(&temp_image);
            }
            else {
                column_sum_vector = sum_image_columns(&power_image);
            }

            if ( column_sum > local_best_sum ) {
                // at this point I will not add 90 to the saved psi if using FFT
                local_best_psi           = psi;
                local_best_sum           = column_sum;
                local_columns_sum_vector = column_sum_vector;
                //wxPrintf("The best psi angle with the highest sum is %f with sum %f \n\n", local_best_psi, local_best_sum);
            }
            power_image.Deallocate( ); // should this be here inside the psi loop or after it?? and reset the power_image to 0 constant ???
        }
        temp_image.Deallocate( );
        // tuning the rotation angle to search with 2 degrees before and after the initial angle
        // The step size will be rotation angle (2) /20

        float tuning_rotation_range = psi_step / 2;
        float tuning_step_size      = psi_step / 20;

        // column wise calculations
        float current_best_psi = local_best_psi; // as the tube_rotation[image_counter] is still zero when having local variables!
        //float current_best_psi       = tube_rotation[image_counter];
        float tuning_psi_lower_range = current_best_psi - tuning_rotation_range;
        float tuning_psi_upper_range = current_best_psi + tuning_rotation_range;

        Image tuning_temp_image;
        Image tuning_power_image;

        if ( use_auto_corr ) {
            tuning_temp_image.Allocate(x_dim, y_dim, true);
            tuning_temp_image.SetToConstant(0.0);
        }
        else {
            tuning_temp_image.Allocate(x_dim, y_dim, false);
            tuning_temp_image.SetToConstant(0.0);
        }

        for ( float tuning_psi = tuning_psi_lower_range; tuning_psi <= tuning_psi_upper_range; tuning_psi += tuning_step_size ) {

            if ( use_auto_corr ) {
                tuning_temp_image.CopyFrom(&current_image);
                tuning_temp_image.Rotate2DInPlace(tuning_psi, FLT_MAX);
            }
            else {
                AnglesAndShifts tuning_rotation_angle;
                tuning_rotation_angle.Init(0.0, 0.0, tuning_psi, 0.0, 0.0);

                Image tuning_rotated_image;
                tuning_rotated_image.Allocate(x_dim, y_dim, false);
                tuning_rotated_image.SetToConstant(0.0);

                tuning_temp_image.CopyFrom(&current_image);
                tuning_temp_image.SwapRealSpaceQuadrants( );
                tuning_temp_image.RotateFourier2D(tuning_rotated_image, tuning_rotation_angle);

                // allocate memory for power image
                tuning_power_image.Allocate(x_dim, y_dim, true);
                tuning_power_image.SetToConstant(0.0);

                tuning_rotated_image.ComputeAmplitudeSpectrumFull2D(&tuning_power_image);
                tuning_rotated_image.Deallocate( );
            }
            float tuning_column_sum;
            float tuning_image_average;
            float tuning_image_sd;
            float tuning_image_threshold;
            Image tuning_binary_mask;

            if ( use_auto_corr ) {
                // Filter the autocorrelation image so that any pixel values above average are kept to ensure correct rotation calculation
                // should we make that in the advanced options something to be optimized by user

                tuning_image_average = tuning_temp_image.ReturnAverageOfRealValues( );
                tuning_binary_mask.Allocate(x_dim, y_dim, true);
                tuning_binary_mask.SetToConstant(0.0);
                tuning_binary_mask.CopyFrom(&tuning_temp_image);
                tuning_binary_mask.Binarise(tuning_image_average);
                //binary_mask.QuickAndDirtyWriteSlice("binary_mask.mrc", 1);

                float filter_edge = 40.0;
                tuning_temp_image.ApplyMask(tuning_binary_mask, cosine_edge / pixel_size, outside_weight, pixel_size / filter_radius, pixel_size / filter_edge, outside_value, use_outside_value);
                tuning_binary_mask.Deallocate( );
                // Use a threshold to find the line corresponding to tube axis
                tuning_column_sum = max_abs_column_sum(&tuning_temp_image);
            }
            else {
                // Use a threshold to find the line corresponding to tube axis angle of rotation
                tuning_image_average = tuning_power_image.ReturnAverageOfRealValues( );
                tuning_image_sd      = sqrt(tuning_power_image.ReturnVarianceOfRealValues( ));
                // Threshold value is 2 sd away from mean to eliminate any outliers
                tuning_image_threshold = tuning_image_average + (2 * tuning_image_sd);
                tuning_binary_mask.CopyFrom(&tuning_power_image);
                tuning_binary_mask.Binarise(tuning_image_threshold);
                //binary_mask.QuickAndDirtyWriteSlice("binary_tuning_image.mrc", 1);
                float filter_edge = 40.0;
                tuning_power_image.ApplyMask(tuning_binary_mask, cosine_edge / pixel_size, outside_weight, pixel_size / filter_radius, pixel_size / filter_edge, outside_value, use_outside_value);
                tuning_column_sum = max_abs_column_sum(&tuning_power_image);
            }

            std::vector<float> tuning_column_sum_vector;

            if ( use_auto_corr ) {
                tuning_column_sum_vector = sum_image_columns(&tuning_temp_image);
            }
            else {
                tuning_column_sum_vector = sum_image_columns(&tuning_power_image);
            }
            if ( tuning_column_sum > local_best_sum ) {
                // at this point I will not add 90 to the saved psi if using FFT
                local_best_psi           = tuning_psi;
                local_best_sum           = tuning_column_sum;
                local_columns_sum_vector = tuning_column_sum_vector;
                //wxPrintf("Tuning best psi angle with the highest sum is %f with sum %f \n\n", local_best_psi, local_best_sum);
            }
            tuning_power_image.Deallocate( ); // should this be here or outside the tuning psi loop?? and reset the image to zero inside the loop?
        }
        tuning_temp_image.Deallocate( );

        //wxPrintf("The best psi angle with the highest sum is %f with sum %f \n\n", tube_rotation[image_counter], best_sum[image_counter]);
        //wxPrintf("The best psi angle with the highest sum from the column sum is %f with sum %f \n\n", tube_rotation[image_counter], best_sum_column[image_counter]);
        tube_rotation[image_counter]   = local_best_psi;
        best_sum_column[image_counter] = local_best_sum;
        all_columns_sum[image_counter] = local_columns_sum_vector;

        // finding the x shift for each image to ensure the sum image is centered
        final_image.Allocate(x_dim, y_dim, true);
        final_image.SetToConstant(0.0);
        // find  the x,y shift
        // ReadSlice requires omp critical to avoid parallel reads, which may lead to the wrong slice being read
        if ( use_memory ) {
            final_image.CopyFrom(&image_stack_filtered_masked[image_counter]);
        }
        else {
#pragma omp critical
            final_image.ReadSlice(&my_input_images, image_counter + 1);
            final_image.Normalize( );
            if ( outer_mask_radius != 0 ) {
                final_image.CircleMask(outer_mask_radius);
            }

            final_image.ForwardFFT( );
            final_image.ZeroCentralPixel( );
            final_image.GaussianLowPassFilter((pixel_size * 2) / low_pass_resolution);
        }
        final_image.BackwardFFT( );

        // to apply Gaussian filter you need to be in Fourier space
        // Nyquist frequency value is the pixel size * 2
        // if we want to apply a gaussian pass filter that will make the image at 150 angestrom to be well smoothened and get better peaks
        // It should be pixel_size/resolution limit ?? or pixel_size * 2
        if ( use_auto_corr ) {
            final_image.Rotate2DInPlace(local_best_psi, FLT_MAX);
        }
        else {
            final_image.Rotate2DInPlace(local_best_psi + 90.0, FLT_MAX);
        }
        all_columns_sum[image_counter] = sum_image_columns(&final_image);
        // calculate the required x shift to center the tubes
        //auto [peak_one_column_sum, peak_two_column_sum] = find_column_sum_peaks(&final_image, min_tube_diameter);
        // New sum peaks will find the outer edges of tubes not the inner ones
        auto [peak_one_column_sum, peak_two_column_sum] = FindOuterTubeEdges(all_columns_sum[image_counter], min_tube_diameter, max_tube_diameter);

        // The next line not needed
        float tube_center_column_sum          = std::abs(peak_one_column_sum - peak_two_column_sum) / 2;
        float distance_from_center_column_sum = -((peak_one_column_sum + peak_two_column_sum) / 2 - center_peak_index);

        x_shift_column[image_counter] = distance_from_center_column_sum; // the x-shift needed

        //final_image.PhaseShift(x_shift_column[image_counter], 0.0, 0.0);

        final_image.Deallocate( );

        if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
            my_progress->Update(image_counter + 1);
    }

    //final_image.Deallocate( );
    delete my_progress;
    // save column sum output file
    //save_all_columns_sum_to_file(all_columns_sum, "column_sums_output.txt");

    // If class_average is set to No or false, this means the input stack is normal images not class average images
    // So will have low contrast that's why we need to add an extra cross-correlation step to improve finding psi

    std::vector<float> best_correlation_score(number_of_input_images, -FLT_MAX);
    std::vector<float> best_psi_value(number_of_input_images, 0.0f);
    std::vector<float> best_x_shift_value(number_of_input_images, 0.0f);

    if ( ! class_average ) {

        Image added_image; // This is the sum image based on the initial rotation angle calculated from auto-correlation
                // Initialize the sum images based on the number specified by the user
        Image sum_images;

        sum_images.Allocate(x_dim, y_dim, true); //allocate in real space
        sum_images.SetToConstant(0.0);

        wxPrintf("\nCreating Initial Sum Images...\n\n");
        ProgressBar* sum_progress = new ProgressBar(number_of_input_images);

        for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
            if ( use_memory ) {
                added_image.CopyFrom(&image_stack_filtered_masked[image_counter]); // no need for counter + 1 anymore
                added_image.BackwardFFT( );
            }
            else {
                added_image.ReadSlice(&my_input_images, image_counter + 1);
            }
            added_image.Normalize( );
            //added_image.QuickAndDirtyWriteSlice("added_image_after_normalization.mrc", image_counter +1);
            current_ctf.Init(ctf_parameters_stack[image_counter].acceleration_voltage, ctf_parameters_stack[image_counter].spherical_aberration, ctf_parameters_stack[image_counter].amplitude_contrast, ctf_parameters_stack[image_counter].defocus_1, ctf_parameters_stack[image_counter].defocus_2, ctf_parameters_stack[image_counter].astigmatism_angle, ctf_parameters_stack[image_counter].lowest_frequency_for_fitting, ctf_parameters_stack[image_counter].highest_frequency_for_fitting, ctf_parameters_stack[image_counter].astigmatism_tolerance, ctf_parameters_stack[image_counter].pixel_size, ctf_parameters_stack[image_counter].additional_phase_shift);

            added_image.ForwardFFT( );
            added_image.ZeroCentralPixel( );

            ApplyCTFAndReturnCTFSumOfSquares(added_image, current_ctf, absolute, apply_beam_tilt, apply_envelope, (*CTFSumOfSquares));

            added_image.BackwardFFT( );
            // Here will adjust the angles based on the method used to find tube rotation
            if ( use_auto_corr ) {
                added_image.Rotate2DInPlace(tube_rotation[image_counter], FLT_MAX);
            }
            else {
                added_image.Rotate2DInPlace(tube_rotation[image_counter] + 90.0, FLT_MAX);
            }

            added_image.PhaseShift(x_shift_column[image_counter], 0.0, 0.0);

            sum_images.AddImage(&added_image);

            if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
                sum_progress->Update(image_counter + 1);
        }

        delete sum_progress;

        // divide the sum image by CTF sum of squares and centering it
        // Vertically summing the image to ensure cross-correlation doesn't correlate by mistake to a wrong area if input images are pre-aligned
        sum_image_direction(&sum_images, 2);
        sum_images.ForwardFFT( );
        divide_by_ctf_sum_of_squares(sum_images, (*CTFSumOfSquares));
        sum_images.BackwardFFT( );
        // shift sum image to the center after padding based on tube peaks
        // This shift is necessary at this point to ensure the cross-correlation shift is centered correctly later
        std::vector<float> column_sum                   = sum_image_columns(&sum_images);
        auto [peak_one_column_sum, peak_two_column_sum] = FindOuterTubeEdges(column_sum, min_tube_diameter, max_tube_diameter);

        float tube_center_column_sum          = std::abs(peak_one_column_sum - peak_two_column_sum) / 2;
        float distance_from_center_column_sum = -((peak_one_column_sum + peak_two_column_sum) / 2 - center_peak_index);
        //to center the sum image
        sum_images.PhaseShift(-distance_from_center_column_sum, 0.0, 0.0); // the x-shift needed to center the sum image

        // delete CTFSumOfSquares;

        // create average_images that are CTF corrected and then average_image copies directly from here
        Image* average_images;
        average_images = new Image[number_of_input_images];
        Image temporary_average_image;
        for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
            temporary_average_image.CopyFrom(&sum_images);

            temporary_average_image.Normalize( );
            temporary_average_image.ForwardFFT( );
            current_ctf.Init(ctf_parameters_stack[image_counter].acceleration_voltage, ctf_parameters_stack[image_counter].spherical_aberration, ctf_parameters_stack[image_counter].amplitude_contrast, ctf_parameters_stack[image_counter].defocus_1, ctf_parameters_stack[image_counter].defocus_2, ctf_parameters_stack[image_counter].astigmatism_angle, ctf_parameters_stack[image_counter].lowest_frequency_for_fitting, ctf_parameters_stack[image_counter].highest_frequency_for_fitting, ctf_parameters_stack[image_counter].astigmatism_tolerance, ctf_parameters_stack[image_counter].pixel_size, ctf_parameters_stack[image_counter].additional_phase_shift);

            temporary_average_image.ApplyCTF(current_ctf);
            temporary_average_image.BackwardFFT( );
            average_images[image_counter].CopyFrom(&temporary_average_image);
            //temporary_average_image.QuickAndDirtyWriteSlice("average_image_with_ctf.mrc", image_counter + 1);
        }

        /// Getting the correct rotation and shift using cross-correlation
        float inner_radius_for_peak_search;
        float outer_radius_for_peak_search;
        inner_radius_for_peak_search = 0.0; // inner radius should be set to 0
        outer_radius_for_peak_search = x_dim / 2;

        Image my_image;
        Image my_image_copy;
        Image my_image_tuned;
        Image average_image;
        Image tuning_average_image;
        float tuned_rotation_range = psi_step; ///2
        float tuned_step_size      = psi_step / 20;
        Image fine_tuning_average_image;

        wxPrintf("\nFinding Tube Rotation Using Cross-Correlation...\n\n");
        ProgressBar* my_aln_progress = new ProgressBar(number_of_input_images);

#pragma omp parallel for schedule(dynamic, 1) num_threads(max_threads) default(none) shared(number_of_input_images, my_input_images, inner_radius_for_peak_search, outer_radius_for_peak_search, low_pass_resolution, x_dim, y_dim, use_memory, average_images, \
                                                                                            best_correlation_score, best_psi_value, best_x_shift_value, psi_step, tube_rotation, outer_mask_radius, use_auto_corr, image_stack_filtered_masked,                 \
                                                                                            max_threads, sum_images, tuned_rotation_range, tuned_step_size, my_aln_progress, x_shift_column, current_image, all_columns_sum, class_average,                     \
                                                                                            current_ctf, ctf_parameters_stack, min_tube_diameter, max_tube_diameter, center_peak_index, pixel_size) private(my_image, average_image, tuning_average_image, final_image, fine_tuning_average_image, my_image_copy, my_image_tuned)

        for ( long aln_image_counter = 0; aln_image_counter < number_of_input_images; aln_image_counter++ ) {

            if ( use_memory ) {
                my_image.CopyFrom(&image_stack_filtered_masked[aln_image_counter]);
            }
            else {
#pragma omp critical
                my_image.ReadSlice(&my_input_images, aln_image_counter + 1);

                my_image.Normalize( );

                if ( outer_mask_radius != 0 ) {
                    my_image.CircleMask(outer_mask_radius);
                }
                else if ( outer_mask_radius == 0 ) {
                    my_image.CircleMask(x_dim * 0.45);
                }
                my_image.ForwardFFT( );
                my_image.ZeroCentralPixel( );
                //testing adding a low pass filter on the original image before getting the correct shift from the correlation and how that can affect the centering of the mask at the end

                my_image.GaussianLowPassFilter((pixel_size * 2) / low_pass_resolution); //150
            }
            my_image.BackwardFFT( );
            //my_image.QuickAndDirtyWriteSlice("low_pass_filtered_image_for_comparison.mrc", aln_image_counter + 1);

            // initial angle search will start from the rotation angle we got from the auto-correlation/FT
            // then change the auto-correlation angle to be within 180
            // do another search within +/- 90 degrees of the auto-correlation psi angle or the FT psi angle (+90)
            // only at this point we need to adjust the angle before cross-correlation but later it will be already correct and no further adjustments
            float local_best_corr_score = -FLT_MAX;
            float local_best_psi        = 0.0f;
            float local_best_x_shift    = 0.0f;
            float local_best_y_shift    = 0.0f;

            float current_best_psi;
            if ( use_auto_corr ) {
                current_best_psi = tube_rotation[aln_image_counter];
            }
            else {
                current_best_psi = tube_rotation[aln_image_counter] + 90.0;
            }

            float angle_range   = 180.0;
            float psi_min_angle = current_best_psi - 0.5 * angle_range;
            float psi_max_angle = current_best_psi + 0.5 * angle_range;

            for ( float psi = psi_min_angle; psi < psi_max_angle; psi += psi_step ) {
                // create a new peak to save the cross-correlation peak values
                Peak current_peak;
                average_image.CopyFrom(&average_images[aln_image_counter]);

                //will rotate the original image to be aligned with the sum image to facilitate the shift calculations
                my_image_copy.Allocate(x_dim, y_dim, true);
                my_image_copy.CopyFrom(&my_image);
                my_image_copy.Rotate2DInPlace(psi, FLT_MAX);
                // // make directional sum to eliminate any vertical signal
                // sum_image_direction(&my_image_copy, 2);
                // // calculate the cross correlation of the reference image with the rotated image
                average_image.CalculateCrossCorrelationImageWith(&my_image_copy);

                if ( outer_mask_radius != 0 ) {
                    average_image.CircleMask(outer_mask_radius);
                }
                else if ( outer_mask_radius == 0 ) {
                    average_image.CircleMask(x_dim * 0.45);
                }

                // find the peak from the cross corrlation to get the values
                current_peak = average_image.FindPeakWithParabolaFit(inner_radius_for_peak_search, outer_radius_for_peak_search);

                if ( current_peak.value > local_best_corr_score ) {
                    local_best_corr_score = current_peak.value;
                    local_best_psi        = psi;
                    local_best_x_shift    = current_peak.x;
                    local_best_y_shift    = current_peak.y;
                }
            }
            // end of initial search for the best psi and shift angles
            // Start the tuning loop
            // save the best_psi as the current_best_psi for that image
            // calculate the tuning psi range which is = to the rotation angle and step size = rotation angle / 20
            current_best_psi = local_best_psi;
            //current_best_psi            = best_psi_value[aln_image_counter];
            float tuned_psi_lower_range = current_best_psi - tuned_rotation_range;
            float tuned_psi_upper_range = current_best_psi + tuned_rotation_range;

            // loop over the range of +/- half the rotation angle
            // increment by 1/10 of the tuned rotation angle (rotation angle/2)/10 degrees for tuning
            for ( float tuned_psi = tuned_psi_lower_range; tuned_psi <= tuned_psi_upper_range; tuned_psi += tuned_step_size ) {

                // create a new peak to save the tuned values
                Peak current_tuned_peak;
                tuning_average_image.CopyFrom(&average_images[aln_image_counter]);
                my_image_tuned.Allocate(x_dim, y_dim, true);
                my_image_tuned.CopyFrom(&my_image);
                my_image_tuned.Rotate2DInPlace(tuned_psi, FLT_MAX);
                // // make directional sum to eliminate any vertical signal
                // sum_image_direction(&my_image_tuned, 2);

                // calculate the cross correlation of the reference image with the rotated image
                tuning_average_image.CalculateCrossCorrelationImageWith(&my_image_tuned); //rotated_image
                // Added this as usually tubes are centered so to avoid any extra shift after aligining
                // especially if using a gaussian filter

                if ( outer_mask_radius != 0 ) {
                    tuning_average_image.CircleMask(outer_mask_radius);
                }
                else if ( outer_mask_radius == 0 ) {
                    tuning_average_image.CircleMask(x_dim * 0.45);
                }
                // find the peak from the cross corrlation to get the values from the tuning_average_image
                current_tuned_peak = tuning_average_image.FindPeakWithParabolaFit(inner_radius_for_peak_search, outer_radius_for_peak_search);

                if ( current_tuned_peak.value > local_best_corr_score ) {
                    local_best_corr_score = current_tuned_peak.value;
                    local_best_psi        = tuned_psi;
                    local_best_x_shift    = current_tuned_peak.x;
                    local_best_y_shift    = current_tuned_peak.y;
                }
            }

            best_correlation_score[aln_image_counter] = local_best_corr_score;
            best_psi_value[aln_image_counter]         = local_best_psi;
            best_x_shift_value[aln_image_counter]     = local_best_x_shift;
            // best_y_shift_value[aln_image_counter]     = local_best_y_shift;

            // Updating the tube diameters
            final_image.Allocate(x_dim, y_dim, true);
            final_image.SetToConstant(0.0);

            if ( use_memory ) {
                final_image.CopyFrom(&image_stack_filtered_masked[aln_image_counter]);
                final_image.BackwardFFT( );
            }
            else {
#pragma omp critical
                final_image.ReadSlice(&my_input_images, aln_image_counter + 1);
                final_image.Normalize( );

                if ( outer_mask_radius != 0 ) {
                    final_image.CircleMask(outer_mask_radius);
                }
                else if ( outer_mask_radius == 0 ) {
                    final_image.CircleMask(x_dim * 0.45);
                }
                final_image.ForwardFFT( );
                final_image.ZeroCentralPixel( );
                final_image.GaussianLowPassFilter((pixel_size * 2) / low_pass_resolution);
                final_image.BackwardFFT( );
            }
            // removed the -psi from here as I want to rotate the image to be aligned with Y-axis as the average image to get the correct x-shift
            final_image.Rotate2DInPlace(best_psi_value[aln_image_counter], FLT_MAX);
            final_image.PhaseShift(best_x_shift_value[aln_image_counter], 0.0);
            // find the outer edges peaks
            all_columns_sum[aln_image_counter]              = sum_image_columns(&final_image);
            auto [peak_one_column_sum, peak_two_column_sum] = FindOuterTubeEdges(all_columns_sum[aln_image_counter], min_tube_diameter, max_tube_diameter);

            final_image.Deallocate( );
            my_image_copy.Deallocate( );
            my_image_tuned.Deallocate( );

            if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
                my_aln_progress->Update(aln_image_counter + 1);
        }

        delete my_aln_progress;
        delete[] average_images;
        final_image.Deallocate( );
    }

    wxPrintf("\nWriting aligned images...\n\n");
    if ( ! my_output_images.IsOpen( ) ) {
        my_output_images.OpenFile(output_images.ToStdString( ), true);
        if ( ! my_output_images.IsOpen( ) ) {
            wxPrintf("ERROR: Could not open '%s' for writing\n", output_images.ToStdString( ));
            DEBUG_ABORT;
        }
    }
    my_output_images.my_header.SetNumberOfImages(number_of_input_images);
    my_output_images.my_header.SetDimensionsImage(x_dim, y_dim);
    my_output_images.SetPixelSize(pixel_size);
    my_output_images.WriteHeader( );
    my_output_images.rewrite_header_on_close = true;

    ProgressBar* writing_progress = new ProgressBar(number_of_input_images);

#pragma omp parallel for schedule(static, 1) num_threads(max_threads) default(none) shared(my_input_images, tube_rotation, number_of_input_images, max_threads, my_output_images, x_dim, y_dim, tube_rotation, x_shift_column, use_auto_corr, class_average, best_psi_value, best_x_shift_value, writing_progress) private(image_counter, final_image)

    // #pragma omp for ordered schedule(static, 1)
    for ( image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
        final_image.Allocate(x_dim, y_dim, true);
        final_image.SetToConstant(0.0);
#pragma omp critical
        final_image.ReadSlice(&my_input_images, image_counter + 1);
        float average;
        //average = final_image.ReturnAverageOfRealValuesOnEdges( );
        average = final_image.ReturnAverageOfRealValues( );
        //padd the image so that when rotated it doesn't mess up the outer edges
        final_image.Resize(2 * x_dim, 2 * y_dim, 1, average);
        if ( class_average ) {
            if ( use_auto_corr ) {
                final_image.Rotate2DInPlace(tube_rotation[image_counter], FLT_MAX);
            }
            else {
                final_image.Rotate2DInPlace(tube_rotation[image_counter] + 90.0, FLT_MAX);
            }
            //final_image.PhaseShift(x_shift_column[image_counter], 0);
        }
        else {
            final_image.Rotate2DInPlace(best_psi_value[image_counter], FLT_MAX);
            // // generate the full rotation matrix
            RotationMatrix temp_matrix;
            float          rotated_x, rotated_y, rotated_z;
            temp_matrix.SetToEulerRotation(0.0, 90.0, (90.0 - best_psi_value[image_counter]));
            // assuming no Y shift will be applied to ensure everything is centered
            temp_matrix.RotateCoords((best_x_shift_value[image_counter]), 0.0, 0.0, rotated_x, rotated_y, rotated_z);
            // correct the shift happenning because of extractslice - Note x and y shift needs to be flipped
            // since projection is centered so negative rotation and shift is needed here
            final_image.PhaseShift(-rotated_y, 0);
        }

        final_image.Resize(x_dim, y_dim, 1);

#pragma omp critical
        final_image.WriteSlice(&my_output_images, image_counter + 1);
        final_image.Deallocate( );
        if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
            writing_progress->Update(image_counter + 1);
    }

    return true;
}

/// try the column sum but add a declaration at the begining of the code
std::vector<float> sum_image_columns(Image* current_image) {
    std::vector<float> column_sum(current_image->logical_x_dimension, 0.0);

    long pixel_counter = 0;

    for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
        for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
            long pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
            column_sum[i] += current_image->real_values[pixel_coord_xy];
            pixel_counter++;
        }
        pixel_counter += current_image->padding_jump_value;
    }

    return column_sum;
}

float max_abs_column_sum(Image* current_image) { // Tim's method to calculate the best rotation based on the row sum
    std::vector<float> column_sum(current_image->logical_x_dimension, 0.0);

    long pixel_counter = 0;

    for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
        for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
            long pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
            column_sum[i] += current_image->real_values[pixel_coord_xy];
            pixel_counter++;
        }
        pixel_counter += current_image->padding_jump_value;
    }

    float max_value = *std::max_element(column_sum.begin( ), column_sum.end( ),
                                        [](float a, float b) { return std::abs(a) < std::abs(b); });

    return abs(max_value);
}

// Detects the two strongest outer-edge peaks in a 1D intensity profile.
// Returns indices of the best peak pair (sorted low->high), or an empty vector if none found.
std::pair<int, int> FindOuterTubeEdges(const std::vector<float>& cols, float min_tube_diameter, float max_tube_diameter) {
    int n = cols.size( );
    if ( n < 3 )
        return {-1, -1}; // need at least 3 points to form a peak

    // 1) Normalize the 1D profile
    float              minVal = *std::min_element(cols.begin( ), cols.end( ));
    float              maxVal = *std::max_element(cols.begin( ), cols.end( ));
    std::vector<float> norm(n);
    for ( int i = 0; i < n; ++i )
        norm[i] = cols[i] - minVal;

    float normMax = *std::max_element(norm.begin( ), norm.end( ));
    if ( normMax <= 0.0f )
        return {-1, -1};

    // Inverted profile for negative peaks
    std::vector<float> normInv(n);
    for ( int i = 0; i < n; ++i )
        normInv[i] = normMax - norm[i];

    // 2) Detect peaks
    std::vector<std::pair<int, float>> posPeaks;
    std::vector<std::pair<int, float>> negPeaks;

    for ( int i = 1; i < n - 1; ++i ) {
        if ( norm[i] > norm[i - 1] && norm[i] > norm[i + 1] ) {
            posPeaks.emplace_back(i, norm[i]);
        }
        if ( normInv[i] > normInv[i - 1] && normInv[i] > normInv[i + 1] ) {
            negPeaks.emplace_back(i, normInv[i]);
        }
    }
    // // debugging and printing the scores
    // std::cerr << "posPeaks (idx,val): ";
    // for ( auto& p : posPeaks )
    //     std::cerr << "(" << p.first << "," << p.second << ") ";
    // std::cerr << "\n";
    // std::cerr << "negPeaks (idx,val): ";
    // for ( auto& p : negPeaks )
    //     std::cerr << "(" << p.first << "," << p.second << ") ";
    // std::cerr << "\n";

    // helper function to find the best pair of peaks based on their height and distance between peaks
    auto bestPair = [&](const std::vector<std::pair<int, float>>& peaks)
            -> std::pair<float, std::pair<int, int>> {
        float               bestScore = -std::numeric_limits<float>::infinity( );
        std::pair<int, int> bestIdx   = {-1, -1};

        // Adding gap penalty and out of range penalty so that we would favor more the peaks within the range, but also if nothing was found within range, out of range peaks are saved and returned
        const float IDEAL_GAP           = min_tube_diameter;
        const float GAP_PENALTY         = 0.1f; // e.g. 0.1 points lost per pixel of gap deviation
        const float OUT_OF_RANGE_FACTOR = 10.0f; // scale factor for out-of-range penalty- changed that from 2 to 10 to heavily penalize out of range to favor in range more

        for ( size_t a = 0; a < peaks.size( ); ++a ) {
            for ( size_t b = a + 1; b < peaks.size( ); ++b ) {
                int i   = peaks[a].first;
                int j   = peaks[b].first;
                int gap = j - i;

                float sumAmp = peaks[a].second + peaks[b].second;
                float score  = sumAmp - GAP_PENALTY * std::fabs(gap - IDEAL_GAP);

                // scale penalty by how far out of range the gap is
                if ( gap < min_tube_diameter ) {
                    score -= OUT_OF_RANGE_FACTOR * (min_tube_diameter - gap);
                }
                else if ( gap > max_tube_diameter ) {
                    score -= OUT_OF_RANGE_FACTOR * (gap - max_tube_diameter);
                }

                if ( score > bestScore ) {
                    bestScore = score;
                    bestIdx   = {i, j};
                }
            }
        }

        return std::make_pair(bestScore, bestIdx);
    };

    // 3) Find best pair among positive peaks and among negative peaks.
    auto [scorePos, bestPos] = bestPair(posPeaks);
    auto [scoreNeg, bestNeg] = bestPair(negPeaks);

    std::pair<int, int> bestPairIdx = {-1, -1};

    // 4) If no valid pairs exist at all, return -1
    if ( scorePos == -std::numeric_limits<float>::infinity( ) &&
         scoreNeg == -std::numeric_limits<float>::infinity( ) ) {
        return {-1, -1};
    }

    // 5) keeping the values of the best negative peaks as reference
    bestPairIdx     = bestNeg;
    float bestScore = scoreNeg;

    // find the highest negative peaks within the range of the expected diameter
    // then find the positive peak before the first negative peak and the positive peak after the second negative peak and those should be the outer edges
    if ( bestNeg.first != -1 && bestNeg.second != -1 ) {
        int iNeg = bestNeg.first;
        int jNeg = bestNeg.second;
        if ( iNeg > jNeg )
            std::swap(iNeg, jNeg); // enforce left->right

        // Find last positive BEFORE iNeg
        int   posBefore = -1;
        float ampBefore = 0;
        for ( auto it = posPeaks.rbegin( ); it != posPeaks.rend( ); ++it ) {
            if ( it->first < iNeg ) {
                posBefore = it->first;
                ampBefore = it->second;
                break;
            }
        }

        // Find first positive AFTER jNeg
        int   posAfter = -1;
        float ampAfter = 0;
        for ( auto& p : posPeaks ) {
            if ( p.first > jNeg ) {
                posAfter = p.first;
                ampAfter = p.second;
                break;
            }
        }
        const float IDEAL_GAP           = min_tube_diameter;
        const float GAP_PENALTY         = 0.1f; // e.g. 0.1 points lost per pixel of gap deviation
        const float OUT_OF_RANGE_FACTOR = 10.0f; // scale factor for out-of-range penalty

        // Step 3: Only refine if both positives exist and are ordered
        if ( posAfter != -1 && posBefore != -1 && posAfter < posBefore ) {
            int   gap    = posBefore - posAfter;
            float sumAmp = ampAfter + ampBefore;
            float score  = sumAmp - GAP_PENALTY * std::fabs(gap - IDEAL_GAP);

            if ( gap < min_tube_diameter )
                score -= OUT_OF_RANGE_FACTOR * (min_tube_diameter - gap);
            else if ( gap > max_tube_diameter )
                score -= OUT_OF_RANGE_FACTOR * (gap - max_tube_diameter);

            // Step 4: Replace if adjacency score is better
            if ( score > bestScore ) {
                bestScore   = score;
                bestPairIdx = {posAfter, posBefore};
            }
        }
    }
    // Final: enforce sorted order before returning
    if ( bestPairIdx.first > bestPairIdx.second )
        std::swap(bestPairIdx.first, bestPairIdx.second);

    return std::make_pair(bestPairIdx.first, bestPairIdx.second);
}

void sum_image_direction(Image* current_image, int dim) {
    // image must be in real-space
    Image directional_image_sum;
    directional_image_sum.Allocate(current_image->logical_x_dimension, current_image->logical_y_dimension, true);
    directional_image_sum.SetToConstant(0.0);

    // x-direction
    if ( dim == 1 ) {

        long pixel_coord_y  = 0;
        long pixel_coord_xy = 0;
        long pixel_counter  = 0;

        // sum columns of my_image_sum (NxM) and store in array (1xN)
        for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
            for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
                pixel_coord_y  = current_image->ReturnReal1DAddressFromPhysicalCoord(0, j, 0);
                pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                directional_image_sum.real_values[pixel_coord_y] += current_image->real_values[pixel_coord_xy];
                pixel_counter++;
            }
            pixel_counter += current_image->padding_jump_value;
        }

        // repeat column sum into my_vertical_sum
        pixel_counter = 0;
        for ( int j = 0; j < directional_image_sum.logical_y_dimension; j++ ) {
            for ( int i = 0; i < directional_image_sum.logical_x_dimension; i++ ) {
                pixel_coord_y                                     = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(0, j, 0);
                pixel_coord_xy                                    = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                directional_image_sum.real_values[pixel_coord_xy] = directional_image_sum.real_values[pixel_coord_y];
                pixel_counter++;
            }
            pixel_counter += directional_image_sum.padding_jump_value;
        }

        directional_image_sum.DivideByConstant(directional_image_sum.logical_x_dimension);
    }
    // y-direction
    else {

        long pixel_coord_x  = 0;
        long pixel_coord_xy = 0;
        long pixel_counter  = 0;

        // sum columns of my_image_sum (NxM) and store in array (1xM)
        for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
            for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
                pixel_coord_x  = current_image->ReturnReal1DAddressFromPhysicalCoord(i, 0, 0);
                pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                directional_image_sum.real_values[pixel_coord_x] += current_image->real_values[pixel_coord_xy];
                pixel_counter++;
            }
            pixel_counter += current_image->padding_jump_value;
        }

        // repeat column sum into my_vertical_sum
        pixel_counter = 0;
        for ( int i = 0; i < directional_image_sum.logical_x_dimension; i++ ) {
            for ( int j = 0; j < directional_image_sum.logical_y_dimension; j++ ) {
                pixel_coord_x                                     = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(i, 0, 0);
                pixel_coord_xy                                    = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                directional_image_sum.real_values[pixel_coord_xy] = directional_image_sum.real_values[pixel_coord_x];
                pixel_counter++;
            }
            pixel_counter += directional_image_sum.padding_jump_value;
        }

        directional_image_sum.DivideByConstant(directional_image_sum.logical_y_dimension);
    }

    // copy to current iamge
    current_image->CopyFrom(&directional_image_sum);
    directional_image_sum.Deallocate( );
}

void InitializeCTFSumOfSquares(Image& current_image, std::vector<float>& ctf_sum_of_squares) {
    int number_of_pixels = current_image.real_memory_allocated / 2;

    // Resize to number_of_pixels and set all values to zero
    ctf_sum_of_squares.assign(number_of_pixels, 0.0f);
}

void ApplyCTFAndReturnCTFSumOfSquares(Image& image, CTF ctf_to_apply, bool absolute, bool apply_beam_tilt, bool apply_envelope, std::vector<float>& ctf_sum_of_squares) {
    MyDebugAssertTrue(image.is_in_memory, "Memory not allocated");
    MyDebugAssertTrue(image.is_in_real_space == false, "Image not in Fourier space");
    MyDebugAssertTrue(image.logical_z_dimension == 1, "Volumes not supported");

    //std::vector<float> squared_ctf_values; // Vector to store squared CTF values

    long pixel_counter = 0;

    for ( int j = 0; j <= image.physical_upper_bound_complex_y; j++ ) {
        float y_coord    = image.ReturnFourierLogicalCoordGivenPhysicalCoord_Y(j) * image.fourier_voxel_size_y;
        float y_coord_sq = pow(y_coord, 2.0);

        for ( int i = 0; i <= image.physical_upper_bound_complex_x; i++ ) {
            float x_coord    = i * image.fourier_voxel_size_x;
            float x_coord_sq = pow(x_coord, 2);
            float azimuth;
            // Compute the azimuth
            if ( i == 0 && j == 0 ) {
                azimuth = 0.0;
            }
            else {
                azimuth = atan2(y_coord, x_coord);
            }

            // Compute the square of the frequency
            float frequency_squared = x_coord_sq + y_coord_sq;
            float ctf_value;

            if ( apply_envelope ) {
                ctf_value = ctf_to_apply.EvaluateWithEnvelope(frequency_squared, azimuth);
            }
            else {
                ctf_value = ctf_to_apply.Evaluate(frequency_squared, azimuth);
            }

            if ( absolute ) {
                ctf_value = fabsf(ctf_value);
            }

            // Apply CTF to the image if needed
            image.complex_values[pixel_counter] *= ctf_value;

            if ( apply_beam_tilt && (ctf_to_apply.GetBeamTiltX( ) != 0.0f || ctf_to_apply.GetBeamTiltY( ) != 0.0f) ) {
                image.complex_values[pixel_counter] *= ctf_to_apply.EvaluateBeamTiltPhaseShift(frequency_squared, azimuth);
            }

            // Add the squared CTF value to the input CTF sum of squares vector
            // but check first if the vector is empty to avoid memory problems
            // if (ctf_sum_of_squares.empty()) {
            //     ctf_sum_of_squares[pixel_counter] = powf(ctf_value, 2);
            // } else {
            ctf_sum_of_squares[pixel_counter] += powf(ctf_value, 2);
            // }

            pixel_counter++;
        }
    }

    //return ctf_sum_of_squares;
}

void divide_by_ctf_sum_of_squares(Image& current_image, std::vector<float>& ctf_sum_of_squares) {
    // normalize by sum of squared CTFs (voxel by voxel)
    long pixel_counter = 0;

    for ( int j = 0; j <= current_image.physical_upper_bound_complex_y; j++ ) {
        for ( int i = 0; i <= current_image.physical_upper_bound_complex_x; i++ ) {
            if ( ctf_sum_of_squares[pixel_counter] != 0.0 )
                current_image.complex_values[pixel_counter] /= sqrtf(ctf_sum_of_squares[pixel_counter]);
            pixel_counter++;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Helper function: Find peaks in a 1D dataset
// Returns a vector of pairs: {index, value}
// Sorted by index
std::vector<std::pair<int, float>> FindPeaks(const std::vector<float>& data, float min_dist, float threshold) {
    std::vector<std::pair<int, float>> peaks;
    int                                n = data.size( );
    if ( n < 3 )
        return peaks;

    // 1. Identify local maxima above threshold
    for ( int i = 1; i < n - 1; ++i ) {
        if ( data[i] > data[i - 1] && data[i] > data[i + 1] ) {
            if ( data[i] >= threshold ) {
                peaks.push_back({i, data[i]});
            }
        }
    }

    // 2. Sort by amplitude (descending) to prioritize processing largest peaks
    std::sort(peaks.begin( ), peaks.end( ), [](const std::pair<int, float>& a, const std::pair<int, float>& b) {
        return a.second > b.second;
    });

    // 3. Filter peaks based on minimum distance
    std::vector<std::pair<int, float>> filtered_peaks;
    for ( const auto& p : peaks ) {
        bool keep = true;
        for ( const auto& accepted : filtered_peaks ) {
            if ( std::abs(p.first - accepted.first) < min_dist ) {
                keep = false;
                break;
            }
        }
        if ( keep ) {
            filtered_peaks.push_back(p);
        }
    }

    // 4. Sort final results by index (ascending) for geometric processing
    std::sort(filtered_peaks.begin( ), filtered_peaks.end( ), [](const std::pair<int, float>& a, const std::pair<int, float>& b) {
        return a.first < b.first;
    });

    return filtered_peaks;
}

// Main Algorithm: Find Outer Tube Edges
// use_half_way defaults to true
std::pair<float, float> FindOuterTubeEdges(const std::vector<float>& cols, float min_tube_diameter, float max_tube_diameter, bool use_half_way = true) {
    int n = cols.size( );
    if ( n < 3 )
        return {-1.0f, -1.0f};

    // 1. Smooth the profile
    std::vector<float> smooth_cols   = cols;
    int                smooth_radius = 2;
    for ( int i = smooth_radius; i < n - smooth_radius; ++i ) {
        double sum = 0;
        for ( int k = -smooth_radius; k <= smooth_radius; ++k ) {
            sum += cols[i + k];
        }
        smooth_cols[i] = sum / (2 * smooth_radius + 1);
    }

    // 2. Normalize
    float              minVal = *std::min_element(smooth_cols.begin( ), smooth_cols.end( ));
    std::vector<float> norm(n);
    for ( int i = 0; i < n; ++i )
        norm[i] = smooth_cols[i] - minVal;

    float normMax = *std::max_element(norm.begin( ), norm.end( ));
    if ( normMax <= 0.0f )
        return {-1.0f, -1.0f};

    // Inverted profile for Negative peaks (Inner Walls)
    std::vector<float> normInv(n);
    for ( int i = 0; i < n; ++i )
        normInv[i] = normMax - norm[i];

    // 3. Find All Peaks
    float min_dist  = 5.0f; // Minimum distance between peaks of same type
    float threshold = 0.0f;

    // posPeaks = Candidates for PL and PR (Outer Walls)
    std::vector<std::pair<int, float>> posPeaks = FindPeaks(norm, min_dist, threshold);
    // negPeaks = Candidates for NL and NR (Inner Walls)
    std::vector<std::pair<int, float>> negPeaks = FindPeaks(normInv, min_dist, threshold);

    // 4. Search Pattern: NL -> PL ... PR <- NR
    float bestScore = -std::numeric_limits<float>::infinity( );
    int   best_NL = -1, best_NR = -1;
    int   best_PL = -1, best_PR = -1;

    float center_idx = (float)(n - 1) / 2.0f;

    // Iterate through all possible Left Negative Peaks (NL)
    for ( const auto& pNL : negPeaks ) {
        int   idx_NL = pNL.first;
        float val_NL = pNL.second;

        // Iterate through all possible Right Negative Peaks (NR)
        for ( const auto& pNR : negPeaks ) {
            int   idx_NR = pNR.first;
            float val_NR = pNR.second;

            // Basic geometric constraints
            if ( idx_NR <= idx_NL )
                continue; // Right must be to the right
            float width = idx_NR - idx_NL;
            if ( width < min_tube_diameter || width > max_tube_diameter )
                continue;

            // Find best PL: Highest Positive peak strictly between NL and Center
            int   idx_PL     = -1;
            float max_val_PL = -1.0f;

            for ( const auto& pPos : posPeaks ) {
                if ( pPos.first > idx_NL && pPos.first < (idx_NL + idx_NR) / 2.0f ) {
                    if ( pPos.second > max_val_PL ) {
                        max_val_PL = pPos.second;
                        idx_PL     = pPos.first;
                    }
                }
            }

            // Find best PR: Highest Positive peak strictly between Center and NR
            int   idx_PR     = -1;
            float max_val_PR = -1.0f;

            for ( const auto& pPos : posPeaks ) {
                if ( pPos.first > (idx_NL + idx_NR) / 2.0f && pPos.first < idx_NR ) {
                    if ( pPos.second > max_val_PR ) {
                        max_val_PR = pPos.second;
                        idx_PR     = pPos.first;
                    }
                }
            }

            // Require both Positive peaks to exist for this pattern
            if ( idx_PL != -1 && idx_PR != -1 ) {

                // --- SCORING ---
                float score = 0.0f;

                // 1. Magnitude Score (Sum of all 4 peaks)
                score += (val_NL + val_NR + max_val_PL + max_val_PR);

                // 2. Symmetry Penalty (Tube should be roughly centered)
                float midpoint = (float)(idx_NL + idx_NR) / 2.0f;
                score -= 5.0f * std::abs(midpoint - center_idx) / n;

                // 3. Wall Thickness Consistency Penalty
                float left_wall_w  = idx_PL - idx_NL;
                float right_wall_w = idx_NR - idx_PR;
                score -= 2.0f * std::abs(left_wall_w - right_wall_w);

                if ( score > bestScore ) {
                    bestScore = score;
                    best_NL   = idx_NL;
                    best_NR   = idx_NR;
                    best_PL   = idx_PL;
                    best_PR   = idx_PR;
                }
            }
        }
    }

    // 5. Return Results
    if ( best_NL != -1 && best_NR != -1 && best_PL != -1 && best_PR != -1 ) {
        if ( use_half_way ) {
            // Average of Inner (Neg) and Outer (Pos) wall positions
            float edge_L = (float)(best_NL + best_PL) / 2.0f;
            float edge_R = (float)(best_NR + best_PR) / 2.0f;
            return {edge_L, edge_R};
        }
        else {
            // Return the Inner Walls (Negative peaks)
            return {(float)best_NL, (float)best_NR};
        }
    }

    return {-1.0f, -1.0f};
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////