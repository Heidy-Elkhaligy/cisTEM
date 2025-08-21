#include "../../core/core_headers.h"
#include <iostream>
#include <fstream>
#include <cmath>

class
        apply_RASTR_phi_constraints : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

float angle_within360(float angle);

IMPLEMENT_APP(apply_RASTR_phi_constraints)

// override the DoInteractiveUserInput

void apply_RASTR_phi_constraints::DoInteractiveUserInput( ) {
    float occupancy_threshold = 0.0;
    int   max_threads;

    UserInput* my_input                = new UserInput("apply_RASTR_phi_constraints", 1.00);
    wxString   input_images            = my_input->GetFilenameFromUser("Input images filename", "Filename of input images to be filtered", "my_stack.mrc", true);
    wxString   output_images           = my_input->GetFilenameFromUser("Output filtered images filename", "The filtered images based on phi angle", "filtered_stack.mrc", false);
    wxString   input_star_filename     = my_input->GetFilenameFromUser("Input star file", "The input star file of the input images", "my_parameters.star", true);
    wxString   reference_star_filename = my_input->GetFilenameFromUser("Input reference star file", "The input reference star filename (usually RASTR output star file)", "reference_parameters.star", true);
    wxString   output_star_filename    = my_input->GetFilenameFromUser("Output star file name", "The output star file after filtering the images", "filtered_parameters.star", false);
    float      angular_range           = my_input->GetFloatFromUser("Allowed phi angular range", "The allowed angular range to keep around the input given angle in the input reference star file", "0.0", 0.0);
    bool       classification_results  = my_input->GetYesNoFromUser("Is the input star file from a 3d classification?", "If input star file is from 3d classification, then occupancy threshold needs to be given", "NO");
    if ( classification_results ) {
        occupancy_threshold = my_input->GetFloatFromUser("Occupancy threshold", "All particles above that threshold will be kept", "80.0", 0.0, 100.0);
    }
#ifdef _OPENMP
    max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "when threading, what is the max threads to run", "1", 1);
#else
    max_threads = 1;
#endif

    delete my_input;

    my_current_job.Reset(9);
    my_current_job.ManualSetArguments("tttttfbfi", input_images.ToUTF8( ).data( ), output_images.ToUTF8( ).data( ), input_star_filename.ToUTF8( ).data( ), reference_star_filename.ToUTF8( ).data( ), output_star_filename.ToUTF8( ).data( ), angular_range, classification_results, occupancy_threshold, max_threads); //update_star_file, input_star_filename.ToUTF8( ).data( ),
}

// override the do calculation method which will be what is actually run..

bool apply_RASTR_phi_constraints::DoCalculation( ) {
    wxString input_images            = my_current_job.arguments[0].ReturnStringArgument( );
    wxString output_images           = my_current_job.arguments[1].ReturnStringArgument( );
    wxString input_star_filename     = my_current_job.arguments[2].ReturnStringArgument( );
    wxString reference_star_filename = my_current_job.arguments[3].ReturnStringArgument( );
    wxString output_star_filename    = my_current_job.arguments[4].ReturnStringArgument( );
    float    angular_range           = my_current_job.arguments[5].ReturnFloatArgument( );
    bool     classification_results  = my_current_job.arguments[6].ReturnBoolArgument( );
    float    occupancy_threshold     = my_current_job.arguments[7].ReturnFloatArgument( );
    int      max_threads             = my_current_job.arguments[8].ReturnIntegerArgument( );

    MRCFile          my_input_images(input_images.ToStdString( ), false);
    MRCFile          my_output_images(output_images.ToStdString( ), true);
    cisTEMParameters input_star_file;
    cisTEMParameters ref_star_file;

    if ( (is_running_locally && ! DoesFileExist(input_star_filename.ToStdString( ))) ) {
        SendErrorAndCrash(wxString::Format("Error: Input star file %s not found\n", input_star_filename));
    }
    input_star_file.ReadFromcisTEMStarFile(input_star_filename);

    if ( (is_running_locally && ! DoesFileExist(reference_star_filename.ToStdString( ))) ) {
        SendErrorAndCrash(wxString::Format("Error: Input star file %s not found\n", reference_star_filename));
    }
    ref_star_file.ReadFromcisTEMStarFile(reference_star_filename);

    long  number_of_input_images = my_input_images.ReturnNumberOfSlices( );
    Image my_image;

    cisTEMParameterLine input_parameters;
    cisTEMParameterLine reference_parameters;
    cisTEMParameters    output_params;
    // setup parameters for the output star file
    output_params.parameters_to_write.SetActiveParameters(POSITION_IN_STACK | IMAGE_IS_ACTIVE | PSI | THETA | PHI | X_SHIFT | Y_SHIFT | DEFOCUS_1 | DEFOCUS_2 | DEFOCUS_ANGLE | PHASE_SHIFT | OCCUPANCY | LOGP | SIGMA | SCORE | PIXEL_SIZE | MICROSCOPE_VOLTAGE | MICROSCOPE_CS | AMPLITUDE_CONTRAST | BEAM_TILT_X | BEAM_TILT_Y | IMAGE_SHIFT_X | IMAGE_SHIFT_Y | ASSIGNED_SUBSET);
    output_params.PreallocateMemoryAndBlank(number_of_input_images); //in case all had occupance > threshold

    long         new_counter = 0;
    ProgressBar* my_progress = new ProgressBar(number_of_input_images);

    for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {

        input_parameters     = input_star_file.ReturnLine(image_counter); // the star file numbering is 0 indexed!!
        reference_parameters = ref_star_file.ReturnLine(image_counter);
        float input_phi      = input_parameters.phi; //should I make this absolute to remove the negative value effect??
        // make sure input phi angle is positive and within 0-360
        input_phi                  = angle_within360(input_phi);
        float reference_phi        = reference_parameters.phi;
        float mirror_reference_phi = reference_phi + 180.0;

        if ( classification_results ) {
            if ( (input_parameters.occupancy >= occupancy_threshold) && (input_phi >= reference_phi - angular_range && input_phi <= reference_phi + angular_range) ) { //|| (input_phi >= mirror_reference_phi - angular_range && input_phi <= mirror_reference_phi + angular_range)
                // Only keep the images with psi around 90 or 270 (any other number will be considered misaligned or wrong particle)
                // float input_psi = input_parameters.psi;
                // float psi_range = 5.0;
                // if ( (input_psi >= 90.0 - psi_range && input_psi <= 90.0 + psi_range) || (input_psi >= 270.0 - psi_range && input_psi <= 270.0 + psi_range) ) {
                my_image.ReadSlice(&my_input_images, image_counter + 1);
                // save the image into the new output MRC file
                my_image.WriteSlice(&my_output_images, new_counter + 1);
                // save the parameter information of the image into the new star file
                //wxPrintf("The input phi %f and the reference phi is %f \n", input_phi, reference_phi);
                //wxPrintf("The image counter is %li and the new counter is %li \n", image_counter + 1, new_counter + 1);
                // new_counter here should start at 0
                output_params.all_parameters[new_counter].position_in_stack                  = new_counter + 1;
                output_params.all_parameters[new_counter].image_is_active                    = input_parameters.image_is_active;
                output_params.all_parameters[new_counter].psi                                = input_parameters.psi;
                output_params.all_parameters[new_counter].theta                              = input_parameters.theta;
                output_params.all_parameters[new_counter].phi                                = input_parameters.phi;
                output_params.all_parameters[new_counter].x_shift                            = input_parameters.x_shift;
                output_params.all_parameters[new_counter].y_shift                            = input_parameters.y_shift;
                output_params.all_parameters[new_counter].defocus_1                          = input_parameters.defocus_1;
                output_params.all_parameters[new_counter].defocus_2                          = input_parameters.defocus_2;
                output_params.all_parameters[new_counter].defocus_angle                      = input_parameters.defocus_angle;
                output_params.all_parameters[new_counter].phase_shift                        = input_parameters.phase_shift;
                output_params.all_parameters[new_counter].occupancy                          = input_parameters.occupancy;
                output_params.all_parameters[new_counter].logp                               = input_parameters.logp;
                output_params.all_parameters[new_counter].sigma                              = input_parameters.sigma;
                output_params.all_parameters[new_counter].score                              = input_parameters.score;
                output_params.all_parameters[new_counter].score_change                       = input_parameters.score_change;
                output_params.all_parameters[new_counter].pixel_size                         = input_parameters.pixel_size;
                output_params.all_parameters[new_counter].microscope_voltage_kv              = input_parameters.microscope_voltage_kv;
                output_params.all_parameters[new_counter].microscope_spherical_aberration_mm = input_parameters.microscope_spherical_aberration_mm;
                output_params.all_parameters[new_counter].amplitude_contrast                 = input_parameters.amplitude_contrast;
                output_params.all_parameters[new_counter].beam_tilt_x                        = input_parameters.beam_tilt_x;
                output_params.all_parameters[new_counter].beam_tilt_y                        = input_parameters.beam_tilt_y;
                output_params.all_parameters[new_counter].image_shift_x                      = input_parameters.image_shift_x;
                output_params.all_parameters[new_counter].image_shift_y                      = input_parameters.image_shift_y;
                if ( input_parameters.position_in_stack % 2 == 1 ) {
                    input_parameters.assigned_subset = 1; // Odd particle number
                }
                else {
                    input_parameters.assigned_subset = 2; // Even particle number
                }
                output_params.all_parameters[new_counter].assigned_subset = input_parameters.assigned_subset; // There is no assigned subset and we need to keep the saved assigned subset as is to ensure that no 2 particles are in the same group even after filteration
                //wxPrintf("input parameters assigned subset is %i \n", input_parameters.assigned_subset);

                new_counter++;
            }
        }
        else {
            // bool keep = false;
            if ( (input_phi >= reference_phi - angular_range && input_phi <= reference_phi + angular_range) ) {
                //if ( (input_phi >= reference_phi - angular_range && input_phi <= reference_phi + angular_range) ) { //|| (input_phi >= mirror_reference_phi - angular_range && input_phi <= mirror_reference_phi + angular_range)
                // Only keep the images with psi around 90 or 270 (any other number will be considered misaligned or wrong particle)
                // float input_psi = input_parameters.psi;
                // float psi_range = 5.0;
                // if ( (input_psi >= 90.0 - psi_range && input_psi <= 90.0 + psi_range) || (input_psi >= 270.0 - psi_range && input_psi <= 270.0 + psi_range) ) {
                my_image.ReadSlice(&my_input_images, image_counter + 1);
                // save the image into the new output MRC file
                my_image.WriteSlice(&my_output_images, new_counter + 1);
                // save the parameter information of the image into the new star file
                //wxPrintf("The input phi %f and the reference phi is %f \n", input_phi, reference_phi);
                //wxPrintf("The image counter is %li and the new counter is %li \n", image_counter + 1, new_counter + 1);
                // new_counter here should start at 0
                output_params.all_parameters[new_counter].position_in_stack                  = new_counter + 1;
                output_params.all_parameters[new_counter].image_is_active                    = input_parameters.image_is_active;
                output_params.all_parameters[new_counter].psi                                = input_parameters.psi;
                output_params.all_parameters[new_counter].theta                              = input_parameters.theta;
                output_params.all_parameters[new_counter].phi                                = input_parameters.phi;
                output_params.all_parameters[new_counter].x_shift                            = input_parameters.x_shift;
                output_params.all_parameters[new_counter].y_shift                            = input_parameters.y_shift;
                output_params.all_parameters[new_counter].defocus_1                          = input_parameters.defocus_1;
                output_params.all_parameters[new_counter].defocus_2                          = input_parameters.defocus_2;
                output_params.all_parameters[new_counter].defocus_angle                      = input_parameters.defocus_angle;
                output_params.all_parameters[new_counter].phase_shift                        = input_parameters.phase_shift;
                output_params.all_parameters[new_counter].occupancy                          = input_parameters.occupancy;
                output_params.all_parameters[new_counter].logp                               = input_parameters.logp;
                output_params.all_parameters[new_counter].sigma                              = input_parameters.sigma;
                output_params.all_parameters[new_counter].score                              = input_parameters.score;
                output_params.all_parameters[new_counter].score_change                       = input_parameters.score_change;
                output_params.all_parameters[new_counter].pixel_size                         = input_parameters.pixel_size;
                output_params.all_parameters[new_counter].microscope_voltage_kv              = input_parameters.microscope_voltage_kv;
                output_params.all_parameters[new_counter].microscope_spherical_aberration_mm = input_parameters.microscope_spherical_aberration_mm;
                output_params.all_parameters[new_counter].amplitude_contrast                 = input_parameters.amplitude_contrast;
                output_params.all_parameters[new_counter].beam_tilt_x                        = input_parameters.beam_tilt_x;
                output_params.all_parameters[new_counter].beam_tilt_y                        = input_parameters.beam_tilt_y;
                output_params.all_parameters[new_counter].image_shift_x                      = input_parameters.image_shift_x;
                output_params.all_parameters[new_counter].image_shift_y                      = input_parameters.image_shift_y;
                if ( input_parameters.position_in_stack % 2 == 1 ) {
                    input_parameters.assigned_subset = 1; // Odd particle number
                }
                else {
                    input_parameters.assigned_subset = 2; // Even particle number
                }
                output_params.all_parameters[new_counter].assigned_subset = input_parameters.assigned_subset; // There is no assigned subset and we need to keep the saved assigned subset as is to ensure that no 2 particles are in the same group even after filteration
                //wxPrintf("input parameters assigned subset is %i \n", input_parameters.assigned_subset);

                new_counter++;
            }
            // else if ( (reference_phi == 90.0 || reference_phi == 270.0) && ((input_phi >= 0.0 - angular_range && input_phi <= 0.0 + angular_range)) ) {
            //     keep = true;
            // }
            // else if ( (reference_phi == 180.0) && ((input_phi >= 0.0 - angular_range && input_phi <= 0.0 + angular_range)) ) {
            //     keep = true;
            // }

            // bool keep = isWithinMaskRegion(reference_phi, input_phi, 0.0f, angular_range);
        }

        if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
            my_progress->Update(image_counter + 1);
    }
    delete my_progress;
    // write the output star file for the matched references
    output_params.WriteTocisTEMStarFile(output_star_filename);
    float percent_kept     = (float(new_counter) / float(number_of_input_images)) * 100.0f;
    float percent_filtered = ((float(number_of_input_images) - float(new_counter)) / number_of_input_images) * 100.0f;
    wxPrintf("\n\n%.3f %% (= %li particles) are kept, and %.3f %% (= %li particles) are filtered out successfully\n", percent_kept, new_counter, percent_filtered, number_of_input_images - new_counter);

    return true;
}

// Function to ensure the angle is within the range [0, 360)
float angle_within360(float angle) {
    if ( angle < 0.0 ) {
        angle += 360.0;
        return angle_within360(angle);
    }
    else if ( angle >= 360.0 ) {
        angle -= 360.0;
        return angle_within360(angle);
    }
    else {
        return angle;
    }
}
