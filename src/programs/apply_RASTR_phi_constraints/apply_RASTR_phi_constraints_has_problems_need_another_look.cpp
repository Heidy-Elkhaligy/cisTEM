#include "../../core/core_headers.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip> // for std::setprecision

class
        apply_RASTR_phi_constraints : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

struct Result {
    int   image_index;
    float input_phi;
    float reference_phi;
    float difference;
    float adj_ref_phi;
    float input_psi;
    float ref_psi;
    float adj_ref_psi;
    float input_theta; //no need for ref_theta as it was set to 90
    float ref_theta;
    float adj_ref_theta;
};

float angle_within360(float angle);
float angle_difference(float a, float b);

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
    int   num_of_models          = 4;
    long  number_of_base_images  = number_of_input_images / num_of_models;
    Image my_image;
    //Added new as OMP was causing problems when writing images to a file that is not opened and have set dimensions and header information
    MRCFile removed_images_output("removed_images.mrc", true);
    if ( ! removed_images_output.IsOpen( ) ) {
        removed_images_output.OpenFile("removed_images.mrc", true);
        if ( ! removed_images_output.IsOpen( ) ) {
            wxPrintf("ERROR: Could not open '%s' for writing\n", "removed_images.mrc");
            DEBUG_ABORT;
        }
    }
    removed_images_output.my_header.SetDimensionsImage(my_input_images.ReturnXSize( ), my_input_images.ReturnYSize( ));
    removed_images_output.SetPixelSize(my_input_images.ReturnPixelSize( ));
    removed_images_output.WriteHeader( );

    std::ofstream angles_file("difference_angles_file.txt"); // Open file once

    if ( ! angles_file.is_open( ) ) {
        std::cerr << "Error: Could not open difference_angles_file.txt\n";
        //return;
    }

    angles_file << std::fixed << std::setprecision(2); // Optional: set float precision

    angles_file << "image_index, input_phi, reference_phi, difference, adj_ref_phi, input_psi, reference_psi, adj_ref_psi, input_theta, ref_theta, adj_ref_theta\n";

    std::vector<Result> results; // start empty, will grow dynamically
    results.reserve(number_of_input_images);

    cisTEMParameterLine input_parameters;
    cisTEMParameterLine reference_parameters;
    cisTEMParameters    output_params;

    // setup parameters for the output star file
    output_params.parameters_to_write.SetActiveParameters(POSITION_IN_STACK | IMAGE_IS_ACTIVE | PSI | THETA | PHI | X_SHIFT | Y_SHIFT | DEFOCUS_1 | DEFOCUS_2 | DEFOCUS_ANGLE | PHASE_SHIFT | OCCUPANCY | LOGP | SIGMA | SCORE | PIXEL_SIZE | MICROSCOPE_VOLTAGE | MICROSCOPE_CS | AMPLITUDE_CONTRAST | BEAM_TILT_X | BEAM_TILT_Y | IMAGE_SHIFT_X | IMAGE_SHIFT_Y | ASSIGNED_SUBSET);
    output_params.PreallocateMemoryAndBlank(number_of_input_images); //in case all had occupance > threshold

    cisTEMParameters output_removed_params;
    // setup parameters for the output star file
    output_removed_params.parameters_to_write.SetActiveParameters(POSITION_IN_STACK | IMAGE_IS_ACTIVE | PSI | THETA | PHI | X_SHIFT | Y_SHIFT | DEFOCUS_1 | DEFOCUS_2 | DEFOCUS_ANGLE | PHASE_SHIFT | OCCUPANCY | LOGP | SIGMA | SCORE | PIXEL_SIZE | MICROSCOPE_VOLTAGE | MICROSCOPE_CS | AMPLITUDE_CONTRAST | BEAM_TILT_X | BEAM_TILT_Y | IMAGE_SHIFT_X | IMAGE_SHIFT_Y | ASSIGNED_SUBSET);
    output_removed_params.PreallocateMemoryAndBlank(number_of_input_images); //in case all had occupance > threshold

    cisTEMParameterLine input_parameters_r;
    cisTEMParameterLine reference_parameters_r;

    long         new_counter           = 0;
    long         removed_image_counter = 0;
    ProgressBar* my_progress           = new ProgressBar(number_of_input_images);

    // main loop
    for ( long image_counter = 0; image_counter < number_of_base_images; image_counter++ ) { //number_of_input_images
        input_parameters     = input_star_file.ReturnLine(image_counter); // the star file numbering is 0 indexed!!
        reference_parameters = ref_star_file.ReturnLine(image_counter);

        // input angles normalized into [0,360)
        float input_phi   = angle_within360(input_parameters.phi);
        float input_psi   = angle_within360(input_parameters.psi);
        float input_theta = input_parameters.theta;
        float ref_phi     = angle_within360(reference_parameters.phi); // do I need angle_within360??
        float ref_psi     = reference_parameters.psi;
        float ref_theta   = reference_parameters.theta;

        // // compute adjusted reference euler angles (use adj_ref_phi for comparison)
        // float          adj_ref_phi   = 0.0f;
        // float          adj_ref_theta = 0.0f;
        // float          adj_ref_psi   = 0.0f;
        // RotationMatrix temp_matrix;
        // temp_matrix.SetToEulerRotation(reference_parameters.phi, reference_parameters.theta, reference_parameters.psi);
        // temp_matrix.ConvertToValidEulerAngles(adj_ref_phi, adj_ref_theta, adj_ref_psi); // Is this needed?
        // adj_ref_phi = angle_within360(adj_ref_phi); // make sure adjusted phi is in [0,360)

        // smallest signed difference (wrap-aware)
        //float diff = angle_difference(input_phi, adj_ref_phi);
        float diff = angle_difference(input_phi, ref_phi);

        // For backwards-compatibility with your "removed" debug logic we still check +/-90 and +/-270 offsets and save them to removed list,
        // BUT we DO NOT accept +180 or +90 as valid matches for keeping. We only keep particles whose input_phi is within angular_range of adj_ref_phi.

        if ( classification_results ) {
            // occupancy normalization for particle (support 0..1 or 0..100 values)
            float particle_occupancy = input_parameters.occupancy;
            if ( (particle_occupancy >= occupancy_threshold) && (fabs(diff) <= angular_range) ) {
                // KEEP particle
                my_image.ReadSlice(&my_input_images, image_counter + 1);
                // save the image into the new output MRC file
                my_image.WriteSlice(&my_output_images, new_counter + 1);

                // Print the values for debugging as requested
                wxPrintf("KEPT (image %li): reference_phi(raw) = %f, input_phi = %f, ref_phi = %f\n",
                         image_counter + 1, reference_parameters.phi, input_phi, ref_phi);

                // save the parameter information of the image into the new star file
                // (kept same long assignments as you requested)
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

                new_counter++;
            }
            else {
                // Not kept. For debugging: check +/-90 & +/-270 offsets and save those into removed list, as in original code.
                float diff_90_offset  = angle_difference(input_phi, ref_phi + 90.0f);
                float diff_270_offset = angle_difference(input_phi, ref_phi + 270.0f);

                if ( fabs(diff_90_offset) <= angular_range || fabs(diff_270_offset) <= angular_range ) {
                    my_image.ReadSlice(&my_input_images, image_counter + 1);
                    my_image.WriteSlice(&removed_images_output, removed_image_counter + 1);
                    // calculating the angular difference for debugging output (raw numbers)
                    float angular_difference = input_parameters.phi - reference_parameters.phi;

                    results.push_back({image_counter + 1,
                                       input_parameters.phi,
                                       reference_parameters.phi,
                                       angular_difference,
                                       ref_phi,
                                       input_psi,
                                       reference_parameters.psi,
                                       ref_psi,
                                       input_theta,
                                       reference_parameters.theta,
                                       ref_theta});

                    output_removed_params.all_parameters[removed_image_counter].position_in_stack                  = removed_image_counter + 1;
                    output_removed_params.all_parameters[removed_image_counter].image_is_active                    = input_parameters.image_is_active;
                    output_removed_params.all_parameters[removed_image_counter].psi                                = input_parameters.psi;
                    output_removed_params.all_parameters[removed_image_counter].theta                              = input_parameters.theta;
                    output_removed_params.all_parameters[removed_image_counter].phi                                = input_parameters.phi;
                    output_removed_params.all_parameters[removed_image_counter].x_shift                            = input_parameters.x_shift;
                    output_removed_params.all_parameters[removed_image_counter].y_shift                            = input_parameters.y_shift;
                    output_removed_params.all_parameters[removed_image_counter].defocus_1                          = input_parameters.defocus_1;
                    output_removed_params.all_parameters[removed_image_counter].defocus_2                          = input_parameters.defocus_2;
                    output_removed_params.all_parameters[removed_image_counter].defocus_angle                      = input_parameters.defocus_angle;
                    output_removed_params.all_parameters[removed_image_counter].phase_shift                        = input_parameters.phase_shift;
                    output_removed_params.all_parameters[removed_image_counter].occupancy                          = input_parameters.occupancy;
                    output_removed_params.all_parameters[removed_image_counter].logp                               = input_parameters.logp;
                    output_removed_params.all_parameters[removed_image_counter].sigma                              = input_parameters.sigma;
                    output_removed_params.all_parameters[removed_image_counter].score                              = input_parameters.score;
                    output_removed_params.all_parameters[removed_image_counter].score_change                       = input_parameters.score_change;
                    output_removed_params.all_parameters[removed_image_counter].pixel_size                         = input_parameters.pixel_size;
                    output_removed_params.all_parameters[removed_image_counter].microscope_voltage_kv              = input_parameters.microscope_voltage_kv;
                    output_removed_params.all_parameters[removed_image_counter].microscope_spherical_aberration_mm = input_parameters.microscope_spherical_aberration_mm;
                    output_removed_params.all_parameters[removed_image_counter].amplitude_contrast                 = input_parameters.amplitude_contrast;
                    output_removed_params.all_parameters[removed_image_counter].beam_tilt_x                        = input_parameters.beam_tilt_x;
                    output_removed_params.all_parameters[removed_image_counter].beam_tilt_y                        = input_parameters.beam_tilt_y;
                    output_removed_params.all_parameters[removed_image_counter].image_shift_x                      = input_parameters.image_shift_x;
                    output_removed_params.all_parameters[removed_image_counter].image_shift_y                      = input_parameters.image_shift_y;
                    if ( input_parameters.position_in_stack % 2 == 1 ) {
                        input_parameters.assigned_subset = 1; // Odd particle number
                    }
                    else {
                        input_parameters.assigned_subset = 2; // Even particle number
                    }
                    output_removed_params.all_parameters[removed_image_counter].assigned_subset = input_parameters.assigned_subset;
                    removed_image_counter++;
                }
            }
        }
        else { // not classification_results -> plain branch
            if ( fabs(diff) <= angular_range ) {

                // Keep the first round (r=0) unconditionally. We'll compare subsequent rounds to the first round's input_phi.
                // if true will keep them all for now.
                //if ( true ) {
                // KEEP particle (no mirror / no +180 acceptance)
                my_image.ReadSlice(&my_input_images, image_counter + 1);
                my_image.WriteSlice(&my_output_images, new_counter + 1);

                // Print the values for debugging
                wxPrintf("KEPT (image %li): reference_phi(raw) = %f, input_phi = %f, ref_phi = %f\n",
                         image_counter + 1, reference_parameters.phi, input_phi, ref_phi);

                // save the parameter information of the image into the new star file
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

                new_counter++;

                for ( size_t r = 1; r < num_of_models; r++ ) {
                    input_parameters_r     = input_star_file.ReturnLine((image_counter + (r * number_of_base_images))); // the star file numbering is 0 indexed!!
                    reference_parameters_r = ref_star_file.ReturnLine((image_counter + (r * number_of_base_images)));

                    // input angles normalized into [0,360)
                    float input_phi_r   = angle_within360(input_parameters_r.phi);
                    float input_psi_r   = angle_within360(input_parameters_r.psi);
                    float input_theta_r = input_parameters_r.theta;
                    float ref_phi_r     = angle_within360(reference_parameters_r.phi); // do I need angle_within360??
                    float ref_psi_r     = reference_parameters_r.psi;
                    float ref_theta_r   = reference_parameters_r.theta;

                    // compare with the reference phi
                    // float diff_r = angle_difference(input_phi_r, ref_phi_r);

                    // will consider the first round phi the reference phi that we want the next rounds to be within range of it + expected angular change
                    float input_phi_ref = input_phi + (r * 90.0);
                    float diff_r        = angle_difference(input_phi_r, input_phi_ref);

                    if ( fabs(diff_r) <= angular_range ) {

                        my_image.ReadSlice(&my_input_images, ((image_counter + 1) + (r * number_of_base_images)));
                        my_image.WriteSlice(&my_output_images, new_counter + 1);

                        // Print the values for debugging as requested
                        // wxPrintf("KEPT (image %li): reference_phi(raw) = %f, input_phi = %f, ref_phi = %f\n",
                        //          ((image_counter + 1) + (r * number_of_base_images)), reference_parameters_r.phi, input_phi_r, ref_phi_r);

                        wxPrintf("KEPT (image %li): reference_phi(raw) = %f, input_phi = %f, ref_phi = %f\n",
                                 ((image_counter + 1) + (r * number_of_base_images)), input_parameters.phi, input_phi_r, input_phi_ref);

                        // save the parameter information of the image into the new star file
                        output_params.all_parameters[new_counter].position_in_stack                  = new_counter + 1;
                        output_params.all_parameters[new_counter].image_is_active                    = input_parameters_r.image_is_active;
                        output_params.all_parameters[new_counter].psi                                = input_parameters_r.psi;
                        output_params.all_parameters[new_counter].theta                              = input_parameters_r.theta;
                        output_params.all_parameters[new_counter].phi                                = input_parameters_r.phi;
                        output_params.all_parameters[new_counter].x_shift                            = input_parameters_r.x_shift;
                        output_params.all_parameters[new_counter].y_shift                            = input_parameters_r.y_shift;
                        output_params.all_parameters[new_counter].defocus_1                          = input_parameters_r.defocus_1;
                        output_params.all_parameters[new_counter].defocus_2                          = input_parameters_r.defocus_2;
                        output_params.all_parameters[new_counter].defocus_angle                      = input_parameters_r.defocus_angle;
                        output_params.all_parameters[new_counter].phase_shift                        = input_parameters_r.phase_shift;
                        output_params.all_parameters[new_counter].occupancy                          = input_parameters_r.occupancy;
                        output_params.all_parameters[new_counter].logp                               = input_parameters_r.logp;
                        output_params.all_parameters[new_counter].sigma                              = input_parameters_r.sigma;
                        output_params.all_parameters[new_counter].score                              = input_parameters_r.score;
                        output_params.all_parameters[new_counter].score_change                       = input_parameters_r.score_change;
                        output_params.all_parameters[new_counter].pixel_size                         = input_parameters_r.pixel_size;
                        output_params.all_parameters[new_counter].microscope_voltage_kv              = input_parameters_r.microscope_voltage_kv;
                        output_params.all_parameters[new_counter].microscope_spherical_aberration_mm = input_parameters_r.microscope_spherical_aberration_mm;
                        output_params.all_parameters[new_counter].amplitude_contrast                 = input_parameters_r.amplitude_contrast;
                        output_params.all_parameters[new_counter].beam_tilt_x                        = input_parameters_r.beam_tilt_x;
                        output_params.all_parameters[new_counter].beam_tilt_y                        = input_parameters_r.beam_tilt_y;
                        output_params.all_parameters[new_counter].image_shift_x                      = input_parameters_r.image_shift_x;
                        output_params.all_parameters[new_counter].image_shift_y                      = input_parameters_r.image_shift_y;
                        if ( input_parameters.position_in_stack % 2 == 1 ) { // on purpose left as input_parameters not input_parameters_r as I want all the RASTR particles from the same image to be in the same group
                            input_parameters_r.assigned_subset = 1; // Odd particle number
                        }
                        else {
                            input_parameters_r.assigned_subset = 2; // Even particle number
                        }
                        output_params.all_parameters[new_counter].assigned_subset = input_parameters_r.assigned_subset; // There is no assigned subset and we need to keep the saved assigned subset as is to ensure that no 2 particles are in the same group even after filteration

                        new_counter++;
                    }
                    else { // if the angle in next rounds not within range then save it in the removed particle stack
                        my_image.ReadSlice(&my_input_images, image_counter + 1);
                        my_image.WriteSlice(&removed_images_output, removed_image_counter + 1);
                        // calculating the angular difference for debugging output (raw numbers)
                        float angular_difference = input_parameters.phi - reference_parameters.phi;

                        results.push_back({image_counter + 1,
                                           input_parameters.phi,
                                           reference_parameters.phi,
                                           angular_difference,
                                           ref_phi,
                                           input_psi,
                                           reference_parameters.psi,
                                           ref_psi,
                                           input_theta,
                                           reference_parameters.theta,
                                           ref_theta});

                        output_removed_params.all_parameters[removed_image_counter].position_in_stack                  = removed_image_counter + 1;
                        output_removed_params.all_parameters[removed_image_counter].image_is_active                    = input_parameters.image_is_active;
                        output_removed_params.all_parameters[removed_image_counter].psi                                = input_parameters.psi;
                        output_removed_params.all_parameters[removed_image_counter].theta                              = input_parameters.theta;
                        output_removed_params.all_parameters[removed_image_counter].phi                                = input_parameters.phi;
                        output_removed_params.all_parameters[removed_image_counter].x_shift                            = input_parameters.x_shift;
                        output_removed_params.all_parameters[removed_image_counter].y_shift                            = input_parameters.y_shift;
                        output_removed_params.all_parameters[removed_image_counter].defocus_1                          = input_parameters.defocus_1;
                        output_removed_params.all_parameters[removed_image_counter].defocus_2                          = input_parameters.defocus_2;
                        output_removed_params.all_parameters[removed_image_counter].defocus_angle                      = input_parameters.defocus_angle;
                        output_removed_params.all_parameters[removed_image_counter].phase_shift                        = input_parameters.phase_shift;
                        output_removed_params.all_parameters[removed_image_counter].occupancy                          = input_parameters.occupancy;
                        output_removed_params.all_parameters[removed_image_counter].logp                               = input_parameters.logp;
                        output_removed_params.all_parameters[removed_image_counter].sigma                              = input_parameters.sigma;
                        output_removed_params.all_parameters[removed_image_counter].score                              = input_parameters.score;
                        output_removed_params.all_parameters[removed_image_counter].score_change                       = input_parameters.score_change;
                        output_removed_params.all_parameters[removed_image_counter].pixel_size                         = input_parameters.pixel_size;
                        output_removed_params.all_parameters[removed_image_counter].microscope_voltage_kv              = input_parameters.microscope_voltage_kv;
                        output_removed_params.all_parameters[removed_image_counter].microscope_spherical_aberration_mm = input_parameters.microscope_spherical_aberration_mm;
                        output_removed_params.all_parameters[removed_image_counter].amplitude_contrast                 = input_parameters.amplitude_contrast;
                        output_removed_params.all_parameters[removed_image_counter].beam_tilt_x                        = input_parameters.beam_tilt_x;
                        output_removed_params.all_parameters[removed_image_counter].beam_tilt_y                        = input_parameters.beam_tilt_y;
                        output_removed_params.all_parameters[removed_image_counter].image_shift_x                      = input_parameters.image_shift_x;
                        output_removed_params.all_parameters[removed_image_counter].image_shift_y                      = input_parameters.image_shift_y;
                        if ( input_parameters.position_in_stack % 2 == 1 ) {
                            input_parameters.assigned_subset = 1; // Odd particle number
                        }
                        else {
                            input_parameters.assigned_subset = 2; // Even particle number
                        }
                        output_removed_params.all_parameters[removed_image_counter].assigned_subset = input_parameters.assigned_subset;
                        removed_image_counter++;
                    }
                }
            }
            else {
                // Not kept. For debugging: check +/-90 & +/-270 offsets and save those into removed list, as in original code.
                // float diff_90_offset  = angle_difference(input_phi, ref_phi + 90.0f);
                // float diff_270_offset = angle_difference(input_phi, ref_phi + 270.0f);

                // if ( fabs(diff_90_offset) <= angular_range || fabs(diff_270_offset) <= angular_range ) {
                my_image.ReadSlice(&my_input_images, image_counter + 1);
                my_image.WriteSlice(&removed_images_output, removed_image_counter + 1);
                // calculating the angular difference for debugging output (raw numbers)
                float angular_difference = input_parameters.phi - reference_parameters.phi;

                results.push_back({image_counter + 1,
                                   input_parameters.phi,
                                   reference_parameters.phi,
                                   angular_difference,
                                   ref_phi,
                                   input_psi,
                                   reference_parameters.psi,
                                   ref_psi,
                                   input_theta,
                                   reference_parameters.theta,
                                   ref_theta});

                output_removed_params.all_parameters[removed_image_counter].position_in_stack                  = removed_image_counter + 1;
                output_removed_params.all_parameters[removed_image_counter].image_is_active                    = input_parameters.image_is_active;
                output_removed_params.all_parameters[removed_image_counter].psi                                = input_parameters.psi;
                output_removed_params.all_parameters[removed_image_counter].theta                              = input_parameters.theta;
                output_removed_params.all_parameters[removed_image_counter].phi                                = input_parameters.phi;
                output_removed_params.all_parameters[removed_image_counter].x_shift                            = input_parameters.x_shift;
                output_removed_params.all_parameters[removed_image_counter].y_shift                            = input_parameters.y_shift;
                output_removed_params.all_parameters[removed_image_counter].defocus_1                          = input_parameters.defocus_1;
                output_removed_params.all_parameters[removed_image_counter].defocus_2                          = input_parameters.defocus_2;
                output_removed_params.all_parameters[removed_image_counter].defocus_angle                      = input_parameters.defocus_angle;
                output_removed_params.all_parameters[removed_image_counter].phase_shift                        = input_parameters.phase_shift;
                output_removed_params.all_parameters[removed_image_counter].occupancy                          = input_parameters.occupancy;
                output_removed_params.all_parameters[removed_image_counter].logp                               = input_parameters.logp;
                output_removed_params.all_parameters[removed_image_counter].sigma                              = input_parameters.sigma;
                output_removed_params.all_parameters[removed_image_counter].score                              = input_parameters.score;
                output_removed_params.all_parameters[removed_image_counter].score_change                       = input_parameters.score_change;
                output_removed_params.all_parameters[removed_image_counter].pixel_size                         = input_parameters.pixel_size;
                output_removed_params.all_parameters[removed_image_counter].microscope_voltage_kv              = input_parameters.microscope_voltage_kv;
                output_removed_params.all_parameters[removed_image_counter].microscope_spherical_aberration_mm = input_parameters.microscope_spherical_aberration_mm;
                output_removed_params.all_parameters[removed_image_counter].amplitude_contrast                 = input_parameters.amplitude_contrast;
                output_removed_params.all_parameters[removed_image_counter].beam_tilt_x                        = input_parameters.beam_tilt_x;
                output_removed_params.all_parameters[removed_image_counter].beam_tilt_y                        = input_parameters.beam_tilt_y;
                output_removed_params.all_parameters[removed_image_counter].image_shift_x                      = input_parameters.image_shift_x;
                output_removed_params.all_parameters[removed_image_counter].image_shift_y                      = input_parameters.image_shift_y;
                if ( input_parameters.position_in_stack % 2 == 1 ) {
                    input_parameters.assigned_subset = 1; // Odd particle number
                }
                else {
                    input_parameters.assigned_subset = 2; // Even particle number
                }
                output_removed_params.all_parameters[removed_image_counter].assigned_subset = input_parameters.assigned_subset;
                removed_image_counter++;
                //}
            }
        }

        if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
            my_progress->Update(image_counter + 1);
    }
    delete my_progress;

    removed_images_output.WriteHeader( );
    // write the output star file for the matched references
    output_params.WriteTocisTEMStarFile(output_star_filename);
    output_removed_params.WriteTocisTEMStarFile("removed_particles_stack.star");

    // Write all results
    for ( const auto& r : results ) {
        angles_file << r.image_index << ", " << r.input_phi << ", " << r.reference_phi << ", " << r.difference << ", " << r.adj_ref_phi << ", " << r.input_psi << ", " << r.ref_psi << ", " << r.adj_ref_psi << ", " << r.input_theta << ", " << r.adj_ref_theta << ", " << r.ref_theta << "\n";
    }

    float percent_kept     = (float(new_counter) / float(number_of_input_images)) * 100.0f;
    float percent_filtered = ((float(number_of_input_images) - float(new_counter)) / number_of_input_images) * 100.0f;
    wxPrintf("\n\n%.3f %% (= %li particles) are kept, and %.3f %% (= %li particles) are filtered out successfully\n", percent_kept, new_counter, percent_filtered, number_of_input_images - new_counter);

    return true;
}

// Function to ensure the angle is within the range [0, 360)
float angle_within360(float angle) {
    // iterative implementation using fmodf
    float a = fmodf(angle, 360.0f);
    if ( a < 0.0f )
        a += 360.0f;
    return a;
}

// Compute smallest signed difference between two angles in (-180, 180] to wrap-around boundary (0°/360°)
float angle_difference(float a, float b) {
    float diff = fmodf(a - b + 540.0f, 360.0f) - 180.0f; // +540 to ensure that the angle is always positive before modulus
    return diff;
}
