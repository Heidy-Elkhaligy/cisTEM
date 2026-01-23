#include "../../core/core_headers.h"

class
        randomize_rastr_phi : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

IMPLEMENT_APP(randomize_rastr_phi)

// override the DoInteractiveUserInput
float              angle_within360(float angle);
std::vector<float> GenerateRandomAngles(size_t count);
std::vector<float> GenerateRandomAnglesWithinRange(size_t count, float angle_range);
std::vector<float> GenerateRandomAnglesWithinRangeFrontBack(size_t count, float angle_range);

void randomize_rastr_phi::DoInteractiveUserInput( ) {

    UserInput* my_input = new UserInput("randomize_rastr_phi", 1.0);

    wxString input_filename       = my_input->GetFilenameFromUser("Input image stack", "Filename of the stack of images to normalize", "input_stack1.mrc", true);
    wxString input_star_filename  = my_input->GetFilenameFromUser("Input star file", "The input star file", "my_parameters.star", true);
    wxString output_star_filename = my_input->GetFilenameFromUser("Output star file", "Output star file with the random phi + reference RASTR phi", "output_parameters.star", false);
    int      number_of_models     = my_input->GetIntFromUser("Number of models used to generate that RASTR particle stack", "The number of models used to generate the RASTR particle stack, minimum number is 1 models", "1", 1);
    delete my_input;

    my_current_job.Reset(5);
    my_current_job.ManualSetArguments("ttti", input_filename.ToUTF8( ).data( ), input_star_filename.ToUTF8( ).data( ), output_star_filename.ToUTF8( ).data( ), number_of_models);
}

// override the do calculation method which will be what is actually run..

bool randomize_rastr_phi::DoCalculation( ) {
    wxString input_filename       = my_current_job.arguments[0].ReturnStringArgument( );
    wxString input_star_filename  = my_current_job.arguments[1].ReturnStringArgument( );
    wxString output_star_filename = my_current_job.arguments[2].ReturnStringArgument( );
    int      number_of_models     = my_current_job.arguments[3].ReturnIntegerArgument( );

    MRCFile my_input_file(input_filename.ToStdString( ), false);

    cisTEMParameters input_star_file;
    if ( (is_running_locally && ! DoesFileExist(input_star_filename.ToStdString( ))) ) {
        SendErrorAndCrash(wxString::Format("Error: Input star file %s not found\n", input_star_filename));
    }
    //CisTEM star
    input_star_file.ReadFromcisTEMStarFile(input_star_filename.ToStdString( ));

    long  number_of_input_images    = my_input_file.ReturnNumberOfSlices( );
    long  number_of_original_images = number_of_input_images / number_of_models;
    float phi_step                  = 360.0 / number_of_models;

    //auto random_phi_angles = GenerateRandomAngles(number_of_original_images);
    auto random_phi_angles = GenerateRandomAnglesWithinRange(number_of_original_images, 45.0); // change 45.0 to something user specific

    //auto random_phi_angles = GenerateRandomAnglesWithinRangeFrontBack(number_of_original_images, 45.0);

    wxPrintf("\nNumber of RASTR images is %li and number of actual images is %li.\n\n", number_of_input_images, number_of_input_images / number_of_models);

    cisTEMParameterLine input_parameters;
    cisTEMParameters    output_params;

    output_params.parameters_to_write.SetActiveParameters(POSITION_IN_STACK | IMAGE_IS_ACTIVE | PSI | THETA | PHI | X_SHIFT | Y_SHIFT | DEFOCUS_1 | DEFOCUS_2 | DEFOCUS_ANGLE | PHASE_SHIFT | OCCUPANCY | LOGP | SIGMA | SCORE | PIXEL_SIZE | MICROSCOPE_VOLTAGE | MICROSCOPE_CS | AMPLITUDE_CONTRAST | BEAM_TILT_X | BEAM_TILT_Y | IMAGE_SHIFT_X | IMAGE_SHIFT_Y | ASSIGNED_SUBSET);
    output_params.PreallocateMemoryAndBlank(number_of_input_images); //in case all had occupance > threshold

    // Image my_image;
    // float input_pixel_size = my_input_file.ReturnPixelSize( );

    wxPrintf("\nRandomizing RASTR Phi Angles..\n\n");
    ProgressBar* my_progress = new ProgressBar(number_of_input_images);

    // main loop
    for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {

        long  original_index = image_counter % number_of_original_images;
        float random_phi     = random_phi_angles[original_index];

        input_parameters = input_star_file.ReturnLine(image_counter); // the star file numbering is 0 indexed!!

        // input phi angles are the reference phi used to generate RASTR ex. for 4 models the input phi should be either 0,90,180, or 270
        float input_phi = input_parameters.phi;
        // now add the reference phi to the randomly generated phi
        // this will ensure each image that was reproduced by RASTR will have the same random phi with the extra RASTR extraction angle so they are all connected correctly
        // the angle will be wrapped again to be within 360
        float output_random_phi = angle_within360(input_phi + random_phi);

        //wxPrintf("Image %li input reference phi is %f, random phi is %f, and the final phi is %f, wraped back to 360 %f\n", image_counter + 1, input_phi, random_phi, output_random_phi, angle_within360(output_random_phi));

        // save the parameter information of the image into the new star file
        // (kept same long assignments as you requested)
        output_params.all_parameters[image_counter].position_in_stack                  = image_counter + 1;
        output_params.all_parameters[image_counter].image_is_active                    = input_parameters.image_is_active;
        output_params.all_parameters[image_counter].psi                                = input_parameters.psi;
        output_params.all_parameters[image_counter].theta                              = input_parameters.theta;
        output_params.all_parameters[image_counter].phi                                = output_random_phi;
        output_params.all_parameters[image_counter].x_shift                            = input_parameters.x_shift;
        output_params.all_parameters[image_counter].y_shift                            = input_parameters.y_shift;
        output_params.all_parameters[image_counter].defocus_1                          = input_parameters.defocus_1;
        output_params.all_parameters[image_counter].defocus_2                          = input_parameters.defocus_2;
        output_params.all_parameters[image_counter].defocus_angle                      = input_parameters.defocus_angle;
        output_params.all_parameters[image_counter].phase_shift                        = input_parameters.phase_shift;
        output_params.all_parameters[image_counter].occupancy                          = input_parameters.occupancy;
        output_params.all_parameters[image_counter].logp                               = input_parameters.logp;
        output_params.all_parameters[image_counter].sigma                              = input_parameters.sigma;
        output_params.all_parameters[image_counter].score                              = input_parameters.score;
        output_params.all_parameters[image_counter].score_change                       = input_parameters.score_change;
        output_params.all_parameters[image_counter].pixel_size                         = input_parameters.pixel_size;
        output_params.all_parameters[image_counter].microscope_voltage_kv              = input_parameters.microscope_voltage_kv;
        output_params.all_parameters[image_counter].microscope_spherical_aberration_mm = input_parameters.microscope_spherical_aberration_mm;
        output_params.all_parameters[image_counter].amplitude_contrast                 = input_parameters.amplitude_contrast;
        output_params.all_parameters[image_counter].beam_tilt_x                        = input_parameters.beam_tilt_x;
        output_params.all_parameters[image_counter].beam_tilt_y                        = input_parameters.beam_tilt_y;
        output_params.all_parameters[image_counter].image_shift_x                      = input_parameters.image_shift_x;
        output_params.all_parameters[image_counter].image_shift_y                      = input_parameters.image_shift_y;
        if ( input_parameters.position_in_stack % 2 == 1 ) {
            input_parameters.assigned_subset = 1; // Odd particle number
        }
        else {
            input_parameters.assigned_subset = 2; // Even particle number
        }
        output_params.all_parameters[image_counter].assigned_subset = input_parameters.assigned_subset; // There is no assigned subset and we need to keep the saved assigned subset as is to ensure that no 2 particles are in the same group even after filteration

        my_progress->Update(image_counter + 1);
    }
    delete my_progress;
    output_params.WriteTocisTEMStarFile(output_star_filename);

    wxPrintf("\n\n");

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

std::vector<float> GenerateRandomAngles(size_t count) {
    std::vector<float> values;
    values.reserve(count);

    // Seed with real entropy → different every run
    std::mt19937 rng(std::random_device{ }( ));

    std::uniform_real_distribution<float> dist(0.0f, 360.0f);

    for ( size_t i = 0; i < count; ++i ) {
        values.push_back(dist(rng));
    }

    return values;
}

std::vector<float> GenerateRandomAnglesWithinRange(size_t count, float angle_range) {
    std::vector<float> values;
    values.reserve(count);

    // Seed with real entropy → different every run
    std::mt19937 rng(std::random_device{ }( ));

    // Uniform distribution between +/- range around 0
    float min_angle = 0.0 - angle_range;
    float max_angle = 0.0 + angle_range;

    std::uniform_real_distribution<float> dist(min_angle, max_angle);

    for ( size_t i = 0; i < count; ++i ) {
        values.push_back(dist(rng));
    }

    return values;
}

std::vector<float> GenerateRandomAnglesWithinRangeFrontBack(size_t count, float angle_range) {
    std::vector<float> values;
    values.reserve(count);

    std::mt19937 rng(std::random_device{ }( ));

    size_t half = count / 2;

    float min_angle_front = 0.0 - angle_range;
    float max_angle_front = 0.0 + angle_range;

    float min_angle_back = 180.0 - angle_range;
    float max_angle_back = 180.0 + angle_range;

    std::uniform_real_distribution<float> around_zero(min_angle_front, max_angle_front);
    std::uniform_real_distribution<float> around_180(min_angle_back, max_angle_back);

    // First half
    for ( size_t i = 0; i < half; ++i ) {
        values.push_back(around_zero(rng));
    }

    // Second half
    for ( size_t i = half; i < count; ++i ) {
        values.push_back(around_180(rng));
    }

    // Shuffle so order is random
    std::shuffle(values.begin( ), values.end( ), rng);

    return values;
}