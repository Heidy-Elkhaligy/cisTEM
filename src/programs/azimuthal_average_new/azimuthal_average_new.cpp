#include "../../core/core_headers.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>


// check scaling 

class
        AzimuthalAverageNew : public MyApp {

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


std::vector<float> sum_image_columns(Image* current_image);
float sum_image_columns_float(Image* current_image);
std::pair<int, int>  find_column_sum_peaks(Image* current_image, float min_gap = 0.0);
std::pair<int, int>  find_row_sum_peaks(Image* current_image, float min_gap = 0.0);
std::vector<float> average_image_columns(Image* current_image);


void  normalize_image(Image* input_image, float pixel_size, float mask_falloff);
void invert_mask(Image* mask_file);


// new way
void create_white_sphere_mask(Image* mask_file, int x_sphere_center, int y_sphere_center, int z_spehere_center, float radius);
void create_black_sphere_mask(Image* mask_file, int x_sphere_center, int y_sphere_center, int z_spehere_center, float radius);

//Functions for the average images bins
void InitializeCTFSumOfSquares(int numBins, Image& current_image, std::vector<std::vector<float>>* ctf_sum_of_squares);
void ApplyCTFAndReturnCTFSumOfSquares(Image& image, CTF ctf_to_apply, bool absolute, bool apply_beam_tilt, bool apply_envelope, std::vector<float>& ctf_sum_of_squares);
void divide_by_ctf_sum_of_squares(Image& current_image, std::vector<float>& ctf_sum_of_squares);
void sum_image_direction(Image* current_image, int dim);
void apply_ctf(Image* current_image, CTF ctf_to_apply, float* ctf_sum_of_squares, bool absolute, bool do_fill_sum_of_squares);
//float angle_within360(float angle);
float ReturnAverageOfRealValuesOnVerticalEdges(Image* current_image);
float calculateRMSD(const std::vector<float>& vec1, const std::vector<float>& vec2);
std::vector<float> GetAverageColumnsValues(const std::vector<std::vector<float>>& average_columns_data);

IMPLEMENT_APP(AzimuthalAverageNew)

// override the DoInteractiveUserInput

void AzimuthalAverageNew::DoInteractiveUserInput( ) {
    // intial parameters
    //int         new_z_size = 1;
    //wxString    text_filename;
    float       pixel_size;
    wxString    output_filename;

    // ctf parameters
    wxString input_star_filename;
    float acceleration_voltage;
    float spherical_aberration;
    float amplitude_contrast;
    float defocus_1               = 0.0;
    float defocus_2               = 0.0;
    float astigmatism_angle       = 0.0;
    float additional_phase_shift  = 0.0;
    bool  input_ctf_values_from_star_file;
    bool  phase_flip_only;
    wxString output_average_per_bin_filename;
    wxString output_azimuthal_average_volume_filename;

    // tube searching parameters
    float min_tube_diameter       = 0.0;
    float max_tube_diameter       = 0.0;
    int   bins_count              = 1; // The number of bins to use when clssifying images based on tube diameter, the min number is 1 class
    int   outer_mask_radius       = 0;
    bool  tubes_centered;
    bool  low_pass;

    // RASTR mask parameters
    bool        RASTR = false;
    bool        input_mask = false;
    wxString    input_mask_filename;
    int         x_mask_center = 1;
    int         y_mask_center = 1;
    int         z_mask_center = 1;
    int         sphere_mask_radius = 1;
    int         number_of_models  = 1;
    bool        center_upweighted = false;
    float       cosine_edge       = 10.0;
    float       outside_weight    = 0.0;
    float       filter_radius     = 0.0;
    float       outside_value     = 0.0;
    bool        use_outside_value = false;
    wxString    RASTR_output_filename;
    wxString    RASTR_output_star_filename;

    // SPOT RASTR
    bool        SPOT_RASTR = false;
    wxString    SPOT_RASTR_output_filename;
    wxString    SPOT_RASTR_output_star_filename;

    // expert options
    bool set_expert_options;
    // rotation (degrees)
    float psi_min                 = 0.0;
    float psi_max                 = 180.0;
    float psi_step                = 5.0;
    float fine_tuning_psi_step    = 0.25;
    float padding_factor          = sqrtf(2); //(sqrt of 2)

    // bool set_expert_options;
    int  max_threads;

    UserInput* my_input = new UserInput("AzimuthalAverageNew", 1.0);

    wxString input_filename         = my_input->GetFilenameFromUser("Input image file name", "Filename of input stack", "input_stack.mrc", true);
      // get CTF from user
    pixel_size                      = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
    acceleration_voltage            = my_input->GetFloatFromUser("Acceleration voltage (keV)", "Acceleration voltage, in keV", "300.0", 0.0, 500.0);
    spherical_aberration            = my_input->GetFloatFromUser("Spherical aberration (mm)", "Objective lens spherical aberration", "2.7", 0.0);
    amplitude_contrast              = my_input->GetFloatFromUser("Amplitude contrast", "Fraction of total contrast attributed to amplitude contrast", "0.07", 0.0);

    input_ctf_values_from_star_file = my_input->GetYesNoFromUser("Use a star file to input defocus values?", "If yes, defocus values will be extracted from star file", "NO");

    if ( input_ctf_values_from_star_file == true ) {
        input_star_filename             = my_input->GetFilenameFromUser("Input star file", "The input star file", "my_parameters.star", true);
    }
    else {
        defocus_1                   = my_input->GetFloatFromUser("Underfocus 1 (A)", "In Angstroms, the objective lens underfocus along the first axis", "1.2");
        defocus_2                   = my_input->GetFloatFromUser("Underfocus 2 (A)", "In Angstroms, the objective lens underfocus along the second axis", "1.2");
        astigmatism_angle           = my_input->GetFloatFromUser("Astigmatism angle", "Angle between the first axis and the x axis of the image", "0.0");
        additional_phase_shift      = my_input->GetFloatFromUser("Additional phase shift (rad)", "Additional phase shift relative to undiffracted beam, as introduced for example by a phase plate", "0.0");
    }

    phase_flip_only                 = my_input->GetYesNoFromUser("Phase Flip Only", "If Yes, only phase flipping is performed", "NO");
    
    // tube diameter arguments 
    min_tube_diameter               = my_input->GetFloatFromUser("Minimum tube diameter", "The minimum tube diameter for searching and bining tubes", "30.0", 0.0);
    max_tube_diameter               = my_input->GetFloatFromUser("Maximum tube diameter", "The maximum tube diameter for searching and bining tubes", "60.0", 0.0);
    bins_count                      = my_input->GetIntFromUser("Number of classes for classifying the tube diameters", "The number of classes (bins) to classify the tubes based on specified min. and max. diameter", "1", 1);
    outer_mask_radius               = my_input->GetIntFromUser("Outer mask radius for tube search (pixels)", "Outer mask radius to use when searching for tubes in pixels", "0", 0);
    tubes_centered                  = my_input->GetYesNoFromUser("Are tubes centered?", "If yes, the outer radius of peak search from cross-correlation will be 1/5 x-dimension", "NO");
    low_pass                        = my_input->GetYesNoFromUser("Apply low-pass gaussian filter when aligning?", "If yes, will apply a gaussian low pass filter to the original images before aligining the azimuthal average projection", "NO");

    // output arguments for azimuthal average volume and class projections
    output_average_per_bin_filename          = my_input->GetFilenameFromUser("Output name for the average images per class stack","The output name for the average images generated per class MRC file", "output_average_per_class.mrc", false );
    output_azimuthal_average_volume_filename = my_input->GetFilenameFromUser("Output name for the azimuthal average volume per class stack","The output name for the azimuthal average volume generated per class MRC file", "output_azimuthal_average_volume.mrc", false );

    //RASTR mask parameters
    RASTR                           = my_input->GetYesNoFromUser("Perform RASTR masking to 3D azimuthal average?", "Mask out region of interest", "NO");
    // get mask properties from user
    if ( RASTR == true ) {
        input_mask                   = my_input->GetYesNoFromUser("Use input mask file?", "Do you want to provide an input mask file?", "No");
        if (input_mask == true) {
            input_mask_filename          = my_input->GetFilenameFromUser("Input mask file name", "The mask to be applied to the 3D azimuthal average model", "my_mask.mrc", true);
        }
        x_mask_center                = my_input->GetIntFromUser("X center of the mask (pixels)", "The X center of the provided mask or the created mask by the program, default is 0.75 * box size", "1", 1);
        y_mask_center                = my_input->GetIntFromUser("Y center pf the mask (pixels)", "The Y center of the provided mask or the created mask by the program, default is 0.5 * box size", "1", 1);
        z_mask_center                = my_input->GetIntFromUser("Z center of the mask (pixels)", "The Z center of the provided mask or the created mask by the program, default is 0.5 * box size", "1", 1);
        sphere_mask_radius           = my_input->GetIntFromUser("Sphere mask radius (pixels)", "The radius of the provided mask or the created mask by the program, default is 0.3 * box size", "1", 1);
        number_of_models             = my_input->GetIntFromUser("Number of models to generate", "3D azimuthal average model is rotated in increments of (360/n) degrees", "4", 1);
        center_upweighted            = my_input->GetYesNoFromUser("Center the masked upweighted regions?", "Do you want to center the masked upweighted regions?", "No");
        cosine_edge                  = my_input->GetFloatFromUser("Width of cosine edge (A)", "Width of the smooth edge to add to the mask in Angstroms", "10.0", 0.0);
        outside_weight               = my_input->GetFloatFromUser("Weight of density outside mask", "Factor to multiply density outside of the mask", "0.0", 0.0, 1.0);
        filter_radius                = my_input->GetFloatFromUser("Low-pass filter outside mask (A)", "Low-pass filter to be applied to the density outside the mask", "0.0", 0.0);
        outside_value                = my_input->GetFloatFromUser("Outside mask value", "Value used to set density outside the mask", "0.0", 0.0);
        use_outside_value            = my_input->GetYesNoFromUser("Use outside mask value", "Should the density outside the mask be set to the user-provided value", "No");
        RASTR_output_filename        = my_input->GetFilenameFromUser("Output masked upweighted regions filename", "The output MRC file containing average subtracted masked upweighted regions", "masked_upweighted_regions_output_filename.mrc", false);
        RASTR_output_star_filename   = my_input->GetFilenameFromUser("Output star file with the psi and shift changes", "The output star file, containing the new psi and shift parameters", "my_RASTR_output_parameters.star", false);

    }

    SPOT_RASTR              = my_input->GetYesNoFromUser("Perform SPOT-RASTR tube subtraction?", "If yes, azimuthal average model projections will be subtracted from original images", "NO");
    if (SPOT_RASTR == true) {
        SPOT_RASTR_output_filename       = my_input->GetFilenameFromUser("Output average subtracted stack filename", "The output average subtracted MRC file", "average_subtracted_output_filename.mrc", false);
        SPOT_RASTR_output_star_filename  = my_input->GetFilenameFromUser("Output star file with the psi and shift changes", "The output star file, containing the new psi and shift parameters", "my_SPOT_RASTR_output_parameters.star", false);

    }

    if ((RASTR == true) || (SPOT_RASTR == true)) {
    }
    // expert options parameters
    set_expert_options              = my_input->GetYesNoFromUser("Set Expert Options?", "Set these for more control, hopefully not needed", "NO");

    // set alignment options from user
    if ( set_expert_options == true ) {
        psi_min                          = my_input->GetFloatFromUser("Minimum rotation for initial search (degrees)", "The minimum angle rotation of initial search will be limited to this value.", "0.0", -180.0);
        psi_max                          = my_input->GetFloatFromUser("Maximum rotation for initial search (degrees)", "The maximum angle rotation of initial search will be limited to this value.", "180.0", 180.0);
        psi_step                         = my_input->GetFloatFromUser("Rotation step size (degrees)", "The step size of each rotation will be limited to this value.", "5.0", 0.0);
        fine_tuning_psi_step             = my_input->GetFloatFromUser("Local refinement angle rotation step size (degrees)", "The local refinement search step size will be this value", "0.25", 0.0);
        padding_factor                   = my_input->GetFloatFromUser("Padding factor", "Factor determining how much the average image is padded to improve subtraction, defau;t is sqrt(2)", "1.4", 1.0);
    }
    
    //output_filename                 = my_input->GetFilenameFromUser("Output average subtracted stack filename", "The output average subtracted MRC file", "average_subtracted_output_filename.mrc", false);




#ifdef _OPENMP
    max_threads                    = my_input->GetIntFromUser("Max. threads to use for calculation", "when threading, what is the max threads to run", "1", 1);
#else
    max_threads = 1;
#endif

    delete my_input;

    my_current_job.Reset(48);
    my_current_job.ManualSetArguments("tffffbtffffbffiibbttbbtiiiiibffffbttbttbfffffi", input_filename.ToUTF8( ).data( ),
                                        pixel_size,
                                        acceleration_voltage,
                                        spherical_aberration,
                                        amplitude_contrast,
                                        input_ctf_values_from_star_file,
                                        input_star_filename.ToUTF8( ).data( ),
                                        defocus_1,
                                        defocus_2,
                                        astigmatism_angle,
                                        additional_phase_shift,
                                        phase_flip_only,
                                        min_tube_diameter,
                                        max_tube_diameter,
                                        bins_count,
                                        outer_mask_radius,
                                        tubes_centered,
                                        low_pass,
                                        output_average_per_bin_filename.ToUTF8( ).data( ),
                                        output_azimuthal_average_volume_filename.ToUTF8( ).data( ),
                                        RASTR,
                                        input_mask,
                                        input_mask_filename.ToUTF8( ).data( ),
                                        x_mask_center,
                                        y_mask_center,
                                        z_mask_center,
                                        sphere_mask_radius,
                                        number_of_models,
                                        center_upweighted,
                                        cosine_edge,
                                        outside_weight,
                                        filter_radius,
                                        outside_value,
                                        use_outside_value,
                                        RASTR_output_filename.ToUTF8( ).data( ),
                                        RASTR_output_star_filename.ToUTF8( ).data( ),
                                        SPOT_RASTR,
                                        SPOT_RASTR_output_filename.ToUTF8( ).data( ),
                                        SPOT_RASTR_output_star_filename.ToUTF8( ).data( ),
                                        set_expert_options,
                                        psi_min,
                                        psi_max,
                                        psi_step,
                                        fine_tuning_psi_step,
                                        padding_factor,
                                        max_threads);
}

// override the do calculation method which will be what is actually run..

bool AzimuthalAverageNew::DoCalculation( ) {
    // get the arguments for this job..
    wxString    input_filename                           = my_current_job.arguments[0].ReturnStringArgument( );
    float       pixel_size                               = my_current_job.arguments[1].ReturnFloatArgument( );
    float       acceleration_voltage                     = my_current_job.arguments[2].ReturnFloatArgument( );
    float       spherical_aberration                     = my_current_job.arguments[3].ReturnFloatArgument( );
    float       amplitude_contrast                       = my_current_job.arguments[4].ReturnFloatArgument( );
    bool        input_ctf_values_from_star_file          = my_current_job.arguments[5].ReturnBoolArgument( );
    wxString    input_star_filename                      = my_current_job.arguments[6].ReturnStringArgument( ); 
    float       defocus_1                                = my_current_job.arguments[7].ReturnFloatArgument( );
    float       defocus_2                                = my_current_job.arguments[8].ReturnFloatArgument( );
    float       astigmatism_angle                        = my_current_job.arguments[9].ReturnFloatArgument( );
    float       additional_phase_shift                   = my_current_job.arguments[10].ReturnFloatArgument( );
    bool        phase_flip_only                          = my_current_job.arguments[11].ReturnBoolArgument( );
    float       min_tube_diameter                        = my_current_job.arguments[12].ReturnFloatArgument( );
    float       max_tube_diameter                        = my_current_job.arguments[13].ReturnFloatArgument( );
    int         bins_count                               = my_current_job.arguments[14].ReturnIntegerArgument( );
    int         outer_mask_radius                        = my_current_job.arguments[15].ReturnIntegerArgument( );
    bool        tubes_centered                           = my_current_job.arguments[16].ReturnBoolArgument( );
    bool        low_pass                                 = my_current_job.arguments[17].ReturnBoolArgument( );
    wxString    output_average_per_bin_filename          = my_current_job.arguments[18].ReturnStringArgument( );
    wxString    output_azimuthal_average_volume_filename = my_current_job.arguments[19].ReturnStringArgument( );    
    bool        RASTR                                    = my_current_job.arguments[20].ReturnBoolArgument( );
    bool        input_mask                               = my_current_job.arguments[21].ReturnBoolArgument( );
    wxString    input_mask_filename                      = my_current_job.arguments[22].ReturnStringArgument( );
    int         x_mask_center                            = my_current_job.arguments[23].ReturnIntegerArgument( );
    int         y_mask_center                            = my_current_job.arguments[24].ReturnIntegerArgument( );
    int         z_mask_center                            = my_current_job.arguments[25].ReturnIntegerArgument( );
    int         sphere_mask_radius                       = my_current_job.arguments[26].ReturnIntegerArgument( );
    int         number_of_models                         = my_current_job.arguments[27].ReturnIntegerArgument( );
    bool        center_upweighted                        = my_current_job.arguments[28].ReturnBoolArgument( );
    float       cosine_edge                              = my_current_job.arguments[29].ReturnFloatArgument( );
    float       outside_weight                           = my_current_job.arguments[30].ReturnFloatArgument( );
    float       filter_radius                            = my_current_job.arguments[31].ReturnFloatArgument( );
    float       outside_value                            = my_current_job.arguments[32].ReturnFloatArgument( );
    bool        use_outside_value                        = my_current_job.arguments[33].ReturnBoolArgument( );
    wxString    RASTR_output_filename                    = my_current_job.arguments[34].ReturnStringArgument( );
    wxString    RASTR_output_star_filename               = my_current_job.arguments[35].ReturnStringArgument( );   
    bool        SPOT_RASTR                               = my_current_job.arguments[36].ReturnBoolArgument( );
    wxString    SPOT_RASTR_output_filename               = my_current_job.arguments[37].ReturnStringArgument( );
    wxString    SPOT_RASTR_output_star_filename          = my_current_job.arguments[38].ReturnStringArgument( );   
    bool        set_expert_options                       = my_current_job.arguments[39].ReturnBoolArgument( );
    float       psi_min                                  = my_current_job.arguments[40].ReturnFloatArgument( );
    float       psi_max                                  = my_current_job.arguments[41].ReturnFloatArgument( );
    float       psi_step                                 = my_current_job.arguments[42].ReturnFloatArgument( );
    float       fine_tuning_psi_step                     = my_current_job.arguments[43].ReturnFloatArgument( );
    float       padding_factor                           = my_current_job.arguments[44].ReturnFloatArgument( );
    int         max_threads                              = my_current_job.arguments[45].ReturnIntegerArgument( );

    
    // initiate I/O variables
    MRCFile          my_input_file(input_filename.ToStdString( ), false);  // check all the functions and things done with the MRCFile and also check the wxPrintF statement 
    MRCFile          my_output_sum_image_filename(output_average_per_bin_filename.ToStdString( ), true);  
;
    MRCFile*         my_output_SPOT_RASTR_filename;

    long             number_of_input_images = my_input_file.ReturnNumberOfSlices( );



    // initiate default parameters for the ApplyCTFAndReturnCTFSumOfSquares function
    // (May be change that later to be expert options inputs???)
    bool absolute = false; 
    bool apply_beam_tilt = false; 
    bool apply_envelope = false; 
 
    // CTF object
    CTF              current_ctf;

    // check if input defocus file is given for the CTF calculations
    // before running any code check if the file containing defocus is given to the program 
    ctf_parameters* ctf_parameters_stack = new ctf_parameters[number_of_input_images];

    if ( input_ctf_values_from_star_file == true ) {
        BasicStarFileReader input_star_file;
        //wxString            star_error_text;
        if ( (is_running_locally && ! DoesFileExist(input_star_filename.ToStdString( ))) ) {
            SendErrorAndCrash(wxString::Format("Error: Input star file %s not found\n", input_star_filename));
        }
        input_star_file.ReadFile(input_star_filename.ToStdString( ));
        
        for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
            ctf_parameters_stack[image_counter].acceleration_voltage          = acceleration_voltage;
            ctf_parameters_stack[image_counter].spherical_aberration          = spherical_aberration ;
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

    current_ctf.Init(acceleration_voltage, spherical_aberration, amplitude_contrast, defocus_1, defocus_2, astigmatism_angle, 0.0, 0.5, 0.0, pixel_size, additional_phase_shift);

    // Initialize sum of squares vector of vectors
    // Create a dynamic memory for vector of vectors to save the CTF sum of squares for all images sumed in each bin
    Image     current_image;

    // calculating the ctf_sum_of_squares
    current_image.ReadSlice(&my_input_file, 1);
    std::vector<std::vector<float>>* CTFSumOfSquares = new std::vector<std::vector<float>>();
    CTFSumOfSquares->resize(bins_count);

    std::vector<std::vector<float>>* CTFSumOfSquaresNew = new std::vector<std::vector<float>>();
    CTFSumOfSquaresNew->resize(bins_count);  

    InitializeCTFSumOfSquares(bins_count, current_image, CTFSumOfSquares);
    InitializeCTFSumOfSquares(bins_count, current_image, CTFSumOfSquaresNew);



    // Initialize the sum images based on the number specified by the user
    Image sum_images[bins_count];

    for (int bin_index = 0; bin_index < bins_count; bin_index++) {
        sum_images[bin_index].Allocate(my_input_file.ReturnXSize( ), my_input_file.ReturnYSize( ), true); //allocate in real space
        sum_images[bin_index].SetToConstant(0.0);
    }


    std::vector<std::vector<float>> all_columns_sum(number_of_input_images, std::vector<float>(my_input_file.ReturnXSize( ), 0.0));

    // calculate the bin range
    float bin_range = (max_tube_diameter - min_tube_diameter) / bins_count;

    float tube_rotation[number_of_input_images] = {0.0};
    float best_sum_column[number_of_input_images] = {0.0};
    float x_shift_column[number_of_input_images] = {0.0};
    float y_shift_row[number_of_input_images] = {0.0};
    float all_diameters[number_of_input_images] = {0.0};
    //std::vector<float> all_diameters(number_of_input_images, 0.0f);
    int   diameter_bins[number_of_input_images];
    float center_peak_index = current_image.logical_y_dimension / 2;


    //Initialize mask file
    MRCFile*         my_mask_file;
    Image*           mask;
    if ( RASTR == true & input_mask == true ) {
            if ( ! DoesFileExist(input_mask_filename.ToStdString( )) ) {
                SendError(wxString::Format("Error: Mask %s not found\n", input_mask_filename.ToStdString( )));
                exit(-1);
            }

            my_mask_file = new MRCFile(input_mask_filename.ToStdString( ), false);
            
            //Curve hist;
            mask = new Image;
            mask->Allocate(my_mask_file->ReturnXSize( ), my_mask_file->ReturnYSize( ), my_mask_file->ReturnZSize( ));
            mask->ReadSlices(my_mask_file, 1, my_mask_file->ReturnNumberOfSlices( ) );
            //mask.ComputeHistogramOfRealValuesCurve(&hist);
            //hist.PrintToStandardOut( );
            if ( mask->ReturnMaximumValue() != 1 & mask->ReturnMinimumValue() != 0){
                SendError(wxString::Format("Error: The minimum value in the mask is %f not 0 and the maximum value in the mask is %f not 1",mask->ReturnMinimumValue(), mask->ReturnMaximumValue( )));
                exit(-1);
            }
            
    }

    // if input mask file is not provided make initialize the x,y,z centers of the mask to the default values
    if (input_mask != true) {
        if (x_mask_center == 1) { 
            x_mask_center =  0.75 * current_image.logical_x_dimension;
        }
        if (y_mask_center == 1) {
            y_mask_center =  0.5  * current_image.logical_x_dimension;
        }
        if (z_mask_center == 1) {
            z_mask_center = 0.5  * current_image.logical_x_dimension;
        }
        if (sphere_mask_radius == 1) {
            sphere_mask_radius = 0.3 * current_image.logical_x_dimension; //previously was 0.1875 
        }
    }
    
    // initiate other needed variables
    Image     temp_image;
    long      image_counter;
    Image     final_image;
    Image     final_image_copy;

    wxPrintf("\nCalculating tube rotation using auto-correlation...\n\n");
    ProgressBar* my_progress = new ProgressBar(number_of_input_images); 
    //my_output_sum_images,
    //sum_image_below_min_diameter, sum_image_above_max_diameter,\
//image_counter , output_file
#pragma omp parallel for ordered schedule(dynamic) num_threads(max_threads) default(none) shared(my_input_file, best_sum_column, tube_rotation, all_columns_sum, number_of_input_images, max_threads, \
                                                                      x_shift_column, all_diameters, psi_step, min_tube_diameter, max_tube_diameter, pixel_size, psi_min, psi_max, tubes_centered,  \
                                                                        defocus_1, defocus_2, astigmatism_angle, additional_phase_shift, current_ctf, my_progress, \
                                                                        CTFSumOfSquares, bins_count , bin_range, sum_images, diameter_bins, outer_mask_radius, center_peak_index,\
                                                                        absolute, apply_beam_tilt, apply_envelope, input_ctf_values_from_star_file, ctf_parameters_stack, I ) \
                                                                    private(current_image, temp_image, final_image, final_image_copy, image_counter) 



    for ( image_counter = 0 ; image_counter < number_of_input_images; image_counter++ ) { 
        //wxPrintf("The image counter for slice reading is %li \n", image_counter+1);

        #pragma omp critical
        // read the current image in the stack
        current_image.ReadSlice(&my_input_file, image_counter + 1);
        // Normalize the image using cisTEM Normalize
        current_image.Normalize( );        
        
        if ((tubes_centered == true) && (outer_mask_radius == 0 || outer_mask_radius != 0)) {
            current_image.CircleMask(my_input_file.ReturnXSize( )/4);
        }else if (tubes_centered == false && outer_mask_radius != 0) {
            current_image.CircleMask(outer_mask_radius); 
        } else if (tubes_centered == false && outer_mask_radius == 0) {
            current_image.CircleMask(my_input_file.ReturnXSize( ) * 0.9);
        }             
        // FT the image
        current_image.ForwardFFT( );
        // convert the central pixel to zero (Is that done in real or Fouriier space??)
        current_image.ZeroCentralPixel( );

        // if defocus values and CTF parameters are from input file then use them
        if ( input_ctf_values_from_star_file ) {
                current_ctf.Init(ctf_parameters_stack[image_counter].acceleration_voltage, ctf_parameters_stack[image_counter].spherical_aberration, ctf_parameters_stack[image_counter].amplitude_contrast, ctf_parameters_stack[image_counter].defocus_1, ctf_parameters_stack[image_counter].defocus_2, ctf_parameters_stack[image_counter].astigmatism_angle, ctf_parameters_stack[image_counter].lowest_frequency_for_fitting, ctf_parameters_stack[image_counter].highest_frequency_for_fitting, ctf_parameters_stack[image_counter].astigmatism_tolerance, ctf_parameters_stack[image_counter].pixel_size, ctf_parameters_stack[image_counter].additional_phase_shift);
        }

        // apply CTF to the current image ????
        current_image.ApplyCTF(current_ctf); //Image::ApplyCTF(CTF ctf_to_apply, bool absolute, bool apply_beam_tilt, bool apply_envelope)
       
        for ( long pixel_counter = 0; pixel_counter < current_image.real_memory_allocated / 2; pixel_counter++ ) {

                // calculating the amplitude is not needed by let's see and print its value
                float amplitude = abs(current_image.complex_values[pixel_counter]);
                //wxPrintf("The image number %li pixel counter is %li \n", image_counter+1, pixel_counter);
                //wxPrintf("Retunring the amplitude from the complex values %f \n", amplitude);
                //wxPrintf("Retunring the complex values before changing them are %f, %f \n", current_image.complex_values[pixel_counter].real(), current_image.complex_values[pixel_counter].imag());

                //As the phase will be zero, the real is just the amplitude and the imaginary is 0
                // so we will set the complex number to be equal to apmlitude + 0 
                current_image.complex_values[pixel_counter] = amplitude*amplitude + I*0.0f ; 

        }
        // return the image to real space again and save them to see the correlation
        current_image.BackwardFFT( );
        current_image.SwapRealSpaceQuadrants( );
        // set the image is centered inside the box as true 
        current_image.object_is_centred_in_box = true;

        // Now let's sum the images horizontally (across the columns and see how they will look) and rotate them then sum again after each rotation
        // finding the initial rotation angle when searching within 180 degrees and a step size 4 degrees

        for (float psi = psi_min; psi <= psi_max; psi+=psi_step) {

            temp_image.CopyFrom(&current_image);
            // rotate by the temporary image by the rotation angle
            temp_image.Rotate2DInPlace(psi, FLT_MAX);
            float column_sum = sum_image_columns_float(&temp_image);
            std::vector<float> column_sum_vector = sum_image_columns(&temp_image);

            // column wise calculations    
            if (column_sum > best_sum_column[image_counter]){
                tube_rotation[image_counter] = psi;
                best_sum_column[image_counter] = column_sum;
                all_columns_sum[image_counter] = column_sum_vector;
            }
        }

        // tuning the rotation angle to search with 2 degrees before and after the initial angle
        // The step size will be rotation angle (2) /20 

        float tuning_rotation_range = psi_step / 2;
        float tuning_step_size = psi_step / 20;

        // column wise calculations
        float current_best_psi_column = tube_rotation[image_counter];
        float tuning_psi_lower_range_column = current_best_psi_column - tuning_rotation_range;
        float tuning_psi_upper_range_column = current_best_psi_column + tuning_rotation_range;



        Image tuning_column_temp_image;


        for (float tuning_psi_column = tuning_psi_lower_range_column; tuning_psi_column <= tuning_psi_upper_range_column; tuning_psi_column+=tuning_step_size) {

            tuning_column_temp_image.CopyFrom(&current_image);
            // rotate by the temporary image by the rotation angle
            tuning_column_temp_image.Rotate2DInPlace(tuning_psi_column, FLT_MAX);

            float tuning_column_sum = sum_image_columns_float(&tuning_column_temp_image);

            std::vector<float> tuning_column_sum_vector = sum_image_columns(&tuning_column_temp_image);

            if (tuning_column_sum > best_sum_column[image_counter]){
                tube_rotation[image_counter] = tuning_psi_column;
                best_sum_column[image_counter] = tuning_column_sum;
                all_columns_sum[image_counter] = tuning_column_sum_vector;
            }

        }

        //wxPrintf("The best psi angle with the highest sum is %f with sum %f \n\n", tube_rotation[image_counter], best_sum[image_counter]);
        //wxPrintf("The best psi angle with the highest sum from the column sum is %f with sum %f \n\n", tube_rotation[image_counter], best_sum_column[image_counter]);

        // find  the x,y shift
        // ReadSlice requires omp critical to avoid parallel reads, which may lead to the wrong slice being read
        #pragma omp critical

        final_image.ReadSlice(&my_input_file, image_counter + 1);
        // if mask is needed to find the tubes then apply it here 
        if ((tubes_centered == true) && (outer_mask_radius == 0 || outer_mask_radius != 0)) {
            final_image.CircleMask(my_input_file.ReturnXSize( )/4);
        }else if (tubes_centered == false && outer_mask_radius != 0) {
            final_image.CircleMask(outer_mask_radius); 
        } else if (tubes_centered == false && outer_mask_radius == 0) {
            final_image.CircleMask(my_input_file.ReturnXSize( ) * 0.9);
        }     

        //final_final_image.Normalize( );
        final_image.Rotate2DInPlace(tube_rotation[image_counter], FLT_MAX); // if not 0.0 it will not crop the images into circle after rotation as no mask will be applied 

        //final_image.QuickAndDirtyWriteSlice("rotated_images_based_on_psi_from_column_sum_b4_shift_calculations_before_gaussian_filter.mrc", image_counter + 1);

        final_image.ForwardFFT( ); 
        // to apply Gaussian filter you need to be in Fourier space
        // Nyquist frequency value is the pixel size * 2
        // if we want to apply a gaussian pass filter that will make the image at 150 angestrom to be well smoothened and get better peaks
        final_image.GaussianLowPassFilter((pixel_size*2)/150); // use a sigma between 0-1 for best results as this will remove the high frequency information
        
        final_image.BackwardFFT( );
        // calculate the required x shift to center the tubes
        //auto [peak_one_column_sum, peak_two_column_sum]  = find_column_sum_inner_diameter_peaks(&final_image, min_tube_diameter);
        auto [peak_one_column_sum, peak_two_column_sum]  = find_column_sum_peaks(&final_image , min_tube_diameter);
        // The next line not needed
        float tube_center_column_sum = std::abs(peak_one_column_sum - peak_two_column_sum)/2;
        float distance_from_center_column_sum = -((peak_one_column_sum + peak_two_column_sum)/2 - center_peak_index); 

        x_shift_column[image_counter] = distance_from_center_column_sum; // the x-shift needed
        // find the diameter and save it
        float tube_diameter = std::abs((peak_one_column_sum - peak_two_column_sum)); // * pixel_size
        all_diameters[image_counter] = tube_diameter;

        // save an index for the class assignment of each image based on its diameter
        int which_bin_index_int = (tube_diameter - min_tube_diameter) / bin_range; // it will always round down so 0.9999 > 0
        // save the bin assignment so that I don't need to recalculate the diameter again outside that loop
        diameter_bins[image_counter] = which_bin_index_int; 

        if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
            my_progress->Update(image_counter + 1);
    
    }

    delete my_progress;


    
    Image added_image; // This is the sum image based on the initial rotation angle calculated from auto-correlation
    wxPrintf("\nCreating Initial Sum Images...\n\n");
    ProgressBar* sum_progress = new ProgressBar(number_of_input_images); 

    for ( image_counter = 0 ; image_counter < number_of_input_images; image_counter++ ) { 

        added_image.ReadSlice(&my_input_file, image_counter + 1);
        added_image.Normalize( );
        //added_image.QuickAndDirtyWriteSlice("added_image_after_normalization.mrc", image_counter +1);
        if ( input_ctf_values_from_star_file ) {
                current_ctf.Init(ctf_parameters_stack[image_counter].acceleration_voltage, ctf_parameters_stack[image_counter].spherical_aberration, ctf_parameters_stack[image_counter].amplitude_contrast, ctf_parameters_stack[image_counter].defocus_1, ctf_parameters_stack[image_counter].defocus_2, ctf_parameters_stack[image_counter].astigmatism_angle, ctf_parameters_stack[image_counter].lowest_frequency_for_fitting, ctf_parameters_stack[image_counter].highest_frequency_for_fitting, ctf_parameters_stack[image_counter].astigmatism_tolerance, ctf_parameters_stack[image_counter].pixel_size, ctf_parameters_stack[image_counter].additional_phase_shift);
        }


        added_image.ForwardFFT( );
        added_image.ZeroCentralPixel( );

        for (int bin_index = 0; bin_index < bins_count; bin_index++){
            if (diameter_bins[image_counter] == bin_index) { // This should catch any diameter within the range
                ApplyCTFAndReturnCTFSumOfSquares(added_image, current_ctf,  absolute, apply_beam_tilt, apply_envelope, (*CTFSumOfSquares)[bin_index]);
            // } else if (diameter_bins[image_counter] < 0 ) { // if the which_bin_index is negative this mean the tube diameter is less than the min. so do nothing
            //     //ApplyCTFAndReturnCTFSumOfSquares(added_image, current_ctf,  absolute, apply_beam_tilt, apply_envelope, (*CTFSumOfSquares)[0]); // The indexing in the std::vector<std::vector<float>> start with 0
            // } else if (diameter_bins[image_counter] >= bins_count ) { // if the which_bin_index is larger than or equal the bins_count then the diameter is more than the max. so will do nothing
            //     //ApplyCTFAndReturnCTFSumOfSquares(added_image, current_ctf,  absolute, apply_beam_tilt, apply_envelope, (*CTFSumOfSquares)[bins_count - 1]); //Indexing start with 0   
            }

        }

        added_image.BackwardFFT( );
        //added_image.QuickAndDirtyWriteSlice("added_image_after_applyactf.mrc", image_counter +1);

        added_image.Rotate2DInPlace(tube_rotation[image_counter], FLT_MAX);
        //added_image.QuickAndDirtyWriteSlice("added_image_after_rotation.mrc", image_counter +1);

        added_image.PhaseShift(x_shift_column[image_counter], 0.0, 0.0);
        //added_image.QuickAndDirtyWriteSlice("added_image_after_phaseshift.mrc", image_counter +1);

        // using the dynamic memory allocation to add the current image to the sum images based on the tube diiameter
        for (int bin_index = 0; bin_index < bins_count; bin_index++){          
            if (diameter_bins[image_counter] == bin_index) { // This should catch any diameter within the range
                sum_images[bin_index].AddImage(&added_image);
            // } else if (diameter_bins[image_counter] < 0 ) { // if the which_bin_index is negative this mean the tube diameter is less than the min. so will be considered with the first bin
            //     sum_image_below_min_diameter.AddImage(&added_image);
            // } else if (diameter_bins[image_counter] >= bins_count ) { // if the which_bin_index is larger than or equal the bins_count then the diameter is more than the max. so will be added to the last bin
            //     sum_image_above_max_diameter.AddImage(&added_image);
            }

        }
    
    if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
        sum_progress->Update(image_counter + 1);
    }

    delete sum_progress;

    
    // divide the sum image by CTF sum of squares and centering it 
    for (int bin_index = 0; bin_index < bins_count; bin_index++){          
        sum_images[bin_index].ForwardFFT( );
        divide_by_ctf_sum_of_squares(sum_images[bin_index], (*CTFSumOfSquares)[bin_index]);
        sum_images[bin_index].BackwardFFT( );
        //sum_images[bin_index].QuickAndDirtyWriteSlice("sum_image_after_divide_ctf_sum_of_squares_before_centering.mrc", bin_index + 1 );
        // shift sum image to the center after padding based on tube peaks?
        auto [peak_one_column_sum, peak_two_column_sum]  = find_column_sum_peaks(&sum_images[bin_index],min_tube_diameter);

        float tube_center_column_sum = std::abs(peak_one_column_sum - peak_two_column_sum)/2;
        float distance_from_center_column_sum = -((peak_one_column_sum + peak_two_column_sum)/2 - center_peak_index); 
        //to center the sum image
        sum_images[bin_index].PhaseShift(distance_from_center_column_sum, 0.0, 0.0); // the x-shift needed to center the sum image
 
        // sum_images[bin_index].QuickAndDirtyWriteSlice("sum_image_after_divide_ctf_sum_of_squares.mrc", bin_index + 1 );

    }
    
    delete CTFSumOfSquares;

    /// Getting the correct rotation and shift using cross-correlation
    float   inner_radius_for_peak_search;
    float   outer_radius_for_peak_search;
    inner_radius_for_peak_search = 0.0; // inner radius should be set to 0
    outer_radius_for_peak_search = my_input_file.ReturnXSize( )/2;



    // find the best rotation, shift for the alignment
    float   best_correlation_score[number_of_input_images] = {-FLT_MAX} ; //= {-FLT_MAX}   
    float   best_psi_value[number_of_input_images] = {0}; // saves the angle within 180 degrees, conversion is done later in the code
    float   best_x_shift_value[number_of_input_images] = {0};
    float   best_y_shift_value[number_of_input_images] = {0};
    
    Image my_image;
    Image average_image;
    Image tuning_average_image;
    float tuned_rotation_range = psi_step; ///2
    float tuned_step_size = fine_tuning_psi_step;
    Image fine_tuning_average_image;

    
    wxPrintf("\nAligning Images...\n\n");
    ProgressBar* my_aln_progress = new ProgressBar(number_of_input_images); 

#pragma omp parallel for ordered schedule(dynamic) num_threads(max_threads) default(none) shared(number_of_input_images, my_input_file , inner_radius_for_peak_search, outer_radius_for_peak_search, tubes_centered, \
                                                                   best_correlation_score, best_psi_value, best_x_shift_value, best_y_shift_value, psi_step, tube_rotation, outer_mask_radius, \
                                                                    max_threads, diameter_bins, sum_images, bins_count,tuned_rotation_range, tuned_step_size,  my_aln_progress, x_shift_column, y_shift_row, all_diameters, bin_range, \
                                                                    input_ctf_values_from_star_file, current_ctf, ctf_parameters_stack, min_tube_diameter, center_peak_index, pixel_size, low_pass) \
                                                                   private(my_image, average_image, tuning_average_image, final_image, fine_tuning_average_image)

    for ( long aln_image_counter = 0; aln_image_counter < number_of_input_images; aln_image_counter++ ) {
        #pragma omp critical
        my_image.ReadSlice(&my_input_file, aln_image_counter + 1); 
        my_image.Normalize( );
        my_image.ForwardFFT( );
        //testing adding a low pass filter on the original image before getting the correct shift from the correlation and how that can affect the centering of the mask at the end
        if (low_pass) {
            my_image.GaussianLowPassFilter((pixel_size*2)/150); 
        }
        my_image.BackwardFFT( );
        if ((tubes_centered == true) && (outer_mask_radius == 0 || outer_mask_radius != 0)) {
            my_image.CircleMask(my_input_file.ReturnXSize( )/4);
        }else if (tubes_centered == false && outer_mask_radius != 0) {
            my_image.CircleMask(outer_mask_radius); 
        } else if (tubes_centered == false && outer_mask_radius == 0) {
            my_image.CircleMask(my_input_file.ReturnXSize( ) * 0.9);
        }

        // initial angle search will start from the rotation angle we got from the auto-correlation
        // then change the auto-correlation angle to be within 180 
        // do another search within +/- 90 degrees of the auto-correlation psi angle
        float current_best_psi = tube_rotation[aln_image_counter];  
        // normalize the angle to be within 180 degrees
        //current_best_psi = fmod((current_best_psi + 360.0), 180.0);  
        float angle_range = 180.0;
        float psi_min_angle = current_best_psi - 0.5 * angle_range;
        float psi_max_angle = current_best_psi + 0.5 * angle_range;

        for (float psi = psi_min_angle; psi < psi_max_angle; psi+=psi_step) {
            // create a new peak to save the cross-correlation peak values
            Peak    current_peak;

            for (int bin_index = 0; bin_index < bins_count; bin_index++){          
                if (diameter_bins[aln_image_counter] == bin_index) { // This should catch any diameter within the range
                    average_image.CopyFrom(&sum_images[bin_index]);
                    average_image.Normalize( );
                } else if (diameter_bins[aln_image_counter] < 0 ) { // if the which_bin_index is negative this mean the tube diameter is less than the min. so will be considered with the first bin
                    average_image.CopyFrom(&sum_images[0]);
                    average_image.Normalize( );
                } else if (diameter_bins[aln_image_counter] >= bins_count ) { // if the which_bin_index is larger than or equal the bins_count then the diameter is more than the max. so will be added to the last bin
                    average_image.CopyFrom(&sum_images[bins_count - 1]);
                    average_image.Normalize( );
                }

            }

            if ( input_ctf_values_from_star_file ) {
                current_ctf.Init(ctf_parameters_stack[aln_image_counter].acceleration_voltage, ctf_parameters_stack[aln_image_counter].spherical_aberration, ctf_parameters_stack[aln_image_counter].amplitude_contrast, ctf_parameters_stack[aln_image_counter].defocus_1, ctf_parameters_stack[aln_image_counter].defocus_2, ctf_parameters_stack[aln_image_counter].astigmatism_angle, ctf_parameters_stack[aln_image_counter].lowest_frequency_for_fitting, ctf_parameters_stack[aln_image_counter].highest_frequency_for_fitting, ctf_parameters_stack[aln_image_counter].astigmatism_tolerance, ctf_parameters_stack[aln_image_counter].pixel_size, ctf_parameters_stack[aln_image_counter].additional_phase_shift);
            }
            average_image.ForwardFFT( );    
            average_image.ApplyCTF(current_ctf);
            average_image.BackwardFFT( );

            // rotate by the sum image by the rotation angle
            average_image.Rotate2DInPlace(psi,0.0);

            // calculate the cross correlation of the reference image with the rotated image
            average_image.CalculateCrossCorrelationImageWith(&my_image); 
            
            // find the peak from the cross corrlation to get the values 
            current_peak = average_image.FindPeakWithParabolaFit(inner_radius_for_peak_search, outer_radius_for_peak_search);

            // check if the current cross-correlation score is better than stored or not
            if ( current_peak.value > best_correlation_score[aln_image_counter] ) {
                best_correlation_score[aln_image_counter]  = current_peak.value;
                best_psi_value[aln_image_counter] = psi; 
                best_x_shift_value[aln_image_counter] = current_peak.x;
                best_y_shift_value[aln_image_counter] = current_peak.y;
                
            } 


        } // end of initial search for the best psi and shift angles
        // Start the tuning loop
        // save the best_psi as the current_best_psi for that image
        // calculate the tuning psi range which is = to the rotation angle and step size = rotation angle / 20
        current_best_psi = best_psi_value[aln_image_counter]; 
        float tuned_psi_lower_range = current_best_psi - tuned_rotation_range;
        float tuned_psi_upper_range = current_best_psi + tuned_rotation_range;

        // loop over the range of +/- half the rotation angle 
        // increment by 1/10 of the tuned rotation angle (rotation angle/2)/10 degrees for tuning 
        for (float tuned_psi = tuned_psi_lower_range; tuned_psi <= tuned_psi_upper_range; tuned_psi+=tuned_step_size) {

            // create a new peak to save the tuned values
            Peak    current_tuned_peak;

            // rotate by the class average image by the rotation angle 

            for (int bin_index = 0; bin_index < bins_count; bin_index++){          
                if (diameter_bins[aln_image_counter] == bin_index) { // This should catch any diameter within the range
                    tuning_average_image.CopyFrom(&sum_images[bin_index]);
                    tuning_average_image.Normalize( );
                } else if (diameter_bins[aln_image_counter] < 0 ) { // if the which_bin_index is negative this mean the tube diameter is less than the min. so will be considered with the first bin
                    tuning_average_image.CopyFrom(&sum_images[0]);
                    tuning_average_image.Normalize( );
                } else if (diameter_bins[aln_image_counter] >= bins_count ) { // if the which_bin_index is larger than or equal the bins_count then the diameter is more than the max. so will be added to the last bin
                    tuning_average_image.CopyFrom(&sum_images[bins_count - 1]);
                    tuning_average_image.Normalize( );
                }

            }

            if ( input_ctf_values_from_star_file ) {
                current_ctf.Init(ctf_parameters_stack[aln_image_counter].acceleration_voltage, ctf_parameters_stack[aln_image_counter].spherical_aberration, ctf_parameters_stack[aln_image_counter].amplitude_contrast, ctf_parameters_stack[aln_image_counter].defocus_1, ctf_parameters_stack[aln_image_counter].defocus_2, ctf_parameters_stack[aln_image_counter].astigmatism_angle, ctf_parameters_stack[aln_image_counter].lowest_frequency_for_fitting, ctf_parameters_stack[aln_image_counter].highest_frequency_for_fitting, ctf_parameters_stack[aln_image_counter].astigmatism_tolerance, ctf_parameters_stack[aln_image_counter].pixel_size, ctf_parameters_stack[aln_image_counter].additional_phase_shift);
            }
            tuning_average_image.ForwardFFT( );    
            tuning_average_image.ApplyCTF(current_ctf);
            tuning_average_image.BackwardFFT( );


            tuning_average_image.Rotate2DInPlace(tuned_psi, 0.0);

            // calculate the cross correlation of the reference image with the rotated image
            tuning_average_image.CalculateCrossCorrelationImageWith(&my_image); //rotated_image
            if ((tubes_centered == true) && (outer_mask_radius == 0 || outer_mask_radius != 0)) {
                tuning_average_image.CircleMask(my_input_file.ReturnXSize( )/4);
            }else if (tubes_centered == false && outer_mask_radius != 0) {
                tuning_average_image.CircleMask(outer_mask_radius); 
            } else if (tubes_centered == false && outer_mask_radius == 0) {
                tuning_average_image.CircleMask(my_input_file.ReturnXSize( ) * 0.9);
            }  
            // find the peak from the cross corrlation to get the values from the tuning_average_image
            current_tuned_peak = tuning_average_image.FindPeakWithParabolaFit(inner_radius_for_peak_search, outer_radius_for_peak_search);
            // debug scores
            //wxPrintf("returning peak values anyway, psi = %f, best score = %lf, current score = %f\n", tuned_psi, best_correlation_score[aln_image_counter], current_tuned_peak.value);

            // check if the current cross-correlation score using the tuned angle is better than stored or not
            if ( current_tuned_peak.value > best_correlation_score[aln_image_counter] ) {
                best_correlation_score[aln_image_counter]  =   current_tuned_peak.value;
                best_psi_value[aln_image_counter] =  tuned_psi; //because we are rotating the sum image and will later want to rotate the sum_image 
                best_x_shift_value[aln_image_counter] = current_tuned_peak.x;
                best_y_shift_value[aln_image_counter] = current_tuned_peak.y;
            } 

        }

        // Finding the tube shift and diameter to correctly shift tubes before adding them for the second time and to correctly classify the tubes based on final psi
        #pragma omp critical
        final_image.ReadSlice(&my_input_file, aln_image_counter + 1);
        if ((tubes_centered == true) && (outer_mask_radius == 0 || outer_mask_radius != 0)) {
            final_image.CircleMask(my_input_file.ReturnXSize( )/4);
        }else if (tubes_centered == false && outer_mask_radius != 0) {
            final_image.CircleMask(outer_mask_radius); 
        } else if (tubes_centered == false && outer_mask_radius == 0) {
            final_image.CircleMask(my_input_file.ReturnXSize( ) * 0.9);
        }  
        //final_final_image.Normalize( );
        final_image.Rotate2DInPlace(-best_psi_value[aln_image_counter], FLT_MAX); // if not 0.0 it will not crop the images into circle after rotation as no mask will be applied 
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // I commented the next line as I need to get the correct shift in x an d y dimensions with respect to the center/////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //final_image.PhaseShift(best_x_shift_value[aln_image_counter], best_y_shift_value[aln_image_counter]);
        
        final_image.ForwardFFT( ); 
        // to apply Gaussian filter you need to be in Fourier space
        // Nyquist value is the pixel size * 2
        // if we want to apply a gaussian pass filter that will make the image at 150 angestrom to be well smoothened and get better peaks
        final_image.GaussianLowPassFilter((pixel_size*2)/150); // use a sigma between 0-1 for best results as this will remove the high frequency information

        final_image.BackwardFFT( );

        
        // calculate the required x shift to center the tubes
        //auto [peak_one_column_sum, peak_two_column_sum]  = find_column_sum_inner_diameter_peaks(&final_image, min_tube_diameter);
        auto [peak_one_column_sum, peak_two_column_sum]  = find_column_sum_peaks(&final_image,min_tube_diameter);
        // update the x_shift required to center the column (will be needed later????)
        float tube_center_column_sum = std::abs(peak_one_column_sum - peak_two_column_sum)/2;
        float distance_from_center_column_sum = -((peak_one_column_sum + peak_two_column_sum)/2 - center_peak_index); 

        x_shift_column[aln_image_counter] = distance_from_center_column_sum; // the x-shift needed

        // find the diameter and save it
        float tube_diameter = std::abs((peak_one_column_sum - peak_two_column_sum)); // * pixel_size
        all_diameters[aln_image_counter] = tube_diameter;

        final_image.Rotate2DInPlace(90.0, FLT_MAX);
        // calculate the required y shift to center the tubes
        // we will need to rotate the tube 90 degrees to be aligned with y-axis
        //auto [peak_one_column_sum, peak_two_column_sum]  = find_column_sum_inner_diameter_peaks(&final_image, min_tube_diameter);
        auto [peak_one_row_sum, peak_two_row_sum]  = find_row_sum_peaks(&final_image,min_tube_diameter);
        // update the x_shift required to center the column (will be needed later????)
        float tube_center_row_sum = std::abs(peak_one_row_sum - peak_two_row_sum)/2;
        float distance_from_center_row_sum = -((peak_one_row_sum + peak_two_row_sum)/2 - center_peak_index); 

        y_shift_row[aln_image_counter] = distance_from_center_row_sum; // the x-shift needed

        //wxPrintf("The image number %li x shift is %f, y-shift is %f and best x shift is %f, best y shift is %f \n\n", aln_image_counter, x_shift_column[aln_image_counter], y_shift_row[aln_image_counter], -best_x_shift_value[aln_image_counter], -best_y_shift_value[aln_image_counter]);

        // save an index for the class assignment of each image based on its diameter
        int which_bin_index_int = (tube_diameter - min_tube_diameter) / bin_range; // it will always round down so 0.9999 > 0
        // save the bin assignment so that I don't need to recalculate the diameter again outside that loop
        diameter_bins[aln_image_counter] = which_bin_index_int; 
        
        if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
            my_aln_progress->Update(aln_image_counter + 1);

    }
    delete my_aln_progress;


    Image added_image_after_aln;
    Image sum_images_after_aln[bins_count];

    wxPrintf("\nCreating Final Sum Images...\n\n");
    ProgressBar* update_sum_progress = new ProgressBar(number_of_input_images); 


    for (int bin_index = 0; bin_index < bins_count; bin_index++) {
        sum_images_after_aln[bin_index].Allocate(my_input_file.ReturnXSize( ), my_input_file.ReturnYSize( ), true); //allocate in real space
        sum_images_after_aln[bin_index].SetToConstant(0.0);
    }

    // creating a sum image after getting the best rotation from the alignment
     for ( image_counter = 0 ; image_counter < number_of_input_images; image_counter++ ) { 

            added_image_after_aln.ReadSlice(&my_input_file, image_counter + 1);
            added_image_after_aln.Normalize( );
            //added_image.QuickAndDirtyWriteSlice("added_image_after_normalization.mrc", image_counter +1);
            if ( input_ctf_values_from_star_file ) {
                    current_ctf.Init(ctf_parameters_stack[image_counter].acceleration_voltage, ctf_parameters_stack[image_counter].spherical_aberration, ctf_parameters_stack[image_counter].amplitude_contrast, ctf_parameters_stack[image_counter].defocus_1, ctf_parameters_stack[image_counter].defocus_2, ctf_parameters_stack[image_counter].astigmatism_angle, ctf_parameters_stack[image_counter].lowest_frequency_for_fitting, ctf_parameters_stack[image_counter].highest_frequency_for_fitting, ctf_parameters_stack[image_counter].astigmatism_tolerance, ctf_parameters_stack[image_counter].pixel_size, ctf_parameters_stack[image_counter].additional_phase_shift);
            }

            added_image_after_aln.ForwardFFT( );
            added_image_after_aln.ZeroCentralPixel( );

            for (int bin_index = 0; bin_index < bins_count; bin_index++){
                if (diameter_bins[image_counter] == bin_index) { // This should catch any diameter within the range
                    ApplyCTFAndReturnCTFSumOfSquares(added_image_after_aln, current_ctf,  absolute, apply_beam_tilt, apply_envelope, (*CTFSumOfSquaresNew)[bin_index]);
                // } else if (diameter_bins[image_counter] < 0 ) { // if the which_bin_index is negative this mean the tube diameter is less than the min. so will be considered with the first bin
                //     //ApplyCTFAndReturnCTFSumOfSquares(added_image_after_aln, current_ctf,  absolute, apply_beam_tilt, apply_envelope, (*CTFSumOfSquaresNew)[0]); // The indexing in the std::vector<std::vector<float>> start with 0
                // } else if (diameter_bins[image_counter] >= bins_count ) { // if the which_bin_index is larger than or equal the bins_count then the diameter is more than the max. so will be added to the last bin
                //     //ApplyCTFAndReturnCTFSumOfSquares(added_image_after_aln, current_ctf,  absolute, apply_beam_tilt, apply_envelope, (*CTFSumOfSquaresNew)[bins_count - 1]); //Indexing start with 0 
                    
                }

            }

            added_image_after_aln.BackwardFFT( );
            //added_image.QuickAndDirtyWriteSlice("added_image_after_applyactf.mrc", image_counter +1);

            added_image_after_aln.Rotate2DInPlace(-best_psi_value[image_counter], FLT_MAX);
            // update the x_shift column to get the shift that is needed to center the tube after rotation

            added_image_after_aln.PhaseShift(best_x_shift_value[image_counter], best_y_shift_value[image_counter], 0.0); // before was x_shift_column if you use it then uncomment that part in aln loop
            //added_image_after_aln.QuickAndDirtyWriteSlice("added_image_after_aln_after_rotation_and_shift.mrc", image_counter +1);

            // using the dynamic memory allocation to add the current image to the sum images based on the tube diiameter
            for (int bin_index = 0; bin_index < bins_count; bin_index++){          
                if (diameter_bins[image_counter] == bin_index) { // This should catch any diameter within the range
                    sum_images_after_aln[bin_index].AddImage(&added_image_after_aln);
                // } else if (diameter_bins[image_counter] < 0 ) { // if the which_bin_index is negative this mean the tube diameter is less than the min. so will be considered with the first bin
                //     sum_image_below_min_diameter.AddImage(&added_image_after_aln);
                // } else if (diameter_bins[image_counter] >= bins_count ) { // if the which_bin_index is larger than or equal the bins_count then the diameter is more than the max. so will be added to the last bin
                //     sum_image_above_max_diameter.AddImage(&added_image_after_aln);
                }

            }
            if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
                update_sum_progress->Update(image_counter+ 1);
    }

    delete update_sum_progress;



    for (int bin_index = 0; bin_index < bins_count; bin_index++){          
        sum_images_after_aln[bin_index].ForwardFFT( );
        divide_by_ctf_sum_of_squares(sum_images_after_aln[bin_index], (*CTFSumOfSquaresNew)[bin_index]);
        sum_images_after_aln[bin_index].BackwardFFT( );
        //sum_images_after_aln[bin_index].QuickAndDirtyWriteSlice("sum_image_new_after_divide_ctf_sum_of_squares_before_centering.mrc", bin_index + 1 );

        // shift sum image to the center after padding based on tube peaks?
        auto [peak_one_column_sum, peak_two_column_sum]  = find_column_sum_peaks(&sum_images_after_aln[bin_index],min_tube_diameter);
        float tube_center_column_sum = std::abs(peak_one_column_sum - peak_two_column_sum)/2;
        float distance_from_center_column_sum = -((peak_one_column_sum + peak_two_column_sum)/2 - center_peak_index); 
        sum_images_after_aln[bin_index].PhaseShift(distance_from_center_column_sum, 0.0, 0.0);
        //sum_images_after_aln[bin_index].QuickAndDirtyWriteSlice("sum_image_new_after_divide_ctf_sum_of_squares_after_centering.mrc", bin_index + 1 );
        sum_images_after_aln[bin_index].AddImage(&sum_images[bin_index]);
        // after we add the sum_images_after_aln to sum_images, we divide the sum_images by 2
        // sum_images_after_aln[bin_index].DivideByConstant(2.0);
        sum_images[bin_index].DivideByConstant(2.0);
    }
    delete CTFSumOfSquaresNew;
    sum_images_after_aln->Deallocate( );

    /////////////////////// Calculate the rotational average from the sum images
    /////////////////////// Then generate a 3D azimuthal average
    
    wxPrintf("\nCreating 3D azimuthal average models...\n\n");
    
    // Find the position of the dot if there is an extension in the filename provided by the user
    size_t extension_pos = output_azimuthal_average_volume_filename.find('.');
    if (extension_pos != std::string::npos) { // If dot is found
        // Extract substring before the dot 
        output_azimuthal_average_volume_filename = output_azimuthal_average_volume_filename.substr(0, extension_pos);
    }

    std::string model_file_name = "";
    // create a copy of the sum images to be used to create the azimuthal average model
    Image sum_images_copy[bins_count];
    long model_dimension = current_image.logical_x_dimension;
    ProgressBar* azimuthal_average_model_progress = new ProgressBar(bins_count * model_dimension); 

    for (int bin_index = 0; bin_index < bins_count; bin_index++) {
        // allocate memory for the sum images copy image
        sum_images_copy[bin_index].Allocate(my_input_file.ReturnXSize( ), my_input_file.ReturnYSize( ), true); //allocate in real space
        sum_images_copy[bin_index].SetToConstant(0.0);

        sum_images_copy[bin_index].CopyFrom(&sum_images[bin_index]);
        sum_image_direction(&sum_images_copy[bin_index], 2);
        //no need as moved before the resizing part 
        //sum_images[bin_index].Resize(current_image.logical_x_dimension, current_image.logical_y_dimension, 1);
        sum_images_copy[bin_index].ApplyRampFilter( );
        sum_images_copy[bin_index].AverageRotationally( );
        //sum_images[bin_index].QuickAndDirtyWriteSlice("sum_images_after_average_rotationally_from_aln.mrc", bin_index +1 );
        model_file_name = output_azimuthal_average_volume_filename + "_" + std::to_string(bin_index + 1) +".mrc";
        
        for (long model_counter = 0; model_counter < model_dimension; model_counter++){
            sum_images_copy[bin_index].QuickAndDirtyWriteSlice(model_file_name, model_counter + 1);
            azimuthal_average_model_progress->Update(bins_count* model_dimension + model_counter + 1);
        } 
    }
    sum_images_copy->Deallocate( );
    delete azimuthal_average_model_progress;


// //////////////////////////////// Delete later this is just testing something//////////////////////////
// // creating a sum image after getting the best rotation from the alignment
//     Image test_image; 
//     for ( image_counter = 0 ; image_counter < number_of_input_images; image_counter++ ) { 
//             //test_image.Allocate(my_input_file.ReturnXSize( ), my_input_file.ReturnYSize( ), true); //allocate in real space
//             //test_image.CopyFrom(&sum_images[0]);
//             test_image.ReadSlice(&my_input_file, image_counter + 1);
//             test_image.ForwardFFT( );
//             test_image.GaussianLowPassFilter((pixel_size*2)/150);
//             test_image.BackwardFFT( );
//             test_image.QuickAndDirtyWriteSlice("final_image_low_pass_filter.mrc", image_counter + 1);

//             //test_image.Rotate2DInPlace(best_psi_value[image_counter]);
//             //test_image.PhaseShift(x_shift_column[image_counter], y_shift_row[image_counter]);
//             //test_image.QuickAndDirtyWriteSlice("test_projection_image_shift.mrc", image_counter + 1);
//      }

    // rotationally averaged 3D reconstruction
    Image               model_volume[bins_count];
    ReconstructedVolume input_3d[bins_count];
    Image               projection_volume_3d;
    Image               projection_volume_image;
    Image               padded_projection_volume_image;
    AnglesAndShifts     my_parameters;

    // generate masked 3D AA models
    Image               my_masked_volume[bins_count];
    Image               my_mask;
    ReconstructedVolume masked_3d[bins_count];
    Image               masked_projection_volume_3d;
    Image               masked_projection_volume_image;
    Image               masked_padded_projection_volume_image;

    // mask needed for masking upweighted regions
    // create the mask
    Image                my_white_mask;
    ReconstructedVolume  mask_volume;
    Image*               mask_projection;

    if (RASTR == true) {
        // Allocate memory for the mask file to be read
        my_mask.Allocate(current_image.logical_x_dimension, current_image.logical_y_dimension, current_image.logical_x_dimension, true);
        
        if (input_mask == true) {
            #pragma omp critical
            my_mask.ReadSlices(my_mask_file, 1, my_mask_file->ReturnNumberOfSlices( ) );
            my_mask.BinariseInverse(0.0f); // any pixel of 0.0 or less will be 1.0

        } else {
            my_mask.SetToConstant(0.0);
            create_black_sphere_mask(&my_mask, x_mask_center, y_mask_center, z_mask_center, sphere_mask_radius);
            //my_mask[bin_index].QuickAndDirtyWriteSlices("my_created_black_mask.mrc", 1, current_image.logical_x_dimension);
        }
        // create a padded mask as the model is padded now
        my_mask.Resize(padding_factor * current_image.logical_x_dimension, padding_factor * current_image.logical_y_dimension, padding_factor * current_image.logical_x_dimension, 1.0);
        //my_mask[bin_index].QuickAndDirtyWriteSlices("my_created_black_mask_after_resizing.mrc", 1, padding_factor * current_image.logical_x_dimension);

        // create the mask for masking the upweighted regions
        // Allocate memory for the mask file to be read
        my_white_mask.Allocate(current_image.logical_x_dimension, current_image.logical_y_dimension, current_image.logical_x_dimension, true);
        if (input_mask == true ) {
            my_white_mask.ReadSlices(my_mask_file, 1, my_mask_file->ReturnNumberOfSlices( ));
            my_white_mask.Resize(padding_factor * my_mask_file->ReturnXSize(), padding_factor * my_mask_file->ReturnYSize(), padding_factor * my_mask_file->ReturnZSize());
            my_white_mask.BinariseInverse(0.0f); //any pixel of 0.0 or less will be 1.0
            my_white_mask.BinariseInverse(0.0f); // now we flip them again so the mask itself is 1.0 and background is 0.0

        } else {
            my_white_mask.SetToConstant(0.0);
            create_white_sphere_mask(&my_white_mask, x_mask_center, y_mask_center, z_mask_center, sphere_mask_radius);
            my_white_mask.Resize(padding_factor * current_image.logical_x_dimension, padding_factor * current_image.logical_y_dimension, padding_factor * current_image.logical_x_dimension);

        }

        //my_white_mask.Binarise(1.0);
        //my_white_mask.QuickAndDirtyWriteSlices("my_mask_final_after_invert_mask.mrc", 1, padding_factor* current_image.logical_x_dimension);

        mask_volume.InitWithDimensions(my_white_mask.logical_x_dimension, my_white_mask.logical_y_dimension, my_white_mask.logical_z_dimension, pixel_size);
        mask_volume.density_map->CopyFrom(&my_white_mask);
        //mask_volume.density_map->QuickAndDirtyWriteSlices("my_mask_volume_after_density_map.mrc", 1, current_image.logical_x_dimension);
        float mask_radius = FLT_MAX; //100 - FLT_MAX
        mask_volume.mask_radius = mask_radius;
        mask_volume.PrepareForProjections(0.0, 2.0 * pixel_size);
        mask_projection = new Image;
        mask_projection->CopyFrom(mask_volume.density_map);
        my_white_mask.Deallocate( );

    }


    wxPrintf("\nPreparing Azimuthal Average model for projection...\n\n");
    ProgressBar* prepare_projections_progress = new ProgressBar(bins_count);     
    //Image azimuthal_average_slice;
#pragma omp parallel num_threads(max_threads)  default(none) shared(SPOT_RASTR, RASTR, prepare_projections_progress, current_image, bins_count,sum_images, model_volume, my_masked_volume , my_mask, input_3d, masked_3d,  \
                                                                 pixel_size, padding_factor, x_mask_center, y_mask_center, z_mask_center, sphere_mask_radius, filter_radius, use_outside_value, \
                                                                 outside_value, outside_weight, cosine_edge, number_of_models, input_mask, my_mask_file, my_output_sum_image_filename) \
                                                                   private(projection_volume_3d, projection_volume_image, padded_projection_volume_image, masked_projection_volume_3d, masked_projection_volume_image, \
                                                                   masked_padded_projection_volume_image, my_white_mask, mask_volume, mask_projection, my_parameters )

#pragma omp for ordered schedule(dynamic, 1)
    for (int bin_index = 0; bin_index < bins_count; bin_index++) {
        
        model_volume[bin_index].Allocate(padding_factor * sum_images[bin_index].logical_x_dimension, padding_factor * sum_images[bin_index].logical_y_dimension, padding_factor * sum_images[bin_index].logical_x_dimension, true);
        model_volume[bin_index].SetToConstant(0.0);

        float edge_value = sum_images[bin_index].ReturnAverageOfRealValuesOnEdges( );
        sum_images[bin_index].Resize(model_volume[bin_index].logical_x_dimension, model_volume[bin_index].logical_y_dimension, 1, edge_value);

        // fill the padded version with the sum image
        // then do average rotationally before filling the volume
        sum_image_direction(&sum_images[bin_index], 2);
        sum_images[bin_index].ApplyRampFilter( );
        sum_images[bin_index].AverageRotationally( );

        // fill in the model volume with the azimuthal average slice
        long pixel_coord_xy  = 0;
        long pixel_coord_xyz = 0;
        long volume_counter  = 0;
        for ( int z = 0; z < model_volume[bin_index].logical_z_dimension; z++ ) {
            for ( int y = 0; y < model_volume[bin_index].logical_y_dimension; y++ ) {
                for ( int x = 0; x < model_volume[bin_index].logical_x_dimension; x++ ) {
                    pixel_coord_xy = sum_images[bin_index].ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
                    model_volume[bin_index].real_values[volume_counter] = sum_images[bin_index].real_values[pixel_coord_xy];
                    volume_counter++;
                }
                volume_counter += sum_images[bin_index].padding_jump_value;
            }
        }

        // if we will not apply mask and will subtract the projection of the azimuthal average as is

        if (SPOT_RASTR == true) {

            input_3d[bin_index].InitWithDimensions(model_volume[bin_index].logical_x_dimension, model_volume[bin_index].logical_y_dimension, model_volume[bin_index].logical_z_dimension, pixel_size);
            input_3d[bin_index].density_map->CopyFrom(&model_volume[bin_index]);
            float mask_radius = FLT_MAX; //100 - FLT_MAX
            input_3d[bin_index].mask_radius = mask_radius;
            //input_3d[bin_index].density_map->CorrectSinc(mask_radius / pixel_size);
            //change those lines to see the effect of resolution change on the created projections
            input_3d[bin_index].PrepareForProjections(0.0, 2.0 * pixel_size); // 0.0, 2.0 * pixel_size float low resolution limit and high resolution limit, bool approximate bining = F and apply_bining = T 

            projection_volume_3d.CopyFrom(input_3d[bin_index].density_map);
            
            projection_volume_image.Allocate(current_image.logical_x_dimension , current_image.logical_y_dimension, true);
            padded_projection_volume_image.Allocate(current_image.logical_x_dimension  * padding_factor, current_image.logical_y_dimension * padding_factor, false); // as my volume now is already padded so no need to add extra padding

            my_parameters.Init(90.0, 90.0, 90.0 , 0.0, 0.0);
            projection_volume_3d.ExtractSlice(padded_projection_volume_image, my_parameters); // 
            padded_projection_volume_image.SwapRealSpaceQuadrants( ); // must do this step as image is not centered in the box
            padded_projection_volume_image.BackwardFFT( );
            padded_projection_volume_image.ClipInto(&projection_volume_image);

            projection_volume_image.WriteSlice(&my_output_sum_image_filename, bin_index + 1);

            padded_projection_volume_image.Deallocate( );
            projection_volume_image.Deallocate( );

        } 
        if (RASTR == true) { // if a mask is applied and masked model is what will be projected and subtracted
            
            my_masked_volume[bin_index].CopyFrom(&model_volume[bin_index]);

            float filter_edge = 40.0;
            float mask_volume_in_voxels;
   
            //wxPrintf("\nMasking Volume...\n");

            if ( ! model_volume[bin_index].HasSameDimensionsAs(&my_mask) ) {
                wxPrintf("\nVolume and mask file have different dimensions\n");
                DEBUG_ABORT;
            }

            // multiply mask and 3D AA model
            //invert_mask(&my_mask[bin_index]);
            //my_mask_final.Binarise(1.0); //try this
            //my_mask[bin_index].QuickAndDirtyWriteSlices("original_mask_inverted.mrc", 1, my_mask_file->ReturnNumberOfSlices( )* padding_factor);
            if ( filter_radius == 0.0 )
                filter_radius = pixel_size;
            mask_volume_in_voxels = my_masked_volume[bin_index].ApplyMask(my_mask, cosine_edge / pixel_size, outside_weight, pixel_size / filter_radius, pixel_size / filter_edge, outside_value, use_outside_value);
        


            // create the density map to initiate the projections
            // pad 3D masked volume
            masked_3d[bin_index].InitWithDimensions(my_masked_volume[bin_index].logical_x_dimension, my_masked_volume[bin_index].logical_y_dimension, my_masked_volume[bin_index].logical_z_dimension, pixel_size);
            masked_3d[bin_index].density_map->CopyFrom(&my_masked_volume[bin_index]);
            float mask_radius = FLT_MAX; //100 - FLT_MAX
            masked_3d[bin_index].mask_radius = mask_radius;
            masked_3d[bin_index].PrepareForProjections(0.0, 2.0 * pixel_size);
            masked_projection_volume_3d.CopyFrom(masked_3d[bin_index].density_map);

            //masked_projection_volume_image.Allocate(current_image.logical_x_dimension , current_image.logical_y_dimension, true);
            //masked_padded_projection_volume_image.Allocate(current_image.logical_x_dimension  * padding_factor, current_image.logical_y_dimension * padding_factor, false); // as my volume now is already padded so no need to add extra padding
            float phi;
            for ( long model_counter = 0; model_counter < number_of_models; model_counter++ ) {
                //Allocate memory for the masked and padded masked projections
                masked_projection_volume_image.Allocate(current_image.logical_x_dimension , current_image.logical_y_dimension, true);
                masked_padded_projection_volume_image.Allocate(current_image.logical_x_dimension  * padding_factor, current_image.logical_y_dimension * padding_factor, false); // as my volume now is already padded so no need to add extra padding
                //calculate the phi angle
                phi = model_counter * 360.0 / number_of_models;
                my_parameters.Init(phi, 90.0, 90.0 , 0.0, 0.0);
                masked_projection_volume_3d.ExtractSlice(masked_padded_projection_volume_image, my_parameters); // 
                masked_padded_projection_volume_image.SwapRealSpaceQuadrants( ); // must do this step as image is not centered in the box
                masked_padded_projection_volume_image.BackwardFFT( );
                masked_padded_projection_volume_image.ClipInto(&masked_projection_volume_image);

                masked_projection_volume_image.WriteSlice(&my_output_sum_image_filename, model_counter * bins_count +  bin_index + 1 );

                masked_padded_projection_volume_image.Deallocate( );
                masked_projection_volume_image.Deallocate( );

            }

            my_masked_volume->Deallocate( );
            
        }
        model_volume->Deallocate( );
        prepare_projections_progress->Update(bin_index +1);
    }
    sum_images->Deallocate( );
    my_mask.Deallocate( );
    
    delete prepare_projections_progress;



    Image               temporary_image;
    Image               subtracted_image;
    Image               y_axis_aligned_image;
    Image               projection_3d;
    Image               projection_image;
    Image               padded_projection_image;
    AnglesAndShifts     my_parameters_for_subtraction;
    float*              adjusted_x_shifts;
    float*              adjusted_y_shifts;



    //Image* all_subtracted_images[number_of_input_images];
 
    if (SPOT_RASTR == true) {
        wxPrintf("\nSubtracting Azimuthal Average Projections...\n\n");

        ProgressBar* subtract_progress = new ProgressBar(number_of_input_images); 


        //std::vector<std::vector<float>> column_sum_before_subtraction_per_image(number_of_input_images, std::vector<float>(my_input_file.ReturnXSize( ), 0.0));
        //std::vector<std::vector<float>> column_sum_after_subtraction_per_image(number_of_input_images, std::vector<float>(my_input_file.ReturnXSize( ), 0.0));

        my_output_SPOT_RASTR_filename = new MRCFile(SPOT_RASTR_output_filename.ToStdString( ), true, true);  

//#pragma omp for ordered schedule(static, 1)
//output_subtraction_file, output_aligned_file 
#pragma omp parallel num_threads(max_threads)  default(none) shared(number_of_input_images, my_input_file ,best_psi_value, best_x_shift_value, best_y_shift_value, current_image, subtract_progress, \
                                                                     ctf_parameters_stack, max_threads, diameter_bins, bins_count, my_output_SPOT_RASTR_filename, \
                                                                    input_ctf_values_from_star_file, current_ctf, pixel_size, padding_factor, input_3d, x_mask_center, y_mask_center, z_mask_center, adjusted_x_shifts, adjusted_y_shifts) \
                                                                   private(temporary_image, subtracted_image, y_axis_aligned_image, projection_3d, projection_image, padded_projection_image, my_parameters_for_subtraction )

        adjusted_x_shifts = new float[number_of_input_images];
        adjusted_y_shifts = new float[number_of_input_images];
#pragma omp for ordered schedule(dynamic, 1)
        for ( long subtraction_image_counter = 0; subtraction_image_counter < number_of_input_images; subtraction_image_counter++ ) {
            #pragma omp critical 
            // read the current image in the stack
            temporary_image.ReadSlice(&my_input_file, subtraction_image_counter + 1);
            //temporary_image.Normalize( );
            if ( input_ctf_values_from_star_file ) {
                    current_ctf.Init(ctf_parameters_stack[subtraction_image_counter].acceleration_voltage, ctf_parameters_stack[subtraction_image_counter].spherical_aberration, ctf_parameters_stack[subtraction_image_counter].amplitude_contrast, ctf_parameters_stack[subtraction_image_counter].defocus_1, ctf_parameters_stack[subtraction_image_counter].defocus_2, ctf_parameters_stack[subtraction_image_counter].astigmatism_angle, ctf_parameters_stack[subtraction_image_counter].lowest_frequency_for_fitting, ctf_parameters_stack[subtraction_image_counter].highest_frequency_for_fitting, ctf_parameters_stack[subtraction_image_counter].astigmatism_tolerance, ctf_parameters_stack[subtraction_image_counter].pixel_size, ctf_parameters_stack[subtraction_image_counter].additional_phase_shift);
            }
            projection_image.Allocate(current_image.logical_x_dimension, current_image.logical_y_dimension, true);
            padded_projection_image.Allocate(current_image.logical_x_dimension * padding_factor , current_image.logical_y_dimension * padding_factor, false); // as my volume now is already padded so no need to add extra padding

            //if (SPOT_RASTR == true) {
            for (int bin_index = 0; bin_index < bins_count; bin_index++){          
                if (diameter_bins[subtraction_image_counter] == bin_index) { // This should catch any diameter within the range
                    projection_3d.CopyFrom(input_3d[bin_index].density_map);

                } else if (diameter_bins[subtraction_image_counter] < 0 ) { // if the which_bin_index is negative this mean the tube diameter is less than the min. so will be considered with the first bin
                    projection_3d.CopyFrom(input_3d[0].density_map);

                } else if (diameter_bins[subtraction_image_counter] >= bins_count ) { // if the which_bin_index is larger than or equal the bins_count then the diameter is more than the max. so will be added to the last bin
                    projection_3d.CopyFrom(input_3d[bins_count - 1].density_map);

                }

            }

            // Extracting a slice from a 3D volume 
            // angles and shifts are negative to align projection with original input stack
            my_parameters_for_subtraction.Init(0.0, 90.0, 90.0 + best_psi_value[subtraction_image_counter], 0.0, 0.0);
            projection_3d.ExtractSlice(padded_projection_image, my_parameters_for_subtraction);
            //padded_projection_image.QuickAndDirtyWriteSlice("padded_projection_image_after_extract_slice.mrc", subtraction_image_counter + 1);
            padded_projection_image.ZeroCentralPixel( );
            // in project 3D the following line is added before apply ctf and apply ctf is done where absolute is false and apply beam tilt is true
            padded_projection_image.complex_values[0] = projection_3d.complex_values[0];
            padded_projection_image.ApplyCTF(current_ctf, false, true);
            padded_projection_image.PhaseShift(-best_x_shift_value[subtraction_image_counter], -best_y_shift_value[subtraction_image_counter]);
            padded_projection_image.SwapRealSpaceQuadrants( );
            padded_projection_image.BackwardFFT( );
            // de-pad projection
            padded_projection_image.ClipInto(&projection_image);
            //projection_image.QuickAndDirtyWriteSlice("frames_to_be_subtracted_from3Dvolume.mrc", subtraction_image_counter + 1); 

            float average = ReturnAverageOfRealValuesOnVerticalEdges(&projection_image);
            projection_image.AddConstant(-average);


            int pixel_counter              = 0;
            float sum_of_pixelwise_product = 0.0;
            float sum_of_squares           = 0.0;
            float scale_factor             = 0.0;

            for ( long j = 0; j < projection_image.logical_y_dimension; j++ ) {
                for ( long i = 0; i < projection_image.logical_x_dimension; i++ ) {

                    sum_of_pixelwise_product += projection_image.real_values[pixel_counter] * temporary_image.real_values[pixel_counter];
                    sum_of_squares += projection_image.real_values[pixel_counter] * projection_image.real_values[pixel_counter];
                    pixel_counter++;
                }
                pixel_counter += projection_image.padding_jump_value;
            }

            scale_factor  = sum_of_pixelwise_product / sum_of_squares;
            //wxPrintf("The scale factor of image %li is %f \n", subtraction_image_counter+1, scale_factor);

            // multiply by the scaling factor calculated
            projection_image.MultiplyByConstant((scale_factor));
            projection_image.QuickAndDirtyWriteSlice("SPOT_RASTR_projection_image_to_be_subtracted.mrc", subtraction_image_counter + 1);


            // read the image that will be subtracted as we need to only normalize it before subtraction
            subtracted_image.Allocate(my_input_file.ReturnXSize( ), my_input_file.ReturnYSize( ), true );
            #pragma omp critical
            subtracted_image.ReadSlice(&my_input_file, subtraction_image_counter + 1);
            //subtracted_image.Normalize( );

            //subtract the averaged sum image after rotation and translation from the original image
            subtracted_image.SubtractImage(&projection_image);

            padded_projection_image.Deallocate( );
            projection_image.Deallocate( );
            // save the adjusted shift as the extract slice function applied a rotation matrix to the model which shifted x, y shifts more
            // so saving those extra shifts is needed to center the images later and to save the correct shift in the output star file
            RotationMatrix temp_matrix;
            float rotated_x, rotated_y, rotated_z;
            // // generate the full rotation matrix
            temp_matrix.SetToEulerRotation(-(90.0 + best_psi_value[subtraction_image_counter]), -90.0, 0.0);

            temp_matrix.RotateCoords((x_mask_center - current_image.physical_address_of_box_center_x) , (y_mask_center - current_image.physical_address_of_box_center_y)  , (z_mask_center  - current_image.physical_address_of_box_center_x) , rotated_x, rotated_y, rotated_z);

            // center the masked upweighted regions to the center
            adjusted_x_shifts[subtraction_image_counter] = best_x_shift_value[subtraction_image_counter] - rotated_x;
            adjusted_y_shifts[subtraction_image_counter] = best_y_shift_value[subtraction_image_counter] - rotated_y; 

            // write the subtracted images
            #pragma omp critical
            subtracted_image.WriteSlice(my_output_SPOT_RASTR_filename, subtraction_image_counter + 1);

            subtracted_image.Deallocate( );

            if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )

                subtract_progress->Update(subtraction_image_counter + 1);


        }
        delete subtract_progress;
        delete my_output_SPOT_RASTR_filename;

    }

    cisTEMParameters SPOT_RASTR_output_params;

    SPOT_RASTR_output_params.parameters_to_write.SetActiveParameters(POSITION_IN_STACK | IMAGE_IS_ACTIVE | PSI | THETA | PHI | X_SHIFT | Y_SHIFT | DEFOCUS_1 | DEFOCUS_2 | DEFOCUS_ANGLE | PHASE_SHIFT | OCCUPANCY | LOGP | SIGMA | SCORE | PIXEL_SIZE | MICROSCOPE_VOLTAGE | MICROSCOPE_CS | AMPLITUDE_CONTRAST | BEAM_TILT_X | BEAM_TILT_Y | IMAGE_SHIFT_X | IMAGE_SHIFT_Y);
    
    if (SPOT_RASTR == true) {
        SPOT_RASTR_output_params.PreallocateMemoryAndBlank(number_of_input_images);            
        #pragma omp for ordered schedule(static, 1)
        for ( image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
            SPOT_RASTR_output_params.all_parameters[image_counter].position_in_stack                  = image_counter + 1;
            SPOT_RASTR_output_params.all_parameters[image_counter].psi                                = 90.0 + best_psi_value[image_counter];
            SPOT_RASTR_output_params.all_parameters[image_counter].theta                              = 90.0f;
            SPOT_RASTR_output_params.all_parameters[image_counter].phi                                = 0;
            SPOT_RASTR_output_params.all_parameters[image_counter].x_shift                            = -adjusted_x_shifts[image_counter];
            SPOT_RASTR_output_params.all_parameters[image_counter].y_shift                            = -adjusted_y_shifts[image_counter];
            SPOT_RASTR_output_params.all_parameters[image_counter].defocus_1                          = ctf_parameters_stack[image_counter].defocus_1;
            SPOT_RASTR_output_params.all_parameters[image_counter].defocus_2                          = ctf_parameters_stack[image_counter].defocus_2;
            SPOT_RASTR_output_params.all_parameters[image_counter].defocus_angle                      = ctf_parameters_stack[image_counter].astigmatism_angle;
            SPOT_RASTR_output_params.all_parameters[image_counter].phase_shift                        = ctf_parameters_stack[image_counter].additional_phase_shift;
            SPOT_RASTR_output_params.all_parameters[image_counter].image_is_active                    = 1;
            SPOT_RASTR_output_params.all_parameters[image_counter].occupancy                          = 100.0f;
            SPOT_RASTR_output_params.all_parameters[image_counter].logp                               = -1000.0f;
            SPOT_RASTR_output_params.all_parameters[image_counter].sigma                              = 10.0f;
            SPOT_RASTR_output_params.all_parameters[image_counter].pixel_size                         = pixel_size;
            SPOT_RASTR_output_params.all_parameters[image_counter].microscope_voltage_kv              = ctf_parameters_stack[image_counter].acceleration_voltage;
            SPOT_RASTR_output_params.all_parameters[image_counter].microscope_spherical_aberration_mm = ctf_parameters_stack[image_counter].spherical_aberration;
            SPOT_RASTR_output_params.all_parameters[image_counter].amplitude_contrast                 = ctf_parameters_stack[image_counter].amplitude_contrast;
            SPOT_RASTR_output_params.all_parameters[image_counter].beam_tilt_x                        = 0.0f;
            SPOT_RASTR_output_params.all_parameters[image_counter].beam_tilt_y                        = 0.0f;
            SPOT_RASTR_output_params.all_parameters[image_counter].image_shift_x                      = 0.0f;
            SPOT_RASTR_output_params.all_parameters[image_counter].image_shift_y                      = 0.0f;

        }
        

        SPOT_RASTR_output_params.WriteTocisTEMStarFile(SPOT_RASTR_output_star_filename.ToStdString( ));
    } 




//////////// Continue with RASTR processing
    Image  subtracted_RASTR_image;
    Image  centered_upweighted_image;
    float* RASTR_adjusted_x_shifts;
    float* RASTR_adjusted_y_shifts;


    // Variables needed for RASTR mask file
    //Image               mask_projection_image;
    //Image               padded_mask_projection_image;
    //AnglesAndShifts     mask_parameters;
    // change this to pointer later ???   

    if (RASTR == true) { // if RASTR is true then we will use the masked model for subtraction
        wxPrintf("\nSubtracting Masked Azimuthal Average Projections...\n\n");

        ProgressBar* mask_subtract_progress = new ProgressBar(number_of_input_images * number_of_models); 
        MRCFile         my_output_RASTR_filename(RASTR_output_filename.ToStdString( ), true);  

#pragma omp parallel num_threads(max_threads)  default(none) shared(number_of_input_images, my_input_file ,best_psi_value, best_x_shift_value, best_y_shift_value, number_of_models, \
                                                                    ctf_parameters_stack, max_threads, diameter_bins, bins_count,current_image, mask_subtract_progress, \
                                                                    input_ctf_values_from_star_file, current_ctf, pixel_size, padding_factor, masked_3d, x_mask_center, y_mask_center, z_mask_center, RASTR_adjusted_x_shifts, RASTR_adjusted_y_shifts, \
                                                                    mask_projection, input_mask, filter_radius, outside_weight, cosine_edge, center_upweighted, sphere_mask_radius, my_output_RASTR_filename ) \
                                                                    private( temporary_image, subtracted_RASTR_image, y_axis_aligned_image, projection_3d, projection_image, padded_projection_image, my_parameters_for_subtraction,  centered_upweighted_image ) //mask_projection_image, padded_mask_projection_image,mask_parameters,

        RASTR_adjusted_x_shifts = new float[number_of_input_images * number_of_models];
        RASTR_adjusted_y_shifts = new float[number_of_input_images * number_of_models];

        //std::vector<std::vector<float>> column_sum_before_subtraction_per_image(number_of_input_images * number_of_models, std::vector<float>(my_input_file.ReturnXSize( ), 0.0));
        //std::vector<std::vector<float>> column_sum_after_subtraction_per_image(number_of_input_images * number_of_models, std::vector<float>(my_input_file.ReturnXSize( ), 0.0));
        

        float phi;

        for ( long model_counter = 0; model_counter < number_of_models; model_counter++ ) {
        //my_output_RASTR_filename = new MRCFile
    
#pragma omp for ordered schedule(dynamic, 1) //1`
            for ( long subtraction_image_counter = 0; subtraction_image_counter < number_of_input_images; subtraction_image_counter++ ) {
                long current_counter = number_of_input_images * model_counter + subtraction_image_counter;
                #pragma omp critical 
                // read the current image in the stack
                temporary_image.ReadSlice(&my_input_file, subtraction_image_counter + 1);
                //temporary_image.Normalize( );
                if ( input_ctf_values_from_star_file ) {
                        current_ctf.Init(ctf_parameters_stack[subtraction_image_counter].acceleration_voltage, ctf_parameters_stack[subtraction_image_counter].spherical_aberration, ctf_parameters_stack[subtraction_image_counter].amplitude_contrast, ctf_parameters_stack[subtraction_image_counter].defocus_1, ctf_parameters_stack[subtraction_image_counter].defocus_2, ctf_parameters_stack[subtraction_image_counter].astigmatism_angle, ctf_parameters_stack[subtraction_image_counter].lowest_frequency_for_fitting, ctf_parameters_stack[subtraction_image_counter].highest_frequency_for_fitting, ctf_parameters_stack[subtraction_image_counter].astigmatism_tolerance, ctf_parameters_stack[subtraction_image_counter].pixel_size, ctf_parameters_stack[subtraction_image_counter].additional_phase_shift);
                }
                projection_image.Allocate(current_image.logical_x_dimension, current_image.logical_y_dimension, true);
                padded_projection_image.Allocate(current_image.logical_x_dimension * padding_factor , current_image.logical_y_dimension * padding_factor, false); // as my volume now is already padded so no need to add extra padding

                for (int bin_index = 0; bin_index < bins_count; bin_index++){          
                    if (diameter_bins[subtraction_image_counter] == bin_index) { // This should catch any diameter within the range
                        projection_3d.CopyFrom(masked_3d[bin_index].density_map);

                    } else if (diameter_bins[subtraction_image_counter] < 0 ) { // if the which_bin_index is negative this mean the tube diameter is less than the min. so will be considered with the first bin
                        projection_3d.CopyFrom(masked_3d[0].density_map);

                    } else if (diameter_bins[subtraction_image_counter] >= bins_count ) { // if the which_bin_index is larger than or equal the bins_count then the diameter is more than the max. so will be added to the last bin
                        projection_3d.CopyFrom(masked_3d[bins_count - 1].density_map);

                    }
                }

                // Extracting a slice from a 3D volume 
                // angles and shifts are negative to align projection with original input stack
                // use the phi based on the number of models
                phi = model_counter * 360.0 / number_of_models;
                my_parameters_for_subtraction.Init(phi, 90.0, 90.0 + best_psi_value[subtraction_image_counter], 0.0 , 0.0);
                projection_3d.ExtractSlice(padded_projection_image, my_parameters_for_subtraction);
                //padded_projection_image.QuickAndDirtyWriteSlice("padded_projection_image_after_extract_slice.mrc", subtraction_image_counter + 1);
                padded_projection_image.ZeroCentralPixel( );
                // in project 3D the following line is added before apply ctf and apply ctf is done where absolute is false and apply beam tilt is true
                padded_projection_image.complex_values[0] = projection_3d.complex_values[0];
                padded_projection_image.ApplyCTF(current_ctf, false, true);
                padded_projection_image.PhaseShift(-best_x_shift_value[subtraction_image_counter], -best_y_shift_value[subtraction_image_counter]);
                padded_projection_image.SwapRealSpaceQuadrants( );
                padded_projection_image.BackwardFFT( );
                // de-pad projection
                padded_projection_image.ClipInto(&projection_image);
                //void Image::ClipInto(Image* other_image, float wanted_padding_value, bool fill_with_noise, float wanted_noise_sigma, int wanted_coordinate_of_box_center_x, int wanted_coordinate_of_box_center_y, int wanted_coordinate_of_box_center_z) 
                //projection_image.QuickAndDirtyWriteSlice("frames_to_be_subtracted_from3Dvolume.mrc", subtraction_image_counter + 1); 

                float edge_average = ReturnAverageOfRealValuesOnVerticalEdges(&projection_image);
                projection_image.AddConstant(-edge_average);


                int pixel_counter              = 0;
                float sum_of_pixelwise_product = 0.0;
                float sum_of_squares           = 0.0;
                float scale_factor             = 0.0;

                for ( long j = 0; j < projection_image.logical_y_dimension; j++ ) {
                    for ( long i = 0; i < projection_image.logical_x_dimension; i++ ) {

                        sum_of_pixelwise_product += projection_image.real_values[pixel_counter] * temporary_image.real_values[pixel_counter];
                        sum_of_squares += projection_image.real_values[pixel_counter] * projection_image.real_values[pixel_counter];
                        pixel_counter++;
                    }
                    pixel_counter += projection_image.padding_jump_value;
                }

                scale_factor  = sum_of_pixelwise_product / sum_of_squares;
                //wxPrintf("The scale factor of image %li is %f \n", subtraction_image_counter+1, scale_factor);

                // multiply by the scaling factor calculated
                projection_image.MultiplyByConstant((scale_factor));
                #pragma omp critical
                projection_image.QuickAndDirtyWriteSlice("RASTR_projection_image_to_be_subtracted.mrc", current_counter + 1);


                // read the image that will be subtracted as we need to only normalize it before subtraction
                subtracted_RASTR_image.Allocate(my_input_file.ReturnXSize( ), my_input_file.ReturnYSize( ), true );
                #pragma omp critical
                subtracted_RASTR_image.ReadSlice(&my_input_file, subtraction_image_counter + 1);

                //subtract the averaged sum image after rotation and translation from the original image
                subtracted_RASTR_image.SubtractImage(&projection_image);

                // Deallocate the projection image and padded projection image
                padded_projection_image.Deallocate( );
                projection_image.Deallocate( );
                // save the adjusted shift as the extract slice function applied a rotation matrix to the model which shifted x, y shifts more
                // so saving those extra shifts is needed to center the images later and to save the correct shift in the output star file
                RotationMatrix temp_matrix;
                float rotated_x, rotated_y, rotated_z;
                // // generate the full rotation matrix
                temp_matrix.SetToEulerRotation(-(90.0 + best_psi_value[subtraction_image_counter]), -90.0, -phi);

                temp_matrix.RotateCoords((x_mask_center - current_image.physical_address_of_box_center_x) , (y_mask_center - current_image.physical_address_of_box_center_y)  , (z_mask_center  - current_image.physical_address_of_box_center_x) , rotated_x, rotated_y, rotated_z);

                // center the masked upweighted regions to the center
                RASTR_adjusted_x_shifts[current_counter] = best_x_shift_value[subtraction_image_counter] - rotated_x;
                RASTR_adjusted_y_shifts[current_counter] = best_y_shift_value[subtraction_image_counter] - rotated_y; 

                // float average = subtracted_RASTR_image.ReturnAverageOfRealValues( );
                // mask_parameters.Init(phi, 90.0, 90.0 + best_psi_value[subtraction_image_counter], 0.0 , 0.0 ); //* pixel_size
                
                // mask_projection_image.Allocate(current_image.logical_x_dimension , current_image.logical_y_dimension, false); // false as it is in FS
                // padded_mask_projection_image.Allocate(padding_factor * current_image.logical_x_dimension , padding_factor * current_image.logical_y_dimension, false);

                // mask_projection->ExtractSlice(padded_mask_projection_image, mask_parameters); 

                // padded_mask_projection_image.SwapRealSpaceQuadrants( ); // must do this step as image is not centered in the box
                // padded_mask_projection_image.BackwardFFT( );
                // // shift the projected mask images
                // padded_mask_projection_image.PhaseShift(-best_x_shift_value[subtraction_image_counter], -best_y_shift_value[subtraction_image_counter]);
                //rebinarize based on the new threshold >=0.01 should be 1 and less should be 0

                //padded_mask_projection_image.Binarise(0.01f);
                
                 
                // padded_mask_projection_image.ClipInto(&mask_projection_image);
                // mask_projection_image.QuickAndDirtyWriteSlice("mask_projection_images.mrc", current_counter + 1);
                // subtracted_RASTR_image.QuickAndDirtyWriteSlice("subtracted_RASTR_images.mrc", current_counter + 1);
                // float filter_edge = 40.0;
                // float mask_volume_in_voxels;
                // if ( filter_radius == 0.0 )
                //     filter_radius = pixel_size;
                
                // //multiply the mask by the mask subtracted image with the upweighted regions
                // mask_volume_in_voxels = subtracted_RASTR_image.ApplyMask(mask_projection_image, cosine_edge / pixel_size, outside_weight, pixel_size / filter_radius, pixel_size / filter_edge, average, true);
                // subtracted_RASTR_image.QuickAndDirtyWriteSlice("masked_subtracted_RASTR_images.mrc", current_counter + 1);
                //mask_projection_image.Deallocate( );
                //padded_mask_projection_image.Deallocate( );

                // if the user specified that the upweighted regions should be centered before saving them
                if (center_upweighted == true) {
                    subtracted_RASTR_image.QuickAndDirtyWriteSlice("upweighted_regions_before_centering.mrc", current_counter + 1);
                    subtracted_RASTR_image.PhaseShift(RASTR_adjusted_x_shifts[current_counter]  , RASTR_adjusted_y_shifts[current_counter]); 
                    // subtracted_RASTR_image.WriteSlice(&my_output_RASTR_filename, current_counter + 1); // no clipping to make sure that the image is recentered correctly before clipping!!!
                    centered_upweighted_image.Allocate(sphere_mask_radius * 2, sphere_mask_radius * 2, true);
                    subtracted_RASTR_image.ClipInto(&centered_upweighted_image);

                    #pragma omp critical
                    centered_upweighted_image.WriteSlice(&my_output_RASTR_filename, current_counter + 1);
                    centered_upweighted_image.Deallocate( );
                }
                else {
                    #pragma omp critical
                    subtracted_RASTR_image.WriteSlice(&my_output_RASTR_filename, current_counter + 1);
                    subtracted_RASTR_image.Deallocate( );
                }

                if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
                    mask_subtract_progress->Update(number_of_input_images * model_counter + subtraction_image_counter + 1);

            }

        }
        delete mask_subtract_progress;
        delete mask_projection;
        //delete my_output_RASTR_filename;  

    }

    
    // write the parameters file
    cisTEMParameters RASTR_output_params;

    RASTR_output_params.parameters_to_write.SetActiveParameters(POSITION_IN_STACK | IMAGE_IS_ACTIVE | PSI | THETA | PHI | X_SHIFT | Y_SHIFT | DEFOCUS_1 | DEFOCUS_2 | DEFOCUS_ANGLE | PHASE_SHIFT | OCCUPANCY | LOGP | SIGMA | SCORE | PIXEL_SIZE | MICROSCOPE_VOLTAGE | MICROSCOPE_CS | AMPLITUDE_CONTRAST | BEAM_TILT_X | BEAM_TILT_Y | IMAGE_SHIFT_X | IMAGE_SHIFT_Y);

    if (RASTR == true) {
        RASTR_output_params.PreallocateMemoryAndBlank(number_of_input_images * number_of_models);
            
        #pragma omp for ordered schedule(static, 1)
        for ( long model_counter = 0; model_counter < number_of_models; model_counter++ ) {
            for ( image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
                long current_counter = number_of_input_images * model_counter + image_counter; // This is for the position in the stack only
                RASTR_output_params.all_parameters[current_counter].position_in_stack                       = current_counter + 1;
                RASTR_output_params.all_parameters[current_counter].psi                                     = 90.0 + best_psi_value[image_counter];
                RASTR_output_params.all_parameters[current_counter].theta                                   = 90.0f;
                RASTR_output_params.all_parameters[current_counter].phi                                     = model_counter * 360.0 / number_of_models; 
                if (center_upweighted == true) {
                    RASTR_output_params.all_parameters[current_counter].x_shift                             = current_image.physical_address_of_box_center_x;
                    RASTR_output_params.all_parameters[current_counter].y_shift                             = current_image.physical_address_of_box_center_y;
                } else {
                    RASTR_output_params.all_parameters[current_counter].x_shift                             = -RASTR_adjusted_x_shifts[current_counter];
                    RASTR_output_params.all_parameters[current_counter].y_shift                             = -RASTR_adjusted_y_shifts[current_counter];
                } 
                RASTR_output_params.all_parameters[current_counter].defocus_1                               = ctf_parameters_stack[image_counter].defocus_1;
                RASTR_output_params.all_parameters[current_counter].defocus_2                               = ctf_parameters_stack[image_counter].defocus_2;
                RASTR_output_params.all_parameters[current_counter].defocus_angle                           = ctf_parameters_stack[image_counter].astigmatism_angle;
                RASTR_output_params.all_parameters[current_counter].phase_shift                             = ctf_parameters_stack[image_counter].additional_phase_shift;
                RASTR_output_params.all_parameters[current_counter].image_is_active                         = 1;
                RASTR_output_params.all_parameters[current_counter].occupancy                               = 100.0f;
                RASTR_output_params.all_parameters[current_counter].logp                                    = -1000.0f;
                RASTR_output_params.all_parameters[current_counter].sigma                                   = 10.0f;
                RASTR_output_params.all_parameters[current_counter].pixel_size                              = pixel_size;
                RASTR_output_params.all_parameters[current_counter].microscope_voltage_kv                   = ctf_parameters_stack[image_counter].acceleration_voltage;
                RASTR_output_params.all_parameters[current_counter].microscope_spherical_aberration_mm      = ctf_parameters_stack[image_counter].spherical_aberration;
                RASTR_output_params.all_parameters[current_counter].amplitude_contrast                      = ctf_parameters_stack[image_counter].amplitude_contrast;
                RASTR_output_params.all_parameters[current_counter].beam_tilt_x                             = 0.0f;
                RASTR_output_params.all_parameters[current_counter].beam_tilt_y                             = 0.0f;
                RASTR_output_params.all_parameters[current_counter].image_shift_x                           = 0.0f;
                RASTR_output_params.all_parameters[current_counter].image_shift_y                           = 0.0f;

            }
        }

        RASTR_output_params.WriteTocisTEMStarFile(RASTR_output_star_filename.ToStdString( ));
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


//// try this method but add the declaration at the begining
float sum_image_columns_float(Image* current_image) { // Tim's method to calculate the best rotation based on the row sum 
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

    float max_value = *std::max_element(column_sum.begin(), column_sum.end(), 
        [](float a, float b) {
        return std::abs(a) < std::abs(b);});

    return abs(max_value);
}

std::vector<float> average_image_columns(Image* current_image) {
    std::vector<float> column_average(current_image->logical_x_dimension, 0.0);

    long pixel_counter = 0;

     for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
            for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
                long pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                column_average[i] += current_image->real_values[pixel_coord_xy];
                pixel_counter++;
            }
            pixel_counter += current_image->padding_jump_value;
            column_average[i] /= current_image->logical_x_dimension;
    }
     
    return column_average;
}

//// What about calling it find_column_peaks and return the peaks to make it easier to calculate both the diameter and the shift without adding extra functions?
/////// Find Y shift needed to center the tubes after rotation
std::pair<int, int> find_column_sum_peaks(Image* current_image, float min_gap) { // it takes the rotated image and calculate the row sum 
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
    // find the highest two peaks in the row sum that are having a min gap = min_gap
    std::vector<int> peaks;
    // I changed this to be i starting from 1 and ends at vector size -1 to avoid checking out of bound ranges
    for (long i = 1; i < column_sum.size() - 1; ++i) {
        // if the value in the 1D array at index i is higher than the surrounding two values then it is a peak
        // and should be added to the vector peaks!????
        // note here I am saving the positions i.e indices
        if (column_sum[i] > column_sum[i - 1] && column_sum[i] > column_sum[i + 1]) {
            peaks.push_back(i);
        }
    }

    // Sort peaks based on their values so we get the highest values at the begining
    std::sort(peaks.begin(), peaks.end(),
              [&column_sum](int i, int j) { return column_sum[i] > column_sum[j]; });
    
    // Initialize variables needed to loop over the peaks 
    int current_main_peak = 0;
    int highest_positive_peak = peaks[current_main_peak];
    int second_highest_positive_peak = -1;
    bool found_both_positive_peaks = false;
    int highest_negative_peak = peaks[current_main_peak];
    int second_highest_negative_peak = -1;
    bool found_both_negative_peaks  = false;


    // Loop to find the highest two positive peaks with a minimum gap
    // I added size -1 to avoid accessing an out of bound index
    for (long index = 0; index < peaks.size()&& !found_both_positive_peaks; index++) {
        // if the absolute difference between the highest peak position and the next one (indices) is higher than or equal to the gap
        // then we have found both peaks
        // check if the value of the peak is positive
        // go over the positive peaks search
        //wxPrintf("%li\n", index);
        if (peaks[index] > 0) {
            if (std::abs(peaks[index + 1] - highest_positive_peak) >= min_gap) { // removed the std::abs from this line 
                second_highest_positive_peak = peaks[index + 1];
                if (column_sum[second_highest_positive_peak] >= 0.8 * column_sum[highest_positive_peak]) {
                    found_both_positive_peaks = true;
                } else { // if the absolute differences is < than the gap
                    current_main_peak++; // add one to the current_main_peak index and consider this the first highest peak
                    highest_positive_peak = peaks[current_main_peak];
                    index = 0;  // Reset index
                }
            }
        } else if ( peaks[index] < 0) {
            current_main_peak = 0;
            std::vector<int> reversed_peaks = peaks;
            std::reverse(reversed_peaks.begin(), reversed_peaks.end());
            if (std::abs(peaks[index + 1] - highest_negative_peak) >= min_gap) {
                second_highest_negative_peak = peaks[index + 1];
                if (column_sum[second_highest_negative_peak] >= 0.8 * column_sum[highest_negative_peak]) { // if the second highest peak is at least 80% tall from the highest
                    found_both_negative_peaks = true;
                } else { // if the absolute differences is < than the gap
                    current_main_peak++; // add one to the current_main_peak index and consider this the first highest peak
                    highest_negative_peak = peaks[current_main_peak];
                    index = 0;  // Reset index
                }
            }

        }

        float center_peak_index = current_image->logical_y_dimension / 2;
        if (found_both_positive_peaks && ! found_both_negative_peaks) {
            // save the peak with the minimum index to peak1 
            int peak1 = std::min(highest_positive_peak, second_highest_positive_peak);
            // save the peak with the highest index to peak2
            int peak2 = std::max(highest_positive_peak, second_highest_positive_peak);
            // float tube_center = std::abs(peak1 - peak2)/2;
            // float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
            // return -distance_from_center;
            return std::make_pair(peak1, peak2);
        }
        if (!found_both_positive_peaks && found_both_negative_peaks) {
            // save the peak with the minimum index to peak1 
            int peak1 = std::min(highest_negative_peak, second_highest_negative_peak);
            // save the peak with the highest index to peak2
            int peak2 = std::max(highest_negative_peak, second_highest_negative_peak);
            // float tube_center = std::abs(peak1 - peak2)/2;
            // float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
            // return -distance_from_center;
            return std::make_pair(peak1, peak2);
        }
        if (found_both_positive_peaks && found_both_negative_peaks) {
            // check if the min gap in the highest positive peak is closer to the specified min gap by the user than the negative peak
            int min_positive_gap = std::abs(highest_positive_peak - second_highest_positive_peak);
            int min_negative_gap = std::abs(highest_negative_peak - second_highest_negative_peak);
            if (std::abs(min_positive_gap - min_gap) <= std::abs(min_negative_gap - min_gap)) {
                    // then use the positive peak information
                    // save the peak with the minimum index to peak1 
                    int peak1 = std::min(highest_positive_peak, second_highest_positive_peak);
                    // save the peak with the highest index to peak2
                    int peak2 = std::max(highest_positive_peak, second_highest_positive_peak);
                    // float tube_center = std::abs(peak1 - peak2)/2;
                    // float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
                    // return -distance_from_center;
                    return std::make_pair(peak1, peak2);
            } else {
                    // save the peak with the minimum index to peak1 
                    int peak1 = std::min(highest_negative_peak, second_highest_negative_peak);
                    // save the peak with the highest index to peak2
                    int peak2 = std::max(highest_negative_peak, second_highest_negative_peak);
                    //float tube_center = std::abs(peak1 - peak2)/2;
                    //float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
                    //return -distance_from_center;
                    return std::make_pair(peak1, peak2);
            }

        }

    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::pair<int, int> find_row_sum_peaks(Image* current_image, float min_gap) { // it takes the rotated image and calculate the row sum 
    std::vector<float> row_sum(current_image->logical_x_dimension, 0.0);

    long pixel_counter = 0;

     for ( int i = 0; i < current_image->logical_y_dimension; i++ ) {
            for ( int j = 0; j < current_image->logical_x_dimension; j++ ) {
                long pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(j, i, 0);
                row_sum[i] += current_image->real_values[pixel_coord_xy];
                pixel_counter++;
            }
            pixel_counter += current_image->padding_jump_value;
    }
    // find the highest two peaks in the row sum that are having a min gap = min_gap
    std::vector<int> peaks;
    for (long i = 0; i < row_sum.size(); ++i) {
        // if the value in the 1D array at index i is higher than the surrounding two values then it is a peak
        // and should be added to the vector peaks!????
        // note here I am saving the positions i.e indices
        if (row_sum[i] > row_sum[i - 1] && row_sum[i] > row_sum[i + 1]) {
            peaks.push_back(i);
        }
    }

    // Sort peaks based on their values so we get the highest values at the begining
    std::sort(peaks.begin(), peaks.end(),
              [&row_sum](int i, int j) { return row_sum[i] > row_sum[j]; });
    
    // Initialize variables needed to loop over the peaks 
    int current_main_peak = 0;
    int highest_positive_peak = peaks[current_main_peak];
    int second_highest_positive_peak = -1;
    bool found_both_positive_peaks = false;
    int highest_negative_peak = peaks[current_main_peak];
    int second_highest_negative_peak = -1;
    bool found_both_negative_peaks  = false;


    // Loop to find the highest two positive peaks with a minimum gap
    
    for (long index = 0; index < peaks.size() && !found_both_positive_peaks; index++) {
        // if the absolute difference between the highest peak position and the next one (indices) is higher than or equal to the gap
        // then we have found both peaks
        // check if the value of the peak is positive
        // go over the positive peaks search
        if (peaks[index] > 0) {
            if (std::abs(peaks[index + 1] - highest_positive_peak) >= min_gap) { // removed the std::abs from this line 
                second_highest_positive_peak = peaks[index + 1];
                if (row_sum[second_highest_positive_peak] >= 0.8 * row_sum[highest_positive_peak]) {
                    found_both_positive_peaks = true;
                } else { // if the absolute differences is < than the gap
                    current_main_peak++; // add one to the current_main_peak index and consider this the first highest peak
                    highest_positive_peak = peaks[current_main_peak];
                    index = 0;  // Reset index
                }
            }
        } else if ( peaks[index] < 0) {
            current_main_peak = 0;
            std::vector<int> reversed_peaks = peaks;
            std::reverse(reversed_peaks.begin(), reversed_peaks.end());
            if (std::abs(peaks[index + 1] - highest_negative_peak) >= min_gap) {
                second_highest_negative_peak = peaks[index + 1];
                if (row_sum[second_highest_negative_peak] >= 0.8 * row_sum[highest_negative_peak]) { // if the second highest peak is at least 80% tall from the highest
                    found_both_negative_peaks = true;
                } else { // if the absolute differences is < than the gap
                    current_main_peak++; // add one to the current_main_peak index and consider this the first highest peak
                    highest_negative_peak = peaks[current_main_peak];
                    index = 0;  // Reset index
                }
            }

        }

        float center_peak_index = current_image->logical_y_dimension / 2;
        if (found_both_positive_peaks && ! found_both_negative_peaks) {
            // save the peak with the minimum index to peak1 
            int peak1 = std::min(highest_positive_peak, second_highest_positive_peak);
            // save the peak with the highest index to peak2
            int peak2 = std::max(highest_positive_peak, second_highest_positive_peak);
            // float tube_center = std::abs(peak1 - peak2)/2;
            // float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
            // return -distance_from_center;
            return std::make_pair(peak1, peak2);
        }
        if (!found_both_positive_peaks && found_both_negative_peaks) {
            // save the peak with the minimum index to peak1 
            int peak1 = std::min(highest_negative_peak, second_highest_negative_peak);
            // save the peak with the highest index to peak2
            int peak2 = std::max(highest_negative_peak, second_highest_negative_peak);
            // float tube_center = std::abs(peak1 - peak2)/2;
            // float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
            // return -distance_from_center;
            return std::make_pair(peak1, peak2);
        }
        if (found_both_positive_peaks && found_both_negative_peaks) {
            // check if the min gap in the highest positive peak is closer to the specified min gap by the user than the negative peak
            int min_positive_gap = std::abs(highest_positive_peak - second_highest_positive_peak);
            int min_negative_gap = std::abs(highest_negative_peak - second_highest_negative_peak);
            if (std::abs(min_positive_gap - min_gap) <= std::abs(min_negative_gap - min_gap)) {
                    // then use the positive peak information
                    // save the peak with the minimum index to peak1 
                    int peak1 = std::min(highest_positive_peak, second_highest_positive_peak);
                    // save the peak with the highest index to peak2
                    int peak2 = std::max(highest_positive_peak, second_highest_positive_peak);
                    // float tube_center = std::abs(peak1 - peak2)/2;
                    // float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
                    // return -distance_from_center;
                    return std::make_pair(peak1, peak2);
            } else {
                    // save the peak with the minimum index to peak1 
                    int peak1 = std::min(highest_negative_peak, second_highest_negative_peak);
                    // save the peak with the highest index to peak2
                    int peak2 = std::max(highest_negative_peak, second_highest_negative_peak);
                    //float tube_center = std::abs(peak1 - peak2)/2;
                    //float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
                    //return -distance_from_center;
                    return std::make_pair(peak1, peak2);
            }

        }

    }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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


void InitializeCTFSumOfSquares(int numBins, Image& current_image, std::vector<std::vector<float>>* ctf_sum_of_squares) {
    int number_of_pixels = current_image.real_memory_allocated / 2;
    // Allocate and initialize images for each bin
    for (int bin_counter = 0; bin_counter < numBins; ++bin_counter) {
        // Allocate memory for CTF sum of squares based on the image size
        // it needs to be dynamically allocated as image size will differ from one input to another
        (*ctf_sum_of_squares)[bin_counter] = std::vector<float>(number_of_pixels, 0.0f);
        //ZeroFloatArray(ctf_sum_of_squares[bin_counter], number_of_pixels / 2);
    }

}

void ApplyCTFAndReturnCTFSumOfSquares(Image& image, CTF ctf_to_apply, bool absolute, bool apply_beam_tilt, bool apply_envelope, std::vector<float>& ctf_sum_of_squares) {
    MyDebugAssertTrue(image.is_in_memory, "Memory not allocated");
    MyDebugAssertTrue(image.is_in_real_space == false, "Image not in Fourier space");
    MyDebugAssertTrue(image.logical_z_dimension == 1, "Volumes not supported");

    //std::vector<float> squared_ctf_values; // Vector to store squared CTF values

    long pixel_counter = 0;

    float y_coord_sq;
    float x_coord_sq;

    float y_coord;
    float x_coord;

    float frequency_squared;
    float azimuth;
    float ctf_value;

    for (int j = 0; j <= image.physical_upper_bound_complex_y; j++) {
        y_coord    = image.ReturnFourierLogicalCoordGivenPhysicalCoord_Y(j) * image.fourier_voxel_size_y;
        y_coord_sq = pow(y_coord, 2.0);

        for (int i = 0; i <= image.physical_upper_bound_complex_x; i++) {
            x_coord    = i * image.fourier_voxel_size_x;
            x_coord_sq = pow(x_coord, 2);

            // Compute the azimuth
            if (i == 0 && j == 0) {
                azimuth = 0.0;
            } else {
                azimuth = atan2(y_coord, x_coord);
            }

            // Compute the square of the frequency
            frequency_squared = x_coord_sq + y_coord_sq;

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

            if (apply_beam_tilt && (ctf_to_apply.GetBeamTiltX() != 0.0f || ctf_to_apply.GetBeamTiltY() != 0.0f)) {
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


void apply_ctf(Image* current_image, CTF ctf_to_apply, float* ctf_sum_of_squares, bool absolute, bool do_fill_sum_of_squares) {
    float y_coord_sq;
    float x_coord_sq;

    float y_coord;
    float x_coord;

    float frequency_squared;
    float azimuth;
    float ctf_value;

    long pixel_counter = 0;

    for ( int j = 0; j <= current_image->physical_upper_bound_complex_y; j++ ) {
        y_coord    = current_image->ReturnFourierLogicalCoordGivenPhysicalCoord_Y(j) * current_image->fourier_voxel_size_y;
        y_coord_sq = powf(y_coord, 2.0);

        for ( int i = 0; i <= current_image->physical_upper_bound_complex_x; i++ ) {
            x_coord    = i * current_image->fourier_voxel_size_x;
            x_coord_sq = powf(x_coord, 2.0);

            // Compute the azimuth
            if ( i == 0 && j == 0 ) {
                azimuth = 0.0;
            }
            else {
                azimuth = atan2f(y_coord, x_coord);
            }

            // Compute the square of the frequency
            frequency_squared = x_coord_sq + y_coord_sq;
            ctf_value         = ctf_to_apply.Evaluate(frequency_squared, azimuth);

            // phase-flip
            if ( absolute )
                ctf_value = fabsf(ctf_value);

            current_image->complex_values[pixel_counter] *= ctf_value;
            if ( do_fill_sum_of_squares ) {
                ctf_sum_of_squares[pixel_counter] += powf(ctf_value, 2);
            }
            pixel_counter++;
        }
    }
}

void normalize_image(Image* input_image, float pixel_size, float mask_falloff) {
    // Normalize background variance and average
    float variance;
    float average;

    // subtract mean value from each image pixel to get a zero-mean
    // divide each pixel value by standard deviation to have unit-variance
    variance = input_image->ReturnVarianceOfRealValues(input_image->physical_address_of_box_center_x - (mask_falloff / pixel_size), 0.0, 0.0, 0.0, true);
    average  = input_image->ReturnAverageOfRealValues(input_image->physical_address_of_box_center_x - (mask_falloff / pixel_size), true);

    if ( variance == 0.0f ) {
        input_image->SetToConstant(0.0f);
    }
    else {
        input_image->AddMultiplyConstant(-average, 1.0 / sqrtf(variance));
    }
}

// // Function to ensure the angle is within the range [0, 360)
// float angle_within360(float angle) {
//     if (angle < 0.0) {
//         angle += 360.0;
//         return angle_within360(angle);
//     } else if (angle >= 360.0) {
//         angle -= 360.0;
//         return angle_within360(angle);
//     } else {
//         return angle;
//     }
// }

// calculates the average of real values on the vertical edges
float ReturnAverageOfRealValuesOnVerticalEdges(Image* current_image) {
    double sum;
    long   number_of_pixels;
    int    pixel_counter;
    int    line_counter;
    int    plane_counter;
    long   address;

    sum              = 0.0;
    number_of_pixels = 0;
    address          = 0;

    if ( current_image->logical_z_dimension == 1 ) {
        // Two-dimensional image
        for ( line_counter = 0; line_counter < current_image->logical_y_dimension; line_counter++ ) {
            sum += current_image->real_values[address];
            address += current_image->logical_x_dimension - 1;
            sum += current_image->real_values[address];
            address += current_image->padding_jump_value + 1;
            number_of_pixels += 2;
        }
    }
    else {
        // Three-dimensional volume
        for ( plane_counter = 0; plane_counter < current_image->logical_z_dimension; plane_counter++ ) {
            for ( line_counter = 0; line_counter < current_image->logical_y_dimension; line_counter++ ) {
                if ( line_counter == 0 || line_counter == current_image->logical_y_dimension - 1 ) {
                    // First and last line of that section
                    for ( pixel_counter = 0; pixel_counter < current_image->logical_x_dimension; pixel_counter++ ) {
                        sum += current_image->real_values[address];
                        address++;
                    }
                    address += current_image->padding_jump_value;
                    number_of_pixels += current_image->logical_x_dimension;
                }
                else {
                    // All other lines (only count first and last pixel)
                    sum += current_image->real_values[address];
                    address += current_image->logical_x_dimension - 1;
                    sum += current_image->real_values[address];
                    address += current_image->padding_jump_value + 1;
                    number_of_pixels += 2;
                }
            }
        }
    }

    return sum / float(number_of_pixels);
}

// Function to calculate RMSD
float calculateRMSD(const std::vector<float>& vec1, const std::vector<float>& vec2) {

    float sumSquaredDifferences = 0.0;
    for (size_t i = 0; i < vec1.size(); ++i) {
        float diff = vec1[i] - vec2[i];
        sumSquaredDifferences += diff * diff;
    }

    float meanSquaredDifferences = sumSquaredDifferences / vec1.size();
    return std::sqrt(meanSquaredDifferences);
}

std::vector<float> GetAverageColumnsValues(const std::vector<std::vector<float>>& sum_columns_data) {
    
    size_t numVectors = sum_columns_data.size();
    size_t vectorSize = sum_columns_data[0].size();
    
    std::vector<float> averages(vectorSize, 0.0);
    
    for (size_t i = 0; i < vectorSize; ++i) {
        float sum = 0.0;
        for (size_t j = 0; j < numVectors; ++j) {
            sum += sum_columns_data[j][i];
        }
        averages[i] = sum / numVectors;
    }
    
    return averages;
}

void invert_mask(Image* mask_file) {
    // inverts binarized mask pixel values (i.e., 0 changes to 1 and 1 changes to 0)
    for ( long pixel_counter = 0; pixel_counter < mask_file->real_memory_allocated; pixel_counter++ ) {
        if ( mask_file->real_values[pixel_counter] == 0.0 )
            mask_file->real_values[pixel_counter] = 1.0; 
        else 
            mask_file->real_values[pixel_counter] = 0.0; 
    }
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
                d = sqrtf(dx*dx + dy*dy + dz*dz);
                if (d < radius) {
                    mask_file->real_values[pixel_counter] = 0.0;
                } else {
                    mask_file->real_values[pixel_counter] = 1.0;
                }
                pixel_counter ++;
            }
            pixel_counter += mask_file->padding_jump_value;
        }

    }
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
                d = sqrtf(dx*dx + dy*dy + dz*dz);
                if (d < radius) {
                    mask_file->real_values[pixel_counter] = 1.0;
                } else {
                    mask_file->real_values[pixel_counter] = 0.0;
                }
                pixel_counter ++;
            }
            pixel_counter += mask_file->padding_jump_value;
        }

    }
    //mask_file->QuickAndDirtyWriteSlices("make_white_sphere_mask_inside_function.mrc", 1, mask_file->logical_z_dimension);
}


// //// What about calling it find_column_peaks and return the peaks to make it easier to calculate both the diameter and the shift without adding extra functions?
// /////// Find Y shift needed to center the tubes after rotation
// std::pair<int, int> find_column_peaks(Image* current_image, float min_gap, float pixel_size) { // it takes the rotated image and calculate the row sum 

//     float center_peak_index = current_image.logical_y_dimension / 2;
//     std::vector<float> column_sum_orig(current_image->logical_x_dimension, 0.0);
//     std::vector<float> column_sum_gaussian(current_image->logical_x_dimension, 0.0);

//     // create a copy of the image where we apply a gaussian filter to be easier to find the peaks in the image
//     Image gaussain_image;
//     gaussain_image.CopyFrom(&current_image);

//     gaussain_image.ForwardFFT( );
//     gaussian_image.GaussianLowPassFilter((pixel_size*2)/150);
//     gaussain_image.BackwardFFT( );

//     // get the column sum of the original image values
//     long pixel_counter = 0;

//      for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
//             for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
//                 long pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
//                 column_sum_orig[i] += current_image->real_values[pixel_coord_xy];
//                 pixel_counter++;
//             }
//             pixel_counter += current_image->padding_jump_value;
//     }

//     // get the column sum of the gaussian image values
//     pixel_counter = 0;

//      for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
//             for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
//                 long pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
//                 column_sum_gaussian[i] += current_image->real_values[pixel_coord_xy];
//                 pixel_counter++;
//             }
//             pixel_counter += current_image->padding_jump_value;
//     }

//     // first we need to normalize the values in the vector by subtracting the min value from each element
//     float min_value_orig = std::min_element(column_sum_orig.begin(), column_sum_orig.end());
//     std::vector<float> column_sum_orig_norm_min= column_sum_orig;

//     for (float& elem : column_sum_orig_norm_min) {
//         elem -= min_value_orig;
//     }
//     // Second we need to normalize the values in the vector by subtracting the max value from each element

//     float max_value_orig = std::max_element(column_sum_orig.begin(), column_sum_orig.end());
//     std::vector<float> column_sum_orig_norm_max= column_sum_orig;

//     for (float& elem : column_sum_orig_norm_max) {
//         elem = max_value_orig - elem;
//     }

//     // find the positive peaks 
//     std::vector<int> orig_positive_peaks;
//     // I changed this to be i starting from 1 and ends at vector size -1 to avoid checking out of bound ranges
//     // and as also we don't care about the peaks on the edges
//     for (long i = 1; i < column_sum_orig_norm_min.size() - 1; ++i) {
//         // if the value in the 1D array at index i is higher than the surrounding two values then it is a peak
//         // and should be added to the vector peaks!????
//         // note here I am saving the positions i.e indices
//         if (column_sum_orig_norm_min[i] > column_sum_orig_norm_min[i - 1] && column_sum_orig_norm_min[i] > column_sum_orig_norm_min[i + 1]) {
//             orig_positive_peaks.push_back(i);
//         }
//     }

//     //find the negative peaks that are found in te inverted vector (vec that was normalized by subtracting the maximum value)
//     std::vector<int> orig_negative_peaks;
//     // I changed this to be i starting from 1 and ends at vector size -1 to avoid checking out of bound ranges
//     // and as also we don't care about the peaks on the edges
//     for (long i = 1; i < column_sum_orig_norm_max.size() - 1; ++i) {
//         // if the value in the 1D array at index i is higher than the surrounding two values then it is a peak
//         // and should be added to the vector peaks!????
//         // note here I am saving the positions i.e indices
//         if (column_sum_orig_norm_max[i] > column_sum_orig_norm_max[i - 1] && column_sum_orig_norm_max[i] > column_sum_orig_norm_max[i + 1]) {
//             orig_negative_peaks.push_back(i);
//         }
//     }

//     // Sort peaks based on their values so we get the highest values at the begining
//     std::sort(orig_positive_peaks.begin(), orig_positive_peaks.end(),
//               [&column_sum_orig_norm_min](int i, int j) { return column_sum_orig_norm_min[i] > column_sum_orig_norm_min[j]; });
//     // Sort peaks based on their values so we get the highest values at the begining
//     std::sort(orig_negative_peaks.begin(), orig_negative_peaks.end(),
//               [&column_sum_orig_norm_max](int i, int j) { return column_sum_orig_norm_max[i] > column_sum_orig_norm_max[j]; });
    
//     // now we have the indices of the peaks sorted based on their value in the normalized vectors in a descending orider where 
//     // the peaks that has the largest value in the sum column are listed first then the one with less values are listed later

//     // as we don't know the size of the diff list we can initiate it without a size
//     std::vector<std::pair<float, int>> diff_list; 
//     std::vector<int> valid_negative_peaks;

//     for (int pos_peak : orig_positive_peaks) {
//         if (pos_peak > center_peak_index){
//             for (int neg_peak : orig_negative_peaks) {
//                 if (neg_peak > pos_peak) {
//                     valid_negative_peaks.push_back(neg_peak);
//                 }
//             }

//         }
//         else {
//             for (int neg_peak : orig_negative_peaks) {
//                 if (neg_peak < pos_peak) {
//                     valid_negative_peaks.push_back(neg_peak);
//                 }
//             }
//         }

//         int closest_neg_peak;
//         if(!valid_neg_peaks.empty()) {
//         //// DO SOMETHING
//             closest_neg_peak = std::min_element(valid_negative_peaks.begin(), valid_negative_peaks.end(), [pos_peak](int a, int b) {return std::abs(a - pos_peak) < std::abs(b - pos_peak)});
//         }

//         float diff = std::abs(column_sum_orig_norm_min[pos_peak] - column_sum_orig_norm_min[closest_neg_peak]);
//         diff_list.push_back({diff, closest_neg_peak});
//     }

    










//     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     // find the highest two peaks in the row sum that are having a min gap = min_gap
//     std::vector<int> gaussian_peaks;
//     // I changed this to be i starting from 1 and ends at vector size -1 to avoid checking out of bound ranges
//     for (long i = 1; i < column_sum_gaussian.size() - 1; ++i) {
//         // if the value in the 1D array at index i is higher than the surrounding two values then it is a peak
//         // and should be added to the vector peaks!????
//         // note here I am saving the positions i.e indices
//         if (column_sum_gaussian[i] > column_sum_gaussian[i - 1] && column_sum_gaussian[i] > column_sum_gaussian[i + 1]) {
//             gaussian_peaks.push_back(i);
//         }
//     }


//     // Initialize variables needed to loop over the peaks 
//     int current_main_peak = 0;
//     int highest_positive_peak = peaks[current_main_peak];
//     int second_highest_positive_peak = -1;
//     bool found_both_positive_peaks = false;
//     int highest_negative_peak = peaks[current_main_peak];
//     int second_highest_negative_peak = -1;
//     bool found_both_negative_peaks  = false;


//     // Loop to find the highest two positive peaks with a minimum gap
//     // I added size -1 to avoid accessing an out of bound index
//     for (long index = 0; index < peaks.size()&& !found_both_positive_peaks; index++) {
//         // if the absolute difference between the highest peak position and the next one (indices) is higher than or equal to the gap
//         // then we have found both peaks
//         // check if the value of the peak is positive
//         // go over the positive peaks search
//         wxPrintf("%li\n", index);
//         if (peaks[index] > 0) {
//             if (std::abs(peaks[index + 1] - highest_positive_peak) >= min_gap) { // removed the std::abs from this line 
//                 second_highest_positive_peak = peaks[index + 1];
//                 if (column_sum[second_highest_positive_peak] >= 0.8 * column_sum[highest_positive_peak]) {
//                     found_both_positive_peaks = true;
//                 } else { // if the absolute differences is < than the gap
//                     current_main_peak++; // add one to the current_main_peak index and consider this the first highest peak
//                     highest_positive_peak = peaks[current_main_peak];
//                     index = 0;  // Reset index
//                 }
//             }
//         } else if ( peaks[index] < 0) {
//             current_main_peak = 0;
//             std::vector<int> reversed_peaks = peaks;
//             std::reverse(reversed_peaks.begin(), reversed_peaks.end());
//             if (std::abs(peaks[index + 1] - highest_negative_peak) >= min_gap) {
//                 second_highest_negative_peak = peaks[index + 1];
//                 if (column_sum[second_highest_negative_peak] >= 0.8 * column_sum[highest_negative_peak]) { // if the second highest peak is at least 80% tall from the highest
//                     found_both_negative_peaks = true;
//                 } else { // if the absolute differences is < than the gap
//                     current_main_peak++; // add one to the current_main_peak index and consider this the first highest peak
//                     highest_negative_peak = peaks[current_main_peak];
//                     index = 0;  // Reset index
//                 }
//             }

//         }

//         float center_peak_index = current_image->logical_y_dimension / 2;
//         if (found_both_positive_peaks && ! found_both_negative_peaks) {
//             // save the peak with the minimum index to peak1 
//             int peak1 = std::min(highest_positive_peak, second_highest_positive_peak);
//             // save the peak with the highest index to peak2
//             int peak2 = std::max(highest_positive_peak, second_highest_positive_peak);
//             // float tube_center = std::abs(peak1 - peak2)/2;
//             // float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
//             // return -distance_from_center;
//             return std::make_pair(peak1, peak2);
//         }
//         if (!found_both_positive_peaks && found_both_negative_peaks) {
//             // save the peak with the minimum index to peak1 
//             int peak1 = std::min(highest_negative_peak, second_highest_negative_peak);
//             // save the peak with the highest index to peak2
//             int peak2 = std::max(highest_negative_peak, second_highest_negative_peak);
//             // float tube_center = std::abs(peak1 - peak2)/2;
//             // float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
//             // return -distance_from_center;
//             return std::make_pair(peak1, peak2);
//         }
//         if (found_both_positive_peaks && found_both_negative_peaks) {
//             // check if the min gap in the highest positive peak is closer to the specified min gap by the user than the negative peak
//             int min_positive_gap = std::abs(highest_positive_peak - second_highest_positive_peak);
//             int min_negative_gap = std::abs(highest_negative_peak - second_highest_negative_peak);
//             if (std::abs(min_positive_gap - min_gap) <= std::abs(min_negative_gap - min_gap)) {
//                     // then use the positive peak information
//                     // save the peak with the minimum index to peak1 
//                     int peak1 = std::min(highest_positive_peak, second_highest_positive_peak);
//                     // save the peak with the highest index to peak2
//                     int peak2 = std::max(highest_positive_peak, second_highest_positive_peak);
//                     // float tube_center = std::abs(peak1 - peak2)/2;
//                     // float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
//                     // return -distance_from_center;
//                     return std::make_pair(peak1, peak2);
//             } else {
//                     // save the peak with the minimum index to peak1 
//                     int peak1 = std::min(highest_negative_peak, second_highest_negative_peak);
//                     // save the peak with the highest index to peak2
//                     int peak2 = std::max(highest_negative_peak, second_highest_negative_peak);
//                     //float tube_center = std::abs(peak1 - peak2)/2;
//                     //float distance_from_center = (peak1 + peak2)/2 - center_peak_index; 
//                     //return -distance_from_center;
//                     return std::make_pair(peak1, peak2);
//             }

//         }

//     }

// }