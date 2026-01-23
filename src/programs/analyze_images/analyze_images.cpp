#include "../../core/core_headers.h"

class
        analyze_images : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

IMPLEMENT_APP(analyze_images)

void analyze_images::DoInteractiveUserInput( ) {

    UserInput* my_input = new UserInput("analyze_images", 1.0);

    wxString input_filename       = my_input->GetFilenameFromUser("Input image stack1", "Filename of the stack of images to normalize", "input_stack1.mrc", true);
    wxString input_image_filename = my_input->GetFilenameFromUser("Input image stack2", "Filename of the stack of images to normalize", "input_stack2.mrc", true);
    wxString star_filename        = my_input->GetFilenameFromUser("Input star file", "The input star file", "my_parameters.star", true);

    // wxString input_star_filename  = my_input->GetFilenameFromUser("Input star file", "The input star file", "my_parameters.star", true);
    // wxString output_star_filename = my_input->GetFilenameFromUser("Output star file", "Output star file with the random phi + reference RASTR phi", "output_parameters.star", false);
    // int      number_of_models     = my_input->GetIntFromUser("Number of models used to generate that RASTR particle stack", "The number of models used to generate the RASTR particle stack, minimum number is 1 models", "1", 1);
    delete my_input;

    my_current_job.Reset(3);
    my_current_job.ManualSetArguments("ttt", input_filename.ToUTF8( ).data( ), input_image_filename.ToUTF8( ).data( ), star_filename.ToUTF8( ).data( )); //, output_star_filename.ToUTF8( ).data( ), number_of_models);
}

// override the do calculation method which will be what is actually run..

bool analyze_images::DoCalculation( ) {
    wxString input_filename       = my_current_job.arguments[0].ReturnStringArgument( );
    wxString input_image_filename = my_current_job.arguments[1].ReturnStringArgument( );
    wxString star_filename        = my_current_job.arguments[2].ReturnStringArgument( );

    // wxString input_star_filename  = my_current_job.arguments[1].ReturnStringArgument( );
    // wxString output_star_filename = my_current_job.arguments[2].ReturnStringArgument( );
    // int      number_of_models     = my_current_job.arguments[3].ReturnIntegerArgument( );

    MRCFile my_input1_file(input_filename.ToStdString( ), false);
    MRCFile my_input2_file(input_image_filename.ToStdString( ), false);

    int  x_dim                  = my_input1_file.ReturnXSize( );
    int  y_dim                  = my_input1_file.ReturnYSize( );
    long number_of_input_images = my_input2_file.ReturnNumberOfSlices( );

    Image ref_img;
    //Image ref_phase_img;
    Image exp_img;

    ref_img.Allocate(x_dim, y_dim, true);
    //ref_phase_img.Allocate(x_dim, y_dim, true);
    exp_img.Allocate(x_dim, y_dim, true);
    ref_img.SetToConstant(0.0);
    //ref_phase_img.SetToConstant(0.0);
    exp_img.SetToConstant(0.0);

    // //CisTEM star
    // cisTEMParameters input_star_file;
    // if ( (is_running_locally && ! DoesFileExist(star_filename.ToStdString( ))) ) {
    //     SendErrorAndCrash(wxString::Format("Error: Input star file %s not found\n", star_filename));
    // }
    // input_star_file.ReadFromcisTEMStarFile(star_filename.ToStdString( ));

    // Relion star
    BasicStarFileReader input_star_file;
    if ( (is_running_locally && ! DoesFileExist(star_filename.ToStdString( ))) ) {
        SendErrorAndCrash(wxString::Format("Error: Input star file %s not found\n", star_filename));
    }
    input_star_file.ReadFile(star_filename.ToStdString( ));

    float pixel_size = 2.0; // will change it manually here

    for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {

        wxPrintf("\n Image number is %li\n", image_counter + 1);

        exp_img.ReadSlice(&my_input2_file, image_counter + 1);

        exp_img.ForwardFFT( );
        exp_img.ZeroCentralPixel( );
        exp_img.BackwardFFT( );
        exp_img.Normalize( );

        float psi = input_star_file.ReturnPsi(image_counter);
        float phi = input_star_file.ReturnPhi(image_counter);

        exp_img.Rotate2DInPlace(-psi + 90.0, FLT_MAX);
        wxPrintf("\n Image psi is %f\n", -psi + 90.0); // -psi - 90.0 is not aligned

        exp_img.QuickAndDirtyWriteSlice("exp_img_rotated.mrc", image_counter + 1);

        ref_img.ReadSlice(&my_input1_file, 1);
        //ref_phase_img.ReadSlice(&my_input1_file, 1);

        ref_img.ForwardFFT( );
        ref_img.ZeroCentralPixel( );
        ref_img.BackwardFFT( );
        ref_img.Normalize( );

        ref_img.QuickAndDirtyWriteSlice("ref_image_normalized.mrc", image_counter + 1);

        Peak current_peak;
        //Peak phase_peak;

        // float peak_scalar = 20.0f;

        // ref_phase_img.CalculatePhaseCrossCorrelationImageWith(exp_img, phase_peak, peak_scalar, true);

        // ref_phase_img.QuickAndDirtyWriteSlice("ref_phase_img_cross_corr.mrc", 1);

        // wxPrintf("\n Phase current peak value is %f\n current peak shift in x is %f\n current peak shift in y is %f\n\n", phase_peak.value, phase_peak.x, phase_peak.y);

        // calculate the cross correlation of the reference image with the rotated image
        ref_img.CalculateCrossCorrelationImageWith(&exp_img); //rotated_image

        // // mask the crosscorr image to only certain part to avoid wrong shift
        // float mask_edge = x_dim * 0.05f; // 5% of the pixels will constitute the soft edge

        // ref_img.CosineMask(x_dim * 0.25, mask_edge, false, false, 0.0f);

        float inner_radius_for_peak_search = 0.0; // inner radius should be set to 0
        float outer_radius_for_peak_search = x_dim / 2;
        // find the peak from the cross corrlation to get the values from the ref_img
        current_peak = ref_img.FindPeakWithParabolaFit(inner_radius_for_peak_search, outer_radius_for_peak_search);
        ref_img.QuickAndDirtyWriteSlice("ref_img_cross_corr.mrc", image_counter + 1);

        wxPrintf("\n current peak value is %f\n current peak shift in x is %f\n current peak shift in y is %f\n", current_peak.value, current_peak.x, current_peak.y);

        // adjusted shift based on RASTR and diameTR
        // Δx​=−d*sin(ψ) Δy=−d*cos(ψ) where d is the distance from center x dim when the tuve is rotated which is the one we get when we do the crosscorrelation calculations
        float new_x_shift = -current_peak.x * sinf(deg_2_rad(psi)); //psi * (pi_v<float> / 180.f)
        float new_y_shift = -current_peak.x * cosf(deg_2_rad(psi)); //psi * (pi_v<float> / 180.f)

        wxPrintf("\nDistance is %f new x shift is %f and new Y shift is %f", -current_peak.x, new_x_shift * pixel_size, new_y_shift * pixel_size);

        Image adjusted_img;
        adjusted_img.Allocate(x_dim, y_dim, true);
        adjusted_img.SetToConstant(0.0);
        adjusted_img.ReadSlice(&my_input2_file, image_counter + 1);
        adjusted_img.PhaseShift(new_x_shift, new_y_shift, 0.0);
        adjusted_img.Rotate2DInPlace(10.0, FLT_MAX);

        RotationMatrix temp_matrix;
        float          rotated_x, rotated_y, rotated_z;
        float          adjusted_x_shifts, adjusted_y_shifts, adjusted_z_shift;
        float          rotated2d_x, rotated2d_y;
        temp_matrix.SetToEulerRotation(-psi + 90.0, 90.0, phi);
        //temp_matrix.ReturnTransposed( );
        // assuming no Y shift will be applied to ensure everything is centered
        temp_matrix.RotateCoords(-current_peak.x, 0.0, 0.0, rotated_x, rotated_y, rotated_z);
        // correct the shift happenning because of extractslice - Note x and y shift needs to be flipped
        // since projection is centered so negative rotation and shift is needed here
        adjusted_x_shifts = rotated_x;
        adjusted_y_shifts = rotated_y;
        adjusted_z_shift  = rotated_z;

        //float x_value = 112.0;
        float y_value = 0.0;
        wxPrintf("\n\nRotateCoords Adjusted x shift is %f and adjusted Y shift is %f and adjusted_z_shift is %f\n", adjusted_x_shifts, adjusted_y_shifts, adjusted_z_shift);
        temp_matrix.RotateCoords2D(current_peak.x, y_value, rotated2d_x, rotated2d_y);
        wxPrintf("\n\nRotateCoords2D Adjusted x shift is %f and adjusted Y shift is %f", rotated2d_x * pixel_size, rotated2d_y * pixel_size);

        // adjusted_img.QuickAndDirtyWriteSlice("adjusted_image_shifted_then_rotated.mrc", image_counter + 1);
    }

    // if ( current_tuned_peak.value > local_best_corr_score ) {
    //     local_best_corr_score = current_tuned_peak.value;
    //     local_best_psi        = tuned_psi;
    //     local_best_x_shift    = current_tuned_peak.x;
    //     local_best_y_shift    = current_tuned_peak.y;
    // }

    wxPrintf("\n\n");

    return true;
}
