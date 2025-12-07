#include "../../core/core_headers.h"

class
        create_average_power_spectrum : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

IMPLEMENT_APP(create_average_power_spectrum)

void create_average_power_spectrum::DoInteractiveUserInput( ) {
    wxString input_images;
    wxString output_power_spectrum;
    float    pixel_size;
    float    padding_factor = 1.0;
    bool     crop;
    int      max_threads;

    UserInput* my_input   = new UserInput("create_average_power_spectrum", 1.00);
    input_images          = my_input->GetFilenameFromUser("Input images file name", "Filen name of helical tube stack aligned vertically", "helical_stack.mrc", true);
    output_power_spectrum = my_input->GetFilenameFromUser("Output average power spectrum file name", "The power spectrum of the tube images file name", "helical_stack_ps.mrc", false);
    pixel_size            = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
    padding_factor        = my_input->GetFloatFromUser("Padding factor", "Factor value to be used for padding the input images", "0.0", 1.0);
    crop                  = my_input->GetYesNoFromUser("Crop the average power spectrum to the original image size?", "Crop the average power spectrum to the original image size or keep it padded", "No");
#ifdef _OPENMP
    max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "when threading, what is the max threads to run", "1", 1);
#else
    max_threads = 1;
#endif
    delete my_input;

    my_current_job.Reset(7);
    my_current_job.ManualSetArguments("ttffbi", input_images.ToUTF8( ).data( ), output_power_spectrum.ToUTF8( ).data( ),
                                      pixel_size, padding_factor, crop, max_threads);
}

// override the do calculation method which will be what is actually run..

bool create_average_power_spectrum::DoCalculation( ) {
    wxString input_images          = my_current_job.arguments[0].ReturnStringArgument( );
    wxString output_power_spectrum = my_current_job.arguments[1].ReturnStringArgument( );
    float    pixel_size            = my_current_job.arguments[2].ReturnFloatArgument( );
    float    padding_factor        = my_current_job.arguments[3].ReturnFloatArgument( );
    bool     crop                  = my_current_job.arguments[4].ReturnBoolArgument( );
    int      max_threads           = my_current_job.arguments[5].ReturnIntegerArgument( );

    MRCFile my_input_images(input_images.ToStdString( ), false);
    MRCFile my_output_power_spectrum(output_power_spectrum.ToStdString( ), true);
    long    number_of_input_images = my_input_images.ReturnNumberOfSlices( );
    int     x_dim                  = my_input_images.ReturnXSize( );
    int     y_dim                  = my_input_images.ReturnYSize( );

    Image my_image;
    Image power_spectrum_sum_image;
    Image power_spectrum_image;
    power_spectrum_sum_image.Allocate(x_dim * padding_factor, y_dim * padding_factor, true);
    power_spectrum_sum_image.SetToConstant(0.0);

    wxPrintf("\nComputing average power spectrum image...\n\n");

    ProgressBar* power_progress = new ProgressBar(number_of_input_images);

#pragma omp parallel for schedule(static, 1) num_threads(max_threads) default(none) shared(my_input_images, pixel_size, padding_factor, crop, number_of_input_images, x_dim, y_dim, pixel_size, my_output_power_spectrum, power_spectrum_sum_image, power_progress) private(my_image, power_spectrum_image)

    for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {

        //my_image.Allocate(x_dim, y_dim, true);
        //my_image.SetToConstant(0.0);
#pragma omp critical
        my_image.ReadSlice(&my_input_images, image_counter + 1);

        my_image.Normalize( );
        // float average;
        // // average = my_image.ReturnAverageOfRealValuesOnEdges( );
        // average = my_image.ReturnAverageOfRealValues( );
        // //padd the image so that when rotated it doesn't mess up the outer edges
        // my_image.Resize(padding_factor * my_image.logical_x_dimension, padding_factor * my_image.logical_y_dimension, 1, average);
        my_image.Resize(padding_factor * my_image.logical_x_dimension, padding_factor * my_image.logical_y_dimension, 1, 0.0); //padding with zero??

        //Image power_spectrum_image;
        power_spectrum_image.Allocate(x_dim * padding_factor, y_dim * padding_factor, false);
        power_spectrum_image.SetToConstant(0.0);
        my_image.ForwardFFT( );
        my_image.ZeroCentralPixel( );
        my_image.ComputeAmplitudeSpectrumFull2D(&power_spectrum_image);
        // use thresholding to decrease noise in PS
        float image_average;
        float image_sd;
        float image_threshold;
        Image binary_mask;
        image_average = power_spectrum_image.ReturnAverageOfRealValues( );
        image_sd      = sqrt(power_spectrum_image.ReturnVarianceOfRealValues( ));
        // Threshold value is 2 sd away from mean to eliminate any outliers
        image_threshold = image_average + (3 * image_sd);
        binary_mask.CopyFrom(&power_spectrum_image);
        binary_mask.Binarise(image_threshold);
        //binary_mask.QuickAndDirtyWriteSlice("binary_mask.mrc", 1);
        float cosine_edge       = 10.0;
        float outside_weight    = 0.0;
        float filter_radius     = 0.0;
        float outside_value     = 0.0;
        bool  use_outside_value = false;
        float filter_edge       = 40.0;
        power_spectrum_image.ApplyMask(binary_mask, cosine_edge / pixel_size, outside_weight, pixel_size / filter_radius, pixel_size / filter_edge, outside_value, use_outside_value);
        power_spectrum_sum_image.AddImage(&power_spectrum_image); // Should I add them in Fourier or Real space?
        //my_image.Deallocate( );
        power_spectrum_image.Deallocate( );

        if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
            power_progress->Update(image_counter + 1);
    }
    // set the image is centered inside the box as true
    power_spectrum_sum_image.object_is_centred_in_box = true;
    if ( crop == true ) {
        power_spectrum_sum_image.Resize(x_dim, y_dim, 1);
        power_spectrum_sum_image.WriteSlice(&my_output_power_spectrum, 1);
        my_output_power_spectrum.my_header.SetDimensionsImage(x_dim, y_dim);
        my_output_power_spectrum.SetPixelSize(pixel_size * padding_factor); //the cropping will change the pixel size by the amount of padding
        wxPrintf("\n\n\nThe new pixel size of the cropped power spectrum is %.2f\n\n", pixel_size * padding_factor);
    }
    else {
        power_spectrum_sum_image.WriteSlice(&my_output_power_spectrum, 1);
        my_output_power_spectrum.my_header.SetDimensionsImage(x_dim * padding_factor, y_dim * padding_factor);
        my_output_power_spectrum.SetPixelSize(pixel_size);
    }

    my_output_power_spectrum.WriteHeader( );

    power_spectrum_sum_image.Deallocate( );

    return true;
}
