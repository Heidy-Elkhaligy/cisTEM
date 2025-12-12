#include "../../core/core_headers.h"
#include <iostream>
#include <fstream>

class
        calculate_autocorrelation : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

IMPLEMENT_APP(calculate_autocorrelation)

// override the DoInteractiveUserInput

void calculate_autocorrelation::DoInteractiveUserInput( ) {
    wxString input_images;
    wxString output_images;
    float    pixel_size;
    int      max_threads;

    UserInput* my_input = new UserInput("calculate_autocorrelation", 1.00);
    input_images        = my_input->GetFilenameFromUser("Input images file name", "Filename of classaverage stack", "input_stack.mrc", true);
    output_images       = my_input->GetFilenameFromUser("Output images", "The auto correlation images", "auto_correlation_stack.mrc", false);
    pixel_size          = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);

#ifdef _OPENMP
    max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "when threading, what is the max threads to run", "1", 1);
#else
    max_threads = 1;
#endif

    delete my_input;

    //	my_current_job.Reset(3);
    my_current_job.ManualSetArguments("ttfi", input_images.ToUTF8( ).data( ), output_images.ToUTF8( ).data( ), pixel_size, max_threads);
}

// override the do calculation method which will be what is actually run..

bool calculate_autocorrelation::DoCalculation( ) {
    wxString input_images  = my_current_job.arguments[0].ReturnStringArgument( );
    wxString output_images = my_current_job.arguments[1].ReturnStringArgument( );
    float    pixel_size    = my_current_job.arguments[2].ReturnFloatArgument( );
    int      max_threads   = my_current_job.arguments[3].ReturnIntegerArgument( );

    MRCFile my_input_filename(input_images.ToStdString( ), false);
    MRCFile my_output_filename(output_images.ToStdString( ), true);

    long number_of_input_images = my_input_filename.ReturnNumberOfSlices( );

    // set the header of the output MRC file to avoid OMP problems
    my_output_filename.my_header.SetNumberOfImages(number_of_input_images);
    my_output_filename.my_header.SetDimensionsImage(my_input_filename.ReturnXSize( ), my_input_filename.ReturnYSize( ));
    my_output_filename.SetPixelSize(pixel_size);
    my_output_filename.WriteHeader( );
    my_output_filename.rewrite_header_on_close = true;

    Image current_image;
    long  image_counter;

    wxPrintf("\nCalculating auto correlation images...\n\n");
    ProgressBar* my_progress = new ProgressBar(number_of_input_images);

#pragma omp parallel for schedule(dynamic, 1) num_threads(max_threads) default(none) shared(my_input_filename, my_output_filename, number_of_input_images, my_progress) private(image_counter, current_image)
    for ( image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
//wxPrintf("The image counter for slice reading is %li \n", image_counter + 1);
// read the current image in the stack
// #pragma omp ordered
//         {
#pragma omp critical
        current_image.ReadSlice(&my_input_filename, image_counter + 1);
        // }
        // Normalize the image using cisTEM Normalize
        current_image.Normalize( );
        // FT the image
        current_image.ForwardFFT( );
        // convert the central pixel to zero (Is that done in real or Fouriier space??)
        current_image.ZeroCentralPixel( );

        for ( long pixel_counter = 0; pixel_counter < current_image.real_memory_allocated / 2; pixel_counter++ ) {

            // calculating the amplitude is not needed by let's see and print its value
            float amplitude = abs(current_image.complex_values[pixel_counter]);
            //wxPrintf("The image number %li pixel counter is %li \n", image_counter+1, pixel_counter);
            //wxPrintf("Retunring the amplitude from the complex values %f \n", amplitude);
            //wxPrintf("Retunring the complex values before changing them are %f, %f \n", current_image.complex_values[pixel_counter].real(), current_image.complex_values[pixel_counter].imag());

            //As the phase will be zero, the real is just the amplitude and the imaginary is 0
            // so we will set the complex number to be equal to apmlitude + 0
            current_image.complex_values[pixel_counter] = amplitude * amplitude + I * 0.0f;
        }
        // return the image to real space again and save them to see the correlation
        current_image.BackwardFFT( );
        current_image.SwapRealSpaceQuadrants( );
        // set the image is centered inside the box as true
        current_image.object_is_centred_in_box = true;
// #pragma omp ordered
//         {
#pragma omp critical
        current_image.WriteSlice(&my_output_filename, image_counter + 1);
        // }
        if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
            my_progress->Update(image_counter + 1);
    }
    delete my_progress;
    return true;
}
