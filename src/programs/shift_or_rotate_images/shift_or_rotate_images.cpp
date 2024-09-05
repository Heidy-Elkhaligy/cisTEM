#include "../../core/core_headers.h"


class
        shift_or_rotate_images : public MyApp {

  public:
    void DoInteractiveUserInput( );
    bool DoCalculation( );

  private:
};

IMPLEMENT_APP(shift_or_rotate_images)

// override the DoInteractiveUserInput

void shift_or_rotate_images::DoInteractiveUserInput( ) {

    int first_image_in_range = 1;
    int last_image_in_range  = 1;
    float x_shift            = 0.0;
    float y_shift            = 0.0;
    float z_shift            = 0.0;
    float rotation_angle     = 0.0;

    UserInput* my_input = new UserInput("shift_or_rotate_images", 1.0);

    wxString input_stack                         = my_input->GetFilenameFromUser("Input images stack", "The input image stack in MRC format", "input_image_stack.mrc", true);
    wxString output_extracted_images_stack       = my_input->GetFilenameFromUser("Output extracted image stack", "The output extracted images in MRC format", "output_image.mrc", false);
    bool     extract_range                       = my_input->GetYesNoFromUser("Extract a specific range from the input stack", "Do you want to Extract a specific ramge of images from the input stack?", "Yes");
    if ( extract_range ) {
          first_image_in_range                   = my_input->GetIntFromUser("First image number in your extracted range", "What is the first image in your extracted range?", "1", 1);
          last_image_in_range                    = my_input->GetIntFromUser("Last image number in your extracted range", "What is the last image in your extracted range?", "2", 1);
    }
    bool     shift_images                        = my_input->GetYesNoFromUser("Shift the input stack", "Do you want to shift the images in the input stack?", "No");
    if (shift_images) {
        x_shift                                  = my_input->GetFloatFromUser("Wanted X-shift", "The wanted X-shift", "0.0", 0.0);
        y_shift                                  = my_input->GetFloatFromUser("Wanted Y-shift", "The wanted Y-shift", "0.0", 0.0);
        z_shift                                  = my_input->GetFloatFromUser("Wanted Z-shift", "The wanted Z-shift", "0.0", 0.0);
    }

    bool    rotate_images                        = my_input->GetYesNoFromUser("Rotate the input stack", "Do you want to rotate the images in the input stack?", "No");
    if (rotate_images) {
        rotation_angle                           = my_input->GetFloatFromUser("Wanted rotation angle", "The wanted rotation angle", "0.0", 0.0);
    }


    delete my_input;

    my_current_job.Reset(11); // the number in the reset corresponds to the number of input arguments of the program
    // The first string should reflect the number of arguments, thus if arguments are 6 so string length should be 6 as well 
    my_current_job.ManualSetArguments("ttbiibfffbf", 
                                      input_stack.ToUTF8( ).data( ),
                                      output_extracted_images_stack.ToUTF8( ).data( ),
                                      extract_range,
                                      first_image_in_range,
                                      last_image_in_range,
                                      shift_images,
                                      x_shift,
                                      y_shift,
                                      z_shift,
                                      rotate_images,
                                      rotation_angle);
}

// override the do calculation method which will be what is actually run..

bool shift_or_rotate_images::DoCalculation( ) {
    wxString input_stack                         = my_current_job.arguments[0].ReturnStringArgument( );
    wxString output_extracted_images_stack       = my_current_job.arguments[1].ReturnStringArgument( );
    bool     extract_range                       = my_current_job.arguments[2].ReturnBoolArgument( );
    int      first_image_in_range                = my_current_job.arguments[3].ReturnIntegerArgument( );
    int      last_image_in_range                 = my_current_job.arguments[4].ReturnIntegerArgument( );
    bool     shift_images                        = my_current_job.arguments[5].ReturnBoolArgument( );
    float    x_shift                             = my_current_job.arguments[6].ReturnFloatArgument( );
    float    y_shift                             = my_current_job.arguments[7].ReturnFloatArgument( );
    float    z_shift                             = my_current_job.arguments[8].ReturnFloatArgument( );
    bool     rotate_images                       = my_current_job.arguments[9].ReturnBoolArgument( );
    float    rotation_angle                      = my_current_job.arguments[10].ReturnFloatArgument( );


    MRCFile my_input_stack(input_stack.ToStdString( ), false); //MRCFile(std::string filename, bool overwrite)
    MRCFile my_output_extracted_images_stack(output_extracted_images_stack.ToStdString( ), true);

    // check if MRC files are loaded correctly by inspecting their functions
    my_input_stack.PrintInfo( );

    Image   my_extracted_image;    
    int     output_pos = 1;
    if ( extract_range == true ) { // if extract and rotate/shift
      //int counter = first_image_in_range;
      int images_count = last_image_in_range - first_image_in_range;
      //my_extracted_image.Allocate(my_input_stack.ReturnXSize( ), my_input_stack.ReturnYSize( ), 1);
      for ( int counter = first_image_in_range; counter <= last_image_in_range; counter++){

        my_extracted_image.ReadSlice(&my_input_stack, counter);

        if (shift_images == true or rotate_images == true) {
            my_extracted_image.Rotate2DInPlace(rotation_angle, 0.0);
            my_extracted_image.PhaseShift(x_shift, y_shift, z_shift); //Image::PhaseShift(float wanted_x_shift, float wanted_y_shift, float wanted_z_shift)
        }
        //my_extracted_image.Rotate2DInPlace(102.5, 0.0);
        //my_extracted_image.PhaseShift(16.0, -8.5, 0.0); //Image::PhaseShift(float wanted_x_shift, float wanted_y_shift, float wanted_z_shift)
        //my_extracted_image.QuickAndDirtyWriteSlice("my_extracted_images_stack.mrc", counter);
        my_extracted_image.WriteSlice(&my_output_extracted_images_stack, output_pos);
        output_pos++;
      }  
      //my_extracted_image.WriteSlices(&my_output_extracted_images_stack, first_image_in_range, last_image_in_range);
      // Add header information to the extracted images
      wxPrintf("The number of slices should be %i \n", images_count+1);
      my_output_extracted_images_stack.my_header.SetDimensionsVolume(my_input_stack.ReturnXSize(), my_input_stack.ReturnYSize(), images_count+1);
      my_output_extracted_images_stack.my_header.SetPixelSize(my_input_stack.ReturnPixelSize());
      my_output_extracted_images_stack.WriteHeader( );

    } else { // if only rotate or shift
        // read the images
        for ( int counter = 0; counter < my_input_stack.ReturnNumberOfSlices( ); counter++){
          my_extracted_image.ReadSlice(&my_input_stack, counter);
          my_extracted_image.Rotate2DInPlace(rotation_angle, 0.0);
          my_extracted_image.PhaseShift(x_shift, y_shift, z_shift); //Image::PhaseShift(float wanted_x_shift, float wanted_y_shift, float wanted_z_shift)
          my_extracted_image.WriteSlice(&my_output_extracted_images_stack, output_pos);
          output_pos++;

          //wxPrintf("The number of slices should be %i \n", images_count+1);
          my_output_extracted_images_stack.my_header.SetDimensionsVolume(my_input_stack.ReturnXSize(), my_input_stack.ReturnYSize(), my_input_stack.ReturnNumberOfSlices( ));
          my_output_extracted_images_stack.my_header.SetPixelSize(my_input_stack.ReturnPixelSize());
          my_output_extracted_images_stack.WriteHeader( );

        }

    }
    
    wxPrintf("\n\n");    

    return true;
}
