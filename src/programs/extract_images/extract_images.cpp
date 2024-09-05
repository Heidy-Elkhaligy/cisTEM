#include "../../core/core_headers.h"


class
        ExtractImagesFromMRC : public MyApp {

  public:
    void DoInteractiveUserInput( );
    bool DoCalculation( );

  private:
};

IMPLEMENT_APP(ExtractImagesFromMRC)

// override the DoInteractiveUserInput

void ExtractImagesFromMRC::DoInteractiveUserInput( ) {

    int first_image_in_range = 1;
    int last_image_in_range  = 1;


    UserInput* my_input = new UserInput("ExtractImagesFromMRC", 1.0);

    wxString input_stack                   = my_input->GetFilenameFromUser("Input images stack", "The input image stack in MRC format", "input_image_stack.mrc", true);
    wxString output_extracted_images_stack       = my_input->GetFilenameFromUser("Output extracted image stack", "The output extracted images in MRC format", "output_image.mrc", false);
    bool     extract_range                       = my_input->GetYesNoFromUser("Extract a specific range from the input stack", "Do you want to Extract a specific ramge of images from the input stack?", "Yes");
    if ( extract_range ) {
          first_image_in_range                   = my_input->GetIntFromUser("First image number in your extracted range", "What is the first image in your extracted range?", "1", 1);
          last_image_in_range                    = my_input->GetIntFromUser("Last image number in your extracted range", "What is the last image in your extracted range?", "2", 1);
    }


    delete my_input;

    my_current_job.Reset(5); // the number in the reset corresponds to the number of input arguments of the program
    // The first string should reflect the number of arguments, thus if arguments are 6 so string length should be 6 as well 
    my_current_job.ManualSetArguments("ttbii", 
                                      input_stack.ToUTF8( ).data( ),
                                      output_extracted_images_stack.ToUTF8( ).data( ),
                                      extract_range,
                                      first_image_in_range,
                                      last_image_in_range);
}

// override the do calculation method which will be what is actually run..

bool ExtractImagesFromMRC::DoCalculation( ) {
    wxString input_stack                         = my_current_job.arguments[0].ReturnStringArgument( );
    wxString output_extracted_images_stack       = my_current_job.arguments[1].ReturnStringArgument( );
    bool     extract_range                       = my_current_job.arguments[2].ReturnBoolArgument( );
    int      first_image_in_range                = my_current_job.arguments[3].ReturnIntegerArgument( );
    int      last_image_in_range                 = my_current_job.arguments[4].ReturnIntegerArgument( );


    MRCFile my_input_stack(input_stack.ToStdString( ), false); //MRCFile(std::string filename, bool overwrite)
    MRCFile my_output_extracted_images_stack(output_extracted_images_stack.ToStdString( ), true);

    // check if MRC files are loaded correctly by inspecting their functions
    my_input_stack.PrintInfo( );

    Image   my_extracted_image;    
    int     output_pos = 1;
    if ( extract_range == true ) {
      
      //int counter = first_image_in_range;
      int images_count = last_image_in_range - first_image_in_range;
      //my_extracted_image.Allocate(my_input_stack.ReturnXSize( ), my_input_stack.ReturnYSize( ), 1);
      for ( int counter = first_image_in_range; counter <= last_image_in_range; counter++){
        my_extracted_image.ReadSlice(&my_input_stack, counter);
        //my_extracted_image.Rotate2DInPlace(102.5, 0.0);
        //my_extracted_image.PhaseShift(-30, -45, 0.0); //Image::PhaseShift(float wanted_x_shift, float wanted_y_shift, float wanted_z_shift)
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

    }
    
    wxPrintf("\n\n");    

    return true;
}
