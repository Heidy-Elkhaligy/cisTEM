#include "../../core/core_headers.h"
#include <fstream>
#include <wx/arrimpl.cpp>
#include <omp.h>
class
        align_classaverages_to_reference : public MyApp {

  public:
    void DoInteractiveUserInput( );
    bool DoCalculation( );

  private:
};


IMPLEMENT_APP(align_classaverages_to_reference)

// override the DoInteractiveUserInput

void align_classaverages_to_reference::DoInteractiveUserInput( ) {

    float rotation_angle        ;
    //float fine_tuning_step_size ;
    int   max_threads           ;


    UserInput* my_input = new UserInput("align_classaverages_to_reference", 1.0);

    wxString input_classaverage_stack    = my_input->GetFilenameFromUser("Input class average stack filename", "The input image stack, containing the 2D class averages images", "my_classaverage_stack.mrc", true);
    wxString input_reference_stack       = my_input->GetFilenameFromUser("Input reference stack filename", "The input image stack, containing your reference images", "my_reference_images.mrc", true);
    wxString input_star_filename         = my_input->GetFilenameFromUser("Input reference projections star filename", "The input parameter file, containing your reference projections parameters", "my_reference_parameters.star", true);
    rotation_angle                       = my_input->GetFloatFromUser("Rotation angle", "The rotational angle applied to the 2D class average image during the initial alignment to reference images", "4.0", 4.0);
    //fine_tuning_step_size                = my_input->GetFloatFromUser("Fine tuning step size", "The fine tuning step size for the last alignment", "0.1", 0.01 , 0.1);
    wxString output_aligned_classaverage = my_input->GetFilenameFromUser("Output aligned class average stack filename", "The aligned class average result", "class_average_output.mrc", false);
    wxString output_reference_stack      = my_input->GetFilenameFromUser("Output matching reference images stack filename", "The matching reference images result", "matching_reference_output.mrc", false);
    wxString output_star_file            = my_input->GetFilenameFromUser("Output matching refernce projections Star File", "The star file containing angles of the matching filtered reference projections", "matching_reference_angles.star", false);
    wxString output_alignment_info       = my_input->GetFilenameFromUser("Output alignment rotation and translation filename", "The output rotation and translation information for the class averages","alignment_info_output.txt", false);

#if defined(_OPENMP)
    max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "When threading, what is the max threads to run", "1", 1);
#else
    max_threads = 1;
#endif
    delete my_input;

    my_current_job.Reset(10); 
    my_current_job.ManualSetArguments("tttftttti", 
                                      input_classaverage_stack.ToUTF8( ).data( ),
                                      input_reference_stack.ToUTF8( ).data( ),
                                      input_star_filename.ToUTF8( ).data( ),
                                      rotation_angle,
                                      //fine_tuning_step_size,
                                      output_aligned_classaverage.ToUTF8( ).data( ),
                                      output_reference_stack.ToUTF8( ).data( ),
                                      output_star_file.ToUTF8( ).data( ),
                                      output_alignment_info.ToUTF8( ).data( ), 
                                      max_threads);
}

// override the do calculation method which will be what is actually run..

bool align_classaverages_to_reference::DoCalculation( ) {
    wxString input_classaverage_stack    = my_current_job.arguments[0].ReturnStringArgument( );
    wxString input_reference_stack       = my_current_job.arguments[1].ReturnStringArgument( );
    wxString input_star_filename         = my_current_job.arguments[2].ReturnStringArgument( );
    float    rotation_angle              = my_current_job.arguments[3].ReturnFloatArgument( );
    //float    fine_tuning_step_size       = my_current_job.arguments[4].ReturnFloatArgument( );
    wxString output_aligned_classaverage = my_current_job.arguments[4].ReturnStringArgument( );
    wxString output_reference_stack      = my_current_job.arguments[5].ReturnStringArgument( );
    wxString output_star_file            = my_current_job.arguments[6].ReturnStringArgument( );
    wxString output_alignment_info       = my_current_job.arguments[7].ReturnStringArgument( );
    int      max_threads                 = my_current_job.arguments[8].ReturnIntegerArgument( );



    MRCFile my_input_classaverage_stack(input_classaverage_stack.ToStdString( ), false); 
    MRCFile my_input_reference_stack(input_reference_stack.ToStdString( ), false);
    MRCFile my_aligned_classaverage(output_aligned_classaverage.ToStdString( ), true);
    MRCFile my_output_reference(output_reference_stack.ToStdString( ), true);



    //Check that two files have same size
    if ( my_input_classaverage_stack.ReturnXSize( ) != my_input_reference_stack.ReturnXSize( ) || my_input_classaverage_stack.ReturnYSize( ) != my_input_reference_stack.ReturnYSize( ) ) {
        MyPrintfRed("\n\n Error: Image dimensions are not the same\n\n");
        exit(-1);
    }

 
    Image   my_classaverage_image;
    Image   my_reference_image;
    int     number_of_classaverage_images = my_input_classaverage_stack.ReturnNumberOfSlices( );
    //wxPrintf("Number of class average input are %i \n", number_of_classaverage_images);
    int     number_of_reference_images = my_input_reference_stack.ReturnNumberOfSlices( );
    //wxPrintf("Number of reference input are %i \n", number_of_reference_images);

    float   best_correlation_score[number_of_classaverage_images] = {-FLT_MAX}  ; //= {-FLT_MAX}   
    float   best_psi_value[number_of_classaverage_images] = {0};
    float   best_x_shift_value[number_of_classaverage_images] = {0};
    float   best_y_shift_value[number_of_classaverage_images] = {0};
    int     best_reference_image_number[number_of_classaverage_images] = {0};

    float   inner_radius_for_peak_search;
    float   outer_radius_for_peak_search;
    inner_radius_for_peak_search = 0.0; // inner radius should be set to 0
    outer_radius_for_peak_search = my_input_reference_stack.ReturnXSize( )/2; //outer radius is better to be half box size


    float   correlation_score[number_of_classaverage_images] = {-FLT_MAX} ; //= {-FLT_MAX}   
    float   psi_value[number_of_classaverage_images] = {0};
    float   x_shift_value[number_of_classaverage_images] = {0};
    float   y_shift_value[number_of_classaverage_images] = {0};
    int     reference_image_number[number_of_classaverage_images] = {0};

    wxPrintf("\nAligning Images...\n\n");
    ProgressBar* my_progress = new ProgressBar(number_of_classaverage_images); 

#pragma omp parallel for ordered schedule(dynamic) num_threads(max_threads) default(none) shared(my_input_classaverage_stack, my_input_reference_stack , my_progress, inner_radius_for_peak_search, outer_radius_for_peak_search, \
                                                                   best_correlation_score, best_psi_value, best_x_shift_value, best_y_shift_value, best_reference_image_number, rotation_angle, reference_image_number,\
                                                                   correlation_score, psi_value, x_shift_value, y_shift_value , number_of_classaverage_images, number_of_reference_images, max_threads) \
                                                                   private(my_classaverage_image, my_reference_image)

    for ( long classaverage_image_counter = 0; classaverage_image_counter < number_of_classaverage_images; classaverage_image_counter++ ) {
	    #pragma omp critical
        my_classaverage_image.ReadSlice(&my_input_classaverage_stack, classaverage_image_counter + 1); 
        my_classaverage_image.Normalize( );

        Image   rotated_image;
        Image   reference_image;

        for ( long reference_image_counter = 0; reference_image_counter < number_of_reference_images; reference_image_counter++ ) {

            // read the reference image 
	        #pragma omp critical
            my_reference_image.ReadSlice(&my_input_reference_stack, reference_image_counter + 1);
            //normalize the image
            my_reference_image.Normalize( );

            //1. copy and rotate the class_average image based on the rotation angle specified 
            //2. calculate the cross correlation and find the peak that represent the highest similarity between the two images
            for (float psi = 0.0; psi < 360.0; psi+=rotation_angle) {
                // create a new peak to save the cross-correlation peak values
                Peak    current_peak;

                rotated_image.CopyFrom(&my_classaverage_image);

                // rotate by the class average image by the rotation angle
                rotated_image.Rotate2DInPlace(psi,0.0);

                reference_image.CopyFrom(&my_reference_image);
                // calculate the cross correlation of the reference image with the rotated image
                rotated_image.CalculateCrossCorrelationImageWith(&reference_image);

                // find the peak from the cross corrlation to get the values 
                current_peak = rotated_image.FindPeakWithParabolaFit(inner_radius_for_peak_search, outer_radius_for_peak_search);

                // check if the current cross-correlation score is better than stored or not

                if ( current_peak.value > correlation_score[classaverage_image_counter] ) {
                    correlation_score[classaverage_image_counter]  =   current_peak.value;
                    psi_value[classaverage_image_counter] = psi; 
                    x_shift_value[classaverage_image_counter] = current_peak.x;
                    y_shift_value[classaverage_image_counter] = current_peak.y;
                    reference_image_number[classaverage_image_counter] = reference_image_counter + 1;
                    //wxPrintf("Inside the loop The reference image counter number is %li \n", reference_image_counter + 1);
                    //wxPrintf("Inside the loop The reference image counter number saved in the array is %i \n", reference_image_number[classaverage_image_counter]);

                } 

            } 

        }
        best_correlation_score[classaverage_image_counter]  =   correlation_score[classaverage_image_counter];
        best_psi_value[classaverage_image_counter] = psi_value[classaverage_image_counter]; 
        best_x_shift_value[classaverage_image_counter] = x_shift_value[classaverage_image_counter];
        best_y_shift_value[classaverage_image_counter] = y_shift_value[classaverage_image_counter];
        best_reference_image_number[classaverage_image_counter] = reference_image_number[classaverage_image_counter];
        
        if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )
        	my_progress->Update(classaverage_image_counter + 1);


    }

    delete my_progress;

    // start tuning the image alignment by using +/- 2 degrees (half the specified rotational angle) 
    Image my_tuned_classaverage_image;
    Image my_rotated_tuned_classaverage;
    Image my_similar_reference_image;
    Image temp_similar_reference;
    float tuned_rotatation_range = rotation_angle / 2;
    float tuned_step_size = tuned_rotatation_range / 10;
    wxPrintf("\nTuning the alignment...\n\n");
    ProgressBar* my_tuning_progress = new ProgressBar(number_of_classaverage_images); 

#pragma omp parallel for ordered schedule(dynamic) num_threads(max_threads) default(none) shared(my_input_classaverage_stack, my_input_reference_stack , my_tuning_progress, inner_radius_for_peak_search, outer_radius_for_peak_search, \
                                                                   best_correlation_score, best_psi_value, best_x_shift_value, best_y_shift_value, best_reference_image_number, rotation_angle, \
                                                                   number_of_classaverage_images, tuned_rotatation_range, max_threads, tuned_step_size ) \
                                                                   private(my_tuned_classaverage_image, my_rotated_tuned_classaverage, my_similar_reference_image, temp_similar_reference)


     for ( long image_counter = 0; image_counter < number_of_classaverage_images; image_counter++ ) {
	    #pragma omp critical
        // read the unrotated and untranslated class average image (initial image)
        my_tuned_classaverage_image.ReadSlice(&my_input_classaverage_stack, image_counter + 1); 
        my_tuned_classaverage_image.Normalize( ); 

        // read its similar reference
        int reference_image_number = best_reference_image_number[image_counter];
        //wxPrintf("The tuning reference image number is %i \n", reference_image_number);
        my_similar_reference_image.ReadSlice(&my_input_reference_stack, reference_image_number);
        my_similar_reference_image.Normalize( );

        float current_best_psi = best_psi_value[image_counter];
        float tuned_psi_lower_range = current_best_psi - tuned_rotatation_range;
        float tuned_psi_upper_range = current_best_psi + tuned_rotatation_range;

        
        // loop over the range of +/- half the rotation angle 
        // increment by 1/10 of the tuned rotation angle (rotation angle/2)/10 degrees for tuning 
        for (float tuned_psi = tuned_psi_lower_range; tuned_psi <= tuned_psi_upper_range; tuned_psi+=tuned_step_size) {
            // create a new peak to save the tuned values
            Peak    current_tuned_peak;
            // copy the class average image to tune the rotation angle later
            my_rotated_tuned_classaverage.CopyFrom(&my_tuned_classaverage_image);
            // rotate by the class average image by the rotation angle 
            my_rotated_tuned_classaverage.Rotate2DInPlace(tuned_psi, 0.0); 
            temp_similar_reference.CopyFrom(&my_similar_reference_image);
            // calculate the cross correlation of the reference image with the rotated image
            my_rotated_tuned_classaverage.CalculateCrossCorrelationImageWith(&temp_similar_reference);

            // find the peak from the cross corrlation to get the values 
            current_tuned_peak = my_rotated_tuned_classaverage.FindPeakWithParabolaFit(inner_radius_for_peak_search, outer_radius_for_peak_search);

            // debug scores
            //wxPrintf("returning peak values anyway, psi = %f, best score = %lf, current score = %f\n", tuned_psi, best_correlation_score[image_counter], current_tuned_peak.value);

            // check if the current cross-correlation score using the tuned angle is better than stored or not
            if ( current_tuned_peak.value > best_correlation_score[image_counter] ) {
                best_correlation_score[image_counter]  =   current_tuned_peak.value;
                best_psi_value[image_counter] =  tuned_psi; 
                best_x_shift_value[image_counter] = current_tuned_peak.x;
                best_y_shift_value[image_counter] = current_tuned_peak.y;
               
                // debug scores
                //wxPrintf("\n\nInside rotation loop returning, psi = %f, best score = %f, current score = %f\n\n",  tuned_psi, best_correlation_score[image_counter], current_tuned_peak.value);
            } 


        } 
        my_tuning_progress->Update(image_counter + 1);
    }

    delete my_tuning_progress;

    //start of fine tuning the alignment of images
    Image my_fine_tuned_classaverage_image;
    Image my_rotated_fine_tuned_classaverage;
    Image my_matched_reference_image;
    Image temp_matched_reference;
    float fine_tuning_rotatation_range = rotation_angle / 4;
    float fine_tuning_step_size = fine_tuning_rotatation_range / 10;

    wxPrintf("\nFine tuning the alignment...\n\n");
    ProgressBar* my_fine_tuning_progress = new ProgressBar(number_of_classaverage_images);

    // start tuning the image alignment by using +/- 2 degrees (half the specified rotational angle) 
#pragma omp parallel for ordered schedule(dynamic) num_threads(max_threads) default(none) shared(my_input_classaverage_stack, my_input_reference_stack , my_fine_tuning_progress, inner_radius_for_peak_search, outer_radius_for_peak_search, \
                                                                   best_correlation_score, best_psi_value, best_x_shift_value, best_y_shift_value, best_reference_image_number, rotation_angle, \
                                                                   number_of_classaverage_images, fine_tuning_step_size, max_threads, fine_tuning_rotatation_range ) \
                                                                   private(my_fine_tuned_classaverage_image, my_rotated_fine_tuned_classaverage, my_matched_reference_image, temp_matched_reference)
     for ( long fine_tune_counter = 0; fine_tune_counter < number_of_classaverage_images; fine_tune_counter++ ) {
	    #pragma omp critical
        // read the unrotated and untranslated class average image (initial image)
        my_fine_tuned_classaverage_image.ReadSlice(&my_input_classaverage_stack, fine_tune_counter + 1); 
        my_fine_tuned_classaverage_image.Normalize( ); 
        // read its similar reference
        int reference_image_number = best_reference_image_number[fine_tune_counter];
        my_matched_reference_image.ReadSlice(&my_input_reference_stack, reference_image_number);
        my_matched_reference_image.Normalize( );
         
        float current_best_psi = best_psi_value[fine_tune_counter];
        float fine_tuned_psi_lower_range = current_best_psi - fine_tuning_rotatation_range;
        float fine_tuned_psi_upper_range = current_best_psi + fine_tuning_rotatation_range;

        
        // loop over the +/- 1 degrees around the rotation angle 
        // increment by fine tuning step size (0.01-0.1)
        for (float fine_tuned_psi = fine_tuned_psi_lower_range; fine_tuned_psi <= fine_tuned_psi_upper_range; fine_tuned_psi+=fine_tuning_step_size) {
            // create a new peak to save the tuned values
            Peak    current_fine_tuned_peak;
            // copy the class average image to fine tune the rotation angle later
            my_rotated_fine_tuned_classaverage.CopyFrom(&my_fine_tuned_classaverage_image);

            // rotate by the class average image by the rotation angle 
            my_rotated_fine_tuned_classaverage.Rotate2DInPlace(fine_tuned_psi, 0.0); // 

            temp_matched_reference.CopyFrom(&my_matched_reference_image);
            // calculate the cross correlation of the reference image with the rotated image
            my_rotated_fine_tuned_classaverage.CalculateCrossCorrelationImageWith(&temp_matched_reference);

            // find the peak from the cross corrlation to get the values 
            current_fine_tuned_peak = my_rotated_fine_tuned_classaverage.FindPeakWithParabolaFit(inner_radius_for_peak_search, outer_radius_for_peak_search);

            // debug scores
            //wxPrintf("returning peak values anyway, fine tuned psi = %f, best score = %lf, current score = %f\n", fine_tuned_psi, best_correlation_score[fine_tune_counter], current_fine_tuned_peak.value);

            // check if the current cross-correlation score using the tuned angle is better than stored or not
            if ( current_fine_tuned_peak.value > best_correlation_score[fine_tune_counter] ) {
                best_correlation_score[fine_tune_counter]  =   current_fine_tuned_peak.value;
                best_psi_value[fine_tune_counter] =  fine_tuned_psi; 
                best_x_shift_value[fine_tune_counter] = current_fine_tuned_peak.x;
                best_y_shift_value[fine_tune_counter] = current_fine_tuned_peak.y;
               
                // debug scores
                //wxPrintf("\n\nInside rotation loop returning, fine tuned psi = %f, best score = %f, current score = %f\n\n",  fine_tuned_psi, best_correlation_score[fine_tune_counter], current_fine_tuned_peak.value);
            } 
        } 
        my_fine_tuning_progress->Update(fine_tune_counter + 1);
    }
    delete my_fine_tuning_progress;

    
    Image my_temp_classaverage_image;    
    Image my_temp_reference_image;
    cisTEMParameters input_star_file;
    cisTEMParameters output_params;
    
    // writing out the rotated and translated class averages and their corresponding reference image
    wxPrintf("\nWriting Images...\n\n");
    ProgressBar* my_writing_progress = new ProgressBar(number_of_classaverage_images); 

    // read the projections star file
    input_star_file.ReadFromcisTEMStarFile(input_star_filename);
    cisTEMParameterLine input_parameters;

    // setup parameters for the output star file
    output_params.parameters_to_write.SetActiveParameters(POSITION_IN_STACK | IMAGE_IS_ACTIVE | PSI | THETA | PHI | X_SHIFT | Y_SHIFT | DEFOCUS_1 | DEFOCUS_2 | DEFOCUS_ANGLE | PHASE_SHIFT | OCCUPANCY | LOGP | SIGMA | SCORE | PIXEL_SIZE | MICROSCOPE_VOLTAGE | MICROSCOPE_CS | AMPLITUDE_CONTRAST | BEAM_TILT_X | BEAM_TILT_Y | IMAGE_SHIFT_X | IMAGE_SHIFT_Y | ASSIGNED_SUBSET);
    output_params.PreallocateMemoryAndBlank(number_of_classaverage_images);

    // write the final class average rotation angles and translation values needed to match the matched reference image to output file
    std::ofstream outputFile(output_alignment_info.ToStdString());
    // write a header line to the output file
    outputFile << "position,psi,xshift,yshift,crosscorr_score" << std::endl;

#pragma omp for ordered schedule(static, 1) 

    for ( long write_counter = 0; write_counter < my_input_classaverage_stack.ReturnNumberOfSlices( ); write_counter++ ) {
	        #pragma omp critical
            // writing the rotated and translated class average images
            my_temp_classaverage_image.ReadSlice(&my_input_classaverage_stack, write_counter + 1); 
            my_temp_classaverage_image.Rotate2DInPlace(best_psi_value[write_counter], 0.0);
            my_temp_classaverage_image.PhaseShift(-best_x_shift_value[write_counter], -best_y_shift_value[write_counter], 0.0); 
            my_temp_classaverage_image.WriteSlice(&my_aligned_classaverage, write_counter + 1);
            
            outputFile << write_counter << "," << best_psi_value[write_counter] << "," << -best_x_shift_value[write_counter] << "," << -best_y_shift_value[write_counter] << "," << best_correlation_score[write_counter] << std::endl;
            // writing the matching reference images
            int reference_image_number = best_reference_image_number[write_counter];
	        #pragma omp critical
            my_temp_reference_image.ReadSlice(&my_input_reference_stack, reference_image_number);
            my_temp_reference_image.WriteSlice(&my_output_reference, write_counter + 1);
            //wxPrintf("reference image number is %i \n", reference_image_number);
            input_parameters = input_star_file.ReturnLine(reference_image_number - 1); // the star file numbering is 0 indexed!!
            // write_counter here should start at 0
            output_params.all_parameters[write_counter].position_in_stack                  = write_counter +1; 
            output_params.all_parameters[write_counter].image_is_active                    = input_parameters.image_is_active;
            output_params.all_parameters[write_counter].psi                                = input_parameters.psi;
            output_params.all_parameters[write_counter].theta                              = input_parameters.theta;
            output_params.all_parameters[write_counter].phi                                = input_parameters.phi;
            output_params.all_parameters[write_counter].x_shift                            = input_parameters.x_shift;
            output_params.all_parameters[write_counter].y_shift                            = input_parameters.y_shift;
            output_params.all_parameters[write_counter].defocus_1                          = input_parameters.defocus_1;
            output_params.all_parameters[write_counter].defocus_2                          = input_parameters.defocus_2;
            output_params.all_parameters[write_counter].defocus_angle                      = input_parameters.defocus_angle;
            output_params.all_parameters[write_counter].phase_shift                        = input_parameters.phase_shift;
            output_params.all_parameters[write_counter].occupancy                          = input_parameters.occupancy;
            output_params.all_parameters[write_counter].logp                               = input_parameters.logp;
            output_params.all_parameters[write_counter].sigma                              = input_parameters.sigma;
            output_params.all_parameters[write_counter].score                              = input_parameters.score;
            output_params.all_parameters[write_counter].score_change                       = input_parameters.score_change;
            output_params.all_parameters[write_counter].pixel_size                         = input_parameters.pixel_size;
            output_params.all_parameters[write_counter].microscope_voltage_kv              = input_parameters.microscope_voltage_kv;
            output_params.all_parameters[write_counter].microscope_spherical_aberration_mm = input_parameters.microscope_spherical_aberration_mm;
            output_params.all_parameters[write_counter].amplitude_contrast                 = input_parameters.amplitude_contrast;
            output_params.all_parameters[write_counter].beam_tilt_x                        = input_parameters.beam_tilt_x;
            output_params.all_parameters[write_counter].beam_tilt_y                        = input_parameters.beam_tilt_y;
            output_params.all_parameters[write_counter].image_shift_x                      = input_parameters.image_shift_x;
            output_params.all_parameters[write_counter].image_shift_y                      = input_parameters.image_shift_y;

            my_writing_progress->Update(write_counter + 1);
        }

    outputFile.close();
    // write the output star file for the matched references
    output_params.WriteTocisTEMStarFile(output_star_file);
    delete my_writing_progress;
    wxPrintf("\n\n");    

    return true;
}

