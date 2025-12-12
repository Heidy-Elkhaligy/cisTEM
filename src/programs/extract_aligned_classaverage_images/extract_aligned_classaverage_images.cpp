#include "../../core/core_headers.h"

class
        extract_aligned_classaverage_images : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

// // Made to ensure proper association of the class_average_index and its members
// typedef struct ClassAveragesAndMembers {
//     int         class_average_index;
//     long        number_of_members;
//     long        member_counter = 0;
//     wxArrayLong class_members;
// } ClassAveragesAndMembers;

IMPLEMENT_APP(extract_aligned_classaverage_images)

// override the DoInteractiveUserInput

void extract_aligned_classaverage_images::DoInteractiveUserInput( ) {
    int        max_threads;
    UserInput* my_input = new UserInput("extract_aligned_classaverage_images", 1.00);

    std::string database_filename       = my_input->GetFilenameFromUser("Database filename that contains relevant classification", "cisTEM .db file containing the class average that should be extracted and aligned", "input_database.db", true);
    int         classification_id       = my_input->GetIntFromUser("Input the classification ID that contains the class average for extraction and alignment", "Classification ID of class average run being used for extraction and alignment", "0");
    int         class_number            = my_input->GetIntFromUser("Input class that will be extracted and aligned", "Class ID of the class average that will be extracted and aligned", "1");
    std::string particle_stack_filename = my_input->GetFilenameFromUser("Input particle stack", "The filename for the relevant .mrc file", "input.mrc", true);
    std::string output_filename         = my_input->GetFilenameFromUser("Output filename of the extracted aligned particles", "Filename to save the output aligned images of the given class.", "extracted_aligned_class_images.mrc", false);

    // #ifdef _OPENMP
    //     max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "when threading, what is the max threads to run", "1", 1);
    // #else
    //     max_threads = 1;
    // #endif

    delete my_input;

    my_current_job.Reset(7);
    my_current_job.ManualSetArguments("tiitt", database_filename.c_str( ), classification_id, class_number, particle_stack_filename.c_str( ), output_filename.c_str( )); //update_star_file, input_star_filename.ToUTF8( ).data( ),, max_threads
}

// override the do calculation method which will be what is actually run..

bool extract_aligned_classaverage_images::DoCalculation( ) {
    wxFileName  database_filename       = wxFileName(my_current_job.arguments[0].ReturnStringArgument( ));
    const int   classification_id       = my_current_job.arguments[1].ReturnIntegerArgument( );
    const int   class_number            = my_current_job.arguments[2].ReturnIntegerArgument( );
    std::string particle_stack_filename = my_current_job.arguments[3].ReturnStringArgument( );
    std::string output_filename         = my_current_job.arguments[4].ReturnStringArgument( );

    //int         max_threads             = my_current_job.arguments[4].ReturnIntegerArgument( );

    // Open db and prepare the refinement packages
    Database selected_db = Database( );
    selected_db.Open(database_filename);
    // // should return how many images are in this class
    // int num_of_images = selected_db.ReturnNumberOf2DClassMembers(classification_id, class_number);
    // // I need  to find a way to get the number of the classaverages for that specific classification id using the database
    // wxPrintf("\nNumber of images in class number %i is %i...\n\n", class_number, num_of_images);
    // Load refinement packages for writing .star
    ///ArrayOfRefinementPackages refinement_package_list;
    //RefinementPackage* temp_package;

    //selected_db.BeginAllRefinementPackagesSelect( ); // do I need this?

    Classification* needed_class = selected_db.GetClassificationByID(classification_id); //Retrieves all metadata and per-particle classification results for a given classification ID from the databas
    // RefinementPackage needed_package          = needed_class->refinement_package_asset_id - 1; // Select refinement package using classification we made
    wxArrayLong input_class_members = selected_db.Return2DClassMembers(classification_id, class_number); // retrieves a list of particle positions (indices) that belong to a specific 2D class in a specific classification from the database.
    selected_db.Close(false);
    long number_of_class_members = input_class_members.GetCount( );
    // wxPrintf("\nNumber of images in class number %i is %li...\n\n", class_number, number_of_class_members);
    // for ( unsigned int i = 0; i < number_of_class_members; i++ ) {
    //     wxPrintf("Member %u: %ld\n", i, input_class_members[i]);
    // }
    MRCFile my_particle_stack(particle_stack_filename.c_str( ), false); // no overwrite only reading
    MRCFile my_output_filename(output_filename.c_str( ), true); //overwrite

    Image        class_image;
    ProgressBar* my_progress = new ProgressBar(number_of_class_members);
    // EXTRACT THE CLASS MEMBERS
    for ( long image_counter = 0; image_counter < number_of_class_members; image_counter++ ) {
        // long class_member_id = input_class_members[image_counter];
        // class_image.ReadSlice(&my_particle_stack, class_member_id); //As the counter of the array starts from 0 and readslice starts from 1?? -1
        long raw_id = input_class_members[image_counter]; // 1‑based
        long idx    = raw_id - 1; // 0‑based

        class_image.ReadSlice(&my_particle_stack, raw_id);

        float psi     = needed_class->classification_results[idx].psi;
        float x_shift = needed_class->classification_results[idx].xshift;
        float y_shift = needed_class->classification_results[idx].yshift;
        //float pixel_size = needed_class->classification_results[idx].pixel_size;

        class_image.QuickAndDirtyWriteSlice("original_unaligned_class_images.mrc", image_counter + 1);

        // float psi     = needed_class->classification_results[class_member_id].psi;
        // float x_shift = needed_class->classification_results[class_member_id].xshift;
        // float y_shift = needed_class->classification_results[class_member_id].yshift;
        wxPrintf("  Member Index %li -> Particle ID: %li\n", image_counter, input_class_members[image_counter]);
        wxPrintf("Image %li rotation xshift and yshift are %f, %f, %f \n", image_counter, -psi, -x_shift, -y_shift);
        // trying rotating by psi not -psi to align particles correctly
        class_image.Rotate2DInPlace(-psi, FLT_MAX);
        class_image.QuickAndDirtyWriteSlice("rotated_unaligned_images.mrc", image_counter + 1);
        class_image.PhaseShift(x_shift, y_shift);
        class_image.WriteSlice(&my_output_filename, image_counter + 1);
        my_progress->Update(image_counter + 1);
    }
    wxPrintf("  Particle ID: 1 has the following classification results if 0 idx %f, %f, %f \n", needed_class->classification_results[0].psi, needed_class->classification_results[0].xshift, needed_class->classification_results[0].yshift);
    wxPrintf("  Particle ID: 1 has the following classification results if 1 idx %f, %f, %f \n", needed_class->classification_results[1].psi, needed_class->classification_results[1].xshift, needed_class->classification_results[1].yshift);

    return true;
}
