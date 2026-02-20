#include "../../core/core_headers.h"
#include <iomanip>

class
        save_1d_projection : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

std::vector<float> sum_image_columns(Image* current_image);
void               save_all_columns_sum_to_file(const std::vector<float>& all_columns_sum, const std::string& filename);

IMPLEMENT_APP(save_1d_projection)

// override the DoInteractiveUserInput

void save_1d_projection::DoInteractiveUserInput( ) {

    wxString input_filename;
    wxString output_1d_filename;
    int      image_number;

    UserInput* my_input = new UserInput("rsave_1d_projection", 1.0);
    input_filename      = my_input->GetFilenameFromUser("Input projection mrc filename", "Name of the input projection filename", "input_projection.mrc", true);
    output_1d_filename  = my_input->GetFilenameFromUser("1D projection filename", "Name of the 1D projection filename", "projected_1d_values.txt", false);
    image_number        = my_input->GetIntFromUser("Slice number to extract", "which slice to extract its 1D projection?", "1", 1);

    delete my_input;

    my_current_job.Reset(4); // + 1 the number shown in the string below
    my_current_job.ManualSetArguments("tti", input_filename.ToUTF8( ).data( ), output_1d_filename.ToUTF8( ).data( ), image_number);
}

// override the do calculation method which will be what is actually run..

bool save_1d_projection::DoCalculation( ) {
    // get the arguments for this job..
    wxString input_filename     = my_current_job.arguments[0].ReturnStringArgument( );
    wxString output_1d_filename = my_current_job.arguments[1].ReturnStringArgument( );
    int      image_number       = my_current_job.arguments[2].ReturnIntegerArgument( );

    MRCFile my_input_file(input_filename.ToStdString( ), false); // check all the functions and things done with the MRCFile and also check the wxPrintF statement
    int     x_dim = my_input_file.ReturnXSize( );
    Image   current_image;

    current_image.ReadSlice(&my_input_file, image_number);

    std::vector<float> column_sum_vector;
    column_sum_vector = sum_image_columns(&current_image);
    save_all_columns_sum_to_file(column_sum_vector, output_1d_filename.ToStdString( ));

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

void save_all_columns_sum_to_file(
        const std::vector<float>& all_columns_sum,
        const std::string&        filename) {
    std::ofstream out_file(filename);
    if ( ! out_file.is_open( ) ) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    out_file << std::fixed << std::setprecision(2);

    for ( size_t i = 0; i < all_columns_sum.size( ); ++i ) {
        out_file << all_columns_sum[i];
        if ( i < all_columns_sum.size( ) - 1 )
            out_file << ", ";
    }
    out_file << '\n';
}
