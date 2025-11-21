#include "../../core/core_headers.h"

class
        rotation_calculation : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

IMPLEMENT_APP(rotation_calculation)

// override the DoInteractiveUserInput

void rotation_calculation::DoInteractiveUserInput( ) {

    wxString output_filename;
    int      angular_step = 90;

    UserInput* my_input = new UserInput("rotation_calculation", 1.0);
    output_filename     = my_input->GetFilenameFromUser("Output calculations filename", "Name of the output rotation calculations filename", "rotation_calculations.txt", false);
    angular_step        = my_input->GetIntFromUser("Angular step size", "Angular step size to use in the rotation calculations", "90", 0, 360);
    delete my_input;

    my_current_job.Reset(3); // + 1 the number shown in the string below
    my_current_job.ManualSetArguments("ti", output_filename.ToUTF8( ).data( ), angular_step);
}

// override the do calculation method which will be what is actually run..

bool rotation_calculation::DoCalculation( ) {
    // get the arguments for this job..
    wxString output_filename = my_current_job.arguments[0].ReturnStringArgument( );
    int      angular_step    = my_current_job.arguments[1].ReturnIntegerArgument( );

    int box_size      = 448;
    int box_center    = box_size / 2;
    int x_mask_center = 0.75 * box_size;
    int y_mask_center = 0.5 * box_size;
    int z_mask_center = 0.5 * box_size;
    wxPrintf("\n\nPsi 90 calculations\n\n");

    for ( int phi = 0; phi < 360; phi += angular_step ) {

        RotationMatrix rotation_matrix;
        RotationMatrix new_matrix;
        RotationMatrix inverse_matrix;
        RotationMatrix valid_matrix;
        RotationMatrix valid_inverse_matrix;

        float valid_phi, valid_theta, valid_psi;

        float rotated_x, rotated_y, rotated_z;
        float new_rotated_x, new_rotated_y, new_rotated_z;
        float valid_rotated_x, valid_rotated_y, valid_rotated_z;

        // // generate the full rotation matrix
        rotation_matrix.SetToEulerRotation(-90.0, -90.0, -phi); // removed the negative from all
        new_matrix.SetToEulerRotation(phi, 90.0, 90.0); // removed the negative from all
        inverse_matrix = new_matrix.ReturnTransposed( );

        new_matrix.ConvertToValidEulerAngles(valid_phi, valid_theta, valid_psi);
        valid_matrix.SetToEulerRotation(valid_phi, valid_theta, valid_psi);
        valid_inverse_matrix = valid_matrix.ReturnTransposed( );

        rotation_matrix.RotateCoords((x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), rotated_x, rotated_y, rotated_z);
        inverse_matrix.RotateCoords((x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), new_rotated_x, new_rotated_y, new_rotated_z);
        valid_inverse_matrix.RotateCoords((x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), valid_rotated_x, valid_rotated_y, valid_rotated_z);

        wxPrintf("original mask location in x, y, z at phi %i and psi %f are %i, %i, %i and rotation matrix new location are %f, %f, %f \n", -phi, -(90.0f), (x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), rotated_x, rotated_y, rotated_z);
        wxPrintf("original mask location in x, y, z at phi %i and psi %f are %i, %i, %i and inverse matrix new location are %f, %f, %f \n", phi, 90.0f, (x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), new_rotated_x, new_rotated_y, new_rotated_z);
        wxPrintf("original mask location in x, y, z at phi %f and psi %f are %i, %i, %i and valid inverse matrix new location are %f, %f, %f \n\n", valid_phi, valid_psi, (x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), valid_rotated_x, valid_rotated_y, valid_rotated_z);
    }
    wxPrintf("Psi 270 calculations\n\n");

    for ( int phi = 0; phi < 360; phi += angular_step ) {

        RotationMatrix rotation_matrix;
        RotationMatrix new_matrix;
        RotationMatrix inverse_matrix;
        RotationMatrix valid_matrix;
        RotationMatrix valid_inverse_matrix;

        float valid_phi, valid_theta, valid_psi;

        float rotated_x, rotated_y, rotated_z;
        float new_rotated_x, new_rotated_y, new_rotated_z;
        float valid_rotated_x, valid_rotated_y, valid_rotated_z;

        // // generate the full rotation matrix
        rotation_matrix.SetToEulerRotation(-270.0, -90.0, -phi); // removed the negative from all
        new_matrix.SetToEulerRotation(phi, 90.0, 270.0); // removed the negative from all
        inverse_matrix = new_matrix.ReturnTransposed( );

        new_matrix.ConvertToValidEulerAngles(valid_phi, valid_theta, valid_psi);
        valid_matrix.SetToEulerRotation(valid_phi, valid_theta, valid_psi);
        valid_inverse_matrix = valid_matrix.ReturnTransposed( );

        rotation_matrix.RotateCoords((x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), rotated_x, rotated_y, rotated_z);
        inverse_matrix.RotateCoords((x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), new_rotated_x, new_rotated_y, new_rotated_z);
        valid_inverse_matrix.RotateCoords((x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), valid_rotated_x, valid_rotated_y, valid_rotated_z);

        wxPrintf("original mask location in x, y, z at phi %i and psi %f are %i, %i, %i and rotation matrix new location are %f, %f, %f \n", -phi, -(270.0f), (x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), rotated_x, rotated_y, rotated_z);
        wxPrintf("original mask location in x, y, z at phi %i and psi %f are %i, %i, %i and inverse matrix new location are %f, %f, %f \n", phi, 270.0f, (x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), new_rotated_x, new_rotated_y, new_rotated_z);
        wxPrintf("original mask location in x, y, z at phi %f and psi %f are %i, %i, %i and valid inverse matrix new location are %f, %f, %f \n\n", valid_phi, valid_psi, (x_mask_center - box_center), (y_mask_center - box_center), (z_mask_center - box_center), valid_rotated_x, valid_rotated_y, valid_rotated_z);
    }
    return true;
}
