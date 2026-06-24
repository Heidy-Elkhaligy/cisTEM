#include "../../core/core_headers.h"

class
        SaveOrthViews : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

struct OrthogonalViews {
    Image slice_x, slice_y, slice_z;
    Image proj_x, proj_y, proj_z;
    bool  has_projections;

    OrthogonalViews( ) : has_projections(false) {}
};

OrthogonalViews CreateOrthogonalViewsAndSave(
        Image& volume,
        bool   include_projections,
        float  scale_factor,
        float  mask_radius_in_pixels) {
    //MyDebugAssertTrue(volume.IsCubic( ) == true, "Only Cubic Volumes Supported");

    OrthogonalViews result;
    result.has_projections = include_projections;

    int i, j, k;

    /* allocate slices */
    result.slice_x.Allocate(volume.logical_x_dimension,
                            volume.logical_y_dimension, true);

    result.slice_y.Allocate(volume.logical_x_dimension,
                            volume.logical_y_dimension, true);

    result.slice_z.Allocate(volume.logical_x_dimension,
                            volume.logical_y_dimension, true);

    result.slice_x.SetToConstant(0.0);
    result.slice_y.SetToConstant(0.0);
    result.slice_z.SetToConstant(0.0);

    /* allocate projections if requested */
    if ( include_projections ) {
        result.proj_x.Allocate(volume.logical_x_dimension,
                               volume.logical_y_dimension, true);

        result.proj_y.Allocate(volume.logical_x_dimension,
                               volume.logical_y_dimension, true);

        result.proj_z.Allocate(volume.logical_x_dimension,
                               volume.logical_y_dimension, true);

        result.proj_x.SetToConstant(0.0);
        result.proj_y.SetToConstant(0.0);
        result.proj_z.SetToConstant(0.0);
    }

    long input_counter = 0;

    /* ---- projections ---- */
    if ( include_projections ) {
        for ( k = 0; k < volume.logical_z_dimension; k++ )
            for ( j = 0; j < volume.logical_y_dimension; j++ ) {
                for ( i = 0; i < volume.logical_x_dimension; i++ ) {
                    result.proj_x.real_values[result.proj_x.ReturnReal1DAddressFromPhysicalCoord(i, j, 0)] += volume.real_values[input_counter];

                    result.proj_y.real_values[result.proj_y.ReturnReal1DAddressFromPhysicalCoord(j, k, 0)] += volume.real_values[input_counter];

                    result.proj_z.real_values[result.proj_z.ReturnReal1DAddressFromPhysicalCoord(i, k, 0)] += volume.real_values[input_counter];

                    input_counter++;
                }
                input_counter += volume.padding_jump_value;
            }
    }

    /* ---- central slices ---- */
    long output_counter = 0;

    for ( j = 0; j < result.slice_x.logical_y_dimension; j++ ) {
        for ( i = 0; i < result.slice_x.logical_x_dimension; i++ ) {
            result.slice_x.real_values[output_counter] =
                    volume.ReturnRealPixelFromPhysicalCoord(
                            i, j, volume.physical_address_of_box_center_z);

            result.slice_y.real_values[output_counter] =
                    volume.ReturnRealPixelFromPhysicalCoord(
                            volume.physical_address_of_box_center_x, i, j);

            result.slice_z.real_values[output_counter] =
                    volume.ReturnRealPixelFromPhysicalCoord(
                            i, volume.physical_address_of_box_center_x, j);

            output_counter++;
        }
        output_counter += result.slice_x.padding_jump_value;
    }

    /* ---- scaling ---- */
    if ( scale_factor != 1.0f ) {
        Image* slices[3] = {&result.slice_x, &result.slice_y, &result.slice_z};

        for ( int s = 0; s < 3; s++ ) {
            slices[s]->ForwardFFT( );
            slices[s]->Resize(
                    myroundint(volume.logical_x_dimension * scale_factor),
                    myroundint(volume.logical_y_dimension * scale_factor), 1);
            slices[s]->BackwardFFT( );
            slices[s]->Normalize(1.0);
        }

        if ( include_projections ) {
            Image* projs[3] = {&result.proj_x, &result.proj_y, &result.proj_z};

            for ( int s = 0; s < 3; s++ ) {
                projs[s]->ForwardFFT( );
                projs[s]->Resize(
                        myroundint(volume.logical_x_dimension * scale_factor),
                        myroundint(volume.logical_y_dimension * scale_factor), 1);
                projs[s]->BackwardFFT( );
                projs[s]->Normalize(1.0);
            }
        }
    }

    /* ---- normalize slices together ---- */
    float minv, maxv, tmpmin, tmpmax;

    result.slice_x.GetMinMax(minv, maxv);
    result.slice_y.GetMinMax(tmpmin, tmpmax);
    minv = std::min(minv, tmpmin);
    maxv = std::max(maxv, tmpmax);
    result.slice_z.GetMinMax(tmpmin, tmpmax);
    minv = std::min(minv, tmpmin);
    maxv = std::max(maxv, tmpmax);

    Image* slices[3] = {&result.slice_x, &result.slice_y, &result.slice_z};
    for ( int s = 0; s < 3; s++ ) {
        slices[s]->AddConstant(-minv);
        slices[s]->DivideByConstant(maxv - minv);
    }

    /* ---- normalize projections ---- */
    if ( include_projections ) {
        result.proj_x.GetMinMax(minv, maxv);
        result.proj_y.GetMinMax(tmpmin, tmpmax);
        minv = std::min(minv, tmpmin);
        maxv = std::max(maxv, tmpmax);
        result.proj_z.GetMinMax(tmpmin, tmpmax);
        minv = std::min(minv, tmpmin);
        maxv = std::max(maxv, tmpmax);

        Image* projs[3] = {&result.proj_x, &result.proj_y, &result.proj_z};
        for ( int s = 0; s < 3; s++ ) {
            projs[s]->AddConstant(-minv);
            projs[s]->DivideByConstant(maxv - minv);
        }
    }

    /* ---- masking ---- */
    if ( mask_radius_in_pixels != 0.0f ) {
        float r = mask_radius_in_pixels * scale_factor;

        for ( int s = 0; s < 3; s++ )
            slices[s]->CircleMaskWithValue(
                    r, slices[s]->ReturnAverageOfRealValuesAtRadius(r));

        if ( include_projections ) {
            Image* projs[3] = {&result.proj_x, &result.proj_y, &result.proj_z};
            for ( int s = 0; s < 3; s++ )
                projs[s]->CircleMaskWithValue(
                        r, projs[s]->ReturnAverageOfRealValuesAtRadius(r));
        }
    }

    // /* ---- AUTO SAVE ---- */
    // result.slice_x.QuickAndDirtyWriteSlice("slice_x.mrc", 1);
    // result.slice_y.QuickAndDirtyWriteSlice("slice_y.mrc", 1);
    // result.slice_z.QuickAndDirtyWriteSlice("slice_z.mrc", 1);

    // if ( include_projections ) {
    //     result.proj_x.QuickAndDirtyWriteSlice("projection_x.mrc", 1);
    //     result.proj_y.QuickAndDirtyWriteSlice("projection_y.mrc", 1);
    //     result.proj_z.QuickAndDirtyWriteSlice("projection_z.mrc", 1);
    // }

    return result; // C++11 NRVO/move handles this efficiently
}

IMPLEMENT_APP(SaveOrthViews)

// override the DoInteractiveUserInput

void SaveOrthViews::DoInteractiveUserInput( ) {

    UserInput* my_input = new UserInput("SaveOrthViews", 1.0);

    wxString input_volume = my_input->GetFilenameFromUser("Input image/volume file name", "Name of input image volume", "input.mrc", true);
    wxString output_image = my_input->GetFilenameFromUser("Output orth views image file name", "Name of output image ", "orth.mrc", false);

    delete my_input;

    my_current_job.Reset(2);
    my_current_job.ManualSetArguments("tt", input_volume.ToUTF8( ).data( ), output_image.ToUTF8( ).data( ));
}

// override the do calculation method which will be what is actually run..

bool SaveOrthViews::DoCalculation( ) {

    wxString input_volume = my_current_job.arguments[0].ReturnStringArgument( );
    wxString output_image = my_current_job.arguments[1].ReturnStringArgument( );

    MRCFile input3d_file(input_volume.ToStdString( ), false);
    MRCFile output_file(output_image.ToStdString( ), true);

    Image my_input_volume;
    Image my_orth_views_image;

    wxPrintf("\nMaking orth views image...\n");

    my_input_volume.ReadSlices(&input3d_file, 1, input3d_file.ReturnNumberOfSlices( ));
    /////////////////////////////////////////////////////////////////////////////
    // Old make orthogonal views which saves all images in one large image
    ////////////////////////////////////////////////////////////////////////////////
    // my_orth_views_image.Allocate(my_input_volume.logical_x_dimension * 3, my_input_volume.logical_y_dimension * 2, 1, true);
    // my_input_volume.CreateOrthogonalProjectionsImage(&my_orth_views_image);
    // my_orth_views_image.WriteSlice(&output_file, 1);
    OrthogonalViews views = CreateOrthogonalViewsAndSave(my_input_volume, true, 1.0f, 0.0f);

    // Saving each orthogonal view independently so it can be used later for 1D projections or any analysis
    views.slice_x.WriteSlice(&output_file, 1);
    views.slice_y.WriteSlice(&output_file, 2);
    views.slice_z.WriteSlice(&output_file, 3);

    if ( views.has_projections ) {
        views.proj_x.WriteSlice(&output_file, 4);
        views.proj_y.WriteSlice(&output_file, 5);
        views.proj_z.WriteSlice(&output_file, 6);
    }

    return true;
}
