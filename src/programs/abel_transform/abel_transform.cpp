#include "../../core/core_headers.h"
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>

class
        abel_transform : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

using ComplexD = std::complex<double>;
using ComplexF = std::complex<float>;

static inline std::complex<double>
ReadComplexFromImage2D(Image* img, int x, int y) {
    // Use cisTEM's accessor exactly as defined
    std::complex<float> cp =
            img->ReturnComplexPixelFromLogicalCoord(
                    x, y, 0,
                    std::complex<float>(0.0f, 0.0f) // out-of-bounds value
            );

    return std::complex<double>(cp.real( ), cp.imag( ));
}

static inline void
WriteComplexToImage2D(Image* img, int x, int y, const std::complex<double>& val) {
    long addr = img->ReturnFourier1DAddressFromLogicalCoord(x, y, 0);
    img->complex_values[addr] =
            std::complex<float>(float(val.real( )), float(val.imag( )));
}

static void BilinearAddComplex2D(std::vector<ComplexD>& plane, int NX, int NY,
                                 double x, double y, const ComplexD& val);
void        FourierRotateToCentralSlice(Image* proj2D,
                                        Image* out2D,
                                        int    nAngles,
                                        int    outSize);

// void GenerateRotationalVolumeFrom2D(Image* current_image);
void sum_image_direction(Image* current_image, int dim);

//void AverageRotationallyFourierSlices(Image* img);

void Generate3DVolumeFromFilament2DClass(Image* image, Image* volume);

// Helper: create 1D Gaussian kernel (normalized)
static float* MakeGaussianKernel(float sigma, int& kernel_radius) {
    if ( sigma <= 0.0f ) {
        kernel_radius = 0;
        return nullptr;
    }
    // choose radius = ceil(3*sigma)
    kernel_radius      = std::max(1, int(std::ceil(3.0f * sigma)));
    int         size   = 2 * kernel_radius + 1;
    float*      kernel = new float[size];
    float       sum    = 0.0f;
    const float denom  = 2.0f * sigma * sigma;
    for ( int i = -kernel_radius; i <= kernel_radius; ++i ) {
        float v                   = std::exp(-(i * i) / denom);
        kernel[i + kernel_radius] = v;
        sum += v;
    }
    // normalize
    if ( sum > 0.0f ) {
        for ( int i = 0; i < size; ++i )
            kernel[i] /= sum;
    }
    return kernel;
}

// Helper: convolve volume along Z (separable 1D convolution along axis z)
// Assumes volume logical dims: X x Y x Z. Uses edge clamp for out-of-bounds.
static void ConvolveAlongZ(Image* vol, float* kernel, int kernel_radius) {
    if ( ! kernel || kernel_radius == 0 )
        return;
    long X   = vol->logical_x_dimension;
    long Y   = vol->logical_y_dimension;
    long Z   = vol->logical_z_dimension;
    long pad = vol->padding_jump_value;

    float* temp = new float[X * Y * Z]; // contiguous scratch (we'll index manually)
    // initialize to zero
    for ( long i = 0; i < X * Y * Z; ++i )
        temp[i] = 0.0f;

    for ( long z = 0; z < Z; ++z ) {
        for ( long y = 0; y < Y; ++y ) {
            for ( long x = 0; x < X; ++x ) {
                float accum = 0.0f;
                for ( int k = -kernel_radius; k <= kernel_radius; ++k ) {
                    long zz = z + k;
                    if ( zz < 0 )
                        zz = 0;
                    if ( zz >= Z )
                        zz = Z - 1;
                    long addr = vol->ReturnReal1DAddressFromPhysicalCoord(x, y, zz);
                    accum += vol->real_values[addr] * kernel[k + kernel_radius];
                }
                long idx  = (z * Y + y) * X + x;
                temp[idx] = accum;
            }
        }
    }
    // copy back
    for ( long z = 0; z < Z; ++z ) {
        for ( long y = 0; y < Y; ++y ) {
            for ( long x = 0; x < X; ++x ) {
                long idx               = (z * Y + y) * X + x;
                long addr              = vol->ReturnReal1DAddressFromPhysicalCoord(x, y, z);
                vol->real_values[addr] = temp[idx];
            }
        }
    }
    delete[] temp;
}

// Helper: convolve along X and Y (in-place separable). Edge clamp used.
static void ConvolveAlongX(Image* vol, float* kernel, int kernel_radius) {
    if ( ! kernel || kernel_radius == 0 )
        return;
    long X = vol->logical_x_dimension;
    long Y = vol->logical_y_dimension;
    long Z = vol->logical_z_dimension;

    float* rowbuf = new float[X];

    for ( long z = 0; z < Z; ++z ) {
        for ( long y = 0; y < Y; ++y ) {
            // extract
            for ( long x = 0; x < X; ++x ) {
                long addr = vol->ReturnReal1DAddressFromPhysicalCoord(x, y, z);
                rowbuf[x] = vol->real_values[addr];
            }
            // convolve into tmp buffer then write back
            for ( long x = 0; x < X; ++x ) {
                float accum = 0.0f;
                for ( int k = -kernel_radius; k <= kernel_radius; ++k ) {
                    long xx = x + k;
                    if ( xx < 0 )
                        xx = 0;
                    if ( xx >= X )
                        xx = X - 1;
                    accum += rowbuf[xx] * kernel[k + kernel_radius];
                }
                long addr              = vol->ReturnReal1DAddressFromPhysicalCoord(x, y, z);
                vol->real_values[addr] = accum;
            }
        }
    }

    delete[] rowbuf;
}

static void ConvolveAlongY(Image* vol, float* kernel, int kernel_radius) {
    if ( ! kernel || kernel_radius == 0 )
        return;
    long X = vol->logical_x_dimension;
    long Y = vol->logical_y_dimension;
    long Z = vol->logical_z_dimension;

    float* colbuf = new float[Y];

    for ( long z = 0; z < Z; ++z ) {
        for ( long x = 0; x < X; ++x ) {
            // extract column
            for ( long y = 0; y < Y; ++y ) {
                long addr = vol->ReturnReal1DAddressFromPhysicalCoord(x, y, z);
                colbuf[y] = vol->real_values[addr];
            }
            // convolve
            for ( long y = 0; y < Y; ++y ) {
                float accum = 0.0f;
                for ( int k = -kernel_radius; k <= kernel_radius; ++k ) {
                    long yy = y + k;
                    if ( yy < 0 )
                        yy = 0;
                    if ( yy >= Y )
                        yy = Y - 1;
                    accum += colbuf[yy] * kernel[k + kernel_radius];
                }
                long addr              = vol->ReturnReal1DAddressFromPhysicalCoord(x, y, z);
                vol->real_values[addr] = accum;
            }
        }
    }

    delete[] colbuf;
}

// 1D Gaussian smoothing for a radial profile (in-place)
static void Smooth1DRadialProfile(float* profile, long length, float sigma) {
    if ( sigma <= 0.0f || length <= 1 )
        return;
    int    radius;
    float* kernel = MakeGaussianKernel(sigma, radius);
    if ( ! kernel )
        return;

    float* tmp = new float[length];
    for ( long i = 0; i < length; ++i ) {
        float accum = 0.0f;
        for ( int k = -radius; k <= radius; ++k ) {
            long idx = i + k;
            if ( idx < 0 )
                idx = 0;
            if ( idx >= length )
                idx = length - 1;
            accum += profile[idx] * kernel[k + radius];
        }
        tmp[i] = accum;
    }
    for ( long i = 0; i < length; ++i )
        profile[i] = tmp[i];

    delete[] tmp;
    delete[] kernel;
}

// The main function: build 3D ab initio volume from 2D image with centered/un-centered modes + smoothing
struct AbInitioOptions {
    float radial_profile_sigma    = 0.0f; // gaussian smoothing on radial profiles (centered case)
    float z_smooth_sigma          = 1.0f; // gaussian smoothing along Z after filling (0 = none)
    bool  apply_inplane_smoothing = false;
    float inplane_sigma           = 0.0f; // gaussian smoothing in XY (applied separably X then Y)
};

void BuildAbInitioVolumeFrom2D(Image* input2d, Image* out_volume, bool is_centered, const AbInitioOptions& opts) {
    // keep FFT semantics
    bool input_in_fourier = false;
    if ( ! input2d->is_in_real_space ) {
        input2d->BackwardFFT( );
        input_in_fourier = true;
    }

    long W = input2d->logical_x_dimension;
    long H = input2d->logical_y_dimension;
    long D = H; // depth = number of rows

    // resize output volume
    out_volume->Resize(W, H, D, 0.0f);

    // center coords
    long cx = input2d->physical_address_of_box_center_x;
    long cy = input2d->physical_address_of_box_center_y;
    // number of radial bins (same logic you used before)
    long  number_of_rings = (W + 1) / 2;
    float max_radius_x    = float(cx);
    // construct ring_axis
    float* ring_axis = new float[number_of_rings];
    if ( number_of_rings > 1 ) {
        for ( long r = 0; r < number_of_rings; ++r )
            ring_axis[r] = 0.0f + r * (max_radius_x - 0.0f) / float(number_of_rings - 1);
    }
    else
        ring_axis[0] = 0.0f;
    float bin_step = (number_of_rings > 1) ? (ring_axis[1] - ring_axis[0]) : 1.0f;
    if ( bin_step <= 1e-8f )
        bin_step = 1.0f;
    long edge_bin = number_of_rings - 1;

    // We'll either produce per-row radial profiles (centered) or fill slices directly (uncentered).
    // For centered: allocate profiles[H][number_of_rings]
    float* profiles = nullptr;
    if ( is_centered ) {
        profiles = new float[H * number_of_rings];
        // build profiles
        for ( long y = 0; y < H; ++y ) {
            // zero
            ZeroArray(&profiles[y * number_of_rings], int(number_of_rings));
            // also weights
            float* weights = new float[number_of_rings];
            ZeroArray(weights, int(number_of_rings));

            // average symmetric pixels about center
            bool width_even   = (W % 2 == 0);
            long center_left  = (width_even) ? (cx - 1) : cx;
            long center_right = cx;

            for ( long r = 0; r < number_of_rings; ++r ) {
                long left_x, right_x;
                if ( width_even ) {
                    left_x  = center_left - r;
                    right_x = center_right + r;
                }
                else {
                    left_x  = center_left - r;
                    right_x = center_right + r;
                }
                if ( left_x < 0 )
                    left_x = 0;
                if ( right_x >= W )
                    right_x = W - 1;

                long  addrL = input2d->ReturnReal1DAddressFromPhysicalCoord(left_x, y, 0);
                long  addrR = input2d->ReturnReal1DAddressFromPhysicalCoord(right_x, y, 0);
                float vL    = input2d->real_values[addrL];
                float vR    = input2d->real_values[addrR];

                if ( left_x == right_x ) {
                    profiles[y * number_of_rings + r] += vL;
                    weights[r] += 1.0f;
                }
                else {
                    profiles[y * number_of_rings + r] += vL + vR;
                    weights[r] += 2.0f;
                }
            }
            // normalize
            for ( long r = 0; r < number_of_rings; ++r ) {
                if ( weights[r] != 0.0f )
                    profiles[y * number_of_rings + r] /= weights[r];
                else
                    profiles[y * number_of_rings + r] = 0.0f;
            }

            // optional radial smoothing on profile
            if ( opts.radial_profile_sigma > 0.0f ) {
                Smooth1DRadialProfile(&profiles[y * number_of_rings], number_of_rings, opts.radial_profile_sigma);
            }

            delete[] weights;
        }
    }

    // Fill the volume slice-by-slice
    // For centered: use profiles[z]
    // For uncentered: for each voxel at (x,y,z=row), compute polar coords (r,phi) and sample input row y at
    //  x_sample = cx + sign(cos(phi)) * r  (linear interpolation)
    for ( long z = 0; z < D; ++z ) {
        long src_row = z; // mapping: row -> z slice
        for ( long y = 0; y < H; ++y ) {
            float y_radius_sq = float(cy - y) * float(cy - y);
            for ( long x = 0; x < W; ++x ) {
                float dx     = float(cx - x);
                float dx_sq  = dx * dx;
                float radius = std::sqrt(y_radius_sq + dx_sq);

                long out_addr = out_volume->ReturnReal1DAddressFromPhysicalCoord(x, y, z);

                if ( is_centered ) {
                    long index_of_bin = long((radius - ring_axis[0]) / bin_step);
                    if ( index_of_bin >= edge_bin ) {
                        out_volume->real_values[out_addr] = profiles[src_row * number_of_rings + edge_bin];
                    }
                    else if ( index_of_bin < 0 ) {
                        out_volume->real_values[out_addr] = profiles[src_row * number_of_rings + 0];
                    }
                    else {
                        float diff = (radius - ring_axis[index_of_bin]) / (ring_axis[index_of_bin + 1] - ring_axis[index_of_bin]);
                        if ( diff < 0.0f )
                            diff = 0.0f;
                        if ( diff > 1.0f )
                            diff = 1.0f;
                        float vlow                        = profiles[src_row * number_of_rings + index_of_bin];
                        float vhigh                       = profiles[src_row * number_of_rings + index_of_bin + 1];
                        out_volume->real_values[out_addr] = vlow * (1.0f - diff) + vhigh * diff;
                    }
                }
                else {
                    // uncentered: rotate row as-is around center.
                    // compute angle phi for the vector from center->(x,y)
                    float vx       = float(x - cx);
                    float vy       = float(y - cy);
                    float phi      = std::atan2(vy, vx); // range [-pi, pi]
                    float cosphi   = std::cos(phi);
                    int   sign     = (cosphi >= 0.0f) ? 1 : -1; // choose right (sign=+1) or left (sign=-1) half of the row
                    float sample_x = float(cx) + float(sign) * radius; // continuous x sample
                    // map sample_x to input row coordinates (0..W-1)
                    if ( sample_x <= 0.0f ) {
                        long addr                         = input2d->ReturnReal1DAddressFromPhysicalCoord(0, src_row, 0);
                        out_volume->real_values[out_addr] = input2d->real_values[addr];
                    }
                    else if ( sample_x >= float(W - 1) ) {
                        long addr                         = input2d->ReturnReal1DAddressFromPhysicalCoord(W - 1, src_row, 0);
                        out_volume->real_values[out_addr] = input2d->real_values[addr];
                    }
                    else {
                        long  xi                          = long(std::floor(sample_x));
                        long  xi2                         = xi + 1;
                        float frac                        = sample_x - float(xi);
                        long  addr1                       = input2d->ReturnReal1DAddressFromPhysicalCoord(xi, src_row, 0);
                        long  addr2                       = input2d->ReturnReal1DAddressFromPhysicalCoord(xi2, src_row, 0);
                        float v1                          = input2d->real_values[addr1];
                        float v2                          = input2d->real_values[addr2];
                        out_volume->real_values[out_addr] = v1 * (1.0f - frac) + v2 * frac;
                    }
                }
            } // x
        } // y
    } // z

    // optional smoothing along Z to enforce continuity between slices (recommended)
    if ( opts.z_smooth_sigma > 0.0f ) {
        int    kr;
        float* k = MakeGaussianKernel(opts.z_smooth_sigma, kr);
        if ( k ) {
            ConvolveAlongZ(out_volume, k, kr);
            delete[] k;
        }
    }

    // optional in-plane smoothing (XY) to remove high-frequency ringing
    if ( opts.apply_inplane_smoothing && opts.inplane_sigma > 0.0f ) {
        int    krx, kry;
        float* kx = MakeGaussianKernel(opts.inplane_sigma, krx);
        float* ky = MakeGaussianKernel(opts.inplane_sigma, kry);
        if ( kx ) {
            ConvolveAlongX(out_volume, kx, krx);
            delete[] kx;
        }
        if ( ky ) {
            ConvolveAlongY(out_volume, ky, kry);
            delete[] ky;
        }
    }

    // cleanup
    delete[] ring_axis;
    if ( profiles )
        delete[] profiles;

    if ( input_in_fourier ) {
        out_volume->ForwardFFT( );
    }
}

IMPLEMENT_APP(abel_transform)

void abel_transform::DoInteractiveUserInput( ) {
    wxString input_images;
    wxString output_image;
    float    pixel_size;
    float    padding_factor = 1.0;
    bool     crop;
    int      max_threads;

    UserInput* my_input = new UserInput("abel_transform", 1.00);
    input_images        = my_input->GetFilenameFromUser("Input images file name", "Filen name of helical tube stack aligned vertically", "helical_stack.mrc", true);
    output_image        = my_input->GetFilenameFromUser("Output file name", "The output images file name", "helical_stack.mrc", false);
    pixel_size          = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
    padding_factor      = my_input->GetFloatFromUser("Padding factor", "Factor value to be used for padding the input images", "0.0", 1.0);
    crop                = my_input->GetYesNoFromUser("Crop the output image to the original image size?", "Crop the output image to the original image size or keep it padded", "No");
#ifdef _OPENMP
    max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "when threading, what is the max threads to run", "1", 1);
#else
    max_threads = 1;
#endif
    delete my_input;

    my_current_job.Reset(7);
    my_current_job.ManualSetArguments("ttffbi", input_images.ToUTF8( ).data( ), output_image.ToUTF8( ).data( ),
                                      pixel_size, padding_factor, crop, max_threads);
}

// override the do calculation method which will be what is actually run..

bool abel_transform::DoCalculation( ) {
    wxString input_images   = my_current_job.arguments[0].ReturnStringArgument( );
    wxString output_image   = my_current_job.arguments[1].ReturnStringArgument( );
    float    pixel_size     = my_current_job.arguments[2].ReturnFloatArgument( );
    float    padding_factor = my_current_job.arguments[3].ReturnFloatArgument( );
    bool     crop           = my_current_job.arguments[4].ReturnBoolArgument( );
    int      max_threads    = my_current_job.arguments[5].ReturnIntegerArgument( );

    MRCFile my_input_images(input_images.ToStdString( ), false);
    MRCFile my_output_image(output_image.ToStdString( ), true);
    long    number_of_input_images = my_input_images.ReturnNumberOfSlices( );
    int     x_dim                  = my_input_images.ReturnXSize( );
    int     y_dim                  = my_input_images.ReturnYSize( );

    // Image my_image;
    // my_image.Allocate(x_dim * padding_factor, y_dim * padding_factor, true);
    // my_image.SetToConstant(0.0);
    // my_image.Deallocate( );

    Image current_image;
    current_image.Allocate(x_dim, y_dim, true);
    current_image.SetToConstant(0.0);
    current_image.ReadSlice(&my_input_images, 1);
    float image_average;
    image_average = current_image.ReturnAverageOfRealValues( );

    // mask the 2D image to avoid ring effect outside the azimuthal average volume
    float mask_edge = x_dim * 0.05f; // 5% of the pixels will constitute the soft edge

    current_image.CosineMask(x_dim * 0.4, mask_edge, false, true, image_average);
    current_image.Resize(x_dim * 2, y_dim * 2, 1, image_average);
    current_image.ApplyRampFilter( );
    current_image.AverageRotationally( );
    current_image.Resize(x_dim, y_dim, 1, 0.0f);
    sum_image_direction(&current_image, 2);
    current_image.AverageRotationally( );
    current_image.ApplyRampFilter( );
    sum_image_direction(&current_image, 2);
    current_image.AverageRotationally( );
    // current_image.ForwardFFT( );
    // current_image.ZeroCentralPixel( ); // to remove the central pixel which has the brightest color in FS
    // current_image.BackwardFFT( );
    //current_image.SwapRealSpaceQuadrants( );
    current_image.WriteSlice(&my_output_image, 1);

    Image input_img;
    Image output_vol;

    input_img.Allocate(x_dim, y_dim, true);
    input_img.SetToConstant(0.0);
    output_vol.Allocate(x_dim, y_dim, x_dim, true, true);
    output_vol.SetToConstant(0.0);

    input_img.ReadSlice(&my_input_images, 1);
    //Generate3DVolumeFromFilament2DClass(&input_img, &output_vol);
    // Set options
    // AbInitioOptions opts;
    // opts.radial_profile_sigma    = 0.5f; // light radial smoothing
    // opts.z_smooth_sigma          = 1.0f; // smooth between slices
    // opts.apply_inplane_smoothing = true;
    // opts.inplane_sigma           = 0.5f;

    // // Build ab initio volume
    // BuildAbInitioVolumeFrom2D(
    //         &input_img,
    //         &output_vol,
    //         false, // is_centered
    //         opts);

    // output_vol.QuickAndDirtyWriteSlices("average_volume_new_method.mrc", 1, x_dim);

    // // FourierRotateToCentralSlice(&current_image, nullptr, 360, -1);
    // // current_image.SwapRealSpaceQuadrants( );
    // // current_image.WriteSlice(&my_output_image, 1);

    // // now current_image holds the reconstructed central slice (real-space)

    // Image centralSlice;
    // FourierRotateToCentralSlice(&current_image, &centralSlice, 720, current_image.logical_x_dimension);
    // // centralSlice now contains the reconstructed real-space central slice
    // centralSlice.SwapRealSpaceQuadrants( );
    // centralSlice.WriteSlice(&my_output_image, 1);

    // // Assuming input_image is a single 2D slice (z=1) and centered in its box
    // Image proj, slice;
    // int   w      = current_image.logical_x_dimension;
    // float center = (w - 1.0f) / 2.0f;
    // proj.Allocate(x_dim, y_dim, true);
    // proj.SetToConstant(0.0);
    // // Step 1: get 1D projection in x-direction
    // proj.CopyFrom(&current_image);
    // sum_image_direction(&proj, 2); // now proj(x,y) = column-average of input
    // proj.MultiplyByConstant(float(current_image.logical_y_dimension)); // undo averaging

    // // Step 2: compute derivative dF/dx
    // std::vector<float> F(w);
    // for ( int x = 0; x < w; x++ ) {
    //     // e.g. take the projection value at y=0 (all y are same)
    //     long addr = proj.ReturnReal1DAddressFromPhysicalCoord(x, 0, 0);
    //     F[x]      = proj.real_values[addr];
    // }
    // std::vector<float> dFdx(w);
    // dFdx[0] = F[1] - F[0]; // forward diff at edge
    // for ( int x = 1; x < w - 1; x++ ) {
    //     dFdx[x] = 0.5f * (F[x + 1] - F[x - 1]); // central diff
    // }
    // dFdx[w - 1] = F[w - 1] - F[w - 2];

    // // Step 3: compute f(r) for r = 0..center
    // int                half = int(center + 0.5f);
    // std::vector<float> fval(half + 1, 0.0f);
    // for ( int ir = 0; ir <= half; ir++ ) {
    //     double r        = ir;
    //     double integral = 0.0;
    //     // sum from x = ir+1 to edge to avoid singularity at x = r
    //     for ( int x = ir + 1; x <= half; x++ ) {
    //         double y     = double(x);
    //         double denom = sqrt(y * y - r * r);
    //         if ( denom > 0.0 ) {
    //             integral += dFdx[x] / denom;
    //         }
    //     }
    //     fval[ir] = float(-integral / M_PI);
    // }

    // // Step 4: write to output slice (radially symmetric)
    // slice.Allocate(w, w, true);
    // slice.SetToConstant(0.0f);
    // for ( int y = 0; y < w; y++ ) {
    //     for ( int x = 0; x < w; x++ ) {
    //         double dx = x - center;
    //         double dy = y - center;
    //         double rr = sqrt(dx * dx + dy * dy);
    //         int    ir = int(rr + 0.5);
    //         if ( ir <= half ) {
    //             long addr               = slice.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
    //             slice.real_values[addr] = fval[ir];
    //         }
    //     }
    // }

    // // //sum_image_direction(&current_image, 2);
    // // // current_image.ApplyRampFilter( );
    // // //GenerateRotationalVolumeFrom2D(&current_image);
    // // //AverageRotationallyFourierSlices(&current_image);
    // // FourierSliceRadialAmplitudeSymmetrize(&current_image);
    // // // ApplyInverseAbelGradient(&trial_image);
    // // //current_image.SwapRealSpaceQuadrants( );

    // // // current_image.ForwardFFT( );
    // // // // current_image.ZeroCentralPixel( );
    // // // current_image.AverageRadially( );
    // // // current_image.BackwardFFT( );
    // slice.WriteSlice(&my_output_image, 1);
    // proj.Deallocate( );
    // slice.Deallocate( );
    current_image.Deallocate( );

    return true;
}

// void Image::GenerateRotationalVolumeFrom2D( Image* current_image) {

//     // 1. Ensure Input is in Fourier Space
//     if ( current_image.is_in_real_space ) {
//         current_image->ForwardFFT( );
//     }

//     // 2. Backup 2D Data and Dimensions
//     long input_logical_x = current_image->logical_x_dimension;
//     long input_logical_y = current_image->logical_y_dimension; // Symmetry axis

//     // Calculate dimensions of the 2D complex array
//     // The physical_upper_bound_complex_x is usually (logical_x / 2)
//     long input_complex_x_dim = current_image->physical_upper_bound_complex_x + 1;
//     long input_data_size     = input_complex_x_dim * input_logical_y;

//     std::vector<std::complex<float>> input_2d_data(input_data_size);
//     for ( long i = 0; i < input_data_size; i++ ) {
//         input_2d_data[i] = current_image->complex_values[i];
//     }

//     // 3. Resize to 3D
//     // New Volume: Width=X, Depth=X, Axis=Y(Old Axis)
//     current_image->Resize(input_logical_x, input_logical_x, input_logical_y);

//     // 4. Calculate New Dimensions for Indexing
//     long output_complex_x_dim = current_image->physical_upper_bound_complex_x + 1;
//     long output_complex_y_dim = current_image->physical_upper_bound_complex_y + 1;

//     // 5. Perform Rotational Expansion
//     for ( int kz = 0; kz <= current_image->physical_upper_bound_complex_z; kz++ ) {

//         // Pre-calculate the offset for current_image Z-slice
//         // Z-slice size = (Y_dim * X_complex_dim)
//         long z_slice_offset = kz * output_complex_y_dim * output_complex_x_dim;

//         // Input Lookup: The input row corresponding to current_image height (Z)
//         // Note: We clamp the Z to the input height just in case, though they should match.
//         long input_row_index = kz;
//         if ( input_row_index >= input_logical_y )
//             input_row_index = input_logical_y - 1;
//         long input_row_offset = input_row_index * input_complex_x_dim;

//         for ( int ky = 0; ky <= current_image->physical_upper_bound_complex_y; ky++ ) {

//             // Pre-calculate the offset for current_image Y-row within the Z-slice
//             long y_row_offset = ky * output_complex_x_dim;

//             float ky_coord = ReturnFourierLogicalCoordGivenPhysicalCoord_Y(ky) * fourier_voxel_size_y;
//             float ky_sq    = ky_coord * ky_coord;

//             for ( int kx = 0; kx <= current_image->physical_upper_bound_complex_x; kx++ ) {

//                 // MANUAL INDEX CALCULATION
//                 long output_index = z_slice_offset + y_row_offset + kx;

//                 float kx_coord = ReturnFourierLogicalCoordGivenPhysicalCoord_X(kx) * fourier_voxel_size_x;

//                 // Calculate Radial Frequency (distance from center of slice)
//                 float radius_freq = sqrtf(kx_coord * kx_coord + ky_sq);

//                 // Map to Input Index
//                 long input_col_index = long(radius_freq / fourier_voxel_size_x);

//                 if ( input_col_index < input_complex_x_dim ) {
//                     long input_index                   = input_row_offset + input_col_index;
//                     current_image->complex_values[output_index] = input_2d_data[input_index];
//                 }
//                 else {
//                     current_image->complex_values[output_index] = std::complex<float>(0.0f, 0.0f);
//                 }
//             }
//         }
//     }

//     // 6. Transform to Real Space
//     current_image->BackwardFFT( );
// }

void GenerateRotationalVolumeFrom2D(Image* current_image) {

    // 1. Ensure Input is in Fourier Space
    if ( current_image->is_in_real_space ) {
        current_image->ForwardFFT( );
    }

    // 2. Backup 2D Data
    // We still need to backup because Resize() clears the memory,
    // even if we aren't changing the resolution.
    long input_logical_x     = current_image->logical_x_dimension;
    long input_logical_y     = current_image->logical_y_dimension; // Symmetry axis
    long input_complex_x_dim = current_image->physical_upper_bound_complex_x + 1;
    long input_data_size     = input_complex_x_dim * input_logical_y;

    std::vector<std::complex<float>> input_2d_data(input_data_size);
    for ( long i = 0; i < input_data_size; i++ ) {
        input_2d_data[i] = current_image->complex_values[i];
    }

    // 3. Resize to 3D (Change Shape Only)
    // No padding here. Just changing from 2D (X,Y) to 3D (X,X,Y)
    current_image->Resize(input_logical_x, input_logical_x, input_logical_y);

    // 4. Calculate Dimensions
    long output_complex_x_dim = current_image->physical_upper_bound_complex_x + 1;
    long output_complex_y_dim = current_image->physical_upper_bound_complex_y + 1;

    // 5. Perform Rotational Expansion
    for ( int kz = 0; kz <= current_image->physical_upper_bound_complex_z; kz++ ) {

        long z_slice_offset = kz * output_complex_y_dim * output_complex_x_dim;

        // Input Lookup: Map Z directly to Input Y (Height)
        long input_row_index = kz;
        if ( input_row_index >= input_logical_y )
            input_row_index = input_logical_y - 1;
        long input_row_offset = input_row_index * input_complex_x_dim;

        for ( int ky = 0; ky <= current_image->physical_upper_bound_complex_y; ky++ ) {

            long  y_row_offset = ky * output_complex_x_dim;
            float ky_coord     = current_image->ReturnFourierLogicalCoordGivenPhysicalCoord_Y(ky) * current_image->fourier_voxel_size_y;
            float ky_sq        = ky_coord * ky_coord;

            for ( int kx = 0; kx <= current_image->physical_upper_bound_complex_x; kx++ ) {

                long  output_index = z_slice_offset + y_row_offset + kx;
                float kx_coord     = current_image->ReturnFourierLogicalCoordGivenPhysicalCoord_X(kx) * current_image->fourier_voxel_size_x;

                // Calculate Radius
                float radius_freq = sqrtf(kx_coord * kx_coord + ky_sq);

                // Simple Lookup (Nearest Neighbor)
                // Since we aren't resizing, 1 pixel distance = 1 input index
                long input_col_index = long(radius_freq / current_image->fourier_voxel_size_x);

                if ( input_col_index < input_complex_x_dim ) {
                    current_image->complex_values[output_index] = input_2d_data[input_row_offset + input_col_index];
                }
                else {
                    current_image->complex_values[output_index] = std::complex<float>(0.0f, 0.0f);
                }
            }
        }
    }

    // 6. Inline Radial Averaging (The "Cleaner")
    // current_image step is vital now. Without upsampling, the nearest-neighbor lookup above
    // is rough. current_image smooths it out into perfect circles.

    float max_slice_radius = sqrtf(powf(current_image->physical_upper_bound_complex_x * current_image->fourier_voxel_size_x, 2) +
                                   powf(current_image->physical_upper_bound_complex_y * current_image->fourier_voxel_size_y, 2));
    int   num_bins         = int(max_slice_radius / current_image->fourier_voxel_size_x) + 2;

    std::vector<double> radial_sum(num_bins);
    std::vector<long>   radial_count(num_bins);
    std::vector<float>  radial_average(num_bins);

    for ( int kz = 0; kz <= current_image->physical_upper_bound_complex_z; kz++ ) {

        std::fill(radial_sum.begin( ), radial_sum.end( ), 0.0);
        std::fill(radial_count.begin( ), radial_count.end( ), 0);
        long z_slice_offset = kz * output_complex_y_dim * output_complex_x_dim;

        // 6a. Accumulate
        for ( int ky = 0; ky <= current_image->physical_upper_bound_complex_y; ky++ ) {
            float ky_coord = current_image->ReturnFourierLogicalCoordGivenPhysicalCoord_Y(ky) * current_image->fourier_voxel_size_y;
            float ky_sq    = ky_coord * ky_coord;
            for ( int kx = 0; kx <= current_image->physical_upper_bound_complex_x; kx++ ) {
                float kx_coord = current_image->ReturnFourierLogicalCoordGivenPhysicalCoord_X(kx) * current_image->fourier_voxel_size_x;
                float r        = sqrtf(kx_coord * kx_coord + ky_sq);
                int   bin_idx  = int(r / current_image->fourier_voxel_size_x);
                if ( bin_idx < num_bins ) {
                    long output_index = z_slice_offset + ky * output_complex_x_dim + kx;
                    radial_sum[bin_idx] += std::abs(current_image->complex_values[output_index]);
                    radial_count[bin_idx]++;
                }
            }
        }

        // 6b. Average
        for ( int b = 0; b < num_bins; b++ ) {
            radial_average[b] = (radial_count[b] > 0) ? float(radial_sum[b] / radial_count[b]) : 0.0f;
        }

        // 6c. Write Back (Interpolated)
        for ( int ky = 0; ky <= current_image->physical_upper_bound_complex_y; ky++ ) {
            float ky_coord = current_image->ReturnFourierLogicalCoordGivenPhysicalCoord_Y(ky) * current_image->fourier_voxel_size_y;
            float ky_sq    = ky_coord * ky_coord;
            for ( int kx = 0; kx <= current_image->physical_upper_bound_complex_x; kx++ ) {
                float kx_coord = current_image->ReturnFourierLogicalCoordGivenPhysicalCoord_X(kx) * current_image->fourier_voxel_size_x;
                float r        = sqrtf(kx_coord * kx_coord + ky_sq);

                float precise_bin = r / current_image->fourier_voxel_size_x;
                int   bin_low     = int(precise_bin);
                int   bin_high    = bin_low + 1;
                float t           = precise_bin - float(bin_low);

                float val_low  = (bin_low < num_bins) ? radial_average[bin_low] : 0.0f;
                float val_high = (bin_high < num_bins) ? radial_average[bin_high] : 0.0f;

                float interpolated_amp = val_low * (1.0f - t) + val_high * t;

                long output_index                           = z_slice_offset + ky * output_complex_x_dim + kx;
                current_image->complex_values[output_index] = std::complex<float>(interpolated_amp, 0.0f);
            }
        }
    }

    // 7. Transform to Real Space
    if ( current_image->is_in_real_space ) {
        current_image->BackwardFFT( );
    }
}

void sum_image_direction(Image* current_image, int dim) {
    // image must be in real-space
    Image directional_image_sum;
    directional_image_sum.Allocate(current_image->logical_x_dimension, current_image->logical_y_dimension, true);
    directional_image_sum.SetToConstant(0.0);

    // x-direction
    if ( dim == 1 ) {

        long pixel_counter = 0;

        // sum columns of my_image_sum (NxM) and store in array (1xN)
        for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
            long pixel_coord_y = current_image->ReturnReal1DAddressFromPhysicalCoord(0, j, 0);
            for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
                long pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                directional_image_sum.real_values[pixel_coord_y] += current_image->real_values[pixel_coord_xy];
                pixel_counter++;
            }
            pixel_counter += current_image->padding_jump_value;
        }

        // repeat column sum into my_vertical_sum
        pixel_counter = 0;
        for ( int j = 0; j < directional_image_sum.logical_y_dimension; j++ ) {
            long pixel_coord_y = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(0, j, 0);
            for ( int i = 0; i < directional_image_sum.logical_x_dimension; i++ ) {
                long pixel_coord_xy                               = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                directional_image_sum.real_values[pixel_coord_xy] = directional_image_sum.real_values[pixel_coord_y];
                pixel_counter++;
            }
            pixel_counter += directional_image_sum.padding_jump_value;
        }

        directional_image_sum.DivideByConstant(directional_image_sum.logical_x_dimension);
    }
    // y-direction
    else {

        long pixel_counter = 0;

        // sum columns of my_image_sum (NxM) and store in array (1xM)
        for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
            long pixel_coord_x = current_image->ReturnReal1DAddressFromPhysicalCoord(i, 0, 0);
            for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
                long pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                directional_image_sum.real_values[pixel_coord_x] += current_image->real_values[pixel_coord_xy];
                pixel_counter++;
            }
            pixel_counter += current_image->padding_jump_value;
        }

        // repeat column sum into my_vertical_sum
        pixel_counter = 0;
        for ( int i = 0; i < directional_image_sum.logical_x_dimension; i++ ) {
            long pixel_coord_x = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(i, 0, 0);
            for ( int j = 0; j < directional_image_sum.logical_y_dimension; j++ ) {
                long pixel_coord_xy                               = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                directional_image_sum.real_values[pixel_coord_xy] = directional_image_sum.real_values[pixel_coord_x];
                pixel_counter++;
            }
            pixel_counter += directional_image_sum.padding_jump_value;
        }

        directional_image_sum.DivideByConstant(directional_image_sum.logical_y_dimension);
    }

    // copy to current iamge
    current_image->CopyFrom(&directional_image_sum);
    directional_image_sum.Deallocate( );
}

// void AverageRotationallyFourierSlices(Image* img) {

//     // Ensure Fourier space
//     bool was_real = img->is_in_real_space;
//     if ( was_real ) {
//         img->ForwardFFT( );
//     }

//     const long nx = img->logical_x_dimension;
//     const long ny = img->logical_y_dimension;
//     const long nz = img->logical_z_dimension;

//     const long cx = img->physical_address_of_box_center_x;
//     const long cy = img->physical_address_of_box_center_y;

//     const float dkx = img->fourier_voxel_size_x;
//     const float dky = img->fourier_voxel_size_y;

//     const long max_bins = nx;

//     std::vector<std::complex<float>> ring_sum(max_bins);
//     std::vector<float>               ring_weight(max_bins);
//     std::vector<std::complex<float>> ring_avg(max_bins);

//     // Loop over Z slices
//     for ( long kz = 0; kz < nz; kz++ ) {

//         std::fill(ring_sum.begin( ), ring_sum.end( ), std::complex<float>(0, 0));
//         std::fill(ring_weight.begin( ), ring_weight.end( ), 0.0f);

//         // --- Accumulate ---
//         for ( long ky = 0; ky < ny; ky++ ) {

//             float ky_coord = img->ReturnFourierLogicalCoordGivenPhysicalCoord_Y(ky) * dky;
//             float ky_sq    = ky_coord * ky_coord;

//             for ( long kx = 0; kx <= img->physical_upper_bound_complex_x; kx++ ) {

//                 if ( img->FourierComponentIsExplicitHermitianMate(kx, ky, kz) )
//                     continue;

//                 float kx_coord = kx * dkx;
//                 float radius   = std::sqrt(kx_coord * kx_coord + ky_sq);

//                 float bin_f = radius / dkx;
//                 long  bin   = long(bin_f);

//                 long addr = img->ReturnFourier1DAddressFromPhysicalCoord(kx, ky, kz);

//                 if ( bin + 1 < max_bins ) {
//                     float t = bin_f - float(bin);
//                     ring_sum[bin] += img->complex_values[addr] * (1.0f - t);
//                     ring_sum[bin + 1] += img->complex_values[addr] * t;
//                     ring_weight[bin] += (1.0f - t);
//                     ring_weight[bin + 1] += t;
//                 }
//             }
//         }

//         // --- Normalize ---
//         for ( long b = 0; b < max_bins; b++ ) {
//             if ( ring_weight[b] > 0.0f )
//                 ring_avg[b] = ring_sum[b] / ring_weight[b];
//             else
//                 ring_avg[b] = std::complex<float>(0, 0);
//         }

//         // --- Write back ---
//         for ( long ky = 0; ky < ny; ky++ ) {

//             float ky_coord = img->ReturnFourierLogicalCoordGivenPhysicalCoord_Y(ky) * dky;
//             float ky_sq    = ky_coord * ky_coord;

//             for ( long kx = 0; kx <= img->physical_upper_bound_complex_x; kx++ ) {

//                 if ( img->FourierComponentIsExplicitHermitianMate(kx, ky, kz) )
//                     continue;

//                 float kx_coord = kx * dkx;
//                 float radius   = std::sqrt(kx_coord * kx_coord + ky_sq);

//                 float bin_f = radius / dkx;
//                 long  bin   = long(bin_f);

//                 std::complex<float> val(0, 0);

//                 if ( bin + 1 < max_bins ) {
//                     float t = bin_f - float(bin);
//                     val     = ring_avg[bin] * (1.0f - t) + ring_avg[bin + 1] * t;
//                 }

//                 long addr                 = img->ReturnFourier1DAddressFromPhysicalCoord(kx, ky, kz);
//                 img->complex_values[addr] = val;
//             }
//         }
//     }

//     // Restore original domain
//     if ( was_real ) {
//         img->BackwardFFT( );
//     }
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FourierRotateReconstruct.cpp
// Paste into your project and include "image.h" (or appropriate header) above this file.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// include your image header
// #include "image.h"

// ------------------------------------------------------------------
// Helper: read complex value from Image using ReturnComplexPixelFromLogicalCoord
// ------------------------------------------------------------------
// We assume image.h declares a function like:
//    ComplexPixel ReturnComplexPixelFromLogicalCoord(int x, int y, int z);
// where ComplexPixel contains fields .r and .i (or similar).
// This helper uses that function to obtain a complex value.
// If your ComplexPixel fields are named differently adjust the .r/.i accessors.
//
// If ReturnComplexPixelFromLogicalCoord returns by pointer or uses another naming,
// replace the implementation below appropriately.

// ------------------------------------------------------------------
// Helper: write complex value into Image complex_values (default)
// ------------------------------------------------------------------
// Most Image implementations expose complex_values[] as std::complex<float> *.
// We'll write into that array using the same indexing function used for real arrays.
// If your code exposes a setter (e.g. SetComplexPixelFromLogicalCoord), swap the body below.

// static void WriteComplexToImage2D(Image* img, int x, int y, const ComplexD& val) {
//     long addr = img->ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
//     // Default: assume complex_values is std::complex<float>*
//     img->complex_values[addr] = ComplexF(float(val.real( )), float(val.imag( )));
//     // If your image has a different API, replace the write above (e.g. img->SetComplexPixelFromLogicalCoord(x,y,cp);)
// }

// ------------------------------------------------------------------
// Bilinear deposit into complex 2D plane
// ------------------------------------------------------------------
static void BilinearAddComplex2D(std::vector<ComplexD>& plane, int NX, int NY,
                                 double x, double y, const ComplexD& val) {
    int    x0 = int(floor(x));
    int    y0 = int(floor(y));
    double dx = x - double(x0);
    double dy = y - double(y0);

    double w00 = (1.0 - dx) * (1.0 - dy);
    double w10 = dx * (1.0 - dy);
    double w01 = (1.0 - dx) * dy;
    double w11 = dx * dy;

    auto add_safe = [&](int xi, int yi, double w) {
        if ( xi >= 0 && xi < NX && yi >= 0 && yi < NY && w != 0.0 ) {
            size_t idx = size_t(yi) * size_t(NX) + size_t(xi);
            plane[idx] += val * w;
        }
    };

    add_safe(x0, y0, w00);
    add_safe(x0 + 1, y0, w10);
    add_safe(x0, y0 + 1, w01);
    add_safe(x0 + 1, y0 + 1, w11);
}

// ------------------------------------------------------------------
// Main: rotate 1D spectrum into kz=0 plane and inverse-FFT it
// proj2D: input image created by sum_image_direction (2D, identical rows assumed)
// out2D: optional output (if null, will overwrite proj2D with the reconstructed slice)
// nAngles: number of angular samples to use when rotating in k-space (e.g. 360)
// outSize: desired output plane size (Nx = Ny). If <=0 uses proj2D width.
//
void FourierRotateToCentralSlice(Image* proj2D, Image* out2D = nullptr,
                                 int nAngles = 360, int outSize = -1) {
    if ( ! proj2D )
        return;

    // sizes
    int W = proj2D->logical_x_dimension;
    int H = proj2D->logical_y_dimension;
    if ( outSize <= 0 )
        outSize = W;
    int NX = outSize;
    int NY = outSize;

    // copy projection and forward FFT (we assume ForwardFFT fills complex_values)
    Image projFFT;
    projFFT.CopyFrom(proj2D);
    projFFT.ForwardFFT( ); // must exist in image.cpp

    // extract 1D complex spectrum along x (pick y=0)
    std::vector<ComplexD> Fk(W);
    for ( int x = 0; x < W; ++x ) {
        int center_x = W / 2;
        Fk[x]        = ReadComplexFromImage2D(&projFFT, x - center_x, 0);
    }

    // prepare k-plane
    std::vector<ComplexD> kplane(size_t(NX) * size_t(NY));
    std::fill(kplane.begin( ), kplane.end( ), ComplexD(0.0, 0.0));

    // centers and scaling
    double kcx        = (NX - 1) / 2.0;
    double kcy        = (NY - 1) / 2.0;
    double src_center = (W - 1) / 2.0;
    double src_half   = src_center;
    double dst_half   = std::min(NX, NY) / 2.0;
    double scale_rad  = (src_half > 0.0) ? (dst_half / src_half) : 1.0;

    // rotate and deposit
    for ( int m = 0; m < nAngles; ++m ) {
        double phi    = 2.0 * M_PI * double(m) / double(nAngles);
        double cosphi = cos(phi), sinphi = sin(phi);

        for ( int kidx = 0; kidx < W; ++kidx ) {
            double krho_src = double(kidx) - src_center;
            double krho_dst = krho_src * scale_rad;
            double kx_f     = kcx + krho_dst * cosphi;
            double ky_f     = kcy + krho_dst * sinphi;

            BilinearAddComplex2D(kplane, NX, NY, kx_f, ky_f, Fk[kidx]);
        }
    }

    // normalize
    double invAngles = 1.0 / double(nAngles);
    for ( size_t i = 0; i < kplane.size( ); ++i )
        kplane[i] *= invAngles;

    // copy kplane into an Image complex field
    Image kplaneImage;
    kplaneImage.Allocate(NX, NY, true);
    kplaneImage.SetToConstant(0.0f);

    for ( int y = 0; y < NY; ++y ) {
        for ( int x = 0; x < NX; ++x ) {
            size_t   idx       = size_t(y) * size_t(NX) + size_t(x);
            ComplexD c         = kplane[idx];
            int      logical_x = x - NX / 2;
            int      logical_y = y - NY / 2;

            long addr = kplaneImage.ReturnFourier1DAddressFromLogicalCoord(logical_x, logical_y, 0);
            kplaneImage.complex_values[addr] =
                    std::complex<float>(float(c.real( )), float(c.imag( )));
        }
    }

    kplaneImage.is_in_real_space = false;

    // inverse FFT to get real central slice
    kplaneImage.BackwardFFT( ); // now kplaneImage.real_values holds the reconstructed slice

    // copy into destination
    Image* dest = out2D ? out2D : proj2D;
    dest->Allocate(NX, NY, true);
    dest->SetToConstant(0.0f);

    for ( int y = 0; y < NY; ++y ) {
        for ( int x = 0; x < NX; ++x ) {
            long addrSrc               = kplaneImage.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
            long addrDst               = dest->ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
            dest->real_values[addrDst] = kplaneImage.real_values[addrSrc];
        }
    }

    // clean up temporaries if required
    projFFT.Deallocate( );
    kplaneImage.Deallocate( );
}

///////////////////////////////////////////////////////////////////////////////

// Build a 3D volume from a 2D image by radializing each input row into a z-slice.
// - image: 2D input (logical_x_dimension x logical_y_dimension)
// - volume: output 3D; will be resized to (image->logical_x_dimension, image->logical_y_dimension, image->logical_y_dimension)
// The function preserves the FFT handling semantics from your original routine.
void Generate3DVolumeFromFilament2DClass(Image* image, Image* volume) {
    // Ensure input is real-space; same pattern as original
    bool input_image_in_fourier_space = false;
    if ( ! image->is_in_real_space ) {
        image->BackwardFFT( );
        input_image_in_fourier_space = true;
    }

    // Dimensions
    long W = image->logical_x_dimension; // width (columns)
    long H = image->logical_y_dimension; // height (rows)
    long D = H; // target depth equals number of rows so final volume will be W x H x H

    // Resize / initialize output volume to desired dims (keep padding value behavior consistent).
    // You may want to choose an appropriate fill value; using 0.0f here.
    volume->Resize(W, H, D, 0.0f);

    // Determine number_of_rings from center -> edge along X (this yields W/2 bins for even W)
    long number_of_rings = (W + 1) / 2; // 224 -> 112, 223 -> 112

    // allocate arrays
    auto ring_axis   = std::make_unique<float[]>(number_of_rings);
    auto ring_values = std::make_unique<float[]>(number_of_rings);
    auto ring_weight = std::make_unique<float[]>(number_of_rings);

    // centers (physical address of box center as used in original code)
    long central_x_pixel = image->physical_address_of_box_center_x;
    long central_y_pixel = image->physical_address_of_box_center_y;
    long central_z_pixel = 0L; // not used for building profiles (kept to mirror original)

    // Build ring_axis mapping bins -> radius
    // Here we map bin indices 0..number_of_rings-1 to radii 0 .. max_radius_along_x
    // max_radius_along_x should correspond to the center -> edge distance in x (consistent with row-based profiles).
    float max_radius_x = float(central_x_pixel); // typical case for even W: center sits between pixels, central_x_pixel equals W/2
    if ( number_of_rings > 1 ) {
        for ( long r = 0; r < number_of_rings; ++r ) {
            ring_axis[r] = 0.0f + r * (max_radius_x - 0.0f) / float(number_of_rings - 1);
        }
    }
    else {
        ring_axis[0] = 0.0f;
    }

    // Precompute denom for bin computation (guard against tiny denom)
    float bin_step = (number_of_rings > 1) ? (ring_axis[1] - ring_axis[0]) : 1.0f;
    if ( bin_step <= 0.0f )
        bin_step = 1.0f;

    // Helpers for center-left/right indices when averaging symmetric pixels
    // If width even: center sits between pixels; we average (center_left, center_right) for r=0.
    // If width odd: center is a pixel and r=0 uses that single pixel.
    bool width_even   = (W % 2 == 0);
    long center_left  = (width_even) ? (central_x_pixel - 1) : central_x_pixel;
    long center_right = central_x_pixel;

    // Edge bin index (last valid bin) - matches logic of using last bin as "edge value"
    long edge_bin_index = number_of_rings - 1;

    // iterate over each row in the input image; each row creates one radial profile and fills z slice = row index
    for ( long y = 0; y < H; ++y ) {
        // zero arrays for this row
        ZeroArray(ring_values.get( ), number_of_rings);
        ZeroArray(ring_weight.get( ), number_of_rings);

        // Build radial profile by averaging symmetric X pixels about center
        for ( long r = 0; r < number_of_rings; ++r ) {
            // compute left/right x indices
            long left_x, right_x;
            if ( width_even ) {
                left_x  = center_left - r; // e.g., for r=0 left=111
                right_x = center_right + r; // e.g., for r=0 right=112
            }
            else {
                left_x  = center_left - r;
                right_x = center_right + r;
            }

            // clamp to image bounds to be safe
            if ( left_x < 0 )
                left_x = 0;
            if ( right_x >= W )
                right_x = W - 1;

            // read values (assumes Image::ReturnReal1DAddressFromPhysicalCoord exists)
            long  addr_left  = image->ReturnReal1DAddressFromPhysicalCoord(left_x, y, 0);
            long  addr_right = image->ReturnReal1DAddressFromPhysicalCoord(right_x, y, 0);
            float v_left     = image->real_values[addr_left];
            float v_right    = image->real_values[addr_right];

            if ( left_x == right_x ) {
                // odd-width center pixel
                ring_values[r] += v_left;
                ring_weight[r] += 1.0f;
            }
            else {
                ring_values[r] += (v_left + v_right); // accumulate sum
                ring_weight[r] += 2.0f; // two samples contributed
            }
        }

        // Normalize radial profile (divide by weights)
        for ( long r = 0; r < number_of_rings; ++r ) {
            if ( ring_weight[r] != 0.0f )
                ring_values[r] /= ring_weight[r];
            else
                ring_values[r] = 0.0f; // fallback if no contributors (shouldn't happen for normal images)
        }

        // Fill z-slice z = y in the output volume using the radial profile
        long z = y; // mapping: row y -> z slice z
        for ( long yy = 0; yy < volume->logical_y_dimension; ++yy ) {
            for ( long xx = 0; xx < volume->logical_x_dimension; ++xx ) {
                // compute radius from center in the slice plane (x,y)
                float dx     = float(central_x_pixel - xx);
                float dy     = float(central_y_pixel - yy);
                float radius = sqrtf(dx * dx + dy * dy);

                // find bin
                long index_of_bin = long((radius - ring_axis[0]) / bin_step);

                // address in volume to write
                long out_addr = volume->ReturnReal1DAddressFromPhysicalCoord(xx, yy, z);

                if ( index_of_bin >= edge_bin_index ) {
                    // outside radial profile coverage -> use edge value
                    volume->real_values[out_addr] = ring_values[edge_bin_index];
                }
                else if ( index_of_bin < 0 ) {
                    // below zero radius (shouldn't happen because ring_axis[0]==0) -> use first bin
                    volume->real_values[out_addr] = ring_values[0];
                }
                else {
                    // linear interpolation between adjacent radial bins
                    float difference = (radius - ring_axis[index_of_bin]) / (ring_axis[index_of_bin + 1] - ring_axis[index_of_bin]);
                    // guard in case of any numerical problems
                    if ( difference < 0.0f )
                        difference = 0.0f;
                    if ( difference > 1.0f )
                        difference = 1.0f;
                    volume->real_values[out_addr] = (ring_values[index_of_bin] * (1.0f - difference)) + (ring_values[index_of_bin + 1] * difference);
                }
            }
        }
    } // end for each row

    // If the input image had been in Fourier space, restore (same behavior as your original function)
    if ( input_image_in_fourier_space ) {
        // Note: original code forward-fft'd the 'volume' variable; we replicate that here.
        volume->ForwardFFT( );
    }
}
