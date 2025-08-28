#include "../../core/core_headers.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iomanip> // for std::fixed and std::setprecision

class
        find_tube_diameters : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

std::vector<float>  sum_image_columns(Image* current_image);
void                save_all_columns_sum_to_file(const std::vector<std::vector<float>>& all_columns_sum, const std::string& filename);
std::pair<int, int> findOuterTubeEdges(const std::vector<float>& cols, float min_tube_diameter, float max_tube_diameter);

IMPLEMENT_APP(find_tube_diameters)

// override the DoInteractiveUserInput

void find_tube_diameters::DoInteractiveUserInput( ) {
    wxString input_images;
    float    pixel_size;
    float    min_tube_diameter = 0.0;
    float    max_tube_diameter;
    float    outer_mask_radius   = 0;
    float    low_pass_resolution = 50.0;

    wxString output_peaks_filename;
    wxString output_diameters_filename;

    int max_threads;

    UserInput* my_input       = new UserInput("Find Tube Diameters", 1.00);
    input_images              = my_input->GetFilenameFromUser("Input images file name", "Filename of helical tube stack", "helical_stack.mrc", true);
    pixel_size                = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
    min_tube_diameter         = my_input->GetFloatFromUser("Minimum tube diameter in pixels", "The expected minimum tube diameter, anything below that diameter will be discarded", "0.0", 0.0);
    max_tube_diameter         = my_input->GetFloatFromUser("Maximum tube diameter in pixels", "The expected maximum tube diameter, anything above that diameter will be discarded", "512.0", 0.0);
    outer_mask_radius         = my_input->GetFloatFromUser("Outer mask radius for masking the images during tube alignment (pixels)", "Outer mask radius to use when searching and aligning tubes in pixels, zero mean no masking should be applied", "0", 0);
    low_pass_resolution       = my_input->GetFloatFromUser("Resolution limit for low pass filtering", "Resolution limit for low pass filter. Only this resolution or worse information will be retained", "50.0", 0.0);
    output_peaks_filename     = my_input->GetFilenameFromUser("Output peaks file name", "Filename of the peaks ", "peaks_output.txt", false);
    output_diameters_filename = my_input->GetFilenameFromUser("Output diameters file name", "Filename of the saved diameters ", "diameters_output.txt", false);

#ifdef _OPENMP
    max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "when threading, what is the max threads to run", "1", 1);
#else
    max_threads = 1;
#endif

    delete my_input;

    my_current_job.Reset(10);
    my_current_job.ManualSetArguments("tffffftti", input_images.ToUTF8( ).data( ),
                                      pixel_size,
                                      min_tube_diameter,
                                      max_tube_diameter,
                                      outer_mask_radius,
                                      low_pass_resolution,
                                      output_peaks_filename.ToUTF8( ).data( ),
                                      output_diameters_filename.ToUTF8( ).data( ),
                                      max_threads); //update_star_file, input_star_filename.ToUTF8( ).data( ),
}

// override the do calculation method which will be what is actually run..

bool find_tube_diameters::DoCalculation( ) {
    wxString input_images              = my_current_job.arguments[0].ReturnStringArgument( );
    float    pixel_size                = my_current_job.arguments[1].ReturnFloatArgument( );
    float    min_tube_diameter         = my_current_job.arguments[2].ReturnFloatArgument( );
    float    max_tube_diameter         = my_current_job.arguments[3].ReturnFloatArgument( );
    float    outer_mask_radius         = my_current_job.arguments[4].ReturnFloatArgument( );
    float    low_pass_resolution       = my_current_job.arguments[5].ReturnFloatArgument( );
    wxString output_peaks_filename     = my_current_job.arguments[6].ReturnStringArgument( );
    wxString output_diameters_filename = my_current_job.arguments[7].ReturnStringArgument( );
    int      max_threads               = my_current_job.arguments[8].ReturnIntegerArgument( );

    MRCFile            my_input_images(input_images.ToStdString( ), false);
    long               number_of_input_images = my_input_images.ReturnNumberOfSlices( );
    std::vector<float> x_shift_column(number_of_input_images, 0.0f);
    Image              current_image;
    int                x_dim;
    int                y_dim;
    current_image.ReadSlice(&my_input_images, 1);
    float center_peak_index = current_image.logical_y_dimension / 2;
    x_dim                   = current_image.logical_x_dimension;
    y_dim                   = current_image.logical_y_dimension;

    std::vector<std::vector<float>> all_columns_sum(number_of_input_images, std::vector<float>(x_dim, 0.0)); // why I am saving those values?????
    std::vector<float>              all_diameters(number_of_input_images, 0.0f);

    // Open the diameter file in write mode
    std::ofstream file(output_diameters_filename.ToStdString( ));
    if ( ! file.is_open( ) ) {
        std::cerr << "Error: Could not open diameters_output.txt\n";
        //return;
    }

    file << std::fixed << std::setprecision(2); // Optional: set float precision
    file << "image_index, diameter\n";

    std::ofstream peak_file(output_peaks_filename.ToStdString( )); // Open file once

    if ( ! peak_file.is_open( ) ) {
        std::cerr << "Error: Could not open peaks_output.txt\n";
        //return;
    }

    peak_file << std::fixed << std::setprecision(2); // Optional: set float precision
    peak_file << "image_index, peak_one_value, peak_two_value\n";

    std::vector<std::pair<int, std::pair<int, int>>> results(number_of_input_images);

    MRCFile my_output("current_image.mrc", true);
    if ( ! my_output.IsOpen( ) ) {
        my_output.OpenFile("current_image.mrc", true);
        if ( ! my_output.IsOpen( ) ) {
            wxPrintf("ERROR: Could not open '%s' for writing\n", "current_image.mrc");
            DEBUG_ABORT;
        }
    }
    my_output.my_header.SetNumberOfImages(number_of_input_images);
    my_output.my_header.SetDimensionsImage(x_dim, y_dim);
    my_output.SetPixelSize(pixel_size);
    my_output.WriteHeader( );
    my_output.rewrite_header_on_close = true;

    ProgressBar* my_progress = new ProgressBar(number_of_input_images);

#pragma omp parallel for schedule(dynamic, 1) num_threads(max_threads) default(none) shared(my_input_images, number_of_input_images, low_pass_resolution, max_threads, x_dim, y_dim, all_diameters, results, min_tube_diameter, max_tube_diameter, \
                                                                                            pixel_size, my_progress, outer_mask_radius, min_tube_diameter, center_peak_index, x_shift_column, all_columns_sum, peak_file, my_output) private(current_image)

    for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
        // wxPrintf("Image number is %li\n\n", image_counter + 1);
        current_image.Allocate(x_dim, y_dim, true);
        current_image.SetToConstant(0.0);
#pragma omp critical
        current_image.ReadSlice(&my_input_images, image_counter + 1);
        current_image.Normalize( );
        if ( outer_mask_radius != 0 ) {
            current_image.CircleMask(outer_mask_radius);
        }
#pragma omp critical
        current_image.WriteSlice(&my_output, image_counter + 1);
        current_image.ForwardFFT( );
        current_image.ZeroCentralPixel( );
        current_image.GaussianLowPassFilter((pixel_size * 2) / low_pass_resolution);
        current_image.BackwardFFT( );

#pragma omp critical
        all_columns_sum[image_counter] = sum_image_columns(&current_image);

        // wxPrintf("Column values for image %li are ", image_counter);
        // for ( auto val : all_columns_sum[image_counter] ) {
        //     wxPrintf("%f, ", val);
        // }
        // wxPrintf("\n\n");
        // find the outer edges of the protein or tube
        auto [peak_one_column_sum, peak_two_column_sum] = findOuterTubeEdges(all_columns_sum[image_counter], min_tube_diameter, max_tube_diameter);

        results[image_counter] = {image_counter, {peak_one_column_sum, peak_two_column_sum}};

        float tube_diameter          = std::abs((peak_one_column_sum - peak_two_column_sum));
        all_diameters[image_counter] = tube_diameter;
        current_image.Deallocate( );

        if ( is_running_locally == true && ReturnThreadNumberOfCurrentThread( ) == 0 )

            my_progress->Update(image_counter + 1);
    }
    delete my_progress;

    //save the peaks to the output file
    if ( peak_file.is_open( ) ) {
        for ( auto& r : results ) {
            peak_file << r.first << ", " << r.second.first << ", " << r.second.second << "\n";
        }
        peak_file.close( );
    }

    // Check if the all diameters file is open
    if ( file.is_open( ) ) {
        for ( size_t i = 0; i < all_diameters.size( ); ++i ) {
            file << i << ", " << all_diameters[i] << '\n';
        }
        file.close( );
    }
    // save column sum output file
    // save_all_columns_sum_to_file(all_columns_sum, "column_sums_output.txt");

    return true;
}

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

void save_all_columns_sum_to_file(const std::vector<std::vector<float>>& all_columns_sum, const std::string& filename) {
    std::ofstream out_file(filename);
    if ( ! out_file.is_open( ) ) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        //return;
    }

    out_file << std::fixed << std::setprecision(2); // Set float precision to 2 decimal places

    for ( const auto& row : all_columns_sum ) {
        for ( size_t i = 0; i < row.size( ); ++i ) {
            out_file << row[i];
            if ( i < row.size( ) - 1 )
                out_file << ", ";
        }
        out_file << '\n';
    }

    out_file.close( );
}

// Detects the two strongest outer-edge peaks in a 1D intensity profile.
// Returns indices of the best peak pair (sorted low->high), or an empty vector if none found.
std::pair<int, int> findOuterTubeEdges(const std::vector<float>& cols, float min_tube_diameter, float max_tube_diameter) {
    int n = cols.size( );
    if ( n < 3 )
        return {-1, -1}; // need at least 3 points to form a peak

    // 1) Normalize the 1D profile
    float              minVal = *std::min_element(cols.begin( ), cols.end( ));
    float              maxVal = *std::max_element(cols.begin( ), cols.end( ));
    std::vector<float> norm(n);
    for ( int i = 0; i < n; ++i )
        norm[i] = cols[i] - minVal;

    float normMax = *std::max_element(norm.begin( ), norm.end( ));
    if ( normMax <= 0.0f )
        return {-1, -1};

    // Inverted profile for negative peaks
    std::vector<float> normInv(n);
    for ( int i = 0; i < n; ++i )
        normInv[i] = normMax - norm[i];

    // 2) Detect peaks
    std::vector<std::pair<int, float>> posPeaks;
    std::vector<std::pair<int, float>> negPeaks;

    for ( int i = 1; i < n - 1; ++i ) {
        if ( norm[i] > norm[i - 1] && norm[i] > norm[i + 1] ) {
            posPeaks.emplace_back(i, norm[i]);
        }
        if ( normInv[i] > normInv[i - 1] && normInv[i] > normInv[i + 1] ) {
            negPeaks.emplace_back(i, normInv[i]);
        }
    }
    // // debugging and printing the scores
    // std::cerr << "posPeaks (idx,val): ";
    // for ( auto& p : posPeaks )
    //     std::cerr << "(" << p.first << "," << p.second << ") ";
    // std::cerr << "\n";
    // std::cerr << "negPeaks (idx,val): ";
    // for ( auto& p : negPeaks )
    //     std::cerr << "(" << p.first << "," << p.second << ") ";
    // std::cerr << "\n";

    // helper function to find the best pair of peaks based on their height and distance between peaks
    auto bestPair = [&](const std::vector<std::pair<int, float>>& peaks)
            -> std::pair<float, std::pair<int, int>> {
        float               bestScore = -std::numeric_limits<float>::infinity( );
        std::pair<int, int> bestIdx   = {-1, -1};

        // Adding gap penalty and out of range penalty so that we would favor more the peaks within the range, but also if nothing was found within range, out of range peaks are saved and returned
        const float IDEAL_GAP           = min_tube_diameter;
        const float GAP_PENALTY         = 0.1f; // e.g. 0.1 points lost per pixel of gap deviation
        const float OUT_OF_RANGE_FACTOR = 10.0f; // scale factor for out-of-range penalty- changed that from 2 to 10 to heavily penalize out of range to favor in range more

        for ( size_t a = 0; a < peaks.size( ); ++a ) {
            for ( size_t b = a + 1; b < peaks.size( ); ++b ) {
                int i   = peaks[a].first;
                int j   = peaks[b].first;
                int gap = j - i;

                float sumAmp = peaks[a].second + peaks[b].second;
                float score  = sumAmp - GAP_PENALTY * std::fabs(gap - IDEAL_GAP);

                // scale penalty by how far out of range the gap is
                if ( gap < min_tube_diameter ) {
                    score -= OUT_OF_RANGE_FACTOR * (min_tube_diameter - gap);
                }
                else if ( gap > max_tube_diameter ) {
                    score -= OUT_OF_RANGE_FACTOR * (gap - max_tube_diameter);
                }

                if ( score > bestScore ) {
                    bestScore = score;
                    bestIdx   = {i, j};
                }
            }
        }

        return std::make_pair(bestScore, bestIdx);
    };

    // 3) Find best pair among positive peaks and among negative peaks.
    auto [scorePos, bestPos] = bestPair(posPeaks);
    auto [scoreNeg, bestNeg] = bestPair(negPeaks);

    std::pair<int, int> bestPairIdx = {-1, -1};

    // 4) If no valid pairs exist at all, return -1
    if ( scorePos == -std::numeric_limits<float>::infinity( ) &&
         scoreNeg == -std::numeric_limits<float>::infinity( ) ) {
        return {-1, -1};
    }

    // 5) keeping the values of the best negative peaks as reference
    bestPairIdx     = bestNeg;
    float bestScore = scoreNeg;

    // find the highest negative peaks within the range of the expected diameter
    // then find the positive peak before the first negative peak and the positive peak after the second negative peak and those should be the outer edges
    if ( bestNeg.first != -1 && bestNeg.second != -1 ) {
        int iNeg = bestNeg.first;
        int jNeg = bestNeg.second;
        if ( iNeg > jNeg )
            std::swap(iNeg, jNeg); // enforce left->right

        // Find last positive BEFORE iNeg
        int   posBefore = -1;
        float ampBefore = 0;
        for ( auto it = posPeaks.rbegin( ); it != posPeaks.rend( ); ++it ) {
            if ( it->first < iNeg ) {
                posBefore = it->first;
                ampBefore = it->second;
                break;
            }
        }

        // Find first positive AFTER jNeg
        int   posAfter = -1;
        float ampAfter = 0;
        for ( auto& p : posPeaks ) {
            if ( p.first > jNeg ) {
                posAfter = p.first;
                ampAfter = p.second;
                break;
            }
        }
        const float IDEAL_GAP           = min_tube_diameter;
        const float GAP_PENALTY         = 0.1f; // e.g. 0.1 points lost per pixel of gap deviation
        const float OUT_OF_RANGE_FACTOR = 10.0f; // scale factor for out-of-range penalty

        // Step 3: Only refine if both positives exist and are ordered
        if ( posAfter != -1 && posBefore != -1 && posAfter < posBefore ) {
            int   gap    = posBefore - posAfter;
            float sumAmp = ampAfter + ampBefore;
            float score  = sumAmp - GAP_PENALTY * std::fabs(gap - IDEAL_GAP);

            if ( gap < min_tube_diameter )
                score -= OUT_OF_RANGE_FACTOR * (min_tube_diameter - gap);
            else if ( gap > max_tube_diameter )
                score -= OUT_OF_RANGE_FACTOR * (gap - max_tube_diameter);

            // Step 4: Replace if adjacency score is better
            if ( score > bestScore ) {
                bestScore   = score;
                bestPairIdx = {posAfter, posBefore};
            }
        }
    }
    // Final: enforce sorted order before returning
    if ( bestPairIdx.first > bestPairIdx.second )
        std::swap(bestPairIdx.first, bestPairIdx.second);

    return std::make_pair(bestPairIdx.first, bestPairIdx.second);
}
