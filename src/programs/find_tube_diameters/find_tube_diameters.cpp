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
float               sum_image_columns_float(Image* current_image);
std::pair<int, int> find_column_sum_peaks(Image* current_image, float min_gap = 0.0, bool find_positive_peaks = true, bool find_negative_peaks = true);
std::pair<int, int> find_column_sum_peaks_new(Image* current_image, float min_gap);
void                save_all_columns_sum_to_file(const std::vector<std::vector<float>>& all_columns_sum, const std::string& filename);
void                sum_image_direction(Image* current_image, int dim);
std::pair<int, int> findOuterTubeEdges(const std::vector<float>& cols, float min_tube_diameter, float max_tube_diameter, bool find_positive_peaks = true, bool find_negative_peaks = true);

IMPLEMENT_APP(find_tube_diameters)

// override the DoInteractiveUserInput

void find_tube_diameters::DoInteractiveUserInput( ) {
    wxString input_images;
    float    pixel_size;
    float    min_tube_diameter = 0.0;
    float    max_tube_diameter;
    float    outer_mask_radius   = 0;
    float    low_pass_resolution = 50.0;
    bool     find_positive_peaks = true;
    bool     find_negative_peaks = true;
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
    find_positive_peaks       = my_input->GetYesNoFromUser("Find the positive peaks?", "If yes, will find the highest positive peaks.", "Yes");
    find_negative_peaks       = my_input->GetYesNoFromUser("Find the negative peaks?", "If yes, will find the highest negative peaks.", "Yes");
    output_peaks_filename     = my_input->GetFilenameFromUser("Output peaks file name", "Filename of the peaks ", "peaks_output.txt", false);
    output_diameters_filename = my_input->GetFilenameFromUser("Output diameters file name", "Filename of the saved diameters ", "diameters_output.txt", false);

#ifdef _OPENMP
    max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "when threading, what is the max threads to run", "1", 1);
#else
    max_threads = 1;
#endif

    delete my_input;

    my_current_job.Reset(12);
    my_current_job.ManualSetArguments("tfffffbbtti", input_images.ToUTF8( ).data( ),
                                      pixel_size,
                                      min_tube_diameter,
                                      max_tube_diameter,
                                      outer_mask_radius,
                                      low_pass_resolution,
                                      find_positive_peaks,
                                      find_negative_peaks,
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
    bool     find_positive_peaks       = my_current_job.arguments[6].ReturnBoolArgument( );
    bool     find_negative_peaks       = my_current_job.arguments[7].ReturnBoolArgument( );
    wxString output_peaks_filename     = my_current_job.arguments[8].ReturnStringArgument( );
    wxString output_diameters_filename = my_current_job.arguments[9].ReturnStringArgument( );
    int      max_threads               = my_current_job.arguments[10].ReturnIntegerArgument( );

    MRCFile            my_input_images(input_images.ToStdString( ), false);
    long               number_of_input_images = my_input_images.ReturnNumberOfSlices( );
    std::vector<float> best_sum_column(number_of_input_images, 0.0f);
    std::vector<float> x_shift_column(number_of_input_images, 0.0f);
    Image              current_image;
    Image              final_image;
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

    ProgressBar* my_progress = new ProgressBar(number_of_input_images);
#pragma omp parallel for schedule(dynamic, 1) num_threads(max_threads) default(none) shared(my_input_images, best_sum_column, number_of_input_images, low_pass_resolution, max_threads, x_dim, y_dim, all_diameters, results, find_positive_peaks, find_negative_peaks, min_tube_diameter, max_tube_diameter, \
                                                                                            pixel_size, my_progress, outer_mask_radius, min_tube_diameter, center_peak_index, x_shift_column, all_columns_sum, peak_file) private(current_image, final_image)

    for ( long image_counter = 0; image_counter < number_of_input_images; image_counter++ ) {
        final_image.Allocate(x_dim, y_dim, true);
        final_image.SetToConstant(0.0);
#pragma omp critical
        final_image.ReadSlice(&my_input_images, image_counter + 1);
        final_image.Normalize( );
        if ( outer_mask_radius != 0 ) {
            final_image.CircleMask(outer_mask_radius);
        }

        final_image.ForwardFFT( );
        final_image.ZeroCentralPixel( );
        final_image.GaussianLowPassFilter((pixel_size * 2) / low_pass_resolution);
        final_image.BackwardFFT( );

#pragma omp critical
        all_columns_sum[image_counter] = sum_image_columns(&final_image);
        // calculate the required x shift to center the tubes
        //auto [peak_one_column_sum, peak_two_column_sum] = find_column_sum_peaks(&final_image, min_tube_diameter, find_positive_peaks, find_negative_peaks);
        auto [peak_one_column_sum, peak_two_column_sum] = findOuterTubeEdges(all_columns_sum[image_counter], min_tube_diameter, max_tube_diameter, find_positive_peaks, find_negative_peaks);
        // New sum peaks will find the outer edges of tubes not the inner ones
        //auto [peak_one_column_sum, peak_two_column_sum] = find_column_sum_peaks_new(&final_image, min_tube_diameter);
        results[image_counter] = {image_counter, {peak_one_column_sum, peak_two_column_sum}};

        // The next line not needed
        float tube_center_column_sum          = std::abs(peak_one_column_sum - peak_two_column_sum) / 2;
        float distance_from_center_column_sum = -((peak_one_column_sum + peak_two_column_sum) / 2 - center_peak_index);

        x_shift_column[image_counter] = distance_from_center_column_sum; // the x-shift needed
        float tube_diameter           = std::abs((peak_one_column_sum - peak_two_column_sum)); // * pixel_size
        all_diameters[image_counter]  = tube_diameter;
        //final_image.PhaseShift(x_shift_column[image_counter], 0.0, 0.0);

        final_image.Deallocate( );

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
    save_all_columns_sum_to_file(all_columns_sum, "column_sums_output.txt");

    return true;
}

float sum_image_columns_float(Image* current_image) { // Tim's method to calculate the best rotation based on the row sum
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

    float max_value = *std::max_element(column_sum.begin( ), column_sum.end( ),
                                        [](float a, float b) { return std::abs(a) < std::abs(b); });

    return abs(max_value);
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

std::pair<int, int> find_column_sum_peaks(Image* current_image, float min_gap, bool find_positive_peaks, bool find_negative_peaks) {
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
    // Step 1: Separate positive and negative local peaks
    std::vector<int> positive_peaks;
    std::vector<int> negative_peaks;

    for ( long i = 1; i < column_sum.size( ) - 1; ++i ) {
        // Local maximum
        if ( column_sum[i] > column_sum[i - 1] && column_sum[i] > column_sum[i + 1] ) {
            positive_peaks.push_back(i);
        }
        // Local minimum
        if ( column_sum[i] < column_sum[i - 1] && column_sum[i] < column_sum[i + 1] ) {
            negative_peaks.push_back(i);
        }
    }

    // Step 2: Sort peaks by amplitude
    // Positive peaks: largest first
    std::sort(positive_peaks.begin( ), positive_peaks.end( ),
              [&column_sum](int i, int j) { return column_sum[i] > column_sum[j]; });
    // Negative peaks: smallest (most negative) first
    std::sort(negative_peaks.begin( ), negative_peaks.end( ),
              [&column_sum](int i, int j) { return column_sum[i] < column_sum[j]; });

    // Step 3: Helper function to find the best pair of peaks
    auto find_peak_pair_with_gap = [&](const std::vector<int>&   peaks_vec,
                                       const std::vector<float>& column_sum) -> std::pair<int, int> {
        if ( peaks_vec.size( ) < 2 )
            return {-1, -1}; // Not enough peaks to form a pair

        std::pair<int, int> best_pair  = {-1, -1};
        float               best_score = -1.0f; // Score = sum of peak amplitudes

        // Consider all possible pairs of peaks
        for ( size_t i = 0; i < peaks_vec.size( ) - 1; ++i ) {
            for ( size_t j = i + 1; j < peaks_vec.size( ); ++j ) {
                int gap = std::abs(peaks_vec[j] - peaks_vec[i]);

                // Only consider pairs separated by at least min_gap
                if ( gap >= min_gap ) {
                    // Optional: second peak must be >= 80% amplitude of first
                    //if ( column_sum[peaks_vec[j]] >= 0.8 * column_sum[peaks_vec[i]] ) {
                    // Score by total amplitude (or max amplitude)
                    float score = column_sum[peaks_vec[i]] + column_sum[peaks_vec[j]];
                    if ( score > best_score ) {
                        best_score = score;
                        best_pair  = {std::min(peaks_vec[i], peaks_vec[j]),
                                      std::max(peaks_vec[i], peaks_vec[j])};
                    }
                    //}
                }
            }
        }
        return best_pair;
    };

    // Step 4: Helper function to choose positive or negative peaks
    auto choose_peaks = [&](const std::pair<int, int>& pos_pair,
                            const std::pair<int, int>& neg_pair) -> std::pair<int, int> {
        if ( find_positive_peaks && pos_pair.first != -1 )
            return pos_pair;
        if ( find_negative_peaks && neg_pair.first != -1 )
            return neg_pair;
        return {-1, -1}; // No valid peaks
    };

    // Step 5: Handle empty cases and find best pairs
    if ( positive_peaks.size( ) < 2 && negative_peaks.size( ) < 2 ) {
        return {-1, -1}; // Not enough peaks
    }
    else {
        // Search positive and negative peaks
        auto pos_pair = find_peak_pair_with_gap(positive_peaks, column_sum);
        auto neg_pair = find_peak_pair_with_gap(negative_peaks, column_sum);
        // Choose the pair to return
        return choose_peaks(pos_pair, neg_pair);
    }

    //     // separated the peaks to positive and negative here to make sure we find two positive or two negative
    //     std::vector<int> positive_peaks;
    //     std::vector<int> negative_peaks;

    //     for ( long i = 1; i < column_sum.size( ) - 1; ++i ) {
    //         if ( column_sum[i] > column_sum[i - 1] && column_sum[i] > column_sum[i + 1] ) {
    //             // local maximum
    //             positive_peaks.push_back(i);
    //         }
    //         if ( column_sum[i] < column_sum[i - 1] && column_sum[i] < column_sum[i + 1] ) {
    //             // local minimum
    //             negative_peaks.push_back(i);
    //         }
    //     }

    //     std::sort(positive_peaks.begin( ), positive_peaks.end( ),
    //               [&column_sum](int i, int j) { return column_sum[i] > column_sum[j]; });

    //     std::sort(negative_peaks.begin( ), negative_peaks.end( ),
    //               [&column_sum](int i, int j) { return column_sum[i] < column_sum[j]; });

    //     // helper function 1
    //     auto find_peak_pair_with_gap = [&](const std::vector<int>& peaks_vec, const std::vector<float>& column_sum) -> std::pair<int, int> {
    //         if ( peaks_vec.size( ) < 2 )
    //             return {-1, -1};
    //         // best_pair stores the peak indices of the best candidate so far
    //         // best_gap_diff stores how far that candidate’s gap is from min_gap
    //         std::pair<int, int> best_pair     = {-1, -1};
    //         int                 best_gap_diff = INT_MAX;
    //         /*
    //         Checks your amplitude condition (column_sum[j] >= 0.8 * column_sum[i]).
    //         Checks if the current gap is closer to min_gap than any previous pair (gap_diff < best_gap_diff).
    //         If true, replace best_pair with the current pair.
    //         Also update best_gap_diff, so only better (closer to min_gap) gaps will replace it later.

    //         The loop considers all valid peak pairs.
    //         best_pair always holds the one with the smallest |gap - min_gap|.
    //         After finishing all iterations, the returned best_pair is the one closest to min_gap, not just the first valid pair
    //         */
    //         for ( size_t i = 0; i < peaks_vec.size( ) - 1; ++i ) {
    //             for ( size_t j = i + 1; j < peaks_vec.size( ); ++j ) {
    //                 int gap      = std::abs(peaks_vec[i] - peaks_vec[j]);
    //                 int gap_diff = std::abs(gap - min_gap);
    //                 if ( column_sum[peaks_vec[j]] >= 0.8 * column_sum[peaks_vec[i]] && gap_diff < best_gap_diff ) {
    //                     best_pair     = {std::min(peaks_vec[i], peaks_vec[j]), std::max(peaks_vec[i], peaks_vec[j])};
    //                     best_gap_diff = gap_diff;
    //                 }
    //             }
    //         }
    //         /*
    // The loop considers all valid peak pairs.
    // best_pair always holds the one with the smallest |gap - min_gap|.
    // After finishing all iterations, the returned best_pair is the one closest to min_gap, not just the first valid pair*/
    //         return best_pair;
    //     };
    //     // helper function 2
    //     auto choose_peaks = [&](const std::pair<int, int>& pos_pair, const std::pair<int, int>& neg_pair) -> std::pair<int, int> {
    //         if ( find_positive_peaks && pos_pair.first != -1 )
    //             return pos_pair;
    //         if ( find_negative_peaks && neg_pair.first != -1 )
    //             return neg_pair;
    //         return {-1, -1}; // if no valid peaks are found
    //     };
    //     // if the positive peaks and negative peaks are empty then return -1,-1
    //     if ( positive_peaks.size( ) < 2 && negative_peaks.size( ) < 2 ) {
    //         return std::make_pair(-1, -1); // not enough peaks
    //     }
    //     else {
    //         // Search positive and negative peaks
    //         auto pos_pair = find_peak_pair_with_gap(positive_peaks, column_sum);
    //         auto neg_pair = find_peak_pair_with_gap(negative_peaks, column_sum);
    //         // find the best that has a diameter closer to the minimum gap
    //         return choose_peaks(pos_pair, neg_pair);
    //     }

    //     // Initialize variables needed to loop over the peaks
    //     int  current_main_peak            = 0;
    //     int  highest_positive_peak        = peaks[current_main_peak];
    //     int  second_highest_positive_peak = -1;
    //     bool found_both_positive_peaks    = false;
    //     int  highest_negative_peak        = peaks[current_main_peak];
    //     int  second_highest_negative_peak = -1;
    //     bool found_both_negative_peaks    = false;

    //     // Loop to find the highest two positive peaks with a minimum gap
    //     // I added size -1 to avoid accessing an out of bound index
    //     for ( long index = 0; index < peaks.size( ) && ! found_both_positive_peaks; index++ ) {
    //         // if the absolute difference between the highest peak position and the next one (indices) is higher than or equal to the gap
    //         // then we have found both peaks
    //         // check if the value of the peak is positive
    //         // go over the positive peaks search
    //         //wxPrintf("%li\n", index);
    //         if ( peaks[index] > 0 ) {
    //             if ( std::abs(peaks[index + 1] - highest_positive_peak) >= min_gap ) { // removed the std::abs from this line
    //                 second_highest_positive_peak = peaks[index + 1];
    //                 if ( column_sum[second_highest_positive_peak] >= 0.8 * column_sum[highest_positive_peak] ) {
    //                     found_both_positive_peaks = true;
    //                 }
    //                 else { // if the absolute differences is < than the gap
    //                     current_main_peak++; // add one to the current_main_peak index and consider this the first highest peak
    //                     highest_positive_peak = peaks[current_main_peak];
    //                     index                 = 0; // Reset index
    //                 }
    //             }
    //         }
    //         else if ( peaks[index] < 0 ) {
    //             current_main_peak               = 0;
    //             std::vector<int> reversed_peaks = peaks;
    //             std::reverse(reversed_peaks.begin( ), reversed_peaks.end( ));
    //             if ( std::abs(peaks[index + 1] - highest_negative_peak) >= min_gap ) {
    //                 second_highest_negative_peak = peaks[index + 1];
    //                 if ( column_sum[second_highest_negative_peak] >= 0.8 * column_sum[highest_negative_peak] ) { // if the second highest peak is at least 80% tall from the highest
    //                     found_both_negative_peaks = true;
    //                 }
    //                 else { // if the absolute differences is < than the gap
    //                     current_main_peak++; // add one to the current_main_peak index and consider this the first highest peak
    //                     highest_negative_peak = peaks[current_main_peak];
    //                     index                 = 0; // Reset index
    //                 }
    //             }
    //         }

    //         float center_peak_index = current_image->logical_y_dimension / 2;
    //         if ( found_both_positive_peaks && ! found_both_negative_peaks ) {
    //             // save the peak with the minimum index to peak1
    //             int peak1 = std::min(highest_positive_peak, second_highest_positive_peak);
    //             // save the peak with the highest index to peak2
    //             int peak2 = std::max(highest_positive_peak, second_highest_positive_peak);
    //             // float tube_center = std::abs(peak1 - peak2)/2;
    //             // float distance_from_center = (peak1 + peak2)/2 - center_peak_index;
    //             // return -distance_from_center;
    //             return std::make_pair(peak1, peak2);
    //         }
    //         if ( ! found_both_positive_peaks && found_both_negative_peaks ) {
    //             // save the peak with the minimum index to peak1
    //             int peak1 = std::min(highest_negative_peak, second_highest_negative_peak);
    //             // save the peak with the highest index to peak2
    //             int peak2 = std::max(highest_negative_peak, second_highest_negative_peak);
    //             // float tube_center = std::abs(peak1 - peak2)/2;
    //             // float distance_from_center = (peak1 + peak2)/2 - center_peak_index;
    //             // return -distance_from_center;
    //             return std::make_pair(peak1, peak2);
    //         }
    //         if ( found_both_positive_peaks && found_both_negative_peaks ) {
    //             // check if the min gap in the highest positive peak is closer to the specified min gap by the user than the negative peak
    //             int min_positive_gap = std::abs(highest_positive_peak - second_highest_positive_peak);
    //             int min_negative_gap = std::abs(highest_negative_peak - second_highest_negative_peak);
    //             if ( std::abs(min_positive_gap - min_gap) <= std::abs(min_negative_gap - min_gap) ) {
    //                 // then use the positive peak information
    //                 // save the peak with the minimum index to peak1
    //                 int peak1 = std::min(highest_positive_peak, second_highest_positive_peak);
    //                 // save the peak with the highest index to peak2
    //                 int peak2 = std::max(highest_positive_peak, second_highest_positive_peak);
    //                 // float tube_center = std::abs(peak1 - peak2)/2;
    //                 // float distance_from_center = (peak1 + peak2)/2 - center_peak_index;
    //                 // return -distance_from_center;
    //                 return std::make_pair(peak1, peak2);
    //             }
    //             else {
    //                 // save the peak with the minimum index to peak1
    //                 int peak1 = std::min(highest_negative_peak, second_highest_negative_peak);
    //                 // save the peak with the highest index to peak2
    //                 int peak2 = std::max(highest_negative_peak, second_highest_negative_peak);
    //                 //float tube_center = std::abs(peak1 - peak2)/2;
    //                 //float distance_from_center = (peak1 + peak2)/2 - center_peak_index;
    //                 //return -distance_from_center;
    //                 return std::make_pair(peak1, peak2);
    //             }
    //         }
    //     }
    // }
}

void save_all_columns_sum_to_file(
        const std::vector<std::vector<float>>& all_columns_sum,
        const std::string&                     filename) {
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

std::pair<int, int> find_column_sum_peaks_new(Image* current_image, float min_gap) {
    std::vector<float> column_sum(current_image->logical_x_dimension, 0.0);

    for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
        for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
            long pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
            column_sum[i] += current_image->real_values[pixel_coord_xy];
        }
    }

    // Normalize column sum by subtracting the minimum
    float min_val = *std::min_element(column_sum.begin( ), column_sum.end( ));
    for ( auto& val : column_sum ) {
        val -= min_val;
    }

    // Invert profile to detect "valleys" (negative peaks)
    float              max_val = *std::max_element(column_sum.begin( ), column_sum.end( ));
    std::vector<float> inverted_column_sum(column_sum.size( ));
    for ( size_t i = 0; i < column_sum.size( ); ++i ) {
        inverted_column_sum[i] = max_val - column_sum[i];
    }

    // Find positive peaks
    std::vector<int> positive_peaks;
    for ( size_t i = 1; i < column_sum.size( ) - 1; ++i ) {
        if ( column_sum[i] > column_sum[i - 1] && column_sum[i] > column_sum[i + 1] )
            positive_peaks.push_back(i);
    }

    // Find negative peaks (from inverted profile)
    std::vector<int> negative_peaks;
    for ( size_t i = 1; i < inverted_column_sum.size( ) - 1; ++i ) {
        if ( inverted_column_sum[i] > inverted_column_sum[i - 1] && inverted_column_sum[i] > inverted_column_sum[i + 1] )
            negative_peaks.push_back(i);
    }

    // Pair positive and negative peaks and compute intensity differences
    int                                  center = column_sum.size( ) / 2;
    std::vector<std::pair<float, float>> diff_list; // pair of diff and average index

    for ( int pos_peak : positive_peaks ) {
        std::vector<int> valid_neg_peaks;
        for ( int neg_peak : negative_peaks ) {
            if ( (pos_peak < center && neg_peak < pos_peak) ||
                 (pos_peak > center && neg_peak > pos_peak) ) {
                valid_neg_peaks.push_back(neg_peak);
            }
        }

        if ( ! valid_neg_peaks.empty( ) ) {
            int   closest_neg_peak = *std::min_element(valid_neg_peaks.begin( ), valid_neg_peaks.end( ),
                                                       [pos_peak](int a, int b) {
                                                         return std::abs(a - pos_peak) < std::abs(b - pos_peak);
                                                     });
            float diff             = std::abs(column_sum[pos_peak] - column_sum[closest_neg_peak]);
            float avg_index        = 0.5f * (pos_peak + closest_neg_peak);
            diff_list.push_back(std::make_pair(diff, avg_index));
        }
    }

    // Sort the peak pairs by difference (descending)
    std::sort(diff_list.begin( ), diff_list.end( ),
              [](const std::pair<float, float>& a, const std::pair<float, float>& b) {
                  return a.first > b.first;
              });

    // Select best left and right peaks
    float left_peak = -1, right_peak = -1;
    for ( const auto& [diff, index] : diff_list ) {
        if ( left_peak != -1 && right_peak != -1 )
            break;
        if ( std::abs(index - center) < min_gap / 2 )
            continue;
        if ( index < center && left_peak == -1 )
            left_peak = index;
        else if ( index > center && right_peak == -1 )
            right_peak = index;
    }

    if ( left_peak == -1 )
        left_peak = 0;
    if ( right_peak == -1 )
        right_peak = 0;

    return std::make_pair(static_cast<int>(left_peak), static_cast<int>(right_peak));
}

void sum_image_direction(Image* current_image, int dim) {
    // image must be in real-space
    Image directional_image_sum;
    directional_image_sum.Allocate(current_image->logical_x_dimension, current_image->logical_y_dimension, true);
    directional_image_sum.SetToConstant(0.0);

    // x-direction
    if ( dim == 1 ) {

        long pixel_coord_y  = 0;
        long pixel_coord_xy = 0;
        long pixel_counter  = 0;

        // sum columns of my_image_sum (NxM) and store in array (1xN)
        for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
            for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
                pixel_coord_y  = current_image->ReturnReal1DAddressFromPhysicalCoord(0, j, 0);
                pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                directional_image_sum.real_values[pixel_coord_y] += current_image->real_values[pixel_coord_xy];
                pixel_counter++;
            }
            pixel_counter += current_image->padding_jump_value;
        }

        // repeat column sum into my_vertical_sum
        pixel_counter = 0;
        for ( int j = 0; j < directional_image_sum.logical_y_dimension; j++ ) {
            for ( int i = 0; i < directional_image_sum.logical_x_dimension; i++ ) {
                pixel_coord_y                                     = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(0, j, 0);
                pixel_coord_xy                                    = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                directional_image_sum.real_values[pixel_coord_xy] = directional_image_sum.real_values[pixel_coord_y];
                pixel_counter++;
            }
            pixel_counter += directional_image_sum.padding_jump_value;
        }

        directional_image_sum.DivideByConstant(directional_image_sum.logical_x_dimension);
    }
    // y-direction
    else {

        long pixel_coord_x  = 0;
        long pixel_coord_xy = 0;
        long pixel_counter  = 0;

        // sum columns of my_image_sum (NxM) and store in array (1xM)
        for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
            for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
                pixel_coord_x  = current_image->ReturnReal1DAddressFromPhysicalCoord(i, 0, 0);
                pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                directional_image_sum.real_values[pixel_coord_x] += current_image->real_values[pixel_coord_xy];
                pixel_counter++;
            }
            pixel_counter += current_image->padding_jump_value;
        }

        // repeat column sum into my_vertical_sum
        pixel_counter = 0;
        for ( int i = 0; i < directional_image_sum.logical_x_dimension; i++ ) {
            for ( int j = 0; j < directional_image_sum.logical_y_dimension; j++ ) {
                pixel_coord_x                                     = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(i, 0, 0);
                pixel_coord_xy                                    = directional_image_sum.ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
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

// Detects the two strongest outer-edge peaks in a 1D intensity profile.
// Returns indices of the best peak pair (sorted low->high), or an empty vector if none found.
std::pair<int, int> findOuterTubeEdges(const std::vector<float>& cols, float min_tube_diameter, float max_tube_diameter, bool find_positive_peaks, bool find_negative_peaks) {
    int n = cols.size( );
    if ( n < 3 )
        return {-1, -1}; // need at least 3 points to form a peak

    // 1) Find all local maxima (positive peaks).
    std::vector<std::pair<int, float>> posPeaks;
    posPeaks.reserve(n / 10);
    for ( int i = 1; i < n - 1; ++i ) {
        if ( cols[i] > cols[i - 1] && cols[i] > cols[i + 1] ) {
            posPeaks.emplace_back(i, cols[i]);
        }
    }

    // 2) Find all local minima (negative peaks), storing their absolute amplitudes.
    std::vector<std::pair<int, float>> negPeaks;
    negPeaks.reserve(n / 10);
    for ( int i = 1; i < n - 1; ++i ) {
        if ( cols[i] < cols[i - 1] && cols[i] < cols[i + 1] ) {
            // Treat negative peak by negating value to get positive amplitude
            negPeaks.emplace_back(i, -cols[i]);
        }
    }

    // // Helper lambda to find best pair (highest score) within a list of peaks
    // auto bestPair = [&](const std::vector<std::pair<int, float>>& peaks) -> std::pair<float, std::pair<int, int>> {
    //     float               bestScoreInRange = -std::numeric_limits<float>::infinity( );
    //     std::pair<int, int> bestIdxInRange   = {-1, -1};

    //     float               bestScoreOutOfRange = -std::numeric_limits<float>::infinity( );
    //     std::pair<int, int> bestIdxOutOfRange   = {-1, -1};
    //     float               bestGapError        = std::numeric_limits<float>::infinity( );

    //     const float IDEAL_GAP   = min_tube_diameter;
    //     const float GAP_PENALTY = 0.1;

    //     for ( size_t a = 0; a < peaks.size( ); ++a ) {
    //         for ( size_t b = a + 1; b < peaks.size( ); ++b ) {
    //             int   i      = peaks[a].first;
    //             int   j      = peaks[b].first;
    //             int   gap    = j - i;
    //             float sumAmp = peaks[a].second + peaks[b].second;
    //             float score  = sumAmp - GAP_PENALTY * std::fabs(gap - IDEAL_GAP);

    //             if ( gap >= min_tube_diameter && gap <= max_tube_diameter ) {
    //                 // candidate within range
    //                 if ( score > bestScoreInRange ) {
    //                     bestScoreInRange = score;
    //                     bestIdxInRange   = {i, j};
    //                 }
    //             }
    //             else {
    //                 // candidate out of range, but keep closest
    //                 float gapError = 0.0f;
    //                 if ( gap < min_tube_diameter )
    //                     gapError = min_tube_diameter - gap;
    //                 else
    //                     gapError = gap - max_tube_diameter;

    //                 if ( gapError < bestGapError ||
    //                      (gapError == bestGapError && score > bestScoreOutOfRange) ) {
    //                     bestGapError        = gapError;
    //                     bestScoreOutOfRange = score;
    //                     bestIdxOutOfRange   = {i, j};
    //                 }
    //             }
    //         }
    //     }
    //     if ( bestIdxInRange.first != -1 )
    //         return std::make_pair(bestScoreInRange, bestIdxInRange);
    //     else
    //         return std::make_pair(bestScoreOutOfRange, bestIdxOutOfRange);
    // };

    // // OLD Helper lambda to find best pair (highest score) within a list of peaks that was working as long as we are within the range
    // auto bestPair = [&](const std::vector<std::pair<int, float>>& peaks) {
    //     float               bestScore   = -std::numeric_limits<float>::infinity( );
    //     std::pair<int, int> bestIdx     = {-1, -1}; // We will prefer gap close to minimum tube diameter by subtracting a small penalty for gap deviation.
    //     const float         IDEAL_GAP   = min_tube_diameter;
    //     const float         GAP_PENALTY = 0.1; // e.g. 0.1 points lost per pixel of gap deviation
    //     // Sort peaks by index to ensure left<right
    //     // (Assumes peaks are already in ascending index order by scan loop.)
    //     for ( size_t a = 0; a < peaks.size( ); ++a ) {
    //         for ( size_t b = a + 1; b < peaks.size( ); ++b ) {
    //             int i   = peaks[a].first;
    //             int j   = peaks[b].first;
    //             int gap = j - i;
    //             if ( gap < min_tube_diameter )
    //                 continue; // enforce minimum gap
    //             if ( gap > max_tube_diameter )
    //                 break; // skip excessively large gaps (optional)
    //             // Sum amplitudes
    //             float sumAmp = peaks[a].second + peaks[b].second;
    //             // Apply a mild penalty for deviating from ideal gap=minimum tube diameter
    //             float score = sumAmp - GAP_PENALTY * std::fabs(gap - IDEAL_GAP);
    //             if ( score > bestScore ) {
    //                 bestScore = score;
    //                 bestIdx   = {i, j};
    //             }
    //         }
    //     }
    //     return std::make_pair(bestScore, bestIdx);
    // };

    /* 
    Removed separate in-range/out-of-range tracking.

    Instead, always compute a single score.

    Apply a big extra penalty if the gap is outside [min_tube_diameter, max_tube_diameter].

    Still returns the best scoring pair overall (even if it’s out of range)
    */
    auto bestPair = [&](const std::vector<std::pair<int, float>>& peaks)
            -> std::pair<float, std::pair<int, int>> {
        float               bestScore = -std::numeric_limits<float>::infinity( );
        std::pair<int, int> bestIdx   = {-1, -1};

        const float IDEAL_GAP           = min_tube_diameter;
        const float GAP_PENALTY         = 0.1f; // e.g. 0.1 points lost per pixel of gap deviation
        const float OUT_OF_RANGE_FACTOR = 2.0f; // scale factor for out-of-range penalty

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

    // Case 1: both positive and negative allowed
    if ( find_positive_peaks && find_negative_peaks ) {
        // 4) Choose the pair with larger combined amplitude (score).
        // No valid pair found in either
        if ( bestPos.first == -1 && bestNeg.first == -1 ) {
            return {-1, -1};
        }
        if ( scorePos > scoreNeg ) {
            bestPairIdx = bestPos;
        }
        else {
            bestPairIdx = bestNeg;
        }

        // Return sorted indices
        if ( bestPairIdx.first > bestPairIdx.second ) {
            std::swap(bestPairIdx.first, bestPairIdx.second);
        }
        return std::make_pair(bestPairIdx.first, bestPairIdx.second);
    }
    else if ( find_positive_peaks == true && find_negative_peaks == false ) {
        // Only use the best positive peak
        bestPairIdx = bestPos;
        // Return sorted indices
        if ( bestPairIdx.first > bestPairIdx.second ) {
            std::swap(bestPairIdx.first, bestPairIdx.second);
        }
        return std::make_pair(bestPairIdx.first, bestPairIdx.second);
    }
    else if ( find_positive_peaks == false && find_negative_peaks == true ) {
        // Only use the best negative peak
        bestPairIdx = bestNeg;
        // Return sorted indices
        if ( bestPairIdx.first > bestPairIdx.second ) {
            std::swap(bestPairIdx.first, bestPairIdx.second);
        }
        return std::make_pair(bestPairIdx.first, bestPairIdx.second);
    }
    else {
        std::make_pair(bestPairIdx.first, bestPairIdx.second); // it is -1 and -1 if anyother case happened
    }
}