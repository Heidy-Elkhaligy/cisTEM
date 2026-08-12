// #include <wx/defs.h>
// #include <wx/utils.h>
#include "../../core/gui_core_headers.h"

// #include <wx/filedlg.h>
// #include <wx/msgdlg.h>
// #include <wx/menu.h>
// #include <wx/sizer.h>
// #include <wx/statusbr.h>
// #include <wx/dcbuffer.h>
// #include <wx/spinctrl.h>
// #include <wx/textctrl.h>
// #include <wx/stattext.h>
// #include <wx/checkbox.h>
// #include <wx/progdlg.h>
// #include <vector>
// #include <cmath>
// #include <cfloat>
// #include <algorithm>
// #include <ctime>
// #include <cstdlib>
// #include <limits>
// #include <map>
// #include "../../gui/PlotHistPanel.h"

// // // ----------------------------------------------------------------------------
// // // Local maxima finder (SciPy-style, plateau aware)
// // // ----------------------------------------------------------------------------
// static void local_maxima_1d(
//         const std::vector<float>& x,
//         std::vector<int>&         midpoints,
//         std::vector<int>&         left_edges,
//         std::vector<int>&         right_edges) {
//     midpoints.clear( );
//     left_edges.clear( );
//     right_edges.clear( );

//     const int n = static_cast<int>(x.size( ));
//     if ( n < 3 )
//         return;

//     midpoints.reserve(n / 2);
//     left_edges.reserve(n / 2);
//     right_edges.reserve(n / 2);

//     int       i     = 1;
//     const int i_max = n - 1;

//     while ( i < i_max ) {
//         if ( x[i - 1] < x[i] ) {
//             int i_ahead = i + 1;

//             // Handle flat plateaus
//             while ( i_ahead < i_max && x[i_ahead] == x[i] ) {
//                 ++i_ahead;
//             }

//             // Confirm peak
//             if ( x[i_ahead] < x[i] ) {
//                 int left  = i;
//                 int right = i_ahead - 1;
//                 int mid   = (left + right) / 2;

//                 left_edges.push_back(left);
//                 right_edges.push_back(right);
//                 midpoints.push_back(mid);

//                 i = i_ahead;
//                 continue;
//             }
//         }
//         ++i;
//     }
// }

// // better but not the best
// static void local_maxima_1d(const std::vector<float>& x, std::vector<int>& midpoints, std::vector<int>& left_edges, std::vector<int>& right_edges) {
//     midpoints.clear( );
//     left_edges.clear( );
//     right_edges.clear( );
//     const int n = (int)x.size( );
//     if ( n < 3 )
//         return;

//     const float EPS = 1e-6f;
//     int         i = 1, i_max = n - 1;
//     while ( i < i_max ) {
//         if ( x[i - 1] + EPS < x[i] ) {
//             int i_ahead = i + 1;
//             while ( i_ahead < i_max && std::fabs(x[i_ahead] - x[i]) < EPS )
//                 ++i_ahead;
//             if ( x[i_ahead] + EPS < x[i] ) {
//                 int left  = i;
//                 int right = i_ahead - 1;
//                 int mid   = (left + right) / 2;
//                 // amplitude threshold relative to neighbors (avoid tiny noise peaks)
//                 float neigh_max = std::max(x[left - 1], x[right + 1]);
//                 if ( x[mid] - neigh_max > 1e-3f ) { // tune this threshold
//                     left_edges.push_back(left);
//                     right_edges.push_back(right);
//                     midpoints.push_back(mid);
//                 }
//                 i = i_ahead;
//                 continue;
//             }
//         }
//         ++i;
//     }
// }

float GetPixelSizeOrDefault(MRCFile& file) {
    float pixel_size = file.ReturnPixelSize( );

    if ( pixel_size <= 0.0f || ! std::isfinite(pixel_size) )
        return 1.0f;

    return pixel_size;
}

// Robust plateau-aware local maxima finder with depth filtering (no prominence)
static void local_maxima_1d(
        const std::vector<float>& x,
        std::vector<int>&         midpoints,
        std::vector<int>&         left_edges,
        std::vector<int>&         right_edges,
        float                     min_depth_abs = 0.0f, // absolute depth threshold (disabled if <= 0)
        float                     min_depth_rel = 0.0f, // relative depth (0..1) of local range, used if >0
        int                       min_distance  = 10) // minimal horizontal separation
{
    midpoints.clear( );
    left_edges.clear( );
    right_edges.clear( );

    const int n = (int)x.size( );
    if ( n < 3 )
        return;

    // global_range used only to scale tiny eps; not for depth decision
    auto [min_it, max_it] = std::minmax_element(x.begin( ), x.end( ));
    float global_range    = *max_it - *min_it;
    if ( global_range <= 0.0f )
        return;

    // window for local statistics (use something related to min_distance)
    int w = std::max(5, min_distance / 2);

    // small epsilon scaled to signal magnitude to detect plateaus robustly
    const float EPS = 1e-6f * std::max(1.0f, std::abs(*max_it));

    // Store depths for non-maximum suppression comparison
    std::vector<float> peak_depths;
    peak_depths.reserve(n / 8);

    int i = 1;
    while ( i < n - 1 ) {
        // detect rising edge into plateau/peak
        if ( x[i] > x[i - 1] + EPS ) {
            int j = i + 1;
            // handle plateau (equal values within EPS)
            while ( j < n - 1 && std::fabs(x[j] - x[i]) < EPS )
                ++j;

            // confirm actual peak (next distinct sample is smaller)
            if ( x[j] < x[i] - EPS ) {
                int left  = i;
                int right = j - 1;
                int mid   = (left + right) / 2;

                // compute local base = minimum in window around the peak midpoint
                int win_lo = std::max(0, mid - w);
                int win_hi = std::min(n - 1, mid + w);

                float base      = x[mid];
                float local_max = x[mid];
                for ( int k = win_lo; k <= win_hi; ++k ) {
                    if ( x[k] < base )
                        base = x[k];
                    if ( x[k] > local_max )
                        local_max = x[k];
                }

                // depth = how far above the local minimum this peak stands
                float depth = x[mid] - base;

                // local_range (for relative thresholding)
                float local_range = local_max - base;
                if ( local_range <= 0.0f )
                    local_range = 1.0f; // avoid div-by-zero

                // Decide acceptance:
                // If an absolute min_depth is provided (>0) use it.
                // Else if a relative min_depth_rel (>0) is provided, require depth >= min_depth_rel * local_range.
                // Else accept any detected local maximum (no depth filtering).
                bool accept = false;
                if ( min_depth_abs > 0.0f ) {
                    accept = (depth >= min_depth_abs);
                }
                else if ( min_depth_rel > 0.0f ) {
                    accept = (depth >= min_depth_rel * local_range);
                }
                else {
                    accept = true; // no depth requirement
                }

                if ( accept ) {
                    // Enforce minimum distance: compare using 'depth' metric
                    if ( ! midpoints.empty( ) && mid - midpoints.back( ) < min_distance ) {
                        // replace the previous peak if this one is stronger (deeper)
                        if ( depth > peak_depths.back( ) ) {
                            midpoints.back( )   = mid;
                            left_edges.back( )  = left;
                            right_edges.back( ) = right;
                            peak_depths.back( ) = depth;
                        }
                    }
                    else {
                        midpoints.push_back(mid);
                        left_edges.push_back(left);
                        right_edges.push_back(right);
                        peak_depths.push_back(depth);
                    }
                }

                i = j;
                continue;
            }
        }
        ++i;
    }
}

// ----------------------------------------------------------------------------
// Application Class
// ----------------------------------------------------------------------------
class FindTubeDiametersGuiApp : public wxApp {
  public:
    virtual bool OnInit( );
};

IMPLEMENT_APP(FindTubeDiametersGuiApp)

// ----------------------------------------------------------------------------
// Image Panel Class
// ----------------------------------------------------------------------------
class ImagePanel : public wxPanel {
    wxImage m_image;
    bool    m_has_image;

    std::vector<float>      m_profile;
    std::pair<float, float> m_edges;
    bool                    m_has_graph;
    bool                    m_show_graph;

  public:
    ImagePanel(wxWindow* parent) : wxPanel(parent), m_has_image(false), m_has_graph(false), m_show_graph(true), m_edges({-1.0f, -1.0f}) {
        SetBackgroundStyle(wxBG_STYLE_PAINT);
    }

    void SetImage(Image* image) {
        if ( ! image || image->logical_x_dimension == 0 )
            return;

        int w = image->logical_x_dimension;
        int h = image->logical_y_dimension;

        float min_val = FLT_MAX;
        float max_val = -FLT_MAX;

        for ( int y = 0; y < h; y++ ) {
            for ( int x = 0; x < w; x++ ) {
                float val = image->real_values[image->ReturnReal1DAddressFromPhysicalCoord(x, y, 0)];
                if ( val < min_val )
                    min_val = val;
                if ( val > max_val )
                    max_val = val;
            }
        }

        float range = max_val - min_val;
        if ( range <= 0 )
            range = 1.0f;

        unsigned char* rgb_data = (unsigned char*)malloc(w * h * 3);
        unsigned char* ptr      = rgb_data;

        for ( int y = 0; y < h; y++ ) {
            for ( int x = 0; x < w; x++ ) {
                float         val       = image->real_values[image->ReturnReal1DAddressFromPhysicalCoord(x, y, 0)];
                unsigned char pixel_val = (unsigned char)((val - min_val) / range * 255.0f);
                *ptr++                  = pixel_val;
                *ptr++                  = pixel_val;
                *ptr++                  = pixel_val;
            }
        }

        m_image     = wxImage(w, h, rgb_data);
        m_has_image = true;
        Refresh( );
    }

    void SetGraphData(const std::vector<float>& profile, std::pair<float, float> edges) {
        m_profile   = profile;
        m_edges     = edges;
        m_has_graph = true;
        Refresh( );
    }

    void ShowGraph(bool show) {
        m_show_graph = show;
        Refresh( );
    }

    void OnPaint(wxPaintEvent& evt) {
        wxAutoBufferedPaintDC dc(this);
        dc.Clear( );

        if ( m_has_image && m_image.IsOk( ) ) {
            wxSize sz = GetClientSize( );
            if ( sz.GetWidth( ) <= 0 || sz.GetHeight( ) <= 0 )
                return;

            int axis_height = 30;
            int available_h = sz.GetHeight( ) - axis_height;
            if ( available_h < 0 )
                available_h = sz.GetHeight( );

            float img_aspect   = (float)m_image.GetWidth( ) / m_image.GetHeight( );
            float panel_aspect = (float)sz.GetWidth( ) / available_h;

            int draw_w, draw_h;

            if ( panel_aspect > img_aspect ) {
                draw_h = available_h;
                draw_w = (int)(draw_h * img_aspect);
            }
            else {
                draw_w = sz.GetWidth( );
                draw_h = (int)(draw_w / img_aspect);
            }

            int x_off = (sz.GetWidth( ) - draw_w) / 2;
            int y_off = (available_h - draw_h) / 2;

            wxBitmap bmp(m_image.Scale(draw_w, draw_h, wxIMAGE_QUALITY_NORMAL));
            dc.DrawBitmap(bmp, x_off, y_off);

            dc.SetPen(*wxBLACK_PEN);
            dc.SetFont(wxFont(8, wxFONTFAMILY_SWISS, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL));

            int axis_y = y_off + draw_h + 2;
            dc.DrawLine(x_off, axis_y, x_off + draw_w, axis_y);

            int img_w = m_image.GetWidth( );
            int tick_interval;

            if ( img_w <= 100 ) {
                // "if image is 100 or less it would be spaced every 20"
                tick_interval = 20;
            }
            else if ( img_w < 500 ) {
                // "if > 100 and less than 500 it would be spaced every 40"
                tick_interval = 40;
            }
            else if ( img_w < 1000 ) {
                // "if more then ... 50" (for medium-large images)
                tick_interval = 50;
            }
            else {
                // "if more then ... 100" (for very large images)
                tick_interval = 100;
            }

            float scale_x = (float)draw_w / img_w;

            for ( int i = 0; i < img_w; i += tick_interval ) {
                int screen_x = x_off + (int)(i * scale_x);
                dc.DrawLine(screen_x, axis_y, screen_x, axis_y + 5);

                wxSize text_sz = dc.GetTextExtent(wxString::Format("%d", i));
                dc.DrawText(wxString::Format("%d", i), screen_x - (text_sz.GetWidth( ) / 2), axis_y + 6);
            }

            if ( m_show_graph && m_has_graph && ! m_profile.empty( ) ) {
                float p_min   = *std::min_element(m_profile.begin( ), m_profile.end( ));
                float p_max   = *std::max_element(m_profile.begin( ), m_profile.end( ));
                float p_range = p_max - p_min;
                if ( p_range <= 0 )
                    p_range = 1.0f;

                float plot_h     = draw_h * 0.5f;
                float plot_off_y = draw_h * 0.5f;

                dc.SetPen(wxPen(*wxCYAN, 2));

                for ( size_t i = 0; i < m_profile.size( ) - 1; ++i ) {
                    float y1 = (1.0f - (m_profile[i] - p_min) / p_range) * plot_h + plot_off_y;
                    float y2 = (1.0f - (m_profile[i + 1] - p_min) / p_range) * plot_h + plot_off_y;

                    int x1_draw = x_off + (int)(i * scale_x);
                    int x2_draw = x_off + (int)((i + 1) * scale_x);
                    int y1_draw = y_off + (int)y1;
                    int y2_draw = y_off + (int)y2;

                    dc.DrawLine(x1_draw, y1_draw, x2_draw, y2_draw);
                }

                if ( m_edges.first != -1 && m_edges.second != -1 ) {
                    dc.SetPen(wxPen(*wxRED, 2, wxPENSTYLE_DOT));
                    int x_edge1 = x_off + (int)(m_edges.first * scale_x);
                    int x_edge2 = x_off + (int)(m_edges.second * scale_x);

                    dc.DrawLine(x_edge1, y_off, x_edge1, y_off + draw_h);
                    dc.DrawLine(x_edge2, y_off, x_edge2, y_off + draw_h);
                }
            }
        }
    }

    void OnSize(wxSizeEvent& evt) {
        Refresh( );
        evt.Skip( );
    }

    DECLARE_EVENT_TABLE( )
};

BEGIN_EVENT_TABLE(ImagePanel, wxPanel)
EVT_PAINT(ImagePanel::OnPaint)
EVT_SIZE(ImagePanel::OnSize)
END_EVENT_TABLE( )

// ----------------------------------------------------------------------------
// Main Frame Class
// ----------------------------------------------------------------------------
class FindTubeMainFrame : public wxFrame {
    ImagePanel*   m_image_panel;
    ImagePanel*   m_image_panel_autocorr;
    ImagePanel*   m_image_panel_FT;
    wxString      m_filename;
    wxTextCtrl*   m_pixel_size_ctrl;
    wxTextCtrl*   m_min_diam_ctrl;
    wxTextCtrl*   m_max_diam_ctrl;
    wxTextCtrl*   m_mask_rad_ctrl;
    wxTextCtrl*   m_lp_res_ctrl;
    wxCheckBox*   m_align_autocorr_check;
    wxCheckBox*   m_align_FT_check;
    wxStaticText* m_result_text;
    //wxCheckBox*   m_chk_halfway;
    wxCheckBox*   m_chk_show_profile;
    wxCheckBox*   m_invert_contrast;
    wxTextCtrl*   m_txtSlice;
    wxButton*     m_btnGo;
    wxStaticText* m_lblTotal;

    wxButton* m_btn_next;
    wxButton* m_btn_prev;
    wxButton* m_btn_random;

    int m_current_slice;
    int m_total_slices;

    enum {
        ID_Open = 1,
        ID_PixelSize,
        ID_MinDiam,
        ID_MaxDiam,
        ID_Next,
        ID_Prev,
        ID_Random,
        ID_Calc,
        ID_Hist,
        ID_Halfway,
        ID_ShowProfile,
        ID_lp_res,
        ID_mask,
        ID_InvertContrast,
        ID_SaveParams
    };

  public:
    FindTubeMainFrame( ) : wxFrame(NULL, wxID_ANY, "Tube Diameter Finder", wxDefaultPosition, wxSize(1600, 800)) {
        wxMenu* menuFile = new wxMenu;
        menuFile->Append(ID_Open, "&Open Image...\tCtrl-O");
        menuFile->AppendSeparator( );
        menuFile->Append(wxID_EXIT);

        wxMenu* menuActions = new wxMenu;
        menuActions->Append(ID_Hist, "&Plot");
        menuActions->Enable(ID_Hist, false);

        menuActions->Append(ID_SaveParams, "&Save Parameters...\tCtrl-S"); // new menu item
        Bind(wxEVT_MENU, &FindTubeMainFrame::OnSaveParameters, this, ID_SaveParams);

        wxMenuBar* menuBar = new wxMenuBar;
        menuBar->Append(menuFile, "&File");
        menuBar->Append(menuActions, "&Actions");

        SetMenuBar(menuBar);

        CreateStatusBar( );
        SetStatusText("Open an MRC stack to begin.");

        wxBoxSizer* mainSizer     = new wxBoxSizer(wxVERTICAL);
        wxBoxSizer* imageRowSizer = new wxBoxSizer(wxHORIZONTAL);

        // Left: Autocorr image
        m_image_panel_autocorr = new ImagePanel(this);
        imageRowSizer->Add(m_image_panel_autocorr, 1, wxEXPAND | wxALL, 5);

        m_image_panel = new ImagePanel(this);
        imageRowSizer->Add(m_image_panel, 1, wxEXPAND | wxALL, 5);

        // Right: Fourier Transform image
        m_image_panel_FT = new ImagePanel(this);
        imageRowSizer->Add(m_image_panel_FT, 1, wxEXPAND | wxALL, 5);

        mainSizer->Add(imageRowSizer, 1, wxEXPAND | wxALL, 5);

        // Image Panel labels
        // Create a larger font
        wxFont labelFont = GetFont( );
        labelFont.SetPointSize(14);
        labelFont.SetWeight(wxFONTWEIGHT_BOLD);

        wxBoxSizer* labelSizer = new wxBoxSizer(wxHORIZONTAL);

        // ---- Column 1 ----
        wxStaticText* label1 = new wxStaticText(this, wxID_ANY, "Auto correlation");
        label1->SetFont(labelFont);

        wxBoxSizer* col1 = new wxBoxSizer(wxHORIZONTAL);
        col1->AddStretchSpacer(1);
        col1->Add(label1, 0, wxALIGN_CENTER);
        col1->AddStretchSpacer(1);

        labelSizer->Add(col1, 1, wxEXPAND);

        // ---- Column 2 ----
        wxStaticText* label2 = new wxStaticText(this, wxID_ANY, "Original");
        label2->SetFont(labelFont);

        wxBoxSizer* col2 = new wxBoxSizer(wxHORIZONTAL);
        col2->AddStretchSpacer(1);
        col2->Add(label2, 0, wxALIGN_CENTER);
        col2->AddStretchSpacer(1);

        labelSizer->Add(col2, 1, wxEXPAND);

        // ---- Column 3 ----
        wxStaticText* label3 = new wxStaticText(this, wxID_ANY, "Fourier Transform");
        label3->SetFont(labelFont);

        wxBoxSizer* col3 = new wxBoxSizer(wxHORIZONTAL);
        col3->AddStretchSpacer(1);
        col3->Add(label3, 0, wxALIGN_CENTER);
        col3->AddStretchSpacer(1);

        labelSizer->Add(col3, 1, wxEXPAND);

        // Add below images
        mainSizer->Add(labelSizer, 0, wxEXPAND | wxTOP, 5);

        //mainSizer->Add(labelSizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxBOTTOM, 5);

        this->Layout( );

        SetSizer(mainSizer);

        // --- Result Text ---
        m_result_text = new wxStaticText(this, wxID_ANY, "Load a file to see results.", wxDefaultPosition, wxDefaultSize, wxALIGN_CENTER);
        wxFont font   = m_result_text->GetFont( );
        font.SetWeight(wxFONTWEIGHT_BOLD);
        font.SetPointSize(12);
        m_result_text->SetFont(font);
        mainSizer->Add(m_result_text, 0, wxALIGN_CENTER_HORIZONTAL | wxALL, 10);
        m_result_text->SetWindowStyleFlag(wxALIGN_CENTER);

        SetSizer(mainSizer);

        // --- Combined Navigation Row ---
        wxBoxSizer* navRowSizer = new wxBoxSizer(wxHORIZONTAL);

        // Image number label
        wxStaticText* lblSlice = new wxStaticText(this, wxID_ANY, "Image #: ");
        navRowSizer->Add(lblSlice, 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 5);

        // Image number input
        m_txtSlice = new wxTextCtrl(this, wxID_ANY, "1", wxDefaultPosition, wxSize(50, -1));
        navRowSizer->Add(m_txtSlice, 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 5);

        // "Go" button
        m_btnGo = new wxButton(this, wxID_ANY, "Go"); // you might want a member for the Go button too
        navRowSizer->Add(m_btnGo, 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 10);
        m_btnGo->Bind(wxEVT_BUTTON, &FindTubeMainFrame::OnGoToSlice, this);
        m_txtSlice->Bind(wxEVT_TEXT_ENTER, &FindTubeMainFrame::OnGoToSlice, this);
        m_txtSlice->SetWindowStyleFlag(wxTE_PROCESS_ENTER);

        // Total images label
        m_lblTotal = new wxStaticText(this, wxID_ANY, "/ 0");
        navRowSizer->Add(m_lblTotal, 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 30);

        // Previous, Next, Random buttons
        m_btn_prev   = new wxButton(this, ID_Prev, "< Prev");
        m_btn_next   = new wxButton(this, ID_Next, "Next >");
        m_btn_random = new wxButton(this, ID_Random, "Random");

        m_btn_prev->Enable(false);
        m_btn_next->Enable(false);
        m_btn_random->Enable(false);

        navRowSizer->Add(m_btn_prev, 0, wxALL, 10); //5
        navRowSizer->Add(m_btn_next, 0, wxALL, 10);
        navRowSizer->Add(m_btn_random, 0, wxALL, 10);

        // Add the combined sizer to mainSizer
        mainSizer->Add(navRowSizer, 0, wxALIGN_CENTER_HORIZONTAL | wxBOTTOM, 10);

        // --- Controls Sizer ---
        wxBoxSizer* controlSizer = new wxBoxSizer(wxHORIZONTAL);

        controlSizer->Add(new wxStaticText(this, wxID_ANY, "Pixel Size (A):"), 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 5);
        m_pixel_size_ctrl = new wxTextCtrl(this, ID_PixelSize, "1.0", wxDefaultPosition, wxSize(60, -1), wxTE_PROCESS_ENTER);
        controlSizer->Add(m_pixel_size_ctrl, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 5);

        controlSizer->Add(new wxStaticText(this, wxID_ANY, "Min Diam (px):"), 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 10);
        m_min_diam_ctrl = new wxTextCtrl(this, ID_MinDiam, "100", wxDefaultPosition, wxSize(60, -1), wxTE_PROCESS_ENTER);
        controlSizer->Add(m_min_diam_ctrl, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 5);

        controlSizer->Add(new wxStaticText(this, wxID_ANY, "Max Diam (px):"), 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 10);
        m_max_diam_ctrl = new wxTextCtrl(this, ID_MaxDiam, "300", wxDefaultPosition, wxSize(60, -1), wxTE_PROCESS_ENTER);
        controlSizer->Add(m_max_diam_ctrl, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 5);

        controlSizer->Add(new wxStaticText(this, wxID_ANY, "Mask Radius (px):"), 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 10);
        m_mask_rad_ctrl = new wxTextCtrl(this, ID_mask, "0", wxDefaultPosition, wxSize(60, -1), wxTE_PROCESS_ENTER);
        controlSizer->Add(m_mask_rad_ctrl, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 5);

        controlSizer->Add(new wxStaticText(this, wxID_ANY, "Low Pass Res:"), 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 10);
        m_lp_res_ctrl = new wxTextCtrl(this, ID_lp_res, "150", wxDefaultPosition, wxSize(60, -1), wxTE_PROCESS_ENTER);
        controlSizer->Add(m_lp_res_ctrl, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 5);

        m_pixel_size_ctrl->Bind(wxEVT_TEXT_ENTER, &FindTubeMainFrame::OnParamChange, this);
        m_min_diam_ctrl->Bind(wxEVT_TEXT_ENTER, &FindTubeMainFrame::OnParamChange, this);
        m_max_diam_ctrl->Bind(wxEVT_TEXT_ENTER, &FindTubeMainFrame::OnParamChange, this);
        m_mask_rad_ctrl->Bind(wxEVT_TEXT_ENTER, &FindTubeMainFrame::OnParamChange, this);
        m_lp_res_ctrl->Bind(wxEVT_TEXT_ENTER, &FindTubeMainFrame::OnParamChange, this);

        // m_chk_halfway = new wxCheckBox(this, ID_Halfway, "Use Half-Max");
        // m_chk_halfway->SetValue(false);
        // controlSizer->Add(m_chk_halfway, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 15);
        wxButton* btnUpdate = new wxButton(this, ID_Calc, "Update");
        btnUpdate->Bind(wxEVT_BUTTON, &FindTubeMainFrame::OnParamChange, this);
        controlSizer->Add(btnUpdate, 0, wxALIGN_CENTER_VERTICAL | wxLEFT | wxRIGHT, 10);

        mainSizer->Add(controlSizer, 0, wxALIGN_CENTER_HORIZONTAL | wxALL, 10);

        // --- Parameter Update & Histogram ---
        wxBoxSizer* paramSizer = new wxBoxSizer(wxHORIZONTAL);

        m_align_autocorr_check = new wxCheckBox(this, ID_Halfway, "Align image (auto-corr)");
        m_align_autocorr_check->SetValue(false);
        paramSizer->Add(m_align_autocorr_check, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 15);

        m_align_FT_check = new wxCheckBox(this, ID_Halfway, "Align image (FT)");
        m_align_FT_check->SetValue(false);
        paramSizer->Add(m_align_FT_check, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 15);

        m_align_autocorr_check->Bind(wxEVT_CHECKBOX, &FindTubeMainFrame::OnAlignCheckbox, this);
        m_align_FT_check->Bind(wxEVT_CHECKBOX, &FindTubeMainFrame::OnAlignCheckbox, this);

        m_chk_show_profile = new wxCheckBox(this, ID_ShowProfile, "Show Profile");
        m_chk_show_profile->SetValue(true);
        m_chk_show_profile->Bind(wxEVT_CHECKBOX, &FindTubeMainFrame::OnToggleProfile, this);
        paramSizer->Add(m_chk_show_profile, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 10);

        m_invert_contrast = new wxCheckBox(this, ID_InvertContrast, "Particles are Black");
        m_invert_contrast->SetValue(false);
        paramSizer->Add(m_invert_contrast, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 10);
        // Bind event so that toggling the checkbox recalculates
        m_invert_contrast->Bind(wxEVT_CHECKBOX, &FindTubeMainFrame::OnInvertContrastToggle, this);

        mainSizer->Add(paramSizer, 0, wxALIGN_CENTER_HORIZONTAL | wxBOTTOM, 10);
    }

    void OnOpen(wxCommandEvent& event) {
        wxFileDialog openFileDialog(this, "Open MRC file", "", "", "MRC files (*.mrc;*.mrcs)|*.mrc;*.mrcs", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
        if ( openFileDialog.ShowModal( ) == wxID_CANCEL )
            return;

        m_filename = openFileDialog.GetPath( );
        MRCFile input_file(m_filename.ToStdString( ), false);

        if ( ! input_file.OpenFile(m_filename.ToStdString( ), false) ) {
            wxMessageBox("Error opening file", "Error", wxICON_ERROR);
            return;
        }
        m_pixel_size_ctrl->SetValue(wxString::Format("%.3f", GetPixelSizeOrDefault(input_file)));

        m_total_slices  = input_file.ReturnNumberOfSlices( );
        m_current_slice = 1;

        m_lblTotal->SetLabel(wxString::Format("/ %d", m_total_slices));
        m_txtSlice->SetValue("1");

        if ( m_total_slices > 1 ) {
            m_btn_next->Enable(true);
            m_btn_random->Enable(true);
        }

        GetMenuBar( )->Enable(ID_Hist, true);

        CalculateAndDisplay( );
    }

    void OnGoToSlice(wxCommandEvent& event) {
        if ( m_total_slices <= 0 )
            return;

        long slice;
        if ( ! m_txtSlice->GetValue( ).ToLong(&slice) ) {
            wxMessageBox("Invalid slice number.", "Error", wxICON_ERROR);
            return;
        }

        if ( slice < 1 || slice > m_total_slices ) {
            wxMessageBox(wxString::Format("Slice number must be between 1 and %d.", m_total_slices), "Error", wxICON_ERROR);
            return;
        }

        m_current_slice = static_cast<int>(slice);
        CalculateAndDisplay( );

        // Update navigation buttons if needed
        UpdateNavButtons( );
    }

    void OnNext(wxCommandEvent& event) {
        if ( m_current_slice < m_total_slices ) {
            m_current_slice++;
            m_txtSlice->SetValue(wxString::Format("%d", m_current_slice));
            CalculateAndDisplay( );
        }
        UpdateNavButtons( );
    }

    void OnPrev(wxCommandEvent& event) {
        if ( m_current_slice > 1 ) {
            m_current_slice--;
            m_txtSlice->SetValue(wxString::Format("%d", m_current_slice));
            CalculateAndDisplay( );
        }
        UpdateNavButtons( );
    }

    void OnRandom(wxCommandEvent& event) {
        if ( m_total_slices > 1 ) {
            m_current_slice = (rand( ) % m_total_slices) + 1;
            m_txtSlice->SetValue(wxString::Format("%d", m_current_slice));
            CalculateAndDisplay( );
        }
        UpdateNavButtons( );
    }

    void OnParamChange(wxCommandEvent& event) {
        CalculateAndDisplay( );
    }

    void OnToggleProfile(wxCommandEvent& event) {
        if ( m_image_panel ) {
            m_image_panel->ShowGraph(m_chk_show_profile->IsChecked( ));
        }
    }

    void UpdateNavButtons( ) {
        m_btn_prev->Enable(m_current_slice > 1);
        m_btn_next->Enable(m_current_slice < m_total_slices);
    }

    void OnInvertContrastToggle(wxCommandEvent& event) {
        // Recalculate and update the display
        CalculateAndDisplay( );
    }

    void OnSaveParameters(wxCommandEvent& WXUNUSED(event)) {
        wxFileDialog saveFileDialog(this, _("Save Parameters"), "", "",
                                    "Text files (*.txt)|*.txt|All files (*.*)|*.*",
                                    wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
        if ( saveFileDialog.ShowModal( ) == wxID_CANCEL )
            return;

        wxString filename = saveFileDialog.GetPath( );

        // Get current values
        wxString pixelSize = m_pixel_size_ctrl->GetValue( );
        wxString minDiam   = m_min_diam_ctrl->GetValue( );
        wxString maxDiam   = m_max_diam_ctrl->GetValue( );
        wxString maskRad   = m_mask_rad_ctrl->GetValue( );
        wxString lpRes     = m_lp_res_ctrl->GetValue( );
        // Determine which alignment option is selected
        wxString alignment;
        if ( m_align_autocorr_check->IsChecked( ) ) {
            alignment = "Autocorrelation";
        }
        else if ( m_align_FT_check->IsChecked( ) ) {
            alignment = "Fourier Transform";
        }
        else {
            alignment = "None";
        }

        // Save to file
        wxTextFile file;
        if ( file.Create(filename) || file.Open(filename) ) {
            file.Clear( );
            file.AddLine(wxString::Format("Finding Tube Diameter Program Parameters"));
            file.AddLine(wxString::Format("File Name=%s", m_filename.ToStdString( )));
            file.AddLine(wxString::Format("PixelSize=%s", pixelSize));
            file.AddLine(wxString::Format("MinDiam=%s", minDiam));
            file.AddLine(wxString::Format("MaxDiam=%s", maxDiam));
            file.AddLine(wxString::Format("MaskRadius=%s", maskRad));
            file.AddLine(wxString::Format("LowPassRes=%s", lpRes));
            file.AddLine(wxString::Format("Alignment=%s", alignment));
            file.Write( );
            file.Close( );

            wxMessageBox("Parameters saved successfully.", "Info", wxOK | wxICON_INFORMATION);
        }
        else {
            wxMessageBox("Failed to save parameters.", "Error", wxOK | wxICON_ERROR);
        }
    }

    // struct EdgeCandidate {
    //     float diff; // Contrast strength (amplitude difference)
    //     float midpoint; // Sub-pixel location
    // };

    // std::pair<int, int> FindOuterTubeEdges(const std::vector<float>& cols,
    //                                        float                     min_tube_diameter,
    //                                        float                     max_tube_diameter,
    //                                        bool                      invert_contrast) {
    //     int n = static_cast<int>(cols.size( ));
    //     if ( n < 3 )
    //         return {-1, -1};

    //     // 1. Preprocessing & Normalization
    //     std::vector<float> profile = cols;
    //     if ( invert_contrast ) {
    //         for ( float& val : profile )
    //             val = -val;
    //     }

    //     float              min_val = *std::min_element(profile.begin( ), profile.end( ));
    //     float              max_val = *std::max_element(profile.begin( ), profile.end( ));
    //     std::vector<float> norm(n);
    //     // // normalizing the values to better find prominent peaks
    //     // float range = max_val - min_val;
    //     // if ( range <= 1e-6f )
    //     //     range = 1.0f;

    //     for ( int i = 0; i < n; ++i )
    //         norm[i] = (profile[i] - min_val); /// range

    //     float              max_norm = *std::max_element(norm.begin( ), norm.end( ));
    //     std::vector<float> norm_inv(n);
    //     for ( int i = 0; i < n; ++i )
    //         norm_inv[i] = max_norm - norm[i];
    //     // for ( int i = 0; i < n; ++i )
    //     //     norm_inv[i] = 1.0f - norm[i];

    //     // 2. Detect Peaks
    //     std::vector<int> pos_mids, p_l, p_r;
    //     std::vector<int> neg_mids, n_l, n_r;
    //     local_maxima_1d(norm, pos_mids, p_l, p_r, 0.0f, 0.20f, 10); //Use relative depth (e.g. require peak >= 20% of local swing):
    //     local_maxima_1d(norm_inv, neg_mids, n_l, n_r, 0.0f, 0.20f, 10); // expected minimum separation between distinct local maxima 10

    //     std::cout << "\n=== PROFILE DEBUG ===\n";
    //     for ( int i = 0; i < n; ++i ) {
    //         std::cout << i
    //                   << "  norm=" << norm[i]
    //                   << "  norm_inv=" << norm_inv[i];

    //         if ( std::find(pos_mids.begin( ), pos_mids.end( ), i) != pos_mids.end( ) )
    //             std::cout << "  <-- POS_PEAK";

    //         if ( std::find(neg_mids.begin( ), neg_mids.end( ), i) != neg_mids.end( ) )
    //             std::cout << "  <-- NEG_PEAK";

    //         std::cout << "\n";
    //     }
    //     // Ensure peaks are sorted by their index (left-to-right)
    //     std::sort(pos_mids.begin( ), pos_mids.end( ));
    //     std::sort(neg_mids.begin( ), neg_mids.end( ));

    //     // 3. Classify Candidates into Left (Rising) and Right (Falling) Lists
    //     std::vector<EdgeCandidate> left_candidates;
    //     std::vector<EdgeCandidate> right_candidates;

    //     // for ( int pos : pos_mids ) {
    //     //     // --- Left Wall (Rising: Neg -> Pos) ---
    //     //     int   best_neg_left = -1;
    //     //     float min_dist_left = std::numeric_limits<float>::max( );

    //     //     for ( int neg : neg_mids ) {
    //     //         if ( neg < pos ) {
    //     //             float dist = static_cast<float>(pos - neg);
    //     //             if ( dist < min_dist_left ) {
    //     //                 min_dist_left = dist;
    //     //                 best_neg_left = neg;
    //     //             }
    //     //         }
    //     //     }
    //     //     if ( best_neg_left != -1 ) {
    //     //         float diff = std::abs(norm[pos] - norm[best_neg_left]);
    //     //         float mid  = (pos + best_neg_left) / 2.0f;
    //     //         left_candidates.push_back({diff, mid});
    //     //     }

    //     //     // --- Right Wall (Falling: Pos -> Neg) ---
    //     //     int   best_neg_right = -1;
    //     //     float min_dist_right = std::numeric_limits<float>::max( );

    //     //     for ( int neg : neg_mids ) {
    //     //         if ( neg > pos ) {
    //     //             float dist = static_cast<float>(neg - pos);
    //     //             if ( dist < min_dist_right ) {
    //     //                 min_dist_right = dist;
    //     //                 best_neg_right = neg;
    //     //             }
    //     //         }
    //     //     }
    //     //     if ( best_neg_right != -1 ) {
    //     //         float diff = std::abs(norm[pos] - norm[best_neg_right]);
    //     //         float mid  = (pos + best_neg_right) / 2.0f;
    //     //         right_candidates.push_back({diff, mid});
    //     //     }
    //     // }

    //     // precompute gradient (optionally smoothing beforehand)
    //     std::vector<float> grad(n, 0.0f);
    //     for ( int i = 1; i < n - 1; ++i )
    //         grad[i] = 0.5f * (norm[i + 1] - norm[i - 1]);
    //     grad[0]     = grad[1];
    //     grad[n - 1] = grad[n - 2];

    //     for ( int pos : pos_mids ) {
    //         // best NEG by contrast to the left
    //         int   best_neg_left      = -1;
    //         float best_contrast_left = -1.0f;
    //         for ( int neg : neg_mids ) {
    //             if ( neg < pos ) {
    //                 float contrast = std::abs(norm[pos] - norm[neg]);
    //                 if ( contrast > best_contrast_left ) {
    //                     best_contrast_left = contrast;
    //                     best_neg_left      = neg;
    //                 }
    //             }
    //         }
    //         if ( best_neg_left != -1 ) {
    //             // pick gradient maximum between neg and pos as midpoint
    //             int   lo      = std::min(best_neg_left, pos);
    //             int   hi      = std::max(best_neg_left, pos);
    //             int   best_k  = lo;
    //             float best_ag = std::fabs(grad[lo]);
    //             for ( int k = lo; k <= hi; ++k ) {
    //                 float ag = std::fabs(grad[k]);
    //                 if ( ag > best_ag ) {
    //                     best_ag = ag;
    //                     best_k  = k;
    //                 }
    //             }
    //             left_candidates.push_back({best_contrast_left, (float)best_k});
    //         }

    //         // best NEG by contrast to the right
    //         int   best_neg_right      = -1;
    //         float best_contrast_right = -1.0f;
    //         for ( int neg : neg_mids ) {
    //             if ( neg > pos ) {
    //                 float contrast = std::abs(norm[pos] - norm[neg]);
    //                 if ( contrast > best_contrast_right ) {
    //                     best_contrast_right = contrast;
    //                     best_neg_right      = neg;
    //                 }
    //             }
    //         }
    //         if ( best_neg_right != -1 ) {
    //             int   lo      = std::min(pos, best_neg_right);
    //             int   hi      = std::max(pos, best_neg_right);
    //             int   best_k  = lo;
    //             float best_ag = std::fabs(grad[lo]);
    //             for ( int k = lo; k <= hi; ++k ) {
    //                 float ag = std::fabs(grad[k]);
    //                 if ( ag > best_ag ) {
    //                     best_ag = ag;
    //                     best_k  = k;
    //                 }
    //             }
    //             right_candidates.push_back({best_contrast_right, (float)best_k});
    //         }
    //     }

    //     // 4. Sort candidates by strength to prioritize better edges during debug
    //     auto sort_fn = [](const EdgeCandidate& a, const EdgeCandidate& b) {
    //         return a.diff > b.diff;
    //     };
    //     std::sort(left_candidates.begin( ), left_candidates.end( ), sort_fn);
    //     std::sort(right_candidates.begin( ), right_candidates.end( ), sort_fn);

    //     // 5. Find Best Pair using Specific Scoring Function
    //     float best_score     = -std::numeric_limits<float>::infinity( );
    //     int   best_left_idx  = -1;
    //     int   best_right_idx = -1;

    //     const float IDEAL_GAP           = min_tube_diameter;
    //     const float GAP_PENALTY         = 0.1f;
    //     const float OUT_OF_RANGE_FACTOR = 10.0f;

    //     for ( const auto& l : left_candidates ) {
    //         for ( const auto& r : right_candidates ) {

    //             // Ensure gap is positive and absolute
    //             // float gap    = std::abs(r.midpoint - l.midpoint);
    //             float gap = r.midpoint - l.midpoint;
    //             if ( gap <= 0.0f )
    //                 continue; // right_candidates sorted by midpoint -> can break early
    //             if ( gap < min_tube_diameter * 0.5f || gap > max_tube_diameter * 1.5f )
    //                 continue;
    //             float sumAmp = l.diff + r.diff;

    //             // Base score: Contrast - Deviation from Ideal
    //             float score = sumAmp - GAP_PENALTY * std::fabs(gap - IDEAL_GAP);

    //             // Apply Out-of-Range Penalties
    //             if ( gap < min_tube_diameter ) {
    //                 score -= OUT_OF_RANGE_FACTOR * (min_tube_diameter - gap);
    //             }
    //             else if ( gap > max_tube_diameter ) {
    //                 score -= OUT_OF_RANGE_FACTOR * (gap - max_tube_diameter);
    //             }

    //             // Update Best
    //             if ( score > best_score ) {
    //                 best_score     = score;
    //                 best_left_idx  = static_cast<int>(std::round(l.midpoint));
    //                 best_right_idx = static_cast<int>(std::round(r.midpoint));
    //             }
    //         }
    //     }

    //     // 6. Return Result
    //     // If score is still negative infinity, no valid pairs existed
    //     if ( best_score == -std::numeric_limits<float>::infinity( ) ) {
    //         return {-1, -1};
    //     }

    //     // Optional sanity check: If the best score is still incredibly low due to
    //     // massive penalties, you might want to return {-1, -1} here too.
    //     // e.g. if (best_score < -100.0f) return {-1, -1};

    //     return {best_left_idx, best_right_idx};
    // }

    struct EdgeCandidate {
        float diff; // Contrast strength
        float midpoint; // Pixel (or subpixel) location
    };

    std::pair<int, int> FindOuterTubeEdges(const std::vector<float>& cols,
                                           float                     min_tube_diameter,
                                           float                     max_tube_diameter,
                                           bool                      invert_contrast) {
        int n = static_cast<int>(cols.size( ));
        if ( n < 3 )
            return {-1, -1};

        const int   center     = n / 2;
        const float MAX_OFFSET = 0.45f * n; // reject far-border edges
        //const float CENTER_W   = 0.6f; // center bias weight

        auto spatially_valid = [&](float mid) {
            return std::abs(mid - center) <= MAX_OFFSET;
        };

        // ------------------------------------------------------------
        // 1. Preprocessing
        // ------------------------------------------------------------
        std::vector<float> profile = cols;
        if ( invert_contrast ) {
            for ( float& v : profile )
                v = -v;
        }

        float min_val = *std::min_element(profile.begin( ), profile.end( ));
        float max_val = *std::max_element(profile.begin( ), profile.end( ));

        std::vector<float> norm(n);
        for ( int i = 0; i < n; ++i )
            norm[i] = profile[i] - min_val;

        float              max_norm = *std::max_element(norm.begin( ), norm.end( ));
        std::vector<float> norm_inv(n);
        for ( int i = 0; i < n; ++i )
            norm_inv[i] = max_norm - norm[i];

        // ------------------------------------------------------------
        // 2. Peak Detection
        // ------------------------------------------------------------
        std::vector<int> pos_mids, neg_mids;
        std::vector<int> tmp;

        local_maxima_1d(norm, pos_mids, tmp, tmp, 0.0f, 0.20f, 10);
        local_maxima_1d(norm_inv, neg_mids, tmp, tmp, 0.0f, 0.20f, 10);

        std::sort(pos_mids.begin( ), pos_mids.end( ));
        std::sort(neg_mids.begin( ), neg_mids.end( ));

        // ------------------------------------------------------------
        // DEBUG: profile + peaks
        // ------------------------------------------------------------
        std::cout << "\n=== PROFILE DEBUG ===\n";
        for ( int i = 0; i < n; ++i ) {
            std::cout << i
                      << " norm=" << norm[i]
                      << " inv=" << norm_inv[i];

            if ( std::find(pos_mids.begin( ), pos_mids.end( ), i) != pos_mids.end( ) )
                std::cout << " <-- POS";

            if ( std::find(neg_mids.begin( ), neg_mids.end( ), i) != neg_mids.end( ) )
                std::cout << " <-- NEG";

            std::cout << "\n";
        }

        // ------------------------------------------------------------
        // 3. Gradient (for edge localization)
        // ------------------------------------------------------------
        std::vector<float> grad(n, 0.0f);
        for ( int i = 1; i < n - 1; ++i )
            grad[i] = 0.5f * (norm[i + 1] - norm[i - 1]);
        grad[0]     = grad[1];
        grad[n - 1] = grad[n - 2];

        // ------------------------------------------------------------
        // 4. Build Edge Candidates
        // ------------------------------------------------------------
        std::vector<EdgeCandidate> left_candidates;
        std::vector<EdgeCandidate> right_candidates;

        for ( int pos : pos_mids ) {

            // ---- LEFT EDGE (NEG -> POS) ----
            int   best_neg_left = -1;
            float best_contrast = -1.0f;

            for ( int neg : neg_mids ) {
                if ( neg < pos ) {
                    float c = std::abs(norm[pos] - norm[neg]);
                    if ( c > best_contrast ) {
                        best_contrast = c;
                        best_neg_left = neg;
                    }
                }
            }

            if ( best_neg_left != -1 ) {
                int   lo     = best_neg_left;
                int   hi     = pos;
                int   best_k = lo;
                float best_g = std::fabs(grad[lo]);
                for ( int k = lo; k <= hi; ++k ) {
                    float g = std::fabs(grad[k]);
                    if ( g > best_g ) {
                        best_g = g;
                        best_k = k;
                    }
                }
                left_candidates.push_back({best_contrast, (float)best_k});
            }

            // ---- RIGHT EDGE (POS -> NEG) ----
            int best_neg_right = -1;
            best_contrast      = -1.0f;

            for ( int neg : neg_mids ) {
                if ( neg > pos ) {
                    float c = std::abs(norm[pos] - norm[neg]);
                    if ( c > best_contrast ) {
                        best_contrast  = c;
                        best_neg_right = neg;
                    }
                }
            }

            if ( best_neg_right != -1 ) {
                int   lo     = pos;
                int   hi     = best_neg_right;
                int   best_k = lo;
                float best_g = std::fabs(grad[lo]);
                for ( int k = lo; k <= hi; ++k ) {
                    float g = std::fabs(grad[k]);
                    if ( g > best_g ) {
                        best_g = g;
                        best_k = k;
                    }
                }
                right_candidates.push_back({best_contrast, (float)best_k});
            }
        }

        // // ------------------------------------------------------------
        // // 5A. OPTION C: Center-first selection
        // // ------------------------------------------------------------
        // bool  found_center      = false;
        // float best_center_score = -1e30f;
        // int   best_l = -1, best_r = -1;

        // std::cout << "\n=== CENTER-FIRST SCORING ===\n";

        // for ( const auto& l : left_candidates ) {
        //     if ( ! spatially_valid(l.midpoint) || l.midpoint >= center )
        //         continue;

        //     for ( const auto& r : right_candidates ) {
        //         if ( ! spatially_valid(r.midpoint) || r.midpoint <= center )
        //             continue;

        //         float gap = r.midpoint - l.midpoint;
        //         if ( gap < min_tube_diameter || gap > max_tube_diameter )
        //             continue;

        //         float center_dist =
        //                 std::abs(l.midpoint - center) +
        //                 std::abs(r.midpoint - center);

        //         float score =
        //                 (l.diff + r.diff) - CENTER_W * center_dist;

        //         std::cout << "L=" << l.midpoint
        //                   << " R=" << r.midpoint
        //                   << " gap=" << gap
        //                   << " contrast=" << (l.diff + r.diff)
        //                   << " center_dist=" << center_dist
        //                   << " score=" << score << "\n";

        //         if ( score > best_center_score ) {
        //             best_center_score = score;
        //             best_l            = (int)std::round(l.midpoint);
        //             best_r            = (int)std::round(r.midpoint);
        //             found_center      = true;
        //         }
        //     }
        // }
        //
        // if ( found_center ) {
        //     std::cout << "\n>>> CENTER-FIRST CHOSEN: "
        //               << best_l << ", " << best_r << "\n";
        //     return {best_l, best_r};
        // }

        // // ------------------------------------------------------------
        // // 5B. FALLBACK: Global scoring (your original logic)
        // // ------------------------------------------------------------
        // float best_score = -1e30f;
        // int   best_left = -1, best_right = -1;

        // std::cout << "\n=== FALLBACK GLOBAL SCORING ===\n";

        // for ( const auto& l : left_candidates ) {
        //     for ( const auto& r : right_candidates ) {
        //         float gap = r.midpoint - l.midpoint;
        //         if ( gap <= 0 )
        //             continue;

        //         float score = (l.diff + r.diff) - 0.1f * std::fabs(gap - min_tube_diameter);

        //         std::cout << "L=" << l.midpoint
        //                   << " R=" << r.midpoint
        //                   << " gap=" << gap
        //                   << " score=" << score << "\n";

        //         if ( score > best_score ) {
        //             best_score = score;
        //             best_left  = (int)std::round(l.midpoint);
        //             best_right = (int)std::round(r.midpoint);
        //         }
        //     }
        // }

        // if ( best_left < 0 || best_right < 0 )
        //     return {-1, -1};

        // std::cout << "\n>>> FALLBACK CHOSEN: "
        //           << best_left << ", " << best_right << "\n";

        // return {best_left, best_right};

        // 5. Find Best Pair using Improved Scoring Function (CENTER-SOFT + GAP-GAUSSIAN)

        float best_score     = -std::numeric_limits<float>::infinity( );
        int   best_left_idx  = -1;
        int   best_right_idx = -1;

        const float IDEAL_GAP = 0.5f * (min_tube_diameter + max_tube_diameter);

        // ---- TUNABLE PARAMETERS ----
        const float SIGMA_GAP = 0.25f * (max_tube_diameter - min_tube_diameter);
        // controls how strict the tube diameter must be

        const float CENTER_W = 0.4f; // soft prior only (do NOT exceed ~0.6)

        const float IMAGE_CENTER = 0.5f * (n - 1);

        std::cout << "\n=== FINAL SCORING ===\n";

        for ( const auto& l : left_candidates ) {
            for ( const auto& r : right_candidates ) {

                float gap = r.midpoint - l.midpoint;
                if ( gap <= 0.0f )
                    continue;

                // Hard rejection outside reasonable physical bounds
                if ( gap < min_tube_diameter || gap > max_tube_diameter )
                    continue;

                // ----------------------------------------------------
                // 1) Edge strength term
                // ----------------------------------------------------
                float contrast = l.diff + r.diff;

                // ----------------------------------------------------
                // 2) Gaussian gap likelihood (KEY FIX)
                // ----------------------------------------------------
                float gap_err   = gap - IDEAL_GAP;
                float gap_score = -(gap_err * gap_err) / (2.0f * SIGMA_GAP * SIGMA_GAP);

                // ----------------------------------------------------
                // 3) Soft center prior (does NOT dominate)
                // ----------------------------------------------------
                float mid_center     = 0.5f * (l.midpoint + r.midpoint);
                float center_dist    = std::abs(mid_center - IMAGE_CENTER);
                float center_penalty = CENTER_W * std::pow(center_dist, 1.5); //std:sqrt

                // ----------------------------------------------------
                // Final score
                // ----------------------------------------------------
                float score = contrast + (gap_score * 2) - center_penalty; //increasing the gap score penalty

                std::cout
                        << "L=" << (int)std::round(l.midpoint)
                        << " R=" << (int)std::round(r.midpoint)
                        << " gap=" << gap
                        << " contrast=" << contrast
                        << " gap_score=" << gap_score
                        << " center_dist=" << center_dist
                        << " score=" << score
                        << "\n";

                if ( score > best_score ) {
                    best_score     = score;
                    best_left_idx  = (int)std::round(l.midpoint);
                    best_right_idx = (int)std::round(r.midpoint);
                }
            }
        }

        std::cout << "\n>>> CHOSEN: L=" << best_left_idx
                  << " R=" << best_right_idx
                  << " score=" << best_score << "\n";

        return {best_left_idx, best_right_idx};
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Simple absolute column sum logic
    float GetMaxAbsColumnSum(Image* current_image) {
        std::vector<float> column_sum(current_image->logical_x_dimension, 0.0);

        //long pixel_counter = 0;

        for ( int i = 0; i < current_image->logical_x_dimension; i++ ) {
            for ( int j = 0; j < current_image->logical_y_dimension; j++ ) {
                long pixel_coord_xy = current_image->ReturnReal1DAddressFromPhysicalCoord(i, j, 0);
                column_sum[i] += current_image->real_values[pixel_coord_xy];
                //pixel_counter++;
            }
            //pixel_counter += current_image->padding_jump_value;
        }

        float max_value = *std::max_element(column_sum.begin( ), column_sum.end( ),
                                            [](float a, float b) { return std::abs(a) < std::abs(b); });

        return abs(max_value);
    }

    void OnAlignCheckbox(wxCommandEvent& event) {
        if ( event.GetEventObject( ) == m_align_autocorr_check ) {
            if ( m_align_autocorr_check->IsChecked( ) )
                m_align_FT_check->SetValue(false);
        }
        else if ( event.GetEventObject( ) == m_align_FT_check ) {
            if ( m_align_FT_check->IsChecked( ) )
                m_align_autocorr_check->SetValue(false);
        }
    }

    // Align the image using Auto-Correlation approach
    void AlignImageAutoCorr(Image& image, float pixel_size) {
        int x_dim = image.logical_x_dimension;
        int y_dim = image.logical_y_dimension;

        // 1. Prepare for alignment (Copy, Normalize, Mask, FFT, LowPass)
        Image search_image;
        search_image.Allocate(x_dim, y_dim, true);
        search_image.CopyFrom(&image);
        search_image.Normalize( );

        if ( x_dim > 0 )
            search_image.CosineMask(x_dim * 0.45f, 10.0f);

        search_image.ForwardFFT( );
        search_image.ZeroCentralPixel( );

        float low_pass_res = 20.0f;
        if ( pixel_size > 0 )
            search_image.GaussianLowPassFilter((pixel_size * 2.0f) / low_pass_res);

        // 2. Convert to Auto-Correlation (Complex = Amp^2 + 0i)
        int num_complex = search_image.real_memory_allocated / 2;
        for ( int i = 0; i < num_complex; ++i ) {
            float amp = std::abs(search_image.complex_values[i]);
            // Use simplified complex assignment
            search_image.complex_values[i] = std::complex<float>(amp * amp, 0.0f);
        }
        search_image.BackwardFFT( );
        search_image.SwapRealSpaceQuadrants( );
        search_image.object_is_centred_in_box = true;

        // 3. Coarse Search (-90 to 90, step 5)
        float best_psi = 0.0f;
        float best_sum = -FLT_MAX;

        Image temp_image;
        temp_image.Allocate(x_dim, y_dim, true);

        for ( float psi = -90.0f; psi <= 90.0f; psi += 4.0f ) {
            temp_image.CopyFrom(&search_image);
            temp_image.Rotate2DInPlace(psi, 0.0f); // Rotate the AC image

            float sum = GetMaxAbsColumnSum(&temp_image);
            if ( sum > best_sum ) {
                best_sum = sum;
                best_psi = psi;
            }
        }

        // 4. Fine Search (+/- 2.5 deg, step 0.5)
        float start_psi = best_psi - 2.0f;
        float end_psi   = best_psi + 2.0f;
        for ( float psi = start_psi; psi <= end_psi; psi += 0.5f ) {
            temp_image.CopyFrom(&search_image);
            temp_image.Rotate2DInPlace(psi, 0.0f);

            float sum = GetMaxAbsColumnSum(&temp_image);
            if ( sum > best_sum ) {
                best_sum = sum;
                best_psi = psi;
            }
        }

        // 5. Apply Alignment to Original Image
        // (AC max column sum implies vertical alignment)
        image.Rotate2DInPlace(best_psi, 0.0f);
    }

    // Align the image using Auto-Correlation approach
    void AlignImageFT(Image& image, float pixel_size) {
        int x_dim = image.logical_x_dimension;
        int y_dim = image.logical_y_dimension;

        float cosine_edge       = 10.0;
        float outside_weight    = 0.0;
        float filter_radius     = 0.0;
        float outside_value     = 0.0;
        bool  use_outside_value = false;
        // 1. Prepare for alignment (Copy, Normalize, Mask, FFT, LowPass)
        Image search_image;
        search_image.Allocate(x_dim, y_dim, true);
        search_image.CopyFrom(&image);
        search_image.Normalize( );

        if ( x_dim > 0 )
            search_image.CosineMask(x_dim * 0.45f, 10.0f);

        search_image.ForwardFFT( );
        search_image.ZeroCentralPixel( );

        float low_pass_res = 20.0f;
        if ( pixel_size > 0 )
            search_image.GaussianLowPassFilter((pixel_size * 2.0f) / low_pass_res);

        // // 2. Convert to Auto-Correlation (Complex = Amp^2 + 0i)
        // int num_complex = search_image.real_memory_allocated / 2;
        // for ( int i = 0; i < num_complex; ++i ) {
        //     float amp = std::abs(search_image.complex_values[i]);
        //     // Use simplified complex assignment
        //     search_image.complex_values[i] = std::complex<float>(amp * amp, 0.0f);
        // }
        // search_image.BackwardFFT( );
        // search_image.SwapRealSpaceQuadrants( );
        // search_image.object_is_centred_in_box = true;

        // 3. Coarse Search (-90 to 90, step 5)
        float best_psi = 0.0f;
        float best_sum = -FLT_MAX;

        Image temp_image;
        Image power_image;
        temp_image.Allocate(x_dim, y_dim, false);

        for ( float psi = -90.0f; psi <= 90.0f; psi += 4.0f ) {
            AnglesAndShifts rotation_angle;
            rotation_angle.Init(0.0, 0.0, psi, 0.0, 0.0);

            Image rotated_image;
            rotated_image.Allocate(x_dim, y_dim, false); // the rotated images real values will be the rotated FT
            rotated_image.SetToConstant(0.0);

            temp_image.CopyFrom(&search_image);
            temp_image.SwapRealSpaceQuadrants( );
            temp_image.RotateFourier2D(rotated_image, rotation_angle);

            // allocate memory for power image
            power_image.Allocate(x_dim, y_dim, true);
            power_image.SetToConstant(0.0);

            rotated_image.ComputeAmplitudeSpectrumFull2D(&power_image);
            rotated_image.Deallocate( );
            // Use a threshold to find the line corresponding to tube axis angle of rotation
            float image_average;
            float image_sd;
            float image_threshold;
            Image binary_mask;
            image_average = power_image.ReturnAverageOfRealValues( );
            image_sd      = sqrt(power_image.ReturnVarianceOfRealValues( ));
            // Threshold value is 2 sd away from mean to eliminate any outliers
            image_threshold = image_average + (2 * image_sd);
            binary_mask.CopyFrom(&power_image);
            binary_mask.Binarise(image_threshold);

            float filter_edge = 40.0;
            power_image.ApplyMask(binary_mask, cosine_edge / pixel_size, outside_weight, pixel_size / filter_radius, pixel_size / filter_edge, outside_value, use_outside_value);

            float sum = GetMaxAbsColumnSum(&power_image);
            if ( sum > best_sum ) {
                best_sum = sum;
                best_psi = psi;
            }
        }

        // 4. Fine Search (+/- 2.5 deg, step 0.5)
        float start_psi = best_psi - 2.0f;
        float end_psi   = best_psi + 2.0f;
        for ( float psi = start_psi; psi <= end_psi; psi += 0.5f ) {
            AnglesAndShifts rotation_angle;
            rotation_angle.Init(0.0, 0.0, psi, 0.0, 0.0);

            Image rotated_image;
            rotated_image.Allocate(x_dim, y_dim, false); // the rotated images real values will be the rotated FT
            rotated_image.SetToConstant(0.0);

            temp_image.CopyFrom(&search_image);
            temp_image.SwapRealSpaceQuadrants( );
            temp_image.RotateFourier2D(rotated_image, rotation_angle);

            // allocate memory for power image
            power_image.Allocate(x_dim, y_dim, true);
            power_image.SetToConstant(0.0);

            rotated_image.ComputeAmplitudeSpectrumFull2D(&power_image);
            rotated_image.Deallocate( );
            // Use a threshold to find the line corresponding to tube axis angle of rotation
            float image_average;
            float image_sd;
            float image_threshold;
            Image binary_mask;
            image_average = power_image.ReturnAverageOfRealValues( );
            image_sd      = sqrt(power_image.ReturnVarianceOfRealValues( ));
            // Threshold value is 2 sd away from mean to eliminate any outliers
            image_threshold = image_average + (2 * image_sd);
            binary_mask.CopyFrom(&power_image);
            binary_mask.Binarise(image_threshold);

            float filter_edge = 40.0;
            power_image.ApplyMask(binary_mask, cosine_edge / pixel_size, outside_weight, pixel_size / filter_radius, pixel_size / filter_edge, outside_value, use_outside_value);

            float sum = GetMaxAbsColumnSum(&temp_image);
            if ( sum > best_sum ) {
                best_sum = sum;
                best_psi = psi;
            }
        }

        // 5. Apply Alignment to Original Image
        // (AC max column sum implies vertical alignment)
        image.Rotate2DInPlace(best_psi + 90.0, 0.0f);
    }

    void CalculateAndDisplay( ) {
        if ( m_filename.IsEmpty( ) )
            return;

        MRCFile input_file(m_filename.ToStdString( ), false);
        if ( ! input_file.OpenFile(m_filename.ToStdString( ), false) )
            return;

        double pixel_size        = wxAtof(m_pixel_size_ctrl->GetValue( ));
        double lp_res            = wxAtof(m_lp_res_ctrl->GetValue( ));
        double mask_rad_ang      = wxAtof(m_mask_rad_ctrl->GetValue( ));
        bool   do_align_autocorr = m_align_autocorr_check->IsChecked( );
        bool   do_align_FT       = m_align_FT_check->IsChecked( );
        // bool   do_halfway         = m_chk_halfway->IsChecked( );
        bool do_invert_contrast = m_invert_contrast->IsChecked( );

        Image my_image;
        my_image.ReadSlice(&input_file, m_current_slice);

        bool particlesAreBlack = m_invert_contrast->GetValue( );

        if ( particlesAreBlack ) {
            my_image.InvertRealValues( );
        }

        int w = my_image.logical_x_dimension;
        int h = my_image.logical_y_dimension;

        SetStatusText(wxString::Format("Image %d / %d",
                                       m_current_slice,
                                       m_total_slices),
                      0);

        if ( lp_res > 0.0 ) {
            my_image.ForwardFFT( );
            my_image.GaussianLowPassFilter((float)((pixel_size * 2.0) / lp_res));
            my_image.BackwardFFT( );
        }

        // 1. Align Image if requested
        if ( do_align_autocorr ) {
            AlignImageAutoCorr(my_image, pixel_size);
        }
        if ( do_align_FT ) {
            AlignImageFT(my_image, pixel_size);
        }
        // set the auto corr image
        Image autocorr_image;
        autocorr_image.CopyFrom(&my_image);
        autocorr_image.ForwardFFT( );
        autocorr_image.ZeroCentralPixel( );

        // 2. Convert to Auto-Correlation (Complex = Amp^2 + 0i)
        int num_complex = autocorr_image.real_memory_allocated / 2;
        for ( int i = 0; i < num_complex; ++i ) {
            float amp = std::abs(autocorr_image.complex_values[i]);
            // Use simplified complex assignment
            autocorr_image.complex_values[i] = std::complex<float>(amp * amp, 0.0f);
        }
        autocorr_image.BackwardFFT( );
        autocorr_image.SwapRealSpaceQuadrants( );
        autocorr_image.object_is_centred_in_box = true;
        m_image_panel_autocorr->SetImage(&autocorr_image);

        //set the FT image
        Image FT_image;
        Image power_image;
        FT_image.CopyFrom(&my_image);
        FT_image.ForwardFFT( );
        FT_image.ZeroCentralPixel( );
        FT_image.SwapRealSpaceQuadrants( );
        power_image.Allocate(w, h, true);
        power_image.SetToConstant(0.0);

        FT_image.ComputeAmplitudeSpectrumFull2D(&power_image);

        if ( mask_rad_ang > 0 ) {
            float mask_rad_pix = mask_rad_ang;
            // int   w            = my_image.logical_x_dimension;
            // int   h            = my_image.logical_y_dimension;
            float cx = w / 2.0f;
            float cy = h / 2.0f;
            for ( int y = 0; y < h; y++ ) {
                for ( int x = 0; x < w; x++ ) {
                    float dx = x - cx;
                    float dy = y - cy;
                    if ( sqrt(dx * dx + dy * dy) > mask_rad_pix ) {
                        long addr                  = my_image.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
                        my_image.real_values[addr] = 0.0f;
                    }
                }
            }
        }

        std::vector<float> profile(w, 0.0f);

        for ( int x = 0; x < w; x++ ) {
            double sum = 0;
            for ( int y = 0; y < h; y++ ) {
                long addr = my_image.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
                sum += my_image.real_values[addr];
            }
            profile[x] = (float)sum;
        }

        double min_diam_pix = wxAtof(m_min_diam_ctrl->GetValue( ));
        double max_diam_pix = wxAtof(m_max_diam_ctrl->GetValue( ));

        // Convert to float, already in pixels
        std::pair<float, float> edges = FindOuterTubeEdges(profile, (float)min_diam_pix, (float)max_diam_pix, do_invert_contrast);
        //std::pair<float, float> edges = find_peak_strategy(profile, (float)min_diam_pix);
        //std::pair<float, float> edges = find_peak_strategy_scored(profile, (float)min_diam_pix, (float)max_diam_pix, false);

        m_image_panel->SetImage(&my_image);

        // Update the panels
        m_image_panel_autocorr->SetImage(&autocorr_image);
        m_image_panel_FT->SetImage(&power_image);

        m_image_panel->SetGraphData(profile, edges);

        if ( edges.first != -1 && edges.second != -1 ) {
            float diameter_pix = edges.second - edges.first;
            float diameter_ang = diameter_pix * pixel_size;

            m_result_text->SetLabel(wxString::Format("Tube Diameter: %.2f px (%.2f A) [Left: %.1f, Right: %.1f]",
                                                     std::abs(diameter_pix), std::abs(diameter_ang), edges.first, edges.second));
            m_result_text->GetParent( )->Layout( ); // forces sizer to recalc sizes
        }
        else {
            m_result_text->SetLabel("Could not detect tube within specified range.");
            m_result_text->GetParent( )->Layout( ); // forces sizer to recalc sizes
        }
    }

    // void OnShowHistogram(wxCommandEvent& event) {
    //     if ( m_filename.IsEmpty( ) )
    //         return;

    //     MRCFile input_file(m_filename.ToStdString( ), false);
    //     if ( ! input_file.OpenFile(m_filename.ToStdString( ), false) )
    //         return;

    //     wxProgressDialog progress("Calculating...", "Calculating diameters for all images...", m_total_slices, this, wxPD_APP_MODAL | wxPD_AUTO_HIDE | wxPD_REMAINING_TIME);

    //     double pixel_size   = wxAtof(m_pixel_size_ctrl->GetValue( ));
    //     double min_diam_pix = wxAtof(m_min_diam_ctrl->GetValue( ));
    //     double max_diam_pix = wxAtof(m_max_diam_ctrl->GetValue( ));
    //     double lp_res       = wxAtof(m_lp_res_ctrl->GetValue( ));
    //     double mask_rad_ang = wxAtof(m_mask_rad_ctrl->GetValue( ));
    //     bool   do_align_autocorr     = m_align_autocorr_check->IsChecked( );
    //     // bool   do_halfway         = m_chk_halfway->IsChecked( );
    //     bool do_invert_contrast = m_invert_contrast->IsChecked( );

    //     std::vector<float> diameters_ang;

    //     // Ensure we calculate for ALL slices in the file
    //     for ( int slice = 1; slice <= m_total_slices; ++slice ) {
    //         Image slice_image;
    //         slice_image.ReadSlice(&input_file, slice);

    //         slice_image.Normalize( );

    //         if ( do_invert_contrast ) {
    //             slice_image.InvertRealValues( );
    //         }

    //         if ( do_align_autocorr )
    //             AlignImageAutoCorr(slice_image, pixel_size);

    //         if ( lp_res > 0.0 ) {
    //             slice_image.ForwardFFT( );
    //             slice_image.GaussianLowPassFilter((float)((pixel_size * 2.0) / lp_res));
    //             slice_image.BackwardFFT( );
    //         }

    //         if ( mask_rad_ang > 0 ) {
    //             float mask_rad_pix = mask_rad_ang / pixel_size;
    //             int   w            = slice_image.logical_x_dimension;
    //             int   h            = slice_image.logical_y_dimension;
    //             float cx           = w / 2.0f;
    //             float cy           = h / 2.0f;
    //             for ( int y = 0; y < h; y++ ) {
    //                 for ( int x = 0; x < w; x++ ) {
    //                     float dx = x - cx;
    //                     float dy = y - cy;
    //                     if ( sqrt(dx * dx + dy * dy) > mask_rad_pix ) {
    //                         long addr                     = slice_image.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
    //                         slice_image.real_values[addr] = 0.0f;
    //                     }
    //                 }
    //             }
    //         }

    //         int w = slice_image.logical_x_dimension;
    //         int h = slice_image.logical_y_dimension;

    //         std::vector<float> profile(w, 0.0f);
    //         for ( int x = 0; x < w; x++ ) {
    //             double sum = 0;
    //             for ( int y = 0; y < h; y++ ) {
    //                 long addr = slice_image.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
    //                 sum += slice_image.real_values[addr];
    //             }
    //             profile[x] = (float)sum;
    //         }

    //         std::pair<float, float> edges = FindOuterTubeEdges(profile, (float)min_diam_pix, (float)max_diam_pix, do_invert_contrast);
    //         //std::pair<float, float> edges = find_peak_strategy(profile, (float)min_diam_pix);
    //         //std::pair<float, float> edges = find_peak_strategy_scored(profile, (float)min_diam_pix, (float)max_diam_pix, false);

    //         if ( edges.first != -1 && edges.second != -1 ) {
    //             float d_pix = edges.second - edges.first;
    //             diameters_ang.push_back(d_pix * pixel_size);
    //         }

    //         if ( ! progress.Update(slice) ) {
    //             break;
    //         }
    //     }

    //     if ( ! diameters_ang.empty( ) ) {
    //         // Plot Histogram
    //         wxDialog* dlg = new wxDialog(this, wxID_ANY, "Diameter Distribution", wxDefaultPosition, wxSize(600, 400), wxDEFAULT_DIALOG_STYLE | wxRESIZE_BORDER);

    //         PlotCurvePanel* plot = new PlotCurvePanel(dlg);
    //         Curve           hist_curve;

    //         // Bin the data
    //         float d_min = *std::min_element(diameters_ang.begin( ), diameters_ang.end( ));
    //         float d_max = *std::max_element(diameters_ang.begin( ), diameters_ang.end( ));

    //         // FIXED: Use a fixed bin width of ~2 pixels (in Angstroms)
    //         float bin_width_pix = 2.0f;
    //         float bin_width_ang = bin_width_pix * pixel_size;

    //         if ( d_max == d_min )
    //             d_max += bin_width_ang;

    //         int bins = (int)ceil((d_max - d_min) / bin_width_ang);
    //         if ( bins < 5 )
    //             bins = 5;

    //         // Simple histogram
    //         std::map<int, int> counts;
    //         int                max_count = 0;

    //         for ( float val : diameters_ang ) {
    //             int idx = (int)((val - d_min) / bin_width_ang);
    //             if ( idx >= bins )
    //                 idx = bins - 1;
    //             counts[idx]++;
    //             if ( counts[idx] > max_count )
    //                 max_count = counts[idx];
    //         }

    //         for ( int i = 0; i < bins; ++i ) {
    //             float x = d_min + i * bin_width_ang;
    //             float y = (float)counts[i];
    //             hist_curve.AddPoint(x, y);
    //             // Draw as bar-ish (add point at next step too)
    //             hist_curve.AddPoint(x + bin_width_ang, y);
    //         }

    //         hist_curve.SetupXAxis(d_min, d_max + bin_width_ang, bins);
    //         plot->AddCurve(hist_curve, *wxBLUE);
    //         plot->Draw( ); // Initial Draw

    //         wxBoxSizer* dlgSizer = new wxBoxSizer(wxVERTICAL);
    //         dlgSizer->Add(plot, 1, wxEXPAND | wxALL, 5);
    //         dlgSizer->Add(new wxButton(dlg, wxID_OK, "Close"), 0, wxALIGN_CENTER | wxALL, 5);

    //         dlg->SetSizer(dlgSizer);
    //         dlg->ShowModal( );
    //         dlg->Destroy( );
    //     }
    //     else {
    //         wxMessageBox("No valid diameters found in stack.", "Info", wxICON_INFORMATION);
    //     }
    // }

    void OnShowHistogram(wxCommandEvent& event) {
        if ( m_filename.IsEmpty( ) )
            return;

        MRCFile input_file(m_filename.ToStdString( ), false);
        if ( ! input_file.OpenFile(m_filename.ToStdString( ), false) )
            return;

        wxProgressDialog progress("Calculating...", "Calculating diameters for all images...", m_total_slices, this, wxPD_APP_MODAL | wxPD_AUTO_HIDE | wxPD_REMAINING_TIME);

        double pixel_size         = wxAtof(m_pixel_size_ctrl->GetValue( ));
        double min_diam_pix       = wxAtof(m_min_diam_ctrl->GetValue( ));
        double max_diam_pix       = wxAtof(m_max_diam_ctrl->GetValue( ));
        double lp_res             = wxAtof(m_lp_res_ctrl->GetValue( ));
        double mask_rad_ang       = wxAtof(m_mask_rad_ctrl->GetValue( ));
        bool   do_align_autocorr  = m_align_autocorr_check->IsChecked( );
        bool   do_align_FT        = m_align_FT_check->IsChecked( );
        bool   do_invert_contrast = m_invert_contrast->IsChecked( );

        std::vector<float> diameters_ang;

        // Calculate diameters for all slices
        for ( int slice = 1; slice <= m_total_slices; ++slice ) {
            Image slice_image;
            slice_image.ReadSlice(&input_file, slice);
            slice_image.Normalize( );

            if ( do_invert_contrast )
                slice_image.InvertRealValues( );

            if ( do_align_autocorr )
                AlignImageAutoCorr(slice_image, pixel_size);

            if ( do_align_FT )
                AlignImageFT(slice_image, pixel_size);

            // if ( lp_res > 0.0 ) {
            //     slice_image.ForwardFFT( );
            //     slice_image.GaussianLowPassFilter((float)((pixel_size * 2.0) / lp_res));
            //     slice_image.BackwardFFT( );
            // }

            if ( mask_rad_ang > 0 ) {
                float mask_rad_pix = mask_rad_ang / pixel_size;
                int   w            = slice_image.logical_x_dimension;
                int   h            = slice_image.logical_y_dimension;
                float cx           = w / 2.0f;
                float cy           = h / 2.0f;
                for ( int y = 0; y < h; y++ ) {
                    for ( int x = 0; x < w; x++ ) {
                        float dx = x - cx;
                        float dy = y - cy;
                        if ( sqrt(dx * dx + dy * dy) > mask_rad_pix ) {
                            long addr                     = slice_image.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
                            slice_image.real_values[addr] = 0.0f;
                        }
                    }
                }
            }

            int w = slice_image.logical_x_dimension;
            int h = slice_image.logical_y_dimension;

            std::vector<float> profile(w, 0.0f);
            for ( int x = 0; x < w; x++ ) {
                double sum = 0;
                for ( int y = 0; y < h; y++ ) {
                    long addr = slice_image.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
                    sum += slice_image.real_values[addr];
                }
                profile[x] = (float)sum;
            }

            std::pair<float, float> edges = FindOuterTubeEdges(profile, (float)min_diam_pix, (float)max_diam_pix, do_invert_contrast);

            if ( edges.first != -1 && edges.second != -1 ) {
                float d_pix = edges.second - edges.first;
                diameters_ang.push_back(d_pix * pixel_size);
            }

            if ( ! progress.Update(slice) )
                break;
        }

        if ( diameters_ang.empty( ) ) {
            wxMessageBox("No valid diameters found in stack.", "Info", wxICON_INFORMATION);
            return;
        }

        // --- Create histogram dialog ---
        wxDialog* dlg = new wxDialog(this, wxID_ANY, "Diameter Distribution", wxDefaultPosition, wxSize(600, 400));

        //PlotHistPanel* hist_panel = new PlotHistPanel(dlg);

        // Compute histogram bins
        float d_min         = *std::min_element(diameters_ang.begin( ), diameters_ang.end( ));
        float d_max         = *std::max_element(diameters_ang.begin( ), diameters_ang.end( ));
        float bin_width_ang = 2.0f * pixel_size; // ~2 px in Angstroms
        int   bins          = std::max(5, (int)ceil((d_max - d_min) / bin_width_ang));

        std::vector<int> counts(bins, 0);
        for ( float val : diameters_ang ) {
            int idx = (int)((val - d_min) / bin_width_ang);
            if ( idx >= bins )
                idx = bins - 1;
            counts[idx]++;
        }

        //hist_panel->SetHistogram(counts, d_min, bin_width_ang);
        //hist_panel->SetTitle("Tube Diameter Distribution");

        // --- Dialog layout ---
        // wxBoxSizer* dlgSizer = new wxBoxSizer(wxVERTICAL);
        // dlgSizer->Add(hist_panel, 1, wxEXPAND | wxALL, 5);
        // dlgSizer->Add(new wxButton(dlg, wxID_OK, "Close"), 0, wxALIGN_CENTER | wxALL, 5);
        // dlg->SetSizer(dlgSizer);
        // dlg->Layout( );

        // dlg->ShowModal( );
        // dlg->Destroy( );
    }

    void OnExit(wxCommandEvent& event) {
        Close(true);
    }

    DECLARE_EVENT_TABLE( )
};

BEGIN_EVENT_TABLE(FindTubeMainFrame, wxFrame)
EVT_MENU(ID_Open, FindTubeMainFrame::OnOpen)
EVT_MENU(wxID_EXIT, FindTubeMainFrame::OnExit)
EVT_BUTTON(ID_Next, FindTubeMainFrame::OnNext)
EVT_BUTTON(ID_Prev, FindTubeMainFrame::OnPrev)
EVT_BUTTON(ID_Random, FindTubeMainFrame::OnRandom)
EVT_MENU(ID_Hist, FindTubeMainFrame::OnShowHistogram)
END_EVENT_TABLE( )

bool FindTubeDiametersGuiApp::OnInit( ) {
    FindTubeMainFrame* frame = new FindTubeMainFrame( );
    frame->Show(true);
    return true;
}
