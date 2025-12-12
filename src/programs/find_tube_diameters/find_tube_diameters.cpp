#include <wx/defs.h>
#include <wx/utils.h>
#include "../../core/gui_core_headers.h"
#include <wx/filedlg.h>
#include <wx/msgdlg.h>
#include <wx/menu.h>
#include <wx/sizer.h>
#include <wx/statusbr.h>
#include <wx/dcbuffer.h>
#include <wx/spinctrl.h>
#include <wx/textctrl.h>
#include <wx/stattext.h>
#include <wx/checkbox.h>
#include <wx/progdlg.h>
#include <vector>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <limits>
#include <map>

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

            // Reserve space at bottom for X-Axis labels (30 pixels)
            int axis_height = 30;
            int available_h = sz.GetHeight( ) - axis_height;
            if ( available_h < 0 )
                available_h = sz.GetHeight( ); // Fallback if too small

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

            // --- DRAW X-AXIS RULER ---
            dc.SetPen(*wxBLACK_PEN);
            dc.SetFont(wxFont(8, wxFONTFAMILY_SWISS, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL));

            int axis_y = y_off + draw_h + 2;
            dc.DrawLine(x_off, axis_y, x_off + draw_w, axis_y);

            int img_w         = m_image.GetWidth( );
            int tick_interval = 50;
            if ( img_w > 1000 )
                tick_interval = 100;
            if ( img_w < 100 )
                tick_interval = 10;

            float scale_x = (float)draw_w / img_w;

            for ( int i = 0; i < img_w; i += tick_interval ) {
                int screen_x = x_off + (int)(i * scale_x);
                dc.DrawLine(screen_x, axis_y, screen_x, axis_y + 5);

                wxSize text_sz = dc.GetTextExtent(wxString::Format("%d", i));
                dc.DrawText(wxString::Format("%d", i), screen_x - (text_sz.GetWidth( ) / 2), axis_y + 6);
            }
            // -------------------------

            if ( m_show_graph && m_has_graph && ! m_profile.empty( ) ) {
                float p_min   = *std::min_element(m_profile.begin( ), m_profile.end( ));
                float p_max   = *std::max_element(m_profile.begin( ), m_profile.end( ));
                float p_range = p_max - p_min;
                if ( p_range <= 0 )
                    p_range = 1.0f;

                float plot_h     = draw_h * 0.8f;
                float plot_off_y = draw_h * 0.1f;

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

    wxDECLARE_EVENT_TABLE( );
};

wxBEGIN_EVENT_TABLE(ImagePanel, wxPanel)
        EVT_PAINT(ImagePanel::OnPaint)
                EVT_SIZE(ImagePanel::OnSize)
                        wxEND_EVENT_TABLE( )

        // ----------------------------------------------------------------------------
        // Main Frame Class
        // ----------------------------------------------------------------------------
        class FindTubeMainFrame : public wxFrame {
  public:
    FindTubeMainFrame(const wxString& title, const wxPoint& pos, const wxSize& size);

    // Event Handlers
    void OnOpen(wxCommandEvent& event);
    void OnExit(wxCommandEvent& event);
    void OnParamChange(wxCommandEvent& event);
    void OnNext(wxCommandEvent& event);
    void OnPrev(wxCommandEvent& event);
    void OnRandom(wxCommandEvent& event);
    void OnToggleGraph(wxCommandEvent& event);
    void OnShowHistogram(wxCommandEvent& event);

    // Logic
    void                    LoadCurrentSlice( );
    void                    CalculateAndDisplay( );
    std::pair<float, float> FindOuterTubeEdges(const std::vector<float>& cols, float min_tube_diameter, float max_tube_diameter, bool use_half_way);
    void                    AlignImage(Image& image);
    float                   GetMaxAbsColumnSum(Image* img);

    std::vector<std::pair<int, float>> FindPeaks(const std::vector<float>& data, float min_dist, float threshold);

    // UI Elements
    wxTextCtrl* m_pixel_size_ctrl;
    wxTextCtrl* m_min_diam_ctrl;
    wxTextCtrl* m_max_diam_ctrl;
    wxTextCtrl* m_mask_rad_ctrl;
    wxTextCtrl* m_lp_res_ctrl;
    wxCheckBox* m_show_graph_check;
    wxCheckBox* m_align_check;
    wxCheckBox* m_half_way_check;

    wxStaticText* m_result_text;
    wxStaticText* m_slice_info;
    wxButton*     m_btn_next;
    wxButton*     m_btn_prev;
    wxButton*     m_btn_random;
    wxButton*     m_btn_hist;

    ImagePanel* m_image_panel;
    // PlotCurvePanel* m_hist_panel; // Not stored as member to avoid state issues

    // Data
    wxString m_current_filename;
    bool     m_image_loaded;
    int      m_current_slice;
    int      m_total_slices;

    wxDECLARE_EVENT_TABLE( );
};

enum {
    ID_Open      = 1,
    ID_Calc      = 2,
    ID_Next      = 3,
    ID_Prev      = 4,
    ID_Random    = 5,
    ID_ShowGraph = 6,
    ID_Align     = 7,
    ID_Hist      = 8,
    ID_HalfWay   = 9
};

wxBEGIN_EVENT_TABLE(FindTubeMainFrame, wxFrame)
        EVT_MENU(ID_Open, FindTubeMainFrame::OnOpen)
                EVT_MENU(wxID_EXIT, FindTubeMainFrame::OnExit)
                        EVT_TEXT_ENTER(wxID_ANY, FindTubeMainFrame::OnParamChange)
                                EVT_BUTTON(ID_Next, FindTubeMainFrame::OnNext)
                                        EVT_BUTTON(ID_Prev, FindTubeMainFrame::OnPrev)
                                                EVT_BUTTON(ID_Random, FindTubeMainFrame::OnRandom)
                                                        EVT_CHECKBOX(ID_ShowGraph, FindTubeMainFrame::OnToggleGraph)
                                                                EVT_CHECKBOX(ID_Align, FindTubeMainFrame::OnParamChange)
                                                                        EVT_CHECKBOX(ID_HalfWay, FindTubeMainFrame::OnParamChange)
                                                                                EVT_BUTTON(ID_Hist, FindTubeMainFrame::OnShowHistogram)
                                                                                        wxEND_EVENT_TABLE( )

                                                                                                FindTubeMainFrame::FindTubeMainFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
    : wxFrame(NULL, wxID_ANY, title, pos, size), m_image_loaded(false), m_current_slice(1), m_total_slices(0) {

    wxMenu* menuFile = new wxMenu;
    menuFile->Append(ID_Open, "&Open Image...\tCtrl-O");
    menuFile->AppendSeparator( );
    menuFile->Append(wxID_EXIT);
    wxMenuBar* menuBar = new wxMenuBar;
    menuBar->Append(menuFile, "&File");
    SetMenuBar(menuBar);

    wxPanel*    topPanel = new wxPanel(this);
    wxBoxSizer* topSizer = new wxBoxSizer(wxVERTICAL);

    wxBoxSizer* paramSizer = new wxBoxSizer(wxHORIZONTAL);
    auto        AddControl = [&](const wxString& label, wxTextCtrl*& ctrl, const wxString& defVal) {
        paramSizer->Add(new wxStaticText(topPanel, wxID_ANY, label), 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 10);
        ctrl = new wxTextCtrl(topPanel, wxID_ANY, defVal, wxDefaultPosition, wxSize(60, -1), wxTE_PROCESS_ENTER);
        paramSizer->Add(ctrl, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 5);
    };

    AddControl("Pixel Size (A):", m_pixel_size_ctrl, "1.0");
    AddControl("Min Diam (Pix):", m_min_diam_ctrl, "100.0");
    AddControl("Max Diam (Pix):", m_max_diam_ctrl, "300.0");
    AddControl("Mask Rad (A):", m_mask_rad_ctrl, "0.0");
    AddControl("Low Pass (A):", m_lp_res_ctrl, "50.0");

    m_show_graph_check = new wxCheckBox(topPanel, ID_ShowGraph, "Show Graph");
    m_show_graph_check->SetValue(true);
    paramSizer->Add(m_show_graph_check, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 15);

    m_align_check = new wxCheckBox(topPanel, ID_Align, "Align Images");
    m_align_check->SetValue(false);
    paramSizer->Add(m_align_check, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 15);

    m_half_way_check = new wxCheckBox(topPanel, ID_HalfWay, "Half-Way Peaks");
    m_half_way_check->SetValue(false);
    paramSizer->Add(m_half_way_check, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 15);

    wxButton* btnUpdate = new wxButton(topPanel, ID_Calc, "Update");
    btnUpdate->Bind(wxEVT_BUTTON, &FindTubeMainFrame::OnParamChange, this);
    paramSizer->Add(btnUpdate, 0, wxALIGN_CENTER_VERTICAL | wxLEFT | wxRIGHT, 10);

    m_btn_hist = new wxButton(topPanel, ID_Hist, "Show Histogram");
    m_btn_hist->Enable(false);
    paramSizer->Add(m_btn_hist, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 10);

    wxBoxSizer* navSizer = new wxBoxSizer(wxHORIZONTAL);
    m_slice_info         = new wxStaticText(topPanel, wxID_ANY, "Image: N/A");
    navSizer->Add(m_slice_info, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 10);

    navSizer->AddStretchSpacer(1);

    m_btn_prev = new wxButton(topPanel, ID_Prev, "Prev Image");
    m_btn_prev->Enable(false);
    navSizer->Add(m_btn_prev, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 10);

    m_btn_next = new wxButton(topPanel, ID_Next, "Next Image");
    m_btn_next->Enable(false);
    navSizer->Add(m_btn_next, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 10);

    m_btn_random = new wxButton(topPanel, ID_Random, "Random Image");
    m_btn_random->Enable(false);
    navSizer->Add(m_btn_random, 0, wxALIGN_CENTER_VERTICAL | wxLEFT | wxRIGHT, 10);

    topSizer->Add(paramSizer, 0, wxEXPAND | wxALL, 5);
    topSizer->Add(navSizer, 0, wxEXPAND | wxBOTTOM, 5);

    topPanel->SetSizer(topSizer);

    m_image_panel = new ImagePanel(this);

    m_result_text = new wxStaticText(this, wxID_ANY, "Tube Diameter: N/A");
    wxFont font   = m_result_text->GetFont( );
    font.SetWeight(wxFONTWEIGHT_BOLD);
    font.SetPointSize(12);
    m_result_text->SetFont(font);

    wxBoxSizer* mainSizer = new wxBoxSizer(wxVERTICAL);
    mainSizer->Add(topPanel, 0, wxEXPAND | wxALL, 5);
    mainSizer->Add(m_image_panel, 1, wxEXPAND | wxALL, 5);
    mainSizer->Add(m_result_text, 0, wxALIGN_CENTER | wxALL, 10);

    SetSizer(mainSizer);
    Layout( );
    CreateStatusBar( );

    srand(time(NULL));
}

void FindTubeMainFrame::OnOpen(wxCommandEvent& event) {
    wxFileDialog openFileDialog(this, "Open MRC file", "", "", "MRC files (*.mrc)|*.mrc", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if ( openFileDialog.ShowModal( ) == wxID_CANCEL )
        return;

    m_current_filename = openFileDialog.GetPath( );

    MRCFile input_file(m_current_filename.ToStdString( ), false);
    m_total_slices = input_file.ReturnNumberOfSlices( );

    if ( m_total_slices == 0 )
        return;

    float ps = input_file.ReturnPixelSize( );
    if ( ps > 0 )
        m_pixel_size_ctrl->SetValue(wxString::Format("%.2f", ps));

    m_current_slice = 1;
    m_image_loaded  = true;
    m_btn_prev->Enable(true);
    m_btn_next->Enable(true);
    m_btn_random->Enable(true);
    m_btn_hist->Enable(true);

    CalculateAndDisplay( );
}

void FindTubeMainFrame::LoadCurrentSlice( ) {
    if ( ! m_image_loaded )
        return;
    m_slice_info->SetLabel(wxString::Format("Image: %d / %d", m_current_slice, m_total_slices));
}

void FindTubeMainFrame::OnParamChange(wxCommandEvent& event) {
    CalculateAndDisplay( );
}

void FindTubeMainFrame::OnNext(wxCommandEvent& event) {
    if ( ! m_image_loaded )
        return;
    m_current_slice++;
    if ( m_current_slice > m_total_slices )
        m_current_slice = 1;
    CalculateAndDisplay( );
}

void FindTubeMainFrame::OnPrev(wxCommandEvent& event) {
    if ( ! m_image_loaded )
        return;
    m_current_slice--;
    if ( m_current_slice < 1 )
        m_current_slice = m_total_slices;
    CalculateAndDisplay( );
}

void FindTubeMainFrame::OnRandom(wxCommandEvent& event) {
    if ( ! m_image_loaded )
        return;
    if ( m_total_slices > 1 ) {
        m_current_slice = (rand( ) % m_total_slices) + 1;
    }
    CalculateAndDisplay( );
}

void FindTubeMainFrame::OnToggleGraph(wxCommandEvent& event) {
    if ( m_image_panel ) {
        m_image_panel->ShowGraph(event.IsChecked( ));
    }
}

void FindTubeMainFrame::OnShowHistogram(wxCommandEvent& event) {
    if ( ! m_image_loaded )
        return;

    wxProgressDialog progress("Generating Histogram",
                              "Calculating diameters for all images...",
                              m_total_slices,
                              this,
                              wxPD_AUTO_HIDE | wxPD_APP_MODAL | wxPD_REMAINING_TIME);

    double pixel_size   = wxAtof(m_pixel_size_ctrl->GetValue( ));
    double lp_res       = wxAtof(m_lp_res_ctrl->GetValue( ));
    double mask_rad_ang = wxAtof(m_mask_rad_ctrl->GetValue( ));
    double min_diam_pix = wxAtof(m_min_diam_ctrl->GetValue( ));
    double max_diam_pix = wxAtof(m_max_diam_ctrl->GetValue( ));
    bool   do_align     = m_align_check->IsChecked( );
    bool   do_halfway   = m_half_way_check->IsChecked( );

    if ( pixel_size <= 0 )
        pixel_size = 1.0;

    std::vector<float> diameters_ang;

    MRCFile input_file(m_current_filename.ToStdString( ), false);
    Image   process_image;

    for ( int slice = 1; slice <= m_total_slices; ++slice ) {
        process_image.ReadSlice(&input_file, slice);
        process_image.Normalize( );

        if ( do_align )
            AlignImage(process_image);

        if ( lp_res > 0.0 ) {
            process_image.ForwardFFT( );
            process_image.GaussianLowPassFilter((float)((pixel_size * 2.0) / lp_res));
            process_image.BackwardFFT( );
        }

        if ( mask_rad_ang > 0 ) {
            float mask_rad_pix = mask_rad_ang / pixel_size;
            int   w            = process_image.logical_x_dimension;
            int   h            = process_image.logical_y_dimension;
            float cx           = w / 2.0f;
            float cy           = h / 2.0f;
            for ( int y = 0; y < h; y++ ) {
                for ( int x = 0; x < w; x++ ) {
                    float dx = x - cx;
                    float dy = y - cy;
                    if ( sqrt(dx * dx + dy * dy) > mask_rad_pix ) {
                        long addr                       = process_image.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
                        process_image.real_values[addr] = 0.0f;
                    }
                }
            }
        }

        int                w = process_image.logical_x_dimension;
        int                h = process_image.logical_y_dimension;
        std::vector<float> profile(w, 0.0f);
        for ( int x = 0; x < w; x++ ) {
            double sum = 0;
            for ( int y = 0; y < h; y++ ) {
                long addr = process_image.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
                sum += process_image.real_values[addr];
            }
            profile[x] = (float)sum;
        }

        std::pair<float, float> edges = FindOuterTubeEdges(profile, (float)min_diam_pix, (float)max_diam_pix, do_halfway);

        if ( edges.first != -1 && edges.second != -1 ) {
            float d_pix = edges.second - edges.first;
            diameters_ang.push_back(d_pix * pixel_size);
        }

        if ( ! progress.Update(slice) ) {
            break;
        }
    }

    wxDialog*       dlg  = new wxDialog(this, wxID_ANY, "Diameter Distribution", wxDefaultPosition, wxSize(600, 400), wxDEFAULT_DIALOG_STYLE | wxRESIZE_BORDER);
    PlotCurvePanel* plot = new PlotCurvePanel(dlg);

    // FIXED: Moved declaration outside the if block so it persists during ShowModal
    Curve hist_curve;

    if ( ! diameters_ang.empty( ) ) {
        float d_min = *std::min_element(diameters_ang.begin( ), diameters_ang.end( ));
        float d_max = *std::max_element(diameters_ang.begin( ), diameters_ang.end( ));

        // FIXED: Use a fixed bin width of ~2 pixels (in Angstroms)
        float bin_width_pix = 2.0f;
        float bin_width_ang = bin_width_pix * pixel_size;

        // Handle case where range is zero (all particles identical)
        if ( d_max == d_min )
            d_max += bin_width_ang;

        // Calculate number of bins based on the 2-pixel width
        int bins = (int)ceil((d_max - d_min) / bin_width_ang);
        if ( bins < 1 )
            bins = 1;

        // Add a small buffer to max to catch the edge cases
        hist_curve.SetupXAxis(d_min, d_max + bin_width_ang, bins);

        std::vector<int> counts(bins, 0);
        for ( float val : diameters_ang ) {
            int idx = (int)((val - d_min) / bin_width_ang);
            if ( idx < 0 )
                idx = 0;
            if ( idx >= bins )
                idx = bins - 1;
            counts[idx]++;
        }

        for ( int i = 0; i < bins; i++ ) {
            hist_curve.data_y[i] = counts[i];
            // Center the bin value on X
            hist_curve.data_x[i] = d_min + (i * bin_width_ang) + (bin_width_ang / 2.0f);
        }

        plot->Initialise("Diameter (A)", "Count", false, true);
        plot->AddCurve(hist_curve, *wxBLUE);
    }
    else {
        plot->Initialise("Diameter (A)", "Count", false, true);
    }

    wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
    sizer->Add(plot, 1, wxEXPAND | wxALL, 5);
    dlg->SetSizer(sizer);
    dlg->ShowModal( );
    dlg->Destroy( );
}

void FindTubeMainFrame::OnExit(wxCommandEvent& event) {
    Close(true);
}

std::vector<std::pair<int, float>> FindTubeMainFrame::FindPeaks(const std::vector<float>& data, float min_dist, float threshold) {
    std::vector<std::pair<int, float>> peaks;
    int                                n = data.size( );
    if ( n < 3 )
        return peaks;

    for ( int i = 1; i < n - 1; ++i ) {
        if ( data[i] > data[i - 1] && data[i] > data[i + 1] ) {
            if ( data[i] >= threshold ) {
                peaks.push_back({i, data[i]});
            }
        }
    }

    std::sort(peaks.begin( ), peaks.end( ), [](const std::pair<int, float>& a, const std::pair<int, float>& b) {
        return a.second > b.second;
    });

    std::vector<std::pair<int, float>> filtered_peaks;
    for ( const auto& p : peaks ) {
        bool keep = true;
        for ( const auto& accepted : filtered_peaks ) {
            if ( std::abs(p.first - accepted.first) < min_dist ) {
                keep = false;
                break;
            }
        }
        if ( keep ) {
            filtered_peaks.push_back(p);
        }
    }

    std::sort(filtered_peaks.begin( ), filtered_peaks.end( ), [](const std::pair<int, float>& a, const std::pair<int, float>& b) {
        return a.first < b.first;
    });

    return filtered_peaks;
}

std::pair<float, float> FindTubeMainFrame::FindOuterTubeEdges(const std::vector<float>& cols, float min_tube_diameter, float max_tube_diameter, bool use_half_way) {
    int n = cols.size( );
    if ( n < 3 )
        return {-1.0f, -1.0f};

    // 1. Smooth the profile
    std::vector<float> smooth_cols   = cols;
    int                smooth_radius = 2;
    for ( int i = smooth_radius; i < n - smooth_radius; ++i ) {
        double sum = 0;
        for ( int k = -smooth_radius; k <= smooth_radius; ++k ) {
            sum += cols[i + k];
        }
        smooth_cols[i] = sum / (2 * smooth_radius + 1);
    }

    // 2. Normalize
    float              minVal = *std::min_element(smooth_cols.begin( ), smooth_cols.end( ));
    std::vector<float> norm(n);
    for ( int i = 0; i < n; ++i )
        norm[i] = smooth_cols[i] - minVal;

    float normMax = *std::max_element(norm.begin( ), norm.end( ));
    if ( normMax <= 0.0f )
        return {-1.0f, -1.0f};

    // Inverted profile for Negative peaks
    std::vector<float> normInv(n);
    for ( int i = 0; i < n; ++i )
        normInv[i] = normMax - norm[i];

    // 3. Find All Peaks
    float min_dist  = 5.0f; // Minimum distance between peaks of same type
    float threshold = 0.0f;

    // posPeaks = Candidates for PL and PR
    std::vector<std::pair<int, float>> posPeaks = FindPeaks(norm, min_dist, threshold);
    // negPeaks = Candidates for NL and NR (Inner Walls)
    std::vector<std::pair<int, float>> negPeaks = FindPeaks(normInv, min_dist, threshold);

    // 4. Search Pattern: NL -> PL ... PR <- NR
    float bestScore = -std::numeric_limits<float>::infinity( );
    int   best_NL = -1, best_NR = -1;
    int   best_PL = -1, best_PR = -1;

    float center_idx = (float)(n - 1) / 2.0f;

    // Iterate through all possible Left Negative Peaks (NL)
    for ( const auto& pNL : negPeaks ) {
        int   idx_NL = pNL.first;
        float val_NL = pNL.second;

        // Iterate through all possible Right Negative Peaks (NR)
        for ( const auto& pNR : negPeaks ) {
            int   idx_NR = pNR.first;
            float val_NR = pNR.second;

            // Basic geometric constraints
            if ( idx_NR <= idx_NL )
                continue; // Right must be to the right
            float width = idx_NR - idx_NL;
            if ( width < min_tube_diameter || width > max_tube_diameter )
                continue;

            // Find best PL: Highest Positive peak strictly between NL and Center
            int   idx_PL     = -1;
            float max_val_PL = -1.0f;

            for ( const auto& pPos : posPeaks ) {
                if ( pPos.first > idx_NL && pPos.first < (idx_NL + idx_NR) / 2.0f ) {
                    if ( pPos.second > max_val_PL ) {
                        max_val_PL = pPos.second;
                        idx_PL     = pPos.first;
                    }
                }
            }

            // Find best PR: Highest Positive peak strictly between Center and NR
            int   idx_PR     = -1;
            float max_val_PR = -1.0f;

            for ( const auto& pPos : posPeaks ) {
                if ( pPos.first > (idx_NL + idx_NR) / 2.0f && pPos.first < idx_NR ) {
                    if ( pPos.second > max_val_PR ) {
                        max_val_PR = pPos.second;
                        idx_PR     = pPos.first;
                    }
                }
            }

            // Require both Positive peaks to exist for this pattern
            if ( idx_PL != -1 && idx_PR != -1 ) {

                // --- SCORING ---
                float score = 0.0f;

                // 1. Magnitude Score (Sum of all 4 peaks)
                score += (val_NL + val_NR + max_val_PL + max_val_PR);

                // 2. Symmetry Penalty (Tube should be roughly centered)
                float midpoint = (float)(idx_NL + idx_NR) / 2.0f;
                score -= 5.0f * std::abs(midpoint - center_idx) / n;

                // 3. Wall Thickness Consistency Penalty
                float left_wall_w  = idx_PL - idx_NL;
                float right_wall_w = idx_NR - idx_PR;
                score -= 2.0f * std::abs(left_wall_w - right_wall_w);

                if ( score > bestScore ) {
                    bestScore = score;
                    best_NL   = idx_NL;
                    best_NR   = idx_NR;
                    best_PL   = idx_PL;
                    best_PR   = idx_PR;
                }
            }
        }
    }

    // 5. Return Results
    if ( best_NL != -1 && best_NR != -1 && best_PL != -1 && best_PR != -1 ) {
        if ( use_half_way ) {
            // Average of Inner (Neg) and Outer (Pos) wall positions
            float edge_L = (float)(best_NL + best_PL) / 2.0f;
            float edge_R = (float)(best_NR + best_PR) / 2.0f;
            return {edge_L, edge_R};
        }
        else {
            // Return the Inner Walls (Negative peaks)
            return {(float)best_NL, (float)best_NR};
        }
    }

    return {-1.0f, -1.0f};
}

// Simple absolute column sum logic
float FindTubeMainFrame::GetMaxAbsColumnSum(Image* img) {
    float max_sum = -FLT_MAX;
    int   w       = img->logical_x_dimension;
    int   h       = img->logical_y_dimension;

    for ( int x = 0; x < w; x++ ) {
        float sum = 0.0f;
        for ( int y = 0; y < h; y++ ) {
            long addr = img->ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
            sum += img->real_values[addr];
        }
        if ( std::abs(sum) > max_sum )
            max_sum = std::abs(sum);
    }
    return max_sum;
}

// Basic Alignment logic: Rotate -90 to +90 to find vertical alignment
void FindTubeMainFrame::AlignImage(Image& image) {
    float best_psi = 0.0f;
    float best_sum = -FLT_MAX;

    Image temp_image;
    temp_image.Allocate(image.logical_x_dimension, image.logical_y_dimension, true);

    // Coarse Search: -90 to 90 degrees, step 5
    for ( float psi = -90.0f; psi <= 90.0f; psi += 5.0f ) {
        temp_image.CopyFrom(&image);
        // Padding with FLT_MAX usually means background/average padding in Rotate2DInPlace
        temp_image.Rotate2DInPlace(psi, 0.0);
        float current_sum = GetMaxAbsColumnSum(&temp_image);

        if ( current_sum > best_sum ) {
            best_sum = current_sum;
            best_psi = psi;
        }
    }

    // Fine Search: +/- 5 degrees, step 1
    float coarse_psi = best_psi;
    for ( float psi = coarse_psi - 5.0f; psi <= coarse_psi + 5.0f; psi += 1.0f ) {
        temp_image.CopyFrom(&image);
        temp_image.Rotate2DInPlace(psi, 0.0);
        float current_sum = GetMaxAbsColumnSum(&temp_image);

        if ( current_sum > best_sum ) {
            best_sum = current_sum;
            best_psi = psi;
        }
    }

    // Apply best rotation to original image
    image.Rotate2DInPlace(best_psi, 0.0);
}

void FindTubeMainFrame::CalculateAndDisplay( ) {
    if ( ! m_image_loaded )
        return;

    LoadCurrentSlice( );

    SetStatusText("Calculating...");

    double pixel_size   = wxAtof(m_pixel_size_ctrl->GetValue( ));
    double lp_res       = wxAtof(m_lp_res_ctrl->GetValue( ));
    double mask_rad_ang = wxAtof(m_mask_rad_ctrl->GetValue( ));
    bool   do_align     = m_align_check->IsChecked( );
    bool   do_halfway   = m_half_way_check->IsChecked( );

    if ( pixel_size <= 0 )
        pixel_size = 1.0;

    Image   process_image;
    MRCFile input_file(m_current_filename.ToStdString( ), false);
    process_image.ReadSlice(&input_file, m_current_slice);
    process_image.Normalize( );

    // 1. Align Image if requested
    if ( do_align ) {
        AlignImage(process_image);
    }

    if ( lp_res > 0.0 ) {
        process_image.ForwardFFT( );
        process_image.GaussianLowPassFilter((float)((pixel_size * 2.0) / lp_res));
        process_image.BackwardFFT( );
    }

    if ( mask_rad_ang > 0 ) {
        float mask_rad_pix = mask_rad_ang / pixel_size;
        int   w            = process_image.logical_x_dimension;
        int   h            = process_image.logical_y_dimension;
        float cx           = w / 2.0f;
        float cy           = h / 2.0f;
        for ( int y = 0; y < h; y++ ) {
            for ( int x = 0; x < w; x++ ) {
                float dx = x - cx;
                float dy = y - cy;
                if ( sqrt(dx * dx + dy * dy) > mask_rad_pix ) {
                    long addr                       = process_image.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
                    process_image.real_values[addr] = 0.0f;
                }
            }
        }
    }

    m_image_panel->SetImage(&process_image);

    int w = process_image.logical_x_dimension;
    int h = process_image.logical_y_dimension;

    std::vector<float> profile(w, 0.0f);

    for ( int x = 0; x < w; x++ ) {
        double sum = 0;
        for ( int y = 0; y < h; y++ ) {
            long addr = process_image.ReturnReal1DAddressFromPhysicalCoord(x, y, 0);
            sum += process_image.real_values[addr];
        }
        profile[x] = (float)sum;
    }

    double min_diam_pix = wxAtof(m_min_diam_ctrl->GetValue( ));
    double max_diam_pix = wxAtof(m_max_diam_ctrl->GetValue( ));

    // Convert to float, already in pixels
    std::pair<float, float> edges = FindOuterTubeEdges(profile, (float)min_diam_pix, (float)max_diam_pix, do_halfway);

    m_image_panel->SetGraphData(profile, edges);

    if ( edges.first != -1 && edges.second != -1 ) {
        float diameter_pix = edges.second - edges.first;
        float diameter_ang = diameter_pix * pixel_size;

        m_result_text->SetLabel(wxString::Format("Tube Diameter: %.2f px (%.2f A) [Left: %.1f, Right: %.1f]",
                                                 diameter_pix, diameter_ang, edges.first, edges.second));
    }
    else {
        m_result_text->SetLabel("Tube Diameter: Not Found");
    }

    SetStatusText("Done.");
}

bool FindTubeDiametersGuiApp::OnInit( ) {
    FindTubeMainFrame* frame = new FindTubeMainFrame("Find Tube Diameters", wxPoint(50, 50), wxSize(1200, 800));
    frame->Show(true);
    return true;
}