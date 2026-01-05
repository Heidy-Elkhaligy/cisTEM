#include <wx/defs.h>
#include <wx/utils.h>
#include "../../core/gui_core_headers.h"
#include <wx/filedlg.h>
#include <wx/msgdlg.h>
#include <wx/menu.h>
#include <wx/sizer.h>
#include <wx/statusbr.h>
#include <wx/dcbuffer.h>
#include <cfloat>

// ----------------------------------------------------------------------------
// Application Class
// ----------------------------------------------------------------------------
class AutocorrelationGuiApp : public wxApp {
  public:
    virtual bool OnInit( );
};

IMPLEMENT_APP(AutocorrelationGuiApp)

// ----------------------------------------------------------------------------
// Image Panel Class (Simple Viewer)
// ----------------------------------------------------------------------------
class ImagePanel : public wxPanel {
    wxImage m_image;
    bool    m_has_image;

  public:
    ImagePanel(wxWindow* parent) : wxPanel(parent), m_has_image(false) {
        SetBackgroundStyle(wxBG_STYLE_PAINT);
    }

    void SetImage(Image* image) {
        if ( ! image || image->logical_x_dimension == 0 )
            return;

        int w = image->logical_x_dimension;
        int h = image->logical_y_dimension;

        // 1. Find Min/Max for scaling
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

        // 2. Create RGB data
        unsigned char* rgb_data = (unsigned char*)malloc(w * h * 3);
        unsigned char* ptr      = rgb_data;

        for ( int y = 0; y < h; y++ ) {
            for ( int x = 0; x < w; x++ ) {
                float val = image->real_values[image->ReturnReal1DAddressFromPhysicalCoord(x, y, 0)];
                // Linear scale to 0-255
                unsigned char pixel_val = (unsigned char)((val - min_val) / range * 255.0f);
                *ptr++                  = pixel_val; // R
                *ptr++                  = pixel_val; // G
                *ptr++                  = pixel_val; // B
            }
        }

        // 3. Store as wxImage (takes ownership of rgb_data)
        m_image     = wxImage(w, h, rgb_data);
        m_has_image = true;
        Refresh( );
    }

    void OnPaint(wxPaintEvent& evt) {
        wxAutoBufferedPaintDC dc(this);
        dc.Clear( );

        if ( m_has_image && m_image.IsOk( ) ) {
            wxSize sz = GetClientSize( );
            if ( sz.GetWidth( ) <= 0 || sz.GetHeight( ) <= 0 )
                return;

            // Maintain aspect ratio
            float img_aspect   = (float)m_image.GetWidth( ) / m_image.GetHeight( );
            float panel_aspect = (float)sz.GetWidth( ) / sz.GetHeight( );

            int draw_w, draw_h;

            if ( panel_aspect > img_aspect ) {
                // Panel is wider relative to image -> limit by height
                draw_h = sz.GetHeight( );
                draw_w = (int)(draw_h * img_aspect);
            }
            else {
                // Panel is taller relative to image -> limit by width
                draw_w = sz.GetWidth( );
                draw_h = (int)(draw_w / img_aspect);
            }

            // Center the image
            int x_off = (sz.GetWidth( ) - draw_w) / 2;
            int y_off = (sz.GetHeight( ) - draw_h) / 2;

            // Scale and Draw
            // Use High Quality if image is smallish to look nice, otherwise Normal for speed
            wxBitmap bmp(m_image.Scale(draw_w, draw_h, wxIMAGE_QUALITY_NORMAL));
            dc.DrawBitmap(bmp, x_off, y_off);
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
        class AutocorrelationMainFrame : public wxFrame {
  public:
    AutocorrelationMainFrame(const wxString& title, const wxPoint& pos, const wxSize& size);

    // Event Handlers
    void OnOpen(wxCommandEvent& event);
    void OnSaveGraph(wxCommandEvent& event);
    void OnSaveImage(wxCommandEvent& event); // Handler for saving the MRC image
    void OnExit(wxCommandEvent& event);

    // UI Elements
    PlotCurvePanel* plot_panel;
    ImagePanel*     image_panel; // New panel for displaying the image

    // Data
    Image calculated_image; // Member to store the autocorrelation image

    // Helper to calculate autocorrelation on a single image
    void CalculateAndPlot(wxString filename);

    wxDECLARE_EVENT_TABLE( );
};

// ----------------------------------------------------------------------------
// Event Table
// ----------------------------------------------------------------------------
enum {
    ID_Open = 1,
    ID_SaveGraph,
    ID_SaveImage
};

wxBEGIN_EVENT_TABLE(AutocorrelationMainFrame, wxFrame)
        EVT_MENU(ID_Open, AutocorrelationMainFrame::OnOpen)
                EVT_MENU(ID_SaveGraph, AutocorrelationMainFrame::OnSaveGraph)
                        EVT_MENU(ID_SaveImage, AutocorrelationMainFrame::OnSaveImage)
                                EVT_MENU(wxID_EXIT, AutocorrelationMainFrame::OnExit)
                                        wxEND_EVENT_TABLE( )

        // ----------------------------------------------------------------------------
        // Main Frame Implementation
        // ----------------------------------------------------------------------------
        AutocorrelationMainFrame::AutocorrelationMainFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
    : wxFrame(NULL, wxID_ANY, title, pos, size) {

    // 1. Setup Menu Bar
    wxMenu* menuFile = new wxMenu;
    menuFile->Append(ID_Open, "&Open MRC File...\tCtrl-O", "Open an image stack to calculate autocorrelation");
    menuFile->Append(ID_SaveGraph, "&Save Graph as PNG...", "Save the current plot as a PNG image");
    menuFile->Append(ID_SaveImage, "Save Autocorrelation &Image...", "Save the calculated autocorrelation image as an MRC file");
    menuFile->AppendSeparator( );
    menuFile->Append(wxID_EXIT);

    wxMenuBar* menuBar = new wxMenuBar;
    menuBar->Append(menuFile, "&File");
    SetMenuBar(menuBar);

    // 2. Setup Panels (Side by Side)
    image_panel = new ImagePanel((wxWindow*)this);
    plot_panel  = new PlotCurvePanel((wxWindow*)this);

    // 3. Layout
    wxBoxSizer* mainSizer = new wxBoxSizer(wxHORIZONTAL);
    // Add Image Panel (Left, 50% width)
    mainSizer->Add(image_panel, 1, wxEXPAND | wxALL, 5);
    // Add Plot Panel (Right, 50% width)
    mainSizer->Add(plot_panel, 1, wxEXPAND | wxALL, 5);

    this->SetSizer(mainSizer);
    this->Layout( );

    // 4. Status Bar
    CreateStatusBar( );
    SetStatusText("Ready. Go to File -> Open to load an image.");

    // Initialise the plot with empty axes
    plot_panel->Initialise("Pixel Position", "Autocorrelation Amplitude", false, true);
    plot_panel->Draw( );
}

void AutocorrelationMainFrame::OnOpen(wxCommandEvent& event) {
    wxFileDialog openFileDialog(this, "Open MRC file", "", "", "MRC files (*.mrc)|*.mrc", wxFD_OPEN | wxFD_FILE_MUST_EXIST);

    if ( openFileDialog.ShowModal( ) == wxID_CANCEL )
        return;

    wxString filename = openFileDialog.GetPath( );
    CalculateAndPlot(filename);
}

void AutocorrelationMainFrame::CalculateAndPlot(wxString filename) {
    SetStatusText("Calculating autocorrelation...");

    // --- LOGIC FROM calculate_autocorrelation.cpp ADAPTED FOR GUI ---

    // 1. Open the file
    MRCFile input_file(filename.ToStdString( ), false);

    if ( input_file.ReturnNumberOfSlices( ) == 0 ) {
        wxMessageBox("Error: File contains no images.", "Error", wxICON_ERROR);
        return;
    }

    // 2. Read the first image (Slice 1) into the member variable
    // We reuse the member variable 'calculated_image' so it can be saved later
    calculated_image.ReadSlice(&input_file, 1);

    // 3. Process the image (The Autocorrelation Sequence)
    // Normalize
    calculated_image.Normalize( );
    // FFT
    calculated_image.ForwardFFT( );
    // Zero Central Pixel
    calculated_image.ZeroCentralPixel( );

    // Compute Power Spectrum (Amplitude Squared) in Fourier Space
    // Note: iterating up to real_memory_allocated/2 for complex array
    for ( long pixel_counter = 0; pixel_counter < calculated_image.real_memory_allocated / 2; pixel_counter++ ) {
        float amplitude = abs(calculated_image.complex_values[pixel_counter]);
        // Set phase to 0, real part = amplitude^2
        calculated_image.complex_values[pixel_counter] = amplitude * amplitude + I * 0.0f;
    }

    // IFFT to get back to Real Space
    calculated_image.BackwardFFT( );
    // Swap Quadrants to center the correlation peak
    calculated_image.SwapRealSpaceQuadrants( );

    // 4. Update the Image Panel
    image_panel->SetImage(&calculated_image);

    // --- EXTRACTION OF 1D PROFILE ---

    int x_size   = calculated_image.logical_x_dimension;
    int y_size   = calculated_image.logical_y_dimension;
    int center_y = y_size / 2; // The central row index

    // Create a curve to store the data
    Curve autocorrelation_curve;

    // Setup X Axis (0 to width)
    // Assuming Pixel size is 1.0 for the graph, or read from file if needed
    float pixel_size = input_file.ReturnPixelSize( );
    autocorrelation_curve.SetupXAxis(0.0, x_size * pixel_size, x_size);

    // Extract the central row (varying X, fixed Y)
    for ( int x = 0; x < x_size; x++ ) {
        // Calculate index for the pixel at (x, center_y) using the correct cisTEM function
        long  address     = calculated_image.ReturnReal1DAddressFromPhysicalCoord(x, center_y, 0);
        float pixel_value = calculated_image.real_values[address];

        // Add to curve.
        autocorrelation_curve.data_y[x] = pixel_value;
    }

    // 5. Update the Plot Panel
    plot_panel->Clear( ); // Remove old curves
    plot_panel->Initialise("Position (A)", "Autocorrelation (Arbitrary)", false, true);
    plot_panel->AddCurve(autocorrelation_curve, *wxBLUE);
    plot_panel->Draw( );

    SetStatusText(wxString::Format("Plotted central row of %s", filename));
}

void AutocorrelationMainFrame::OnSaveGraph(wxCommandEvent& event) {
    wxFileDialog saveFileDialog(this, "Save PNG", "", "graph.png", "PNG files (*.png)|*.png", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);

    if ( saveFileDialog.ShowModal( ) == wxID_CANCEL )
        return;

    // Capture the content of the plot_panel
    wxWindowDC clientDC(plot_panel);
    wxBitmap   bitmap(plot_panel->GetSize( ).GetWidth( ), plot_panel->GetSize( ).GetHeight( ));

    wxMemoryDC memDC;
    memDC.SelectObject(bitmap);

    // Blit the window DC to the Memory DC (Bitmap)
    memDC.Blit(0, 0, plot_panel->GetSize( ).GetWidth( ), plot_panel->GetSize( ).GetHeight( ), &clientDC, 0, 0);
    memDC.SelectObject(wxNullBitmap); // Release the bitmap

    if ( bitmap.SaveFile(saveFileDialog.GetPath( ), wxBITMAP_TYPE_PNG) ) {
        SetStatusText("Graph saved successfully.");
    }
    else {
        wxMessageBox("Failed to save image.", "Error", wxICON_ERROR);
    }
}

void AutocorrelationMainFrame::OnSaveImage(wxCommandEvent& event) {
    if ( calculated_image.logical_x_dimension == 0 ) {
        wxMessageBox("No image calculated yet. Please Open a file first.", "Error", wxICON_ERROR);
        return;
    }

    wxFileDialog saveFileDialog(this, "Save Autocorrelation MRC", "", "autocorrelation.mrc", "MRC files (*.mrc)|*.mrc", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);

    if ( saveFileDialog.ShowModal( ) == wxID_CANCEL )
        return;

    std::string output_filename = saveFileDialog.GetPath( ).ToStdString( );

    // Create an MRC file for writing (second argument true = overwrite/create)
    MRCFile output_file(output_filename, true);

    // Write the stored image to the file
    calculated_image.WriteSlice(&output_file, 1);

    SetStatusText("Autocorrelation image saved.");
}

void AutocorrelationMainFrame::OnExit(wxCommandEvent& event) {
    Close(true);
}

// ----------------------------------------------------------------------------
// App Initialization
// ----------------------------------------------------------------------------
bool AutocorrelationGuiApp::OnInit( ) {
    AutocorrelationMainFrame* frame = new AutocorrelationMainFrame("Autocorrelation 1D Profile", wxPoint(50, 50), wxSize(1000, 600));
    frame->Show(true);
    return true;
}