#include "wx/wx.h"
#include "wx/laywin.h"
#include "wx/filedlg.h"
#include "wx/statbmp.h"

#include "mainFrame.h"

int winNumber = 1;

BEGIN_EVENT_TABLE(MainFrame, wxMDIParentFrame)
    EVT_MENU(VIEWER_ABOUT, MainFrame::OnAbout)
    EVT_SIZE(MainFrame::OnSize)
    EVT_MENU(VIEWER_QUIT, MainFrame::OnQuit)
    EVT_MENU(VIEWER_LOAD, MainFrame::OnLoad)
    EVT_MOUSE_EVENTS(MainFrame::OnMouseEvent)
END_EVENT_TABLE()


// Define my frame constructor
MainFrame::MainFrame(wxWindow *parent, const wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size,
    const long style):
  wxMDIParentFrame(parent, id, title, pos, size, style)
{
	// A window to the left of the client window
	wxSashLayoutWindow* win = new wxSashLayoutWindow(this, ID_WINDOW_LEFT,
                               wxDefaultPosition, wxSize(200, 30),
                               wxNO_BORDER | wxSW_3D | wxCLIP_CHILDREN);
	win->SetDefaultSize(wxSize(200, 1020));
	win->SetOrientation(wxLAYOUT_VERTICAL);
	win->SetAlignment(wxLAYOUT_LEFT);
	win->SetBackgroundColour(wxColour(0, 255, 0));
	win->SetSashVisible(wxSASH_RIGHT, false);
	
	win->SetExtraBorderSize(10);
	
	m_textWindow = new wxTextCtrl(win, wxID_ANY, wxEmptyString, 
			  wxDefaultPosition, wxDefaultSize,
			  wxTE_MULTILINE|wxSUNKEN_BORDER);
	m_textWindow->SetValue(_T("Stats"));
	m_leftWindow = win;
	
	// Main Window
	win = new wxSashLayoutWindow(this, ID_WINDOW_RIGHT, 
		  wxDefaultPosition, wxSize(200, 30),
		  wxNO_BORDER | wxSW_3D | wxCLIP_CHILDREN);
	win->SetDefaultSize(wxSize(1020, 765));
	win->SetOrientation(wxLAYOUT_VERTICAL);
	win->SetAlignment(wxLAYOUT_LEFT);
	win->SetBackgroundColour(wxColour(0,0,0));
	m_mainWindow = win;

  // navigation window with three sub windows for gl widgets 
    win = new wxSashLayoutWindow(m_mainWindow, ID_WINDOW_NAV_MAIN, 
  		  wxDefaultPosition, wxSize(200, 30),
  		  wxNO_BORDER | wxSW_3D | wxCLIP_CHILDREN);
    win->SetDefaultSize(wxSize(255, 765));
    win->SetOrientation(wxLAYOUT_VERTICAL);
    win->SetAlignment(wxLAYOUT_LEFT);
    win->SetBackgroundColour(wxColour(0, 0, 0));
    m_navWindow = win;

    // main window right side, holds the big gl widget
    win = new wxSashLayoutWindow(m_mainWindow, ID_WINDOW_RIGHT1, 
  		  wxDefaultPosition, wxSize(200, 30),
  		  wxNO_BORDER | wxSW_3D | wxCLIP_CHILDREN);
    win->SetDefaultSize(wxSize(765, 765));
    win->SetOrientation(wxLAYOUT_HORIZONTAL);
    win->SetAlignment(wxLAYOUT_RIGHT);
    win->SetBackgroundColour(wxColour(180, 180, 180));
    m_rightWindow = win;

    win = new wxSashLayoutWindow(m_navWindow, ID_WINDOW_NAV1, 
  		  wxDefaultPosition, wxSize(200, 30),
  		  wxRAISED_BORDER | wxSW_3D | wxCLIP_CHILDREN);
	win->SetDefaultSize(wxSize(255, 510));
    win->SetOrientation(wxLAYOUT_HORIZONTAL);
    win->SetAlignment(wxLAYOUT_TOP);
    win->SetBackgroundColour(wxColour(0, 0, 0));
    m_navWindow1 = win;
        
    win = new wxSashLayoutWindow(m_navWindow, ID_WINDOW_NAV2, 
  		  wxDefaultPosition, wxSize(200, 30),
  		  wxRAISED_BORDER | wxSW_3D | wxCLIP_CHILDREN);
    win->SetDefaultSize(wxSize(255, 255));
    win->SetOrientation(wxLAYOUT_HORIZONTAL);
    win->SetAlignment(wxLAYOUT_BOTTOM);
    win->SetBackgroundColour(wxColour(0, 0, 0));
    m_navWindow2 = win;

    win = new wxSashLayoutWindow(m_navWindow1, ID_WINDOW_NAV_X, 
  		  wxDefaultPosition, wxSize(200, 30),
  		  wxRAISED_BORDER | wxSW_3D | wxCLIP_CHILDREN);
    win->SetDefaultSize(wxSize(255, 255));
    win->SetOrientation(wxLAYOUT_HORIZONTAL);
    win->SetAlignment(wxLAYOUT_TOP);
    win->SetBackgroundColour(wxColour(0, 0, 0));
    m_topNavWindow = win;
    
    win = new wxSashLayoutWindow(m_navWindow1, ID_WINDOW_NAV_Y, 
  		  wxDefaultPosition, wxSize(200, 30),
  		  wxRAISED_BORDER | wxSW_3D | wxCLIP_CHILDREN);
    win->SetDefaultSize(wxSize(255, 255));
    win->SetOrientation(wxLAYOUT_HORIZONTAL);
    win->SetAlignment(wxLAYOUT_TOP);
    win->SetBackgroundColour(wxColour(0, 0, 0));
    m_middleNavWindow = win;

    win = new wxSashLayoutWindow(m_navWindow2, ID_WINDOW_NAV_Z, 
     		  wxDefaultPosition, wxSize(200, 30),
     		  wxRAISED_BORDER | wxSW_3D | wxCLIP_CHILDREN);
       win->SetDefaultSize(wxSize(255, 255));
       win->SetOrientation(wxLAYOUT_HORIZONTAL);
       win->SetAlignment(wxLAYOUT_TOP);
       win->SetBackgroundColour(wxColour(0, 0, 0));
    m_bottomNavWindow = win;
       
    win = new wxSashLayoutWindow(m_navWindow2, ID_WINDOW_NAV3, 
      		  wxDefaultPosition, wxSize(200, 30),
       		  wxNO_BORDER | wxSW_3D | wxCLIP_CHILDREN);
          win->SetDefaultSize(wxSize(255, 10));
          win->SetOrientation(wxLAYOUT_HORIZONTAL);
          win->SetAlignment(wxLAYOUT_TOP);
          win->SetBackgroundColour(wxColour(255, 255, 255));
    
      
    m_gl1 = new NavigationCanvas(m_topNavWindow, ID_GL_NAV_X, wxDefaultPosition,
    	        wxDefaultSize, 0, _T("MyGLCanvas"));
    m_gl2 = new NavigationCanvas(m_middleNavWindow, ID_GL_NAV_Y, wxDefaultPosition,
        	        wxDefaultSize, 0, _T("MyGLCanvas"));
    m_gl3 = new NavigationCanvas(m_bottomNavWindow, ID_GL_NAV_Z, wxDefaultPosition,
       	        wxDefaultSize, 0, _T("MyGLCanvas"));
    
    m_xclick = 0;
    m_yclick = 0;
    m_zclick = 0;
    
}

void MainFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
	Close(true);
}

void MainFrame::OnLoad(wxCommandEvent& WXUNUSED(event))
{
	wxString caption = wxT("Choose a file");
	wxString wildcard = wxT("Header files (*.hea)|*.hea|*.*|*.*");
	wxString defaultDir = wxEmptyString;
	wxString defaultFilename = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxOPEN);
	if (dialog.ShowModal() == wxID_OK)
	{
		wxString path = dialog.GetPath();
		m_dataset = new TheDataset();
		if (!m_dataset->load(path)) 
		{
			wxMessageBox(wxT("Fehler"),  wxT(""), wxNO_DEFAULT|wxYES_NO|wxCANCEL|wxICON_INFORMATION, NULL);
		}
		else 
		{ 
			
			m_gl1->setDataset(m_dataset, 0);
			m_gl2->setDataset(m_dataset, 1);
			m_gl3->setDataset(m_dataset, 2);
			
		}
		m_textWindow->SetValue(m_dataset->getInfoString());
		m_xclick = m_dataset->getColumns()/2;
		m_yclick = m_dataset->getRows()/2;
		m_zclick = m_dataset->getFrames()/2;
	}
}

void MainFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
      (void)wxMessageBox(_T("wxWidgets 2.0 Sash Demo\nAuthor: Julian Smart (c) 1998"), _T("About Sash Demo"));
}

void MainFrame::OnSashDrag(wxSashEvent& event)
{
    if (event.GetDragStatus() == wxSASH_STATUS_OUT_OF_RANGE)
        return;

    switch (event.GetId())
    {
    		case ID_WINDOW_LEFT:
            {
    	  if (m_leftWindow != 0)
    		  m_leftWindow->SetDefaultSize(wxSize(event.GetDragRect().width, GetSize().GetHeight()));
    	  if (m_rightWindow != 0)
    	    m_rightWindow->SetDefaultSize(wxSize(GetSize().GetWidth() - event.GetDragRect().width, GetSize().GetHeight()));
                break;
            }
            case ID_WINDOW_RIGHT:
            {
    	  if (m_rightWindow != 0)
      	    m_rightWindow->SetDefaultSize(wxSize(event.GetDragRect().width, GetSize().GetHeight()));
    	  if (m_leftWindow != 0)
    	    m_leftWindow->SetDefaultSize(wxSize(GetSize().GetWidth() - event.GetDragRect().width, GetSize().GetHeight()));
                break;
            }
    }

#if wxUSE_MDI_ARCHITECTURE
    wxLayoutAlgorithm layout;
    layout.LayoutMDIFrame(this);
#endif // wxUSE_MDI_ARCHITECTURE

    // Leaves bits of itself behind sometimes
    GetClientWindow()->Refresh();
}

void MainFrame::OnMouseEvent(wxMouseEvent& event)
{
	wxBell();
	wxPoint pt = event.GetPosition();
	m_yclick = pt.x;
	m_zclick = pt.y;
	switch (event.GetId())
	    {
	    	case ID_GL_NAV_X:
	    		m_yclick = pt.x;
	    		m_zclick = pt.y;
	    		break;
	    	case ID_GL_NAV_Y:
	    		m_xclick = pt.x;
	    		m_zclick = pt.y;
	    		break;
	    	case ID_GL_NAV_Z:
	    		m_xclick = pt.x;
	    		m_yclick = pt.y;
	    		break;
	    }
	
	wxString click = wxString::Format(wxT("\n\n\nx click: %d\ny click: %d\nz click: %d\n"),m_xclick, m_yclick, m_zclick); 
	if (m_dataset != NULL)
		m_textWindow->SetValue(m_dataset->getInfoString() + click);
	
	this->Refresh();
}

void MainFrame::OnSize(wxSizeEvent& WXUNUSED(event))
{
#if wxUSE_MDI_ARCHITECTURE
    wxLayoutAlgorithm layout;
    layout.LayoutMDIFrame(this);
#endif // wxUSE_MDI_ARCHITECTURE
    GetClientWindow()->Refresh();
    this->Refresh();
}
