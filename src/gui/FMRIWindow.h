/////////////////////////////////////////////////////////////////////////////
// Name:            FMRIWindow.h
// Author:          ---
// Creation Date:   ---
//
// Description: This is the implementation file for the RT-fMRI exploration method.
/////////////////////////////////////////////////////////////////////////////


#ifndef FMRIWINDOW_H_
#define FMRIWINDOW_H_

#include "MainCanvas.h"
#include "MyListCtrl.h"

#include "../misc/Algorithms/Helper.h"

#include <wx/scrolwin.h>
#include <wx/statline.h>

class MainFrame;
class wxToggleButton;

class FMRIWindow: public wxScrolledWindow
{
public:
    FMRIWindow(){};
    FMRIWindow( wxWindow *pParent, MainFrame *pMf, wxWindowID id, const wxPoint &pos, const wxSize &size );

    ~FMRIWindow(){};
    void OnPaint( wxPaintEvent &event );
    void OnSize( wxSizeEvent &event );
    wxSizer* getWindowSizer();

public:
	void OnSelectFMRI                     ( wxCommandEvent& event );

private:
    MainFrame           *m_pMainFrame;
	wxButton            *m_pBtnSelectFMRI;
    
private:
    wxSizer *m_pFMRISizer;
    FMRIWindow( wxWindow *pParent, wxWindowID id, const wxPoint &pos, const wxSize &size );
    DECLARE_DYNAMIC_CLASS( FMRIWindow )
    DECLARE_EVENT_TABLE()
};

#endif /*FMRIWINDOW*/
