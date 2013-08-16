#include "FMRIWindow.h"

#include "MainFrame.h"
#include "SceneManager.h"
#include "SelectionBox.h"
#include "SelectionEllipsoid.h"
#include "../main.h"
#include "../dataset/Anatomy.h"
#include "../dataset/Fibers.h"
#include "../misc/IsoSurface/CIsoSurface.h"
#include "../misc/IsoSurface/TriangleMesh.h"

#include <wx/checkbox.h>
#include <wx/grid.h>
#include <wx/tglbtn.h>
#include <wx/treectrl.h>


IMPLEMENT_DYNAMIC_CLASS( FMRIWindow, wxScrolledWindow )

BEGIN_EVENT_TABLE( FMRIWindow, wxScrolledWindow )
END_EVENT_TABLE()


FMRIWindow::FMRIWindow( wxWindow *pParent, MainFrame *pMf, wxWindowID id, const wxPoint &pos, const wxSize &size)
:   wxScrolledWindow( pParent, id, pos, size, wxBORDER_NONE, _T("fMRI resting-state networks") ),
    m_pMainFrame( pMf )
{
    SetBackgroundColour( *wxLIGHT_GREY );
    m_pFMRISizer = new wxBoxSizer( wxVERTICAL );
    SetSizer( m_pFMRISizer );
    SetAutoLayout( true );

    m_pBtnSelectFMRI = new wxButton( this, wxID_ANY,wxT("Load resting-state"), wxPoint(30,0), wxSize(115, -1) );
	pMf->Connect( m_pBtnSelectFMRI->GetId(), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MainFrame::onLoad) );
    m_pBtnSelectFMRI->SetBackgroundColour(wxColour( 255, 147, 147 ));

	wxBoxSizer *pBoxRow1 = new wxBoxSizer( wxHORIZONTAL );
	pBoxRow1->Add( m_pBtnSelectFMRI, 0, wxALIGN_CENTER_HORIZONTAL | wxALL, 1 );
	m_pFMRISizer->Add( pBoxRow1, 0, wxFIXED_MINSIZE | wxALL, 2 );

}

void FMRIWindow::OnSize( wxSizeEvent &WXUNUSED(event) )
{
	
}

void FMRIWindow::OnPaint( wxPaintEvent &WXUNUSED(event) )
{
    wxPaintDC dc( this );
}

wxSizer* FMRIWindow::getWindowSizer()
{
    return m_pFMRISizer;
}

void FMRIWindow::OnSelectFMRI( wxCommandEvent& WXUNUSED(event) )
{

}
