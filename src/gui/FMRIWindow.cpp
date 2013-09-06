#include "FMRIWindow.h"

#include "MainFrame.h"
#include "SceneManager.h"
#include "SelectionBox.h"
#include "SelectionEllipsoid.h"
#include "../main.h"
#include "../dataset/Anatomy.h"
#include "../dataset/Fibers.h"
#include "../dataset/RestingStateNetwork.h"
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

    m_pBtnSelectFMRI = new wxButton( this, wxID_ANY,wxT("Load resting-state"), wxPoint(30,0), wxSize(150, -1) );
	pMf->Connect( m_pBtnSelectFMRI->GetId(), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MainFrame::onLoadAsRestingState) );
    m_pBtnSelectFMRI->SetBackgroundColour(wxColour( 255, 147, 147 ));

	wxBoxSizer *pBoxRow1 = new wxBoxSizer( wxHORIZONTAL );
	pBoxRow1->Add( m_pBtnSelectFMRI, 0, wxALIGN_CENTER_HORIZONTAL | wxALL, 1 );
	m_pFMRISizer->Add( pBoxRow1, 0, wxFIXED_MINSIZE | wxALL, 2 );

	m_pTextRest = new wxStaticText( this, wxID_ANY, wxT("Volume"), wxPoint(0,30), wxSize(60, -1), wxALIGN_CENTER );
    m_pSliderRest = new MySlider( this, wxID_ANY, 0, 0, 107, wxPoint(60,30), wxSize(130, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
    //m_pSliderRest->SetValue( 0 );
    Connect( m_pSliderRest->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(FMRIWindow::OnSliderRestMoved) );
    m_pTxtRestBox = new wxTextCtrl( this, wxID_ANY, wxT("0"), wxPoint(190,30), wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

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

void FMRIWindow::SetSelectButton()
{
	DatasetIndex indx = DatasetManager::getInstance()->m_pRestingStateNetwork->getIndex();
	Anatomy* pNewAnatomy = (Anatomy *)DatasetManager::getInstance()->getDataset( indx );
	m_pBtnSelectFMRI->SetLabel( pNewAnatomy->getName() );
    m_pBtnSelectFMRI->SetBackgroundColour(wxNullColour);
	
	//Set slider max value according to number of timelaps
	m_pSliderRest->SetMax((int)DatasetManager::getInstance()->m_pRestingStateNetwork->getBands()-1);
}

void FMRIWindow::OnSliderRestMoved( wxCommandEvent& WXUNUSED(event) )
{
	int sliderValue = m_pSliderRest->GetValue();
    m_pTxtRestBox->SetValue( wxString::Format( wxT( "%i"), sliderValue ) );
	DatasetManager::getInstance()->m_pRestingStateNetwork->SetTextureFromSlider( sliderValue );
}
