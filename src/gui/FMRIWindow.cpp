#include "FMRIWindow.h"

#include "MainFrame.h"
#include "SceneManager.h"
#include "SelectionBox.h"
#include "SelectionEllipsoid.h"
#include "../main.h"
#include "../dataset/Anatomy.h"
#include "../dataset/Fibers.h"
#include "../dataset/RestingStateNetwork.h"
#include "../dataset/RTFMRIHelper.h"
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
    m_pMainFrame( pMf ),
	showRawData( true )
{
    SetBackgroundColour( *wxLIGHT_GREY );
    m_pFMRISizer = new wxBoxSizer( wxVERTICAL );
    SetSizer( m_pFMRISizer );
    SetAutoLayout( true );

    m_pBtnSelectFMRI = new wxButton( this, wxID_ANY,wxT("Load resting-state"), wxPoint(30,0), wxSize(230, -1) );
	pMf->Connect( m_pBtnSelectFMRI->GetId(), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MainFrame::onLoadAsRestingState) );
    m_pBtnSelectFMRI->SetBackgroundColour(wxColour( 255, 147, 147 ));

	wxBoxSizer *pBoxRow1 = new wxBoxSizer( wxHORIZONTAL );
	pBoxRow1->Add( m_pBtnSelectFMRI, 0, wxALIGN_CENTER_HORIZONTAL | wxALL, 1 );
	m_pFMRISizer->Add( pBoxRow1, 0, wxFIXED_MINSIZE | wxALL, 2 );

	m_pTextDisplayMode = new wxStaticText( this, wxID_ANY, wxT( "Display:" ), wxPoint(0,30), wxSize(200, -1) );
    m_pRadShowRawData = new wxRadioButton( this,  wxID_ANY, wxT( "Raw Data" ), wxPoint(30,60), wxSize(160, -1) );
	m_pRadShowRawData->Disable();
	m_pRadShowNetwork = new wxRadioButton( this,  wxID_ANY, wxT( "Network" ), wxPoint(30,90), wxSize(160, -1) );
	m_pRadShowNetwork->Disable();
	Connect( m_pRadShowRawData->GetId(), wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler( FMRIWindow::onSwitchViewRaw ) );
	Connect( m_pRadShowNetwork->GetId(), wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler( FMRIWindow::onSwitchViewNet ) );

	wxBoxSizer *pBoxRow2 = new wxBoxSizer( wxVERTICAL );
	pBoxRow2->Add( m_pTextDisplayMode, 0, wxALIGN_CENTER_VERTICAL | wxALL, 1 );
	pBoxRow2->Add( m_pRadShowRawData, 0, wxALIGN_CENTER, 1 );
	pBoxRow2->Add( m_pRadShowNetwork, 0, wxALIGN_CENTER, 1 );
	m_pFMRISizer->Add( pBoxRow2, 0, wxFIXED_MINSIZE | wxALL, 2 );

	m_pTextVolumeId = new wxStaticText( this, wxID_ANY, wxT("Volume"), wxPoint(0,120), wxSize(60, -1), wxALIGN_CENTER );
	m_pSliderRest = new MySlider( this, wxID_ANY, 0, 0, 107, wxPoint(60,120), wxSize(130, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
	m_pSliderRest->Disable();
	Connect( m_pSliderRest->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(FMRIWindow::OnSliderRestMoved) );
    m_pTxtRestBox = new wxTextCtrl( this, wxID_ANY, wxT("0"), wxPoint(190,120), wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

	wxBoxSizer *pBoxRow3 = new wxBoxSizer( wxHORIZONTAL );
    pBoxRow3->Add( m_pTextVolumeId, 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 1 );
    pBoxRow3->Add( m_pSliderRest,   0, wxALIGN_LEFT | wxEXPAND | wxALL, 1);
	pBoxRow3->Add( m_pTxtRestBox,   0, wxALIGN_LEFT | wxALL, 1);
	m_pFMRISizer->Add( pBoxRow3, 0, wxFIXED_MINSIZE | wxEXPAND, 0 );

	m_pBtnStart = new wxToggleButton( this, wxID_ANY,wxT("Start correlation"), wxPoint(0,150), wxSize(230, 50) );
    Connect( m_pBtnStart->GetId(), wxEVT_COMMAND_TOGGLEBUTTON_CLICKED, wxCommandEventHandler(FMRIWindow::OnStartRTFMRI) );
    m_pBtnStart->Enable(false);

	wxBoxSizer *pBoxRow4 = new wxBoxSizer( wxHORIZONTAL );
	pBoxRow4->Add( m_pBtnStart, 0, wxALIGN_CENTER_HORIZONTAL | wxALL, 1 );
	m_pFMRISizer->Add( pBoxRow4, 0, wxFIXED_MINSIZE | wxALL, 2 );

	m_pTextCorrThreshold = new wxStaticText( this, wxID_ANY, wxT("Z-Threshold"), wxPoint(0,180), wxSize(60, -1), wxALIGN_CENTER );
	m_pSliderCorrThreshold = new MySlider( this, wxID_ANY, 0, 1, 150, wxPoint(60,180), wxSize(130, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
	m_pSliderCorrThreshold->SetValue( 30 );
	Connect( m_pSliderCorrThreshold->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(FMRIWindow::OnSliderCorrThreshMoved) );
    m_pTxtCorrThreshBox = new wxTextCtrl( this, wxID_ANY, wxT("3.0"), wxPoint(190,180), wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

	wxBoxSizer *pBoxRow5 = new wxBoxSizer( wxHORIZONTAL );
    pBoxRow5->Add( m_pTextCorrThreshold, 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 1 );
    pBoxRow5->Add( m_pSliderCorrThreshold,   0, wxALIGN_LEFT | wxEXPAND | wxALL, 1);
	pBoxRow5->Add( m_pTxtCorrThreshBox,   0, wxALIGN_LEFT | wxALL, 1);
	m_pFMRISizer->Add( pBoxRow5, 0, wxFIXED_MINSIZE | wxEXPAND, 0 );

	m_pTextColorMap = new wxStaticText( this, wxID_ANY, wxT("Max Z-Color range"), wxPoint(0,210), wxSize(60, -1), wxALIGN_CENTER );
	m_pSliderColorMap = new MySlider( this, wxID_ANY, 0, 1, 100, wxPoint(60,210), wxSize(130, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
	m_pSliderColorMap->SetValue( 50 );
	Connect( m_pSliderColorMap->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(FMRIWindow::OnSliderColorMoved) );
    m_pTxtColorMapBox = new wxTextCtrl( this, wxID_ANY, wxT("5.0"), wxPoint(190,210), wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

	wxBoxSizer *pBoxRow6 = new wxBoxSizer( wxHORIZONTAL );
    pBoxRow6->Add( m_pTextColorMap, 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 1 );
    pBoxRow6->Add( m_pSliderColorMap,   0, wxALIGN_LEFT | wxEXPAND | wxALL, 1);
	pBoxRow6->Add( m_pTxtColorMapBox,   0, wxALIGN_LEFT | wxALL, 1);
	m_pFMRISizer->Add( pBoxRow6, 0, wxFIXED_MINSIZE | wxEXPAND, 0 );
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
	
	//m_pSliderRest->Enable();
	//Set slider max value according to number of timelaps
	m_pSliderRest->SetMax((int)DatasetManager::getInstance()->m_pRestingStateNetwork->getBands()-1);

	m_pRadShowRawData->Enable();
	m_pRadShowNetwork->Enable();
	m_pRadShowNetwork->SetValue(true);	
}

void FMRIWindow::onSwitchViewRaw( wxCommandEvent& WXUNUSED(event) )
{
	showRawData = true;
	m_pSliderRest->Enable();
	m_pTextVolumeId->Enable();
	m_pTxtRestBox->Enable();

	int sliderValue = m_pSliderRest->GetValue();
    m_pTxtRestBox->SetValue( wxString::Format( wxT( "%i"), sliderValue ) );
	DatasetManager::getInstance()->m_pRestingStateNetwork->SetTextureFromSlider( sliderValue );
}

void FMRIWindow::onSwitchViewNet( wxCommandEvent& WXUNUSED(event) )
{
	showRawData = false;
	m_pSliderRest->Disable();
	m_pTextVolumeId->Disable();
	m_pTxtRestBox->Disable();
	DatasetManager::getInstance()->m_pRestingStateNetwork->SetTextureFromNetwork();
}
void FMRIWindow::OnSliderRestMoved( wxCommandEvent& WXUNUSED(event) )
{
	int sliderValue = m_pSliderRest->GetValue();
    m_pTxtRestBox->SetValue( wxString::Format( wxT( "%i"), sliderValue ) );
	DatasetManager::getInstance()->m_pRestingStateNetwork->SetTextureFromSlider( sliderValue );
}

void FMRIWindow::OnSliderCorrThreshMoved(  wxCommandEvent& WXUNUSED(event) )
{
	float sliderValue = m_pSliderCorrThreshold->GetValue() / 10.0f;
    m_pTxtCorrThreshBox->SetValue( wxString::Format( wxT( "%.2f"), sliderValue ) );
	DatasetManager::getInstance()->m_pRestingStateNetwork->SetCorrThreshold( sliderValue );
	RTFMRIHelper::getInstance()->setRTFMRIDirty( true );
}

void FMRIWindow::OnSliderColorMoved(  wxCommandEvent& WXUNUSED(event) )
{
	float sliderValue = m_pSliderColorMap->GetValue() / 10.0f;
	m_pTxtColorMapBox->SetValue( wxString::Format( wxT( "%.2f"), sliderValue ) );
	DatasetManager::getInstance()->m_pRestingStateNetwork->SetColorSliderValue( sliderValue );
	RTFMRIHelper::getInstance()->setRTFMRIDirty( true );
}

void FMRIWindow::OnStartRTFMRI( wxCommandEvent& WXUNUSED(event) )
{
	RTFMRIHelper::getInstance()->toggleRTFMRIReady();
    RTFMRIHelper::getInstance()->setRTFMRIDirty( true );

    if( !RTFMRIHelper::getInstance()->isRTFMRIReady() )
    {
        //m_pMainFrame->m_pMainGL->m_pRealTimeFibers->clearFibersRTT();
        //m_pMainFrame->m_pMainGL->m_pRealTimeFibers->clearColorsRTT();
        RTFMRIHelper::getInstance()->setRTFMRIDirty( false );
        m_pBtnStart->SetLabel(wxT("Start correlation"));
    }
    else
    {
        m_pBtnStart->SetLabel(wxT("Stop correlation"));
	}
}

