#ifndef MYLISTCTRL_H_
#define MYLISTCTRL_H_

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "wx/treectrl.h"


class MainFrame;

class MyTreeCtrl: public wxTreeCtrl 
{
public:
    MyTreeCtrl(MainFrame *pMainFrame, wxWindow *pParent, const wxWindowID id, const wxPoint& pos, const wxSize& size, long style );
    int getSelectedType();

private:
    MainFrame *m_mainFrame;
    void OnChar(wxKeyEvent& event);
    void OnRightClick(wxMouseEvent& event);
    void OnToggleAndNot(wxCommandEvent& event);
    void OnDeleteBox(wxCommandEvent& event);


    DECLARE_EVENT_TABLE()
};

class MySlider: public wxSlider 
{
public:
    MySlider( wxWindow *parent, const wxWindowID id, int value , int minValue, int maxValue, const wxPoint& pos, const wxSize& size, long style ) :
        wxSlider (parent, id, value, minValue, maxValue, pos, size, style)
    {
        m_leftIsDown = false;
    };

    bool leftDown() {return m_leftIsDown;};
private:
    void OnMouseEvent(wxMouseEvent& event);

    bool m_leftIsDown;

    DECLARE_EVENT_TABLE()
};

#endif /*MYLISTCTRL_H_*/
