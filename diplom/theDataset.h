#ifndef THEDATASET_H_
#define THEDATASET_H_

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "wx/textfile.h"
#include "wx/file.h"
#include "wx/image.h"
#include "wx/wfstream.h"
#include "wx/datstrm.h"
#include "wx/txtstrm.h"

#include <vector>
#include "ArcBall.h"

#include "datasetInfo.h"
#include "selectionBox.h"

#include "mainFrame.h"
#include "theScene.h"
#include "point.h"

class TheDataset
{
private:

public:
	static bool load(int index = 0, wxString filename = wxT(""));
	static void finishLoading(DatasetInfo*);
	static bool loadSettings(wxString);
	static void save(wxString);

	static std::vector<std::vector<SelectionBox*> > getSelectionBoxes();
	static void updateAllSelectionBoxes();
	static void updateTreeDims();
	static void updateTreeDS(int);

	/*
	 * Called from MainFrame when a kdTree thread signals it's finished
	 */
	static void treeFinished();

	/*
	 * Helper functions
	 */
	static void printTime();
	static void printwxT(wxString);
	/*
	 * Check for GL error
	 */
	static bool GLError();
	static void printGLError(wxString function = wxT(""));


	static int rows;
	static int columns;
	static int frames;
	static bool useVBO;
	static unsigned int countFibers;
	static Matrix4fT m_transform;
	static wxString lastError;
	static GLenum lastGLError;
	static bool anatomy_loaded;
	static bool fibers_loaded;
	static bool surface_loaded;
	static bool surface_isDirty;
	static wxString lastPath;
	static int threadsActive;

	static MainFrame* mainFrame;
	static Point* m_lastSelectedPoint;
	static TheScene* m_scene;



};

#define ID_KDTREE_FINISHED	50

#endif /*THEDATASET_H_*/
