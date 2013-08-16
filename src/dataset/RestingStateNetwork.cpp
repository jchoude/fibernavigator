/////////////////////////////////////////////////////////////////////////////
// Name:            RestingState.cpp
// Author:          Maxime Chamberland
/////////////////////////////////////////////////////////////////////////////

#include "RestingStateNetwork.h"

#include "DatasetManager.h"
#include "../Logger.h"
#include "../gfx/ShaderHelper.h"
#include "../gui/MyListCtrl.h"
#include "../gui/SceneManager.h"
#include "../misc/nifti/nifti1_io.h"

#include <GL/glew.h>
#include <wx/math.h>
#include <wx/xml/xml.h>

#include <algorithm>
#include <fstream>
#include <limits>
#include <vector>
using std::vector;

#if defined(__WXMAC__) || defined(__WXMSW__)
#ifndef isnan
inline bool isnan(double x) {
    return x != x;
}
#endif
#endif

///////////////////////////////////////////
RestingStateNetwork::RestingStateNetwork():
m_dataType( 16 ),
m_bands( 108 )
{
	m_rows = DatasetManager::getInstance()->getRows();
	m_columns = DatasetManager::getInstance()->getColumns();
	m_frames =  DatasetManager::getInstance()->getFrames();

}

//////////////////////////////////////////////////////////////////////////
RestingStateNetwork::~RestingStateNetwork()
{
    Logger::getInstance()->print( wxT( "RestingStateNetwork destructor called but nothing to do." ), LOGLEVEL_DEBUG );
}

//////////////////////////////////////////////////////////////////////////
bool RestingStateNetwork::load( nifti_image *pHeader, nifti_image *pBody )
{
    int datasetSize = pHeader->dim[1] * pHeader->dim[2] * pHeader->dim[3];
	m_rows = pHeader->dim[1];
	m_columns = pHeader->dim[2];
	m_frames = pHeader->dim[3];
	m_bands = pHeader->dim[4];
    
    m_fileFloatData.assign( datasetSize * m_bands, 0.0f);
    float* pData = (float*)pBody->data;

    for( int i( 0 ); i < datasetSize; ++i )
    {
        for( int j( 0 ); j < m_bands; ++j )
        {
            if(!isnan(pData[j * datasetSize + i]))
                m_fileFloatData[i * m_bands + j] = pData[j * datasetSize + i];
        }
    }
    
    createStructure( m_fileFloatData );

	Logger::getInstance()->print( wxT( "Resting-state network initialized" ), LOGLEVEL_MESSAGE );
    return true;
}


//////////////////////////////////////////////////////////////////////////
bool RestingStateNetwork::createStructure  ( std::vector< float > &i_fileFloatData )
{
	int size = m_rows * m_columns * m_frames;
    m_signal.resize( size );

    vector< float >::iterator it;
    int i = 0;

    //Fetching the directions
    for( it = i_fileFloatData.begin(), i = 0; it != i_fileFloatData.end(); it += m_bands, ++i )
    { 
        m_signal[i].insert( m_signal[i].end(), it, it + m_bands );
    }

    return true;
}

