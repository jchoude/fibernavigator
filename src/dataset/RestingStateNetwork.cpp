/////////////////////////////////////////////////////////////////////////////
// Name:            RestingState.cpp
// Author:          Maxime Chamberland
/////////////////////////////////////////////////////////////////////////////

#include "RestingStateNetwork.h"

#include "DatasetManager.h"
#include "AnatomyHelper.h"
#include "RTFMRIHelper.h"
#include "../Logger.h"
#include "../gfx/ShaderHelper.h"
#include "../gfx/TheScene.h"
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

#include "../gfx/Image.h"
#include "../gfx/BitmapHandling.h"
#include "../gfx/TextureHandling.h"
#include "../main.h"

//#include <cuda.h>
//#include <cuda_runtime.h>
//#include <device_launch_parameters.h>

//#define N (67445)
//#define M (100)

//__global__ void cuCorrelation(float *buf)
//{
//	int i = threadIdx.x + blockIdx.x * blockDim.x;
//	if(i < N)
//		buf[i] = buf[i]+0.2f;
//}

//////////////////////////////////////////
//Constructor
//////////////////////////////////////////
RestingStateNetwork::RestingStateNetwork():
m_zMin( 999.0f ),
m_zMax( 0.0f ),
m_alpha( 0.5f),
m_pointSize( 5.0f ),
m_isRealTimeOn( false ),
m_dataType( 16 ),
m_bands( 108 ),
m_corrThreshold( 1.96f ),
m_colorSliderValue( 5.0f ),
m_normalize( true ),
m_boxMoving( false ),
m_originL(0,0,0),
m_origin(0,0,0)
{
	m_rowsL = DatasetManager::getInstance()->getRows();
	m_columnsL = DatasetManager::getInstance()->getColumns();
	m_framesL =  DatasetManager::getInstance()->getFrames();

	m_xL = DatasetManager::getInstance()->getVoxelX();
	m_yL = DatasetManager::getInstance()->getVoxelY();
	m_zL =  DatasetManager::getInstance()->getVoxelZ();

	m_datasetSizeL = m_rowsL * m_columnsL * m_framesL;

	FMatrix &t = DatasetManager::getInstance()->getNiftiTransform();
	m_originL.x = floor(abs(t(0,3)) / m_xL);
	m_originL.y = floor(abs(t(1,3)) / m_yL);
	m_originL.z = floor(abs(t(2,3)) / m_zL);

}

//////////////////////////////////////////
//Destructor
//////////////////////////////////////////
RestingStateNetwork::~RestingStateNetwork()
{
    Logger::getInstance()->print( wxT( "RestingStateNetwork destructor called but nothing to do." ), LOGLEVEL_DEBUG );
}

//////////////////////////////////////////
//Load
//////////////////////////////////////////
bool RestingStateNetwork::load( nifti_image *pHeader, nifti_image *pBody )
{
    m_datasetSize = pHeader->dim[1] * pHeader->dim[2] * pHeader->dim[3];
	m_rows = pHeader->dim[1];
	m_columns = pHeader->dim[2];
	m_frames = pHeader->dim[3];
	m_bands = pHeader->dim[4];

    m_voxelSizeX = pHeader->dx;
    m_voxelSizeY = pHeader->dy;
    m_voxelSizeZ = pHeader->dz;

	if( pHeader->sform_code > 0 )
    {
		m_origin.x = floor(pHeader->sto_ijk.m[0][3]);
		m_origin.y = floor(pHeader->sto_ijk.m[1][3]);
		m_origin.z = floor(pHeader->sto_ijk.m[2][3]);
    }
    else if( pHeader->qform_code > 0 )
    {
        m_origin.x = floor(pHeader->qto_ijk.m[0][3]);
		m_origin.y = floor(pHeader->qto_ijk.m[1][3]);
        m_origin.z = floor(pHeader->qto_ijk.m[2][3]);
    }
    else
    {
        Logger::getInstance()->print( wxT( "No transformation encoded in the nifti file. Alignement will fail." ), LOGLEVEL_WARNING );
    }
 
	std::vector<short int> fileFloatData( m_datasetSize * m_bands, 0);

	if(pHeader->datatype == 4)
	{
		short int* pData = (short int*)pBody->data;
		//Prepare the data into a 1D vector, side by side
		for( int i( 0 ); i < m_datasetSize; ++i )
		{
			for( int j( 0 ); j < m_bands; ++j )
			{
				fileFloatData[i * m_bands + j] = pData[j * m_datasetSize + i];
			}
		}
	}
	else
	{
		float* pData = (float*)pBody->data;
		//Prepare the data into a 1D vector, side by side
		for( int i( 0 ); i < m_datasetSize; ++i )
		{
			for( int j( 0 ); j < m_bands; ++j )
			{
				fileFloatData[i * m_bands + j] = pData[j * m_datasetSize + i];
			}
		}
		
	}

	//Assign structure to a 2D vector of timelaps
    createStructure( fileFloatData );

	//Load fMRI sprite texture.
	Image<ColorRGB> TmpImage;
	wxString name = wxT ("fMRI.bmp");

	wxString iconPath = MyApp::iconsPath;
	wxString fullname = iconPath.append(name);
	std::string stlstring = std::string(fullname.mb_str());
	
    //Load the color scheme #1 image and send it to the GPU as a texture.
    LoadBmp(stlstring,TmpImage);
    m_lookupTex = LoadTexture(TmpImage);

	//Logger::getInstance()->print( wxT( "Resting-state network initialized" ), LOGLEVEL_MESSAGE );
    return true;
}


//////////////////////////////////////////
//Create structure
//////////////////////////////////////////
bool RestingStateNetwork::createStructure( std::vector< short int > &i_fileFloatData )
{
	int size = m_rows * m_columns * m_frames;
    m_signal.resize( size );
	m_signalNormalized.resize ( size );
    vector< short int >::iterator it;
    int i = 0;

    //Fetching the directions
    for( it = i_fileFloatData.begin(), i = 0; it != i_fileFloatData.end(); it += m_bands, ++i )
    { 
		m_signal[i].insert( m_signal[i].end(), it, it + m_bands );
    }
	
	//Find min/max for normalization
	vector<float> dataMax, dataMin;
	dataMax.assign(size, -std::numeric_limits<float>::infinity());
	dataMin.assign(size, std::numeric_limits<float>::infinity());
    for( int s(0); s < size; ++s )
    {
		for( int b(0); b < m_bands; ++b )
		{
			if (m_signal[s][b] > dataMax[s])
			{
				dataMax[s] = m_signal[s][b];
			}
			if (m_signal[s][b] < dataMin[s])
			{
				dataMin[s] = m_signal[s][b];
			}
		}
    }

	//Min max Rescale
    for( int s(0); s < size; ++s )
    {
		for( int b(0); b < m_bands; ++b )
		{
			if((m_signal[s][b] == 0 && dataMin[s] == 0) || (m_signal[s][b] == 16767 && dataMin[s] == 16767)) //Ensure that we dont divide by 0.
				m_signalNormalized[s].push_back(0);
			else
				m_signalNormalized[s].push_back ((m_signal[s][b] - dataMin[s]) / (dataMax[s] - dataMin[s]));
		}
    }

	m_volumes.resize(m_bands);
	m_meansAndSigmas.resize(size);
	//Transpose signal for easy acces of timelaps
    for( int s(0); s < size; ++s )
    {
		calculateMeanAndSigma(m_signalNormalized[s], m_meansAndSigmas[s]);
    }

	//Create texture made of 1st timelaps
	data.assign(size, 0.0f);

    return true;
}

//////////////////////////////////////////////////////////////////////////
//Get 3D indexes to fill the 1x1x1 texture from a 3x3x3 x,y,z point 
//////////////////////////////////////////////////////////////////////////
vector<int> RestingStateNetwork::get3DIndexes(int x, int y, int z)
{
	std::vector<int> indexes;

	for( int padx = 0; padx < 4; padx++)
	{
		for( int pady = 0; pady < 4; pady++)
		{
			for( int padz = 0; padz < 4; padz++)
			{
				int i = (z * floor(float(m_framesL/m_frames)) + padz) * m_columnsL * m_rowsL + (y *floor(float(m_rowsL/m_rows)) + pady) *m_columnsL + (x * floor(float(m_columnsL/m_columns)) + padx);
				indexes.push_back( i );
			}
		}
	}

	return indexes;
}

//////////////////////////////////////////
//Set raw data texture from sliderValue
//////////////////////////////////////////
void RestingStateNetwork::SetTextureFromSlider(int sliderValue)
{
	std::vector<float> vol(m_datasetSizeL* 3, 0.0f);
	std::vector<int> indexes;
	for( float x = 0; x < m_columns; x++)
	{
		for( float y = 0; y < m_rows; y++)
		{
			for( float z = 0; z < m_frames; z++)
			{
				//from fspace to t1space
				int i = z * m_columns * m_rows + y *m_columns + x;
				int zz = ((z - m_origin.z) * m_voxelSizeZ / m_zL) + m_originL.z;
				int yy = ((y - m_origin.y) * m_voxelSizeY / m_yL) + m_originL.y;
				int xx = ((x - m_origin.x) * m_voxelSizeX / m_xL) + m_originL.x;
				int s = zz * m_columnsL * m_rowsL + yy * m_columnsL + xx ; // O
				
				vol[s*3] = m_signalNormalized[i][sliderValue];
				vol[s*3 + 1] = m_signalNormalized[i][sliderValue];
				vol[s*3 + 2] = m_signalNormalized[i][sliderValue];

				//Patch arround for 1x1x1
				//if(m_framesL != m_frames)
				//{
				//	indexes = get3DIndexes(x,y,z);
				//	for(unsigned int s = 0; s < indexes.size(); s++)
				//	{
				//		int id = indexes[s];
				//		vol[id*3] = m_signalNormalized[i][sliderValue];
				//		vol[id*3 + 1] = m_signalNormalized[i][sliderValue];
				//		vol[id*3 + 2] = m_signalNormalized[i][sliderValue];
				//	}
				//}
			}
		}
	}

	Anatomy* pNewAnatomy = (Anatomy *)DatasetManager::getInstance()->getDataset( m_index );
	pNewAnatomy->setFloatDataset(vol);
	pNewAnatomy->generateTexture();
}

//////////////////////////////////////////////////////////////////////////////////////////
//Set texture from Network fmri: NOTE: doesnt work functionally yet, data should be set
//////////////////////////////////////////////////////////////////////////////////////////
void RestingStateNetwork::SetTextureFromNetwork()
{
	Anatomy* pNewAnatomy = (Anatomy *)DatasetManager::getInstance()->getDataset( m_index );
	pNewAnatomy->setFloatDataset(data);
	pNewAnatomy->generateTexture();
}

//////////////////////////////////////////////////////////////////////////////////////////
//Initiate the seed-based algorithm
//////////////////////////////////////////////////////////////////////////////////////////
void RestingStateNetwork::seedBased()
{
	m_3Dpoints.clear();
	m_zMin = 999.0f;
	m_zMax = 0.0f;
	m_boxMoving = true;
	 
	std::vector<float> positions; 

    Vector minCorner, maxCorner, middle;
    SelectionTree::SelectionObjectVector selObjs = SceneManager::getInstance()->getSelectionTree().getAllObjects();

	for( unsigned int b = 0; b < selObjs.size(); b++ )
	{
		minCorner.x = (int)(floor(selObjs[b]->getCenter().x - selObjs[b]->getSize().x * m_xL /  2.0f ) / m_xL );
		minCorner.y = (int)(floor(selObjs[b]->getCenter().y - selObjs[b]->getSize().y * m_yL /  2.0f ) / m_yL );
		minCorner.z = (int)(floor(selObjs[b]->getCenter().z - selObjs[b]->getSize().z * m_zL /  2.0f ) / m_zL );
		maxCorner.x = (int)(floor(selObjs[b]->getCenter().x + selObjs[b]->getSize().x * m_xL /  2.0f ) / m_xL );
		maxCorner.y = (int)(floor(selObjs[b]->getCenter().y + selObjs[b]->getSize().y * m_yL /  2.0f ) / m_yL );
		maxCorner.z = (int)(floor(selObjs[b]->getCenter().z + selObjs[b]->getSize().z * m_zL /  2.0f ) / m_zL );
		
		for( float x = minCorner.x; x <= maxCorner.x; x++)
		{
			for( float y = minCorner.y; y <= maxCorner.y; y++)
			{
				for( float z = minCorner.z; z <= maxCorner.z; z++)
				{
					//Switch to 3x3x3 from t1space
					int zz = ((z - m_originL.z) * m_zL / m_voxelSizeZ) + m_origin.z;
					int yy = ((y - m_originL.y) * m_yL/ m_voxelSizeY) + m_origin.y;
					int xx = ((x - m_originL.x) * m_xL /m_voxelSizeX) + m_origin.x;
					int i = zz * m_columns * m_rows + yy * m_columns + xx ; // O
					positions.push_back( i );
				}
			}
		}
		correlate(positions);
	}
	
	//TODO can be done in rendering directly while looping, change from fspace to t1space
    for(unsigned int s(0); s < m_3Dpoints.size(); ++s )
    {
		m_3Dpoints[s].first.x = ((m_3Dpoints[s].first.x - m_origin.x) * m_voxelSizeX / m_xL) + m_originL.x;
		m_3Dpoints[s].first.y = ((m_3Dpoints[s].first.y - m_origin.y) * m_voxelSizeY / m_yL) + m_originL.y;
		m_3Dpoints[s].first.z = ((m_3Dpoints[s].first.z - m_origin.z) * m_voxelSizeZ / m_zL) + m_originL.z;
    }

	render3D(false);
	RTFMRIHelper::getInstance()->setRTFMRIDirty(false);
}

//////////////////////////////////////////////////////////////////////////////////////////
//Rendering function, for both 3D sprites and textures options.
//////////////////////////////////////////////////////////////////////////////////////////
void RestingStateNetwork::render3D(bool recalculateTexture)
{
	if( m_3Dpoints.size() > 0 )
    {
		std::vector<float> texture(m_datasetSizeL*3, 0.0f);
		
		//Apply ColorMap
		for (unsigned int s = 0; s < m_3Dpoints.size(); s++)
		{
			float R,G,B;
            bool render = true;

			float mid = (m_zMin + m_zMax) / 2.0f;
			float quart = 1.0f* (m_zMin + m_zMax) / 4.0f;
			float trois_quart = 3.0f* (m_zMin + m_zMax) / 4.0f;
			float v = (m_3Dpoints[s].second - m_zMin) / (m_zMax - m_zMin);

			if(m_3Dpoints[s].second < quart)
			{
                R = (m_3Dpoints[s].second - m_zMin) / (quart - m_zMin);
                G = 0.0f;
                B = 0.0f;
                render = false;
			}
			else if(m_3Dpoints[s].second >= quart && m_3Dpoints[s].second < trois_quart)
			{
				R = 1.0f;
				G = (m_3Dpoints[s].second - quart) / (trois_quart - quart);
				B = 0.0f;
			}
			else
			{
				R = 1.0f;
				G = 1.0f;
				B = v;
			}

			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_POINT_SPRITE);
			glPointSize(m_3Dpoints[s].second * m_pointSize + 1.0f);
			glColor4f(R,G,B,(m_3Dpoints[s].second / m_zMax) * m_alpha + 0.1f);

			//glActiveTexture(GL_TEXTURE0);
			//glEnable( GL_TEXTURE_2D );
			//glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
			//glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_BLEND);
			//glBindTexture(GL_TEXTURE_2D, m_lookupTex);

            if(render)
            {
			    glBegin(GL_POINTS);
				    glVertex3f(m_3Dpoints[s].first.x, m_3Dpoints[s].first.y, m_3Dpoints[s].first.z);
			    glEnd();
                render = true;
            }

			//glDisable( GL_TEXTURE_2D );
			glDisable(GL_POINT_SPRITE);
			glDisable(GL_BLEND);

			//check for Unnecesssary computation when static scene
			if(recalculateTexture)
			{
				std::vector<int> indexes;
				//Must be in t1space
				int i = m_3Dpoints[s].first.z * m_columnsL * m_rowsL + m_3Dpoints[s].first.y * m_columnsL + m_3Dpoints[s].first.x ; // O
				
				texture[i*3] = R;
				texture[i*3 + 1] = G;
				texture[i*3 + 2] = B;

				//Patch arround for 1x1x1
				//if(m_framesL != m_frames)
				//{
				//	indexes = get3DIndexes(floor(float(m_3Dpoints[s].first.x * m_columns/m_columnsL)),floor(float(m_3Dpoints[s].first.y * m_rows/m_rowsL)),floor(float(m_3Dpoints[s].first.z * m_frames/m_framesL)));
				//	for(unsigned int u = 0; u < indexes.size(); u++)
				//	{
				//		int id = indexes[u];
				//		texture[id*3] = R;
				//		texture[id*3 + 1] = G;
				//		texture[id*3 + 2] = B;
				//	}
				//}
			}
		}
		//TEXTURE
		if(recalculateTexture)
		{
			Anatomy* pNewAnatomy = (Anatomy *)DatasetManager::getInstance()->getDataset( m_index );
			pNewAnatomy->setFloatDataset(texture);
			pNewAnatomy->generateTexture();
		}
	}

}

//////////////////////////////////////////////////////////////////////////////////////////
//Correlation function given a position, with all other time series
//////////////////////////////////////////////////////////////////////////////////////////
void RestingStateNetwork::correlate(std::vector<float>& positions)
{
	//Mean signal inside box
	std::vector<float> meanSignal;
	for(int i=0; i < m_bands; i++)
	{
		float sum = 0;
		for(unsigned int j=0; j < positions.size(); j++)
		{	
			int idx = positions[j];
			sum += m_signalNormalized[idx][i];
		}
		sum /= positions.size();
		meanSignal.push_back( sum );
	}

	//Get mean and sigma of it
	std::pair<float, float> RefMeanAndSigma;
	calculateMeanAndSigma(meanSignal, RefMeanAndSigma);
	std::vector<float> corrFactors;
	corrFactors.assign(m_datasetSize, 0.0f);
	float corrSum = 0.0f;
	int nb = 0;

	//Correlate with rest of the brain, i.e find corr factors
	for( float x = 0; x < m_columns; x++)
	{
		for( float y = 0; y < m_rows; y++)
		{
			for( float z = 0; z < m_frames; z++)
			{
				int i = z * m_columns * m_rows + y *m_columns + x;
				if(m_meansAndSigmas[i].first != 0)
				{
					float num = 0.0f;
					float denum = 0.0f;
					
					for(int j = 0; j < m_bands; j++)
					{
						num += (meanSignal[j] - RefMeanAndSigma.first) * ( m_signalNormalized[i][j] - m_meansAndSigmas[i].first);
					}
					float value = num / ( RefMeanAndSigma.second * m_meansAndSigmas[i].second);
					value /= (m_bands);
				
					if(value > 0)
					{
						corrSum+=value;
						corrFactors[i] = value;
						nb++;
					}
					else
						corrFactors[i] = -1;
				}
				else
					corrFactors[i] = 0.0f;
			}
		}
	}

	//Find mean and sigma of all corr factors.
	float meanCorr = corrSum / nb;
	float sigma = 0.0f;
	for( float x = 0; x < m_columns; x++)
	{
		for( float y = 0; y < m_rows; y++)
		{
			for( float z = 0; z < m_frames; z++)
			{
				int i = z * m_columns * m_rows + y *m_columns + x;
				if(corrFactors[i] > 0.0f)
				{
					sigma += (corrFactors[i] - meanCorr)*(corrFactors[i] - meanCorr);	
				}		
			}
		}
	}

	//Calculate z-scores, and save them.
	sigma /= nb;
	sigma = sqrt(sigma);
	for( float x = 0; x < m_columns; x++)
	{
		for( float y = 0; y < m_rows; y++)
		{
			for( float z = 0; z < m_frames; z++)
			{
				int i = z * m_columns * m_rows + y *m_columns + x;
				
				if(m_corrThreshold == 0.0f && corrFactors[i] != 0)
					{	
						m_3Dpoints.push_back(std::pair<Vector,float>(Vector(x,y,z),0.0f));
					}

				if(corrFactors[i] > 0)
				{
					float zScore = (corrFactors[i] - meanCorr) / sigma;
					if(zScore < m_zMin && zScore > 0.0f)
						m_zMin = zScore;
					if(zScore > m_zMax)
						m_zMax = zScore;
					if(zScore > m_corrThreshold)
					{
						m_3Dpoints.push_back(std::pair<Vector,float>(Vector(x,y,z),zScore));
					}
				
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//Calculate Mean and Sigma for the signal inside the box
//////////////////////////////////////////////////////////////////////////////////////////
void RestingStateNetwork::calculateMeanAndSigma(std::vector<float> signal, std::pair<float, float>& params)
{
	float mean = 0.0f;
	float sigma = 0.0f;
	
	//mean
	for(unsigned int i=0; i < signal.size(); i++)
	{
		mean+=signal[i];
	}
	mean /= signal.size();

	//sigma
    for(unsigned int i = 0; i < signal.size(); i++)
    {
         sigma += (signal[i] - mean) * (signal[i] - mean) ;
    }
    sigma /= signal.size();

	params.first = mean;
	params.second = sqrt(sigma);
}