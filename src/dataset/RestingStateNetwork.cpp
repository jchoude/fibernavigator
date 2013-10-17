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

///////////////////////////////////////////
RestingStateNetwork::RestingStateNetwork():
m_isRealTimeOn( false ),
m_dataType( 16 ),
m_bands( 108 ),
m_corrThreshold( 3.0f ),
m_colorSliderValue( 5.0f )
{
	m_rows = DatasetManager::getInstance()->getRows();
	m_columns = DatasetManager::getInstance()->getColumns();
	m_frames =  DatasetManager::getInstance()->getFrames();

}

//////////////////////////////////////////////////////////////////////////
RestingStateNetwork::~RestingStateNetwork()
{
    Logger::getInstance()->print( wxT( "RestingStateNetwork destructor called but nothing to do." ), LOGLEVEL_DEBUG );
	//cudaFree(d_data);
}

//////////////////////////////////////////////////////////////////////////
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
    
	std::vector<short int> fileFloatData( m_datasetSize * m_bands, 0);
	//cuData = new short int[m_datasetSize*m_bands];

	if(pHeader->datatype == 4)
	{
		short int* pData = (short int*)pBody->data;
		
		//Prepare the data into a 1D vector, side by side
		for( int i( 0 ); i < m_datasetSize; ++i )
		{
			for( int j( 0 ); j < m_bands; ++j )
			{
				//if(!isnan(pData[j * datasetSize + i]))
					//cuData[i * m_bands + j] = pData[j * m_datasetSize + i];
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
				//if(!isnan(pData[j * datasetSize + i]))
					fileFloatData[i * m_bands + j] = pData[j * m_datasetSize + i];
			}
		}
		
	}
	//std::cout << "Before: " << cuData[2000];
	//cudaMalloc(&d_data, m_datasetSize * m_bands * sizeof(short int));
	//cudaMemcpy(d_data, cuData, m_datasetSize * m_bands * sizeof(short int), cudaMemcpyHostToDevice);
	//Assign structure to a 2D vector of timelaps
    createStructure( fileFloatData );

	Logger::getInstance()->print( wxT( "Resting-state network initialized" ), LOGLEVEL_MESSAGE );
    return true;
}


//////////////////////////////////////////////////////////////////////////
bool RestingStateNetwork::createStructure  ( std::vector< short int > &i_fileFloatData )
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
		for( int b(0); b < m_bands; ++b )
		{
			m_volumes[b].push_back(m_signalNormalized[s][b]);
		}
		calculateMeanAndSigma(m_signalNormalized[s], m_meansAndSigmas[s]);
    }

	//Create texture made of 1st timelaps
	data.assign(size, 0.0f);

    return true;
}

void RestingStateNetwork::SetTextureFromSlider(int sliderValue)
{
	Anatomy* pNewAnatomy = (Anatomy *)DatasetManager::getInstance()->getDataset( m_index );
	pNewAnatomy->setFloatDataset(m_volumes[sliderValue]);
	pNewAnatomy->generateTexture();
}

void RestingStateNetwork::SetTextureFromNetwork()
{
	Anatomy* pNewAnatomy = (Anatomy *)DatasetManager::getInstance()->getDataset( m_index );
	pNewAnatomy->setFloatDataset(data);
	pNewAnatomy->generateTexture();
}

void RestingStateNetwork::seedBased()
{
	m_3Dpoints.clear();
	std::vector<float> texture(m_datasetSize, 0.0f);
	 
    float xVoxel = DatasetManager::getInstance()->getVoxelX();
    float yVoxel = DatasetManager::getInstance()->getVoxelY();
    float zVoxel = DatasetManager::getInstance()->getVoxelZ();

	int columns = DatasetManager::getInstance()->getColumns();
    int rows    = DatasetManager::getInstance()->getRows();
	std::vector<float> positions; 

    Vector minCorner, maxCorner, middle;
    SelectionTree::SelectionObjectVector selObjs = SceneManager::getInstance()->getSelectionTree().getAllObjects();

	for( unsigned int b = 0; b < selObjs.size(); b++ )
	{
		minCorner.x = (int)(floor(selObjs[b]->getCenter().x - selObjs[b]->getSize().x * xVoxel /  2.0f ) / xVoxel );
		minCorner.y = (int)(floor(selObjs[b]->getCenter().y - selObjs[b]->getSize().y * yVoxel /  2.0f ) / yVoxel );
		minCorner.z = (int)(floor(selObjs[b]->getCenter().z - selObjs[b]->getSize().z * zVoxel /  2.0f ) / zVoxel );
		maxCorner.x = (int)(floor(selObjs[b]->getCenter().x + selObjs[b]->getSize().x * xVoxel /  2.0f ) / xVoxel );
		maxCorner.y = (int)(floor(selObjs[b]->getCenter().y + selObjs[b]->getSize().y * yVoxel /  2.0f ) / yVoxel );
		maxCorner.z = (int)(floor(selObjs[b]->getCenter().z + selObjs[b]->getSize().z * zVoxel /  2.0f ) / zVoxel );
		
		for( float x = minCorner.x; x <= maxCorner.x; x++)
		{
			for( float y = minCorner.y; y <= maxCorner.y; y++)
			{
				for( float z = minCorner.z; z <= maxCorner.z; z++)
				{
					positions.push_back( z * columns * rows + y *columns + x );
				}
			}
		}
		correlate(texture, positions);
	}
	
	render3D();
	Anatomy* pNewAnatomy = (Anatomy *)DatasetManager::getInstance()->getDataset( m_index );
	pNewAnatomy->setFloatDataset(texture);
	pNewAnatomy->generateTexture();
	RTFMRIHelper::getInstance()->setRTFMRIDirty(false);
}

void RestingStateNetwork::render3D()
{
	if( m_3Dpoints.size() > 0 )
    {
		for (unsigned int s = 0; s < m_3Dpoints.size(); s++)
		{
			glEnable(GL_POINT_SPRITE);
			glPointSize(m_3Dpoints[s].second*10.0f);
			glColor3f(m_3Dpoints[s].second/m_colorSliderValue,0,0);
			glBegin(GL_POINTS);
				glVertex3f(m_3Dpoints[s].first.x * m_voxelSizeX, m_3Dpoints[s].first.y * m_voxelSizeY, m_3Dpoints[s].first.z * m_voxelSizeZ);
			glEnd();
		}
	}
}
void RestingStateNetwork::correlate(std::vector<float>& texture, std::vector<float>& positions)
{
	 //float data[N]; int count = 0;
	//int N = m_rows * m_columns * m_bands;
	//cuData = new float[m_rows][m_columns];
	//for(int i =0; i<N;i++)
	//	cuData[i] = i+0.1f;

	//std::cout << "Before: " << cuData[1];
	//cudaMalloc(&d_data, N * sizeof(float));
	//cudaMemcpy(d_data, cuData, N * sizeof(float), cudaMemcpyHostToDevice);
	//int block_size = 512;
	//int n_blocks = N/block_size + (N%block_size == 0 ? 0:1);
	//cuCorrelation<<<n_blocks,block_size>>>(d_data);
 //   cudaMemcpy(cuData, d_data, N * sizeof(float), cudaMemcpyDeviceToHost);
 //    
	//std::cout <<"after " << cuData[1];

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
					value /= (m_bands - 1);
				
					corrSum+=value;
					corrFactors[i] = value;
					nb++;
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
				if(corrFactors[i] != 0)
				{
					sigma += sqrt((corrFactors[i] - meanCorr)*(corrFactors[i] - meanCorr));	
				}		
			}
		}
	}

	//Calculate z-scores, and write into texture.
	sigma /= nb;
	for( float x = 0; x < m_columns; x++)
	{
		for( float y = 0; y < m_rows; y++)
		{
			for( float z = 0; z < m_frames; z++)
			{
				int i = z * m_columns * m_rows + y *m_columns + x;
				if(corrFactors[i] != 0)
				{
					float zScore = (corrFactors[i] - meanCorr) / sigma;
					if(zScore > m_corrThreshold)
					{
						texture[i] = zScore/m_colorSliderValue;
						m_3Dpoints.push_back(std::pair<Vector,float>(Vector(x,y,z),zScore));
					}
				}
			}
		}
	}
}

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