/////////////////////////////////////////////////////////////////////////////
// Name:            RestingStateNetwork.h
// Author:          Maxime Chamberland
/////////////////////////////////////////////////////////////////////////////
#ifndef RESTINGSTATENETWORK_H_
#define RESTINGSTATENETWORK_H_

#include "Anatomy.h"
#include "DatasetInfo.h"
#include "DatasetIndex.h"
#include "../misc/nifti/nifti1_io.h"

class RestingStateNetwork
{
public:

    // Constructor/Destructor
    RestingStateNetwork();
    virtual ~RestingStateNetwork();

	bool load( nifti_image *pHeader, nifti_image *pBody );
	void setNetworkInfo( DatasetIndex index ) { m_index = index; }
	void SetTextureFromSlider( int sliderValue );
	void SetTextureFromNetwork();
	
	void seedBased();

	std::vector<std::vector<float> >* getSignal() { return &m_signal; }
	DatasetIndex getIndex()   { return m_index; }
	DatasetIndex getColumns() { return m_columns; }
	DatasetIndex getRows()    { return m_rows; }
	DatasetIndex getFrames()  { return m_frames; }
	DatasetIndex getBands()   { return m_bands; }
	
	std::vector<float> data; //Used for texture mapping

private:
    
    bool createStructure  ( std::vector< float > &i_fileFloatData );

    std::vector<std::vector<float> >   m_signal; //2D containing the data
	std::vector<std::vector<float> >   m_signalNormalized; //2D containing the data normalized
	std::vector<std::vector<float> >   m_volumes; //2D containing the data normalized volume-wise aligned

	bool m_isRealTimeOn;
    int m_dataType;
	int m_rows;
	int m_columns;
	int m_frames;
	int m_bands;
	DatasetIndex m_index;

};

#endif /* RESTINGSTATENETWORK_H_ */
