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
	std::vector<std::vector<float> >* getSignal() { return &m_signal; }

	std::vector<float> data;

private:
    
    bool createStructure  ( std::vector< float > &i_fileFloatData );
  
    std::vector<std::vector<float> >   m_signal;
    std::vector< float > m_fileFloatData;
    int m_dataType;
	int m_rows;
	int m_columns;
	int m_frames;
	int m_bands;
	DatasetIndex m_index;

};

#endif /* RESTINGSTATENETWORK_H_ */
