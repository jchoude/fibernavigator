/////////////////////////////////////////////////////////////////////////////
// Name:            RestingStateNetwork.h
// Author:          Maxime Chamberland
/////////////////////////////////////////////////////////////////////////////
#ifndef RESTINGSTATENETWORK_H_
#define RESTINGSTATENETWORK_H_

#include "DatasetInfo.h"
#include "../misc/nifti/nifti1_io.h"

class RestingStateNetwork
{
public:
    // Constructor/Destructor
    RestingStateNetwork();
    virtual ~RestingStateNetwork();
	bool load( nifti_image *pHeader, nifti_image *pBody );

private:
    
    bool createStructure  ( std::vector< float > &i_fileFloatData );
  
    std::vector<std::vector<float> >   m_signal;
    std::vector< float > m_fileFloatData;
    int m_dataType;
	int m_rows;
	int m_columns;
	int m_frames;
	int m_bands;

};

#endif /* RESTINGSTATENETWORK_H_ */
