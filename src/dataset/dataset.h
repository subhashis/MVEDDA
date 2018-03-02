// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DATASET_H_
#define DATASET_H_

#include <vector>
#include <iostream>
#include <memory>
#include "dataset/distr_array.h"
#include "curvilinear_grid.h"
#include "distributions/variant.h"
#include "grid.h"
namespace edda {

///
/// \brief Holds all information of a dataset, which includes: 1) Geometry and 2) Data array.
///
/// You may use make_Dataset() to create a Dataset with shorter codes.
/// \param T The return type of at_phys() and at_comp()
///
template <typename T>
class Dataset {
protected:
	Grid *pGrid;
	DistrArray *pArray = 0; //will be deprecated
	std::vector<DistrArray *> pVector;
public:
	Dataset(){
		std::cout << "test\n";
	}
	Dataset(Grid *pGrid)
	{
		this->pGrid = pGrid;
	}
	Dataset(Grid *pGrid, DistrArray *pArray) {
		this->pGrid = pGrid;
		this->pArray = pArray;
	}
	Dataset(Grid *pGrid, std::vector<DistrArray *>& pVector)
	{
		this->pGrid = pGrid;
		this->pVector = pVector;
	}
	~Dataset() {
		if (pGrid)
			delete pGrid;
		if (pArray)
			delete pArray;
		pVector.clear();
	}

	Grid *getGrid() { return pGrid; }
	DistrArray *getArray() { return pArray; } //will be deprecated

	int getNumDistrArray(){ return pVector.size(); }
	DistrArray *getArray(int i) {
		if (i < 0 || i >= pVector.size()){
			//TODO: throw a proper exception
		}
		else{
			return pVector[i];
		}
	}


    ///
    /// \brief Get the dimension of the cartesian-grid data.
    ///
    /// Non-cartesian grid data are currently not supported.
    ///
    int *getDimension() {
      CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(pGrid) ;
      if (!cartesianGrid)
        throw NotImplementedException();
      return cartesianGrid->getDimension();
    }

	///
	/// \brief Get the spacing of the cartesian-grid data.
	///
	/// Non-cartesian grid data are currently not supported.
	void getSpacing(float& xspace, float& yspace, float& zspace) {
		CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(pGrid);
		if (!cartesianGrid)
			throw NotImplementedException();
		cartesianGrid->getGridSpacing(0, xspace, yspace, zspace);
		return;
	}

    ///
    /// \brief Return the interpolated data at a given position
    /// \param pos Position to query
    /// \param[out] The query result.  Valid only when the function returns SUCCESS
    /// \return SUCCESS if success, or OUT_OF_BOUND if the position is out of boundary.
    ///
    ReturnStatus at_phys(const VECTOR3 &pos, T &output) const {
        ReturnStatus r;
        switch (pGrid->getInterpType())
        {
        case TRI_LERP:
        {
            PointInfo pinfo;
            pinfo.phyCoord = pos ;
            r = pGrid->phys_to_cell(pinfo);
            if (r != SUCCESS)
                return r;

            std::vector<size_t> vVertices;
            r = pGrid->getCellVertices(pinfo.inCell, vVertices);
            if (r != SUCCESS)
                return r;

            int i;
            std::vector<T> vData(vVertices.size());
            for (i=0; i<vVertices.size(); i++)
            {
              getSample(vVertices[i], vData[i]);
            }

            output = triLerp(vData[0], vData[1], vData[2], vData[3],
                             vData[4], vData[5], vData[6], vData[7],
                             pinfo.interpolant.getData());
            break;
        }
        default:
            throw NotImplementedException();
            break;
        }

        return SUCCESS;
    }

    ///
    /// \brief Return the interpolated data at a given position
    /// \param pos Position to query
    /// \param[out] The query result.  Valid only when the function returns SUCCESS
    /// \return SUCCESS if success, or OUT_OF_BOUND if the position is out of boundary.
    ///
    ReturnStatus at_phys_new(const VECTOR3 &pos, std::vector<T> &output) const {
        ReturnStatus r;
        switch (pGrid->getInterpType())
        {
        case TRI_LERP:
        {
            PointInfo pinfo;
            pinfo.phyCoord = pos ;
            r = pGrid->phys_to_cell(pinfo);
            if (r != SUCCESS)
                return r;

            std::vector<size_t> vVertices;
            r = pGrid->getCellVertices(pinfo.inCell, vVertices);
            if (r != SUCCESS)
                return r;

            int i,d;
            //std::vector<T> vData(vVertices.size());
            std::vector<std::vector<T>> vData(pVector.size(), std::vector<T>(vVertices.size()));
            for(d=0; d<pVector.size(); d++)\
            {
              for (i=0; i<vVertices.size(); i++)
              {
                //getSample(vVertices[i], vData[i]);
                getSampleNew(vVertices[i], vData[d][i], d);
              }
            }

            for(int d=0; d<pVector.size(); d++)
            {
              T result = triLerp(vData[d][0], vData[d][1], vData[d][2], vData[d][3],
                               vData[d][4], vData[d][5], vData[d][6], vData[d][7],
                               pinfo.interpolant.getData());
              output.push_back(result);
            }
            break;
        }
        default:
            throw NotImplementedException();
            break;
        }

        return SUCCESS;
    }


    ///
    /// \brief Return the data in the computational space.
    /// \param i,j,k The coordinates in the computational space.
    /// \return The query result.
    /// Only for structured grids.
    /// If <i,j,k> is out of boundary, the function will throw OutOfBoundException.
    ///
    const T at_comp(int i, int j, int k) const {
      ReturnStatus r;
      CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(pGrid) ;

      if (cartesianGrid) {
        size_t idx;
        r = cartesianGrid->getIndex(i,j,k, idx);
        if (r!=SUCCESS)
          throw OutOfBoundException();

        //return boost::any_cast<const T&>( pArray->getItem(idx) );
        T data;
        getSample(idx, data);
        return data;
      }

      // TODO for other grid types
      throw NotImplementedException();
    }

    const std::vector<T> at_comp_new(int i, int j, int k) const {
      ReturnStatus r;
      CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(pGrid) ;

      if (cartesianGrid) 
      {
        size_t idx;
        r = cartesianGrid->getIndex(i,j,k, idx);
        if (r!=SUCCESS)
          throw OutOfBoundException();

        //return boost::any_cast<const T&>( pArray->getItem(idx) );
        std::vector<T> sampleVector;
        for(int i=0; i<pVector.size(); i++)
        {
          T data;
          getSampleNew(idx, data, i);
          sampleVector.push_back(data);
        }
        
        return sampleVector;
      }
      else
      {
        // TODO for other grid types
        throw NotImplementedException();
      }

      
    }

    ///
    /// set distribution
    ///
    ReturnStatus set_distr(int i, int j, int k, int distId, dist::Variant newDist) const {
      ReturnStatus r;
      CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(pGrid) ;

      if (cartesianGrid) {
        size_t idx;
        r = cartesianGrid->getIndex(i,j,k, idx);
        if (r!=SUCCESS)
          throw OutOfBoundException();

        dist::Variant oldDist = pVector[distId]->getDistr(idx);
        std::cout << "oldDistributionType = " << getName(oldDist) << std::endl;
        std::cout << "newDistributionType = " << getName(newDist) << std::endl;
        
      }
      else{
        // TODO for other grid types
        throw NotImplementedException();
      }
    }

    ///
    /// \brief Return the distribution in the computational space.
    /// \param i,j,k The coordinates in the computational space.
    /// \return The query result.
    /// Only for structured grids.
    /// If <i,j,k> is out of boundary, the function will throw OutOfBoundException.
    ///
    const dist::Variant at_comp_distr(int i, int j, int k) const {
      ReturnStatus r;
      CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(pGrid) ;

      if (cartesianGrid) {
        size_t idx;
        r = cartesianGrid->getIndex(i,j,k, idx);
        if (r!=SUCCESS)
          throw OutOfBoundException();

        return pArray->getDistr(idx);
      }

      // TODO for other grid types
      throw NotImplementedException();
    }

    const std::vector<dist::Variant> at_comp_distr_new(int i, int j, int k) const {
      ReturnStatus r;
      CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(pGrid) ;

      if (cartesianGrid) {
        size_t idx;
        r = cartesianGrid->getIndex(i,j,k, idx);
        if (r!=SUCCESS)
          throw OutOfBoundException();

        std::vector<dist::Variant> distrList;
        for(int i=0; i<pVector.size(); i++)
        {
          distrList.push_back(pVector[i]->getDistr(idx));
        }

        return distrList;
      }

      // TODO for other grid types
      throw NotImplementedException();
    }

    ///
    /// \brief Return the distribution vector in the computational space.
    /// \param i,j,k The coordinates in the computational space.
    /// \return The query result, in std::vector.
    /// Only for structured grids.
    /// If <i,j,k> is out of boundary, the function will throw OutOfBoundException.
    ///
    const std::vector<dist::Variant>  at_comp_distr_vector(int i, int j, int k) const {
      ReturnStatus r;
      CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(pGrid) ;

      if (cartesianGrid) {
        size_t idx;
        r = cartesianGrid->getIndex(i,j,k, idx);
        if (r!=SUCCESS)
          throw OutOfBoundException();

        return pArray->getDistrVector(idx);
      }

      // TODO for other grid types
      throw NotImplementedException();
    }

protected:
    void getSample(int idx, Real &data) const {
      data = pArray->getScalar(idx) ;
    }

    void getSampleNew(int idx, Real &data, size_t distrId) const {
      data = pVector[distrId]->getScalar(idx) ;
    }

    template <int N>
    void getSample(int idx, Vector<Real, N> &data) const {
      std::vector<Real> varvec =  pArray->getVector(idx) ;
      assert(N == varvec.size());
      for (int c=0; c<N; c++) {
         data[c] = varvec[c];
      }
    }

};


///
/// \brief Create a shared pointer of Dataset
///
template <typename T>  // Return type of at_phys
inline std::shared_ptr< Dataset<T> >
make_Dataset(Grid *pGrid, DistrArray *pArray)
{
    return std::shared_ptr< Dataset<T> >(
                new Dataset<T> (pGrid, pArray) );
}

} // namespace edda

#endif // STRUCTURED_GRID_DATASET_H_
