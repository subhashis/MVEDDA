/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Streamlines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EDDA_GRID_H_
#define EDDA_GRID_H_

#include <cassert>
#include <vector>
#include <boost/any.hpp>
#include "common.h"
#include "core/interpolator.h"
#include "core/vector_matrix.h"
#include "edda_export.h"
//#include "Cell.h"

namespace edda{

enum CellType {
	TRIANGLE,
	CUBE,
	POLYGON,
	TETRAHEDRON,
};

enum InterpType {
    TRI_LERP,
    WEIGHTED_SUM
};



//////////////////////////////////////////////////////////////////////////
/// about the advection point
//////////////////////////////////////////////////////////////////////////
typedef struct PointInfo
{
    VECTOR3 phyCoord;
    VECTOR3 interpolant;    // interpolation coefficients
    int fromCell;   // advecting result from which cell, mainly used for unstructured grid
    int inCell;     // in which cell

    PointInfo()
    {
        fromCell = inCell = -1;
    }

    PointInfo(const VECTOR3& pcoord, const VECTOR3& coeff, int fCell, int iCell)
    {
        phyCoord = pcoord;
        interpolant = coeff;
        fromCell = fCell;
        inCell = iCell;
    }
}PointInfo;

//////////////////////////////////////////////////////////////////////////
///
/// \brief base class for grid
///
//////////////////////////////////////////////////////////////////////////
class EDDA_EXPORT Grid
{
public:
    Grid() {}
    virtual ~Grid() {}
	// physical coordinate of vertex verIdx
    virtual ReturnStatus at_vertex(int verIdx, VECTOR3& pos) = 0;
    // Return interpolation result
    //virtual ReturnStatus at_phys(VECTOR3& pos, AbstractDataArray *array, boost::any output) = 0;
	// get vertex list of a cell
    virtual ReturnStatus getCellVertices(int cellId, std::vector<size_t>& vVertices) = 0;
	// get the cell id and also interpolating coefficients for the given physical position
    virtual ReturnStatus phys_to_cell(PointInfo& pInfo) = 0;
    // interpolation: replaced by templated interpolator
    //virtual ReturnStatus interpolate(std::vector<boost::any>& vData, VECTOR3 coeff, boost::any &output) = 0;
	// the volume of cell
    virtual double cellVolume(int cellId) = 0;
	// type of cell
    virtual CellType getCellType() = 0;
    // how to interpolate
    virtual InterpType getInterpType() = 0;
	// get min and maximal boundary
    virtual void boundary(VECTOR3& minB, VECTOR3& maxB) = 0;
	// set bounding box
    virtual void setBoundary(VECTOR3& minB, VECTOR3& maxB) = 0;
	// get grid spacing in x,y,z dimensions
    virtual void getGridSpacing(int cellId, float& xspace, float& yspace, float& zspace) = 0;
	// boundary intersection
    virtual void boundaryIntersection(VECTOR3& intersectP, VECTOR3& startP,
                      VECTOR3& endP,float* stepSize, float oldStepSize) = 0;
	// whether the point is in the bounding box
    virtual bool isInBBox(VECTOR3& pos) = 0;
	// whether the point is in the bounding box not counting ghost cells
    virtual bool isInRealBBox(VECTOR3& p) = 0;
    virtual bool isInRealBBox(VECTOR3& pos, float t) = 0;

protected:
	// reset parameters
	virtual void Reset(void) = 0;
	// compute bounding box
    virtual void computeBBox(void) = 0;
	// whether in a cell
	virtual bool isInCell(PointInfo& pInfo, const int cellId) = 0;
};

//////////////////////////////////////////////////////////////////////////
///
/// \brief Cartesian Grid (Regular and Irregular)
///
//////////////////////////////////////////////////////////////////////////
class EDDA_EXPORT CartesianGrid : public Grid
{
public:
  // constructor and destructor
  CartesianGrid(int xdim, int ydim, int zdim);
  CartesianGrid();
  ~CartesianGrid();

  // get dimensions
  virtual int *getDimension() ;
  // physical coordinate of vertex verIdx
  virtual ReturnStatus at_vertex(int verIdx, VECTOR3& pos) =0;
  // Return interpolation result
  virtual ReturnStatus getIndex(int i, int j, int k, size_t &idx);
  // get vertex list of a cell
  virtual ReturnStatus getCellVertices(int cellId, std::vector<size_t>& vVertices) =0;
  // get the cell id and also interpolating coefficients for the given physical position
  virtual ReturnStatus phys_to_cell(PointInfo& pInfo) =0;

  // interpolation
  //template <class T>
  //virtual void interpolate(std::vector<boost::any>& vData, VECTOR3 coeff, boost::any &output) = 0;
  // the volume of cell
  virtual double cellVolume(int cellId) = 0;
  // type of cell
  virtual CellType getCellType() = 0;
  // how to interpolate
  virtual InterpType getInterpType() { return TRI_LERP; }
  // get min and maximal boundary
  virtual void boundary(VECTOR3& minB, VECTOR3& maxB) = 0;
  // set bounding box (includes ghost cells)
  virtual void setBoundary(VECTOR3& minB, VECTOR3& maxB) = 0;
  // set bounding box (does not include ghost cells)
  virtual void setRealBoundary(VECTOR4& minB, VECTOR4& maxB) = 0;
  // get grid spacing in x,y,z dimensions
  virtual void getGridSpacing(int cellId, float& xspace, float& yspace, float& zspace) = 0;
  // boundary intersection
  virtual void boundaryIntersection(VECTOR3& intersectP, VECTOR3& startP,
                                    VECTOR3& endP,float* stepSize, float oldStepSize) = 0;
protected:
  // reset parameters
  void Reset(void);
  // dimension related
  inline int xdim(void) { return m_nDimension[0];}
  inline int ydim(void) { return m_nDimension[1];}
  inline int zdim(void) { return m_nDimension[2];}
  inline int xcelldim(void) {return (m_nDimension[0] - 1);}
  inline int ycelldim(void) {return (m_nDimension[1] - 1);}
  inline int zcelldim(void) {return (m_nDimension[2] - 1);}

  int m_nDimension[3];				// dimension

  // min and maximal boundary (includes ghost cells)
  VECTOR3 m_vMinBound, m_vMaxBound;

  // min and maximal boundary (does not include ghost cells)
  VECTOR4 m_vMinRealBound, m_vMaxRealBound;
};

//////////////////////////////////////////////////////////////////////////
///
/// \brief regular cartesian grid
///
//////////////////////////////////////////////////////////////////////////
// map coordinates in computational space to physical space
#define UCGridPhy2Comp(x, y, f) (((x) - (y))*(f))

class EDDA_EXPORT RegularCartesianGrid : public CartesianGrid
{
private:
  float mappingFactorX;				// mapping from physical space to computational space
  float mappingFactorY;
  float mappingFactorZ;
  float oneOvermappingFactorX;
  float oneOvermappingFactorY;
  float oneOvermappingFactorZ;
  float gridSpacing;			        // the minimal grid spacing of all dimensions

public:
  RegularCartesianGrid(int xdim, int ydim, int zdim);
  RegularCartesianGrid();
  ~RegularCartesianGrid();
  // physical coordinate of vertex verIdx
  ReturnStatus at_vertex(int verIdx, VECTOR3& pos);
  // Return interpolation result
  //ReturnStatus at_phys(VECTOR3& pos, AbstractDataArray *array, boost::any output) ;
  // get vertex list of a cell
  ReturnStatus getCellVertices(int cellId, std::vector<size_t>& vVertices);
  // get the cell id and also interpolating coefficients for the given physical position
  ReturnStatus phys_to_cell(PointInfo& pInfo);
  // the volume of cell
  double cellVolume(int cellId);
  // cell type
  CellType getCellType() {return CUBE;}
  // set bounding box (includes ghost cells)
  void setBoundary(VECTOR3& minB, VECTOR3& maxB);
  // set bounding box (does not include ghost cells)
  void setRealBoundary(VECTOR4& minB, VECTOR4& maxB);
  // get min and maximal boundary
  void boundary(VECTOR3& minB, VECTOR3& maxB);
  // get grid spacing in x,y,z dimensions
  void getGridSpacing(int cellId, float& xspace, float& yspace, float& zspace)
  { xspace = oneOvermappingFactorX; yspace = oneOvermappingFactorY; zspace = oneOvermappingFactorZ; }
  void boundaryIntersection(VECTOR3&, VECTOR3&, VECTOR3&, float*, float);
  // whether the point is in the bounding box
  bool isInBBox(VECTOR3& pos);
  // whether the point is in the bounding box (not counting ghost cells)
  bool isInRealBBox(VECTOR3& pos);
  bool isInRealBBox(VECTOR3& pos, float t);


protected:
  void reset(void);
  // compute bounding box
  void computeBBox(void);
  // whether in a cell
  bool isInCell(PointInfo& pInfo, const int cellId);
};

/* 
//Comment the following code out since they have not been implemented
// 
//////////////////////////////////////////////////////////////////////////
//
//	irregular cartesian grid
//
//////////////////////////////////////////////////////////////////////////
class IrregularCartesianGrid : public CartesianGrid
{
private:
	float* m_pXSpacing;			// space array for x, y, z dimension
	float* m_pYSpacing;
	float* m_pZSpacing;

public:
	IrregularCartesianGrid(int xdim, int ydim, int zdim);
	IrregularCartesianGrid();
	~IrregularCartesianGrid();

protected:
	void Reset(void);
};

//////////////////////////////////////////////////////////////////////////
//
// curvilinear grid
//
//////////////////////////////////////////////////////////////////////////
class CurvilinearGrid : public Grid
{
private:
  int m_nDimension[3];				// dimension

public:
  // constructor and deconstructor
  CurvilinearGrid(int xdim, int ydim, int zdim);
  CurvilinearGrid();
  ~CurvilinearGrid();
};

*/ 


#if 0
//////////////////////////////////////////////////////////////////////////
//
// irregular grid
//
//////////////////////////////////////////////////////////////////////////
class IrregularGrid : public Grid
{
private:
	int m_nNodeNum;						// number of nodes
	int m_nTetraNum;					// number of tetras
	CVertex* m_pVertexGeom;				// geometry of all vertices
	CTetra* m_pTetra;					// tetra
	TetraInfo* m_pTetraInfo;			// pre-computed tetra information
	TVertex* m_pVertexTopo;				// vertex topology
	VECTOR3 m_vMinBound, m_vMaxBound;	// min and maximal boundary
	bool m_bTetraInfoInit;				// whether the tetra information is pre-computed

public:
	// constructor and deconstructor
	IrregularGrid();
	IrregularGrid(int nodeNum, int tetraNum, CVertex* pVertexGeom, CTetra* pTetra, TVertex* pVertexTopo);
	~IrregularGrid();

	// from virtual functions
    virtual void reset(void);
    virtual void getDimension(int& xdim, int& ydim, int& zdim) ;
    virtual ReturnStatus at_vertex(int verIdx, VECTOR3& pos);
    //virtual bool at_phys(VECTOR3& pos);
    virtual ReturnStatus getCellVertices(int cellId, std::vector<int>& vVertices);
    virtual ReturnStatus phys_to_cell(PointInfo& pInfo);
    //virtual void interpolate(std::vector<boost::any>& vData, VECTOR3 coeff, boost::any &output);
    virtual double cellVolume(int cellId);
    virtual bool isInCell(PointInfo& pInfo, const int cellId);
    virtual CellType getCellType() {return TETRAHEDRON;}

    virtual void computeBBox(void);
    virtual void setBoundary(VECTOR3& minB, VECTOR3& maxB);
    virtual void boundary(VECTOR3& minB, VECTOR3& maxB);
    virtual bool isInBBox(VECTOR3& pos);
    virtual bool isInRealBBox(VECTOR3& pos);
    virtual bool isInRealBBox(VECTOR3& pos, float t);

	// irregular specific functions
    virtual void SetTetraInfoInit(bool bInit);
    virtual bool GetTetraInfoInit(void);
    virtual int nextTetra(PointInfo& pInfo, int tetraId);
    virtual void PreGetP2NMatrix(MATRIX3& m, int cellId);
    virtual bool Physical2NaturalCoord(VECTOR3& nCoord, VECTOR3& pCoord, int cellId);

    virtual void GetGridSpacing(int cellId, float& xspace, float& yspace, float& zspace) {assert(false);}
	virtual void BoundaryIntersection(VECTOR3& intersectP, VECTOR3& startP, 
                      VECTOR3& endP,float* stepSize, float oldStepSize){assert(false);}

};

#endif
float getStepSize(VECTOR3& p, VECTOR3& p1, VECTOR3& p2, float oldStepSize);

}  // namesapce edda
#endif
