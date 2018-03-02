// C++ matrix inversion (boost::ublas)
// Another recurrent question is "how do I invert a matrix in [my library of choice]?. Neither boost::ublas nor MTL documentation answers this satisfactorily.
// If you're using boost::ublas, you can find the InvertMatrix funcion here: http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
// And a simple example below:

#include "common.h"
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace edda
{
  // for linear algebra
  typedef boost::numeric::ublas::vector<Real> ublas_vector;
  typedef boost::numeric::ublas::matrix<Real> ublas_matrix;
  using namespace boost::numeric;

 /* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool invert_matrix(const ublas::matrix<T>& input, ublas::matrix<T>& inverse)
{
	typedef ublas::permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	ublas::matrix<T> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = ublas::lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
        //inverse.assign(ublas::identity_matrix<T> (A.size1()));
        inverse = ublas::identity_matrix<T> (A.size1());

	// backsubstitute to get the inverse
	ublas::lu_substitute(A, pm, inverse);

	return true;
}

inline
int determinant_sign(const boost::numeric::ublas::permutation_matrix<std::size_t>& pm)
{
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}

inline
double determinant(const boost::numeric::ublas::matrix<Real>& m ) {
  ublas_matrix mm = m;
    boost::numeric::ublas::permutation_matrix<std::size_t> pm(m.size1());
    double det = 1.0;
    if( boost::numeric::ublas::lu_factorize(mm,pm) ) {
        det = 0.0;
    } else {
        for(int i = 0; i < m.size1(); i++)
            det *= m(i,i); // multiply by elements on diagonal
        det = det * determinant_sign( pm );
    }
    return det;
}

} //edda
