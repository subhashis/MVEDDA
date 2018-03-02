
#include <core/vector_matrix.h>

namespace edda {
	// set matrix to identity matrix
	void MATRIX3::Identity()
	{
		for (int i = 0; i < Dimension(); i++) {
			mat[i].Zero();
			mat[i][i] = 1.0;
		};
	}


	// determinant of the matrix
	float MATRIX3::det() {
		float d01d12md11d02 = mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2];
		float d01d22md21d02 = mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2];
		float d11d22md21d12 = mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2];

		float det = mat[0][0] * d11d22md21d12 - mat[1][0] * d01d22md21d02 +
			mat[2][0] * d01d12md11d02;

		return det;
	}

	//added by lijie
	MATRIX3 MATRIX3::transpose()
	{
		MATRIX3 inv;
		for (int i = 0; i<3; i++)
		{
			for (int j = 0; j<3; j++)
			{
				inv[j][i] = mat[i][j];
			}
		}
		return inv;
	}

	//added by lijie
	// Jimmy: process in double for better precision
	int MATRIX3::inverse(MATRIX3& m)
	{
		double mat[3][3] =
		{ { this->mat[0][0], this->mat[0][1], this->mat[0][2] },
		{ this->mat[1][0], this->mat[1][1], this->mat[1][2] },
		{ this->mat[2][0], this->mat[2][1], this->mat[2][2] } };

		// 9 floating-point ops
		double d01d12md11d02 = mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2];
		double d01d22md21d02 = mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2];
		double d11d22md21d12 = mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2];

		// 5 floating-point ops
		double det =
			mat[0][0] * d11d22md21d12 -
			mat[1][0] * d01d22md21d02 +
			mat[2][0] * d01d12md11d02;

		if (0.0 == det) return 0;

		double det_inv = 1.0 / det;

		// 19 floating-point ops
		m[0][0] = d11d22md21d12 *det_inv;
		m[0][1] = -d01d22md21d02 *det_inv;
		m[0][2] = d01d12md11d02 *det_inv;
		m[1][0] = (mat[2][0] * mat[1][2] - mat[1][0] * mat[2][2])*det_inv;
		m[1][1] = (mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2])*det_inv;
		m[1][2] = (mat[1][0] * mat[0][2] - mat[0][0] * mat[1][2])*det_inv;
		m[2][0] = (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1])*det_inv;
		m[2][1] = (mat[2][0] * mat[0][1] - mat[0][0] * mat[2][1])*det_inv;
		m[2][2] = (mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1])*det_inv;

		return 1;

	}


	// return m0 * v0
	VECTOR3 operator *(const MATRIX3 & m0, const VECTOR3 & v0)
	{
		VECTOR3 result;

		result[0] = dot(m0(0), v0);
		result[1] = dot(m0(1), v0);
		result[2] = dot(m0(2), v0);

		return(result);
	}

}