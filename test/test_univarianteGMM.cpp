#include "distributions/joint_GMM.h"
#include "distributions/joint_gaussian.h"
#include "vector"
#include "bmp_image.h"
#include "dataset/distr_array.h"
#include "io/edda_vtk_reader.h"
#include "io/edda_vtk_writer.h"
#include "io/edda_reader.h"
#include "io/edda_writer.h"
#include "distributions/gmm.h"
#include "dataset/distr_array.h"
#include "distributions/estimate_gmm.h"

using namespace edda;
using namespace std;
using namespace edda::dist;

float getVoxel(float* data, int d0, int d1, int d2, int D0, int D1, int D2)
{
	int idx = d2*D1*D0 + d1*D0 + d0;
	return data[idx];
}

void setVoxel(float* data, int d0, int d1, int d2, int D0, int D1, int D2, float value)
{
	int idx = d2*D1*D0 + d1*D0 + d0;
	data[idx] = value;
}

//using image to test joint Gaussian GMM
int main()
{
	int Dim0 = 100;
	int Dim1 = 100;
	int Dim2 = 100;
	
	//modify this to the local path
	//this sampled data is in the "sampled_data" folder of the source folder
	FILE* fpr = fopen("C:\\GravityLabDataSet\\Edda\\edda\\sample_data\\subIsabel.raw", "rb");//load the raw test file
	FILE* fpw = fopen("resampled_subIsabel.raw", "wb");//output resampled data

	float* rbuffer = (float*)malloc(sizeof(float)* Dim0 * Dim1 * Dim2);
	float* wbuffer = (float*)malloc(sizeof(float)* Dim0 * Dim1 * Dim2);

	fread(rbuffer, sizeof(float), Dim0 * Dim1 * Dim2, fpr);
	
	//down-sampled block size
	int blockSize = 10;
	int nGmmCompMin = 3;//number of Gaussian compoenents
	int nGmmCompMax = 7;//number of Gaussian compoenents

	////number of row and col after down-sampleing
	int dsD0 = Dim0 / blockSize;
	int dsD1 = Dim1 / blockSize;
	int dsD2 = Dim2 / blockSize;
	
	//GMM array
	shared_ary<GMM> array(new GMM[dsD0*dsD1*dsD2], dsD0*dsD1*dsD2);
	thrust::default_random_engine rng;//random engine for getJointSample()

	double* trainingData = (double*)malloc(sizeof(double)* blockSize* blockSize* blockSize);

	//loop: go through each block
	for (int dsd0 = 0; dsd0 < dsD0; dsd0++){
		for (int dsd1 = 0; dsd1 < dsD1; dsd1++){
			for (int dsd2 = 0; dsd2 < dsD2; dsd2++){
			printf("%d %d %d\n", dsd0, dsd1, dsd2);

			//prepare training array of this local block
			int cnt = 0;
			for (int d0 = dsd0*blockSize; d0<(dsd0 + 1)*blockSize; d0++) {
				for (int d1 = dsd1*blockSize; d1<(dsd1 + 1)*blockSize; d1++) {
					for (int d2 = dsd2*blockSize; d2 < (dsd2 + 1)*blockSize; d2++) {
						trainingData[cnt] = getVoxel(rbuffer, d0, d1, d2, Dim0, Dim1, Dim2);
						cnt++;
					}
				}
			}

			//EM in Edda
			int distAryIdx = dsd0*dsD1*dsD2 + dsd1*dsD2 + dsd2;
			int nGmmComp = (int)(nGmmCompMax - nGmmCompMin) * (rand() / (float)RAND_MAX) + nGmmCompMin;
			array[distAryIdx] = eddaComputeGMM(trainingData, blockSize*blockSize*blockSize, nGmmComp);

			//resample this block and write to image		
			for (int i = 0; i < blockSize; i++){
				for (int j = 0; j < blockSize; j++){
					for (int k = 0; k < blockSize; k++){
						float sample = getSample(array[distAryIdx], rng);

						int rawd0 = dsd0*blockSize + i;//u and v in the original image resolution
						int rawd1 = dsd1*blockSize + j;
						int rawd2 = dsd2*blockSize + k;

						setVoxel(wbuffer, rawd0, rawd1, rawd2, Dim0, Dim1, Dim2, sample);
					}
				}
			}
			}
		}
	}

	fwrite(wbuffer, sizeof(float), Dim0 * Dim1 * Dim2, fpw);

	free(rbuffer);
	free(wbuffer);
	fclose(fpr);
	fclose(fpw);
}

