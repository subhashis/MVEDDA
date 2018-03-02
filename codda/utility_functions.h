#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>



#include <Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

#include <boost/math/distributions/normal.hpp> // for normal_distribution
  using boost::math::normal; // typedef provides default type is double.
#include <limits>
  using std::numeric_limits;
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/assignment.hpp> 

#include "distributions/joint_gaussian.h"

#include <vinecopulib.hpp>

#include <opencv2/opencv.hpp>
#include <opencv2/ml/ml.hpp>
#include <opencv2/legacy/legacy.hpp>


using namespace vinecopulib;
using namespace std;
using namespace edda;
using namespace cv;


/*
  generating samples from multivariate standard normal distribution.
*/
namespace Eigen {
namespace internal {
template<typename Scalar> 
struct scalar_normal_dist_op 
{
  static boost::mt19937 rng;    // The uniform pseudo-random algorithm
  mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

  EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

  template<typename Index>
  inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
};

template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;

template<typename Scalar>
struct functor_traits<scalar_normal_dist_op<Scalar> >
{ enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };
} // end namespace internal
} // end namespace Eigen

/*
  Draw nn samples from a size-dimensional normal distribution
  with a specified mean and covariance
*/


void genSNMVsample(int size, int nn, float  *arr, float  *sampleArray )
{
  //int size = 3; // Dimensionality (rows)
  //int nn=5;     // How many samples (columns) to draw
  Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
  Eigen::internal::scalar_normal_dist_op<double>::rng.seed(1); // Seed the rng

  // Define mean and covariance of the distribution
  Eigen::VectorXd mean(size);       
  Eigen::MatrixXd covar(size,size);

  cout << "covar:" << endl;
  for(int i=0; i<size; i++)
  {
    mean[i] = 0;
    for(int j=0; j<size; j++)
    {
      covar(i,j) = double(arr[i*size + j]);
      cout << "<" <<covar(i,j) << "," << arr[i*size + j] << ">\t" ;
    }
    cout << endl;
  } 

  Eigen::MatrixXd normTransform(size,size);

  Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);

  // We can only use the cholesky decomposition if 
  // the covariance matrix is symmetric, pos-definite.
  // But a covariance matrix might be pos-semi-definite.
  // In that case, we'll go to an EigenSolver
  if (cholSolver.info()==Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    // Use eigen solver
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
    normTransform = eigenSolver.eigenvectors() 
                   * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }

  Eigen::MatrixXd samples = (normTransform 
                           * Eigen::MatrixXd::NullaryExpr(size,nn,randN)).colwise() 
                           + mean;

  //std::cout << "Mean\n" << mean << std::endl;
  //std::cout << "Covar\n" << covar << std::endl;
  //std::cout << "Samples\n" << samples << std::endl;

  //double  *sampleArray = new double [nn*size];

  for(int i=0; i<nn; i++)
  {
    for(int j=0; j<size; j++)
    {
      sampleArray[i*size + j] = float(samples(j,i));
      //cout << sampleArray[i*size + j];

    }
    
  }

}

bool isNaN(float x) { 
  return x != x;
}


void genSNMVsample_EDDA(int size, int nsample, float  *arr, float  *sampleArray )
{
  int nVar = size;
  int nComp = 1;

  thrust::default_random_engine rng;

  
  ublas_vector mean = ublas_vector(nVar,0);
  ublas_matrix cov = ublas::zero_matrix<Real>(nVar,nVar);

  for(int i=0; i<nVar; i++)
  {
    for(int j=0; j<nVar; j++)
    {
      cov(i,j) = arr[i*size + j];
    }
    
  }

  dist::JointGaussian sn_distr(mean,cov);

  /*std::vector<Real> curSample = sn_distr.getJointSample(rng);

  cout << "\nmean0=" << curSample[0];
  cout << "\nmean1=" << curSample[1];
  cout << "\nmean2=" << curSample[2];
  */


  for(int s=0; s<nsample; s++)
  {
    std::vector<Real> curSample = sn_distr.getJointSample(rng);
    
    
    for(int v=0; v<curSample.size(); v++)
    {
      sampleArray[s*size + v] = curSample[v];
    }




  }

  


}

float unif_rand_func()
{
    std::random_device rd;    

    // Use Mersenne twister engine to generate pseudo-random numbers.
    std::mt19937 engine(rd());

    // "Filter" MT engine's output to generate pseudo-random integer values,
    // **uniformly distributed** on the closed interval [0, 99].  
    // (Note that the range is [inclusive, inclusive].)
    std::uniform_real_distribution<double> dist(0, 1);

    return dist(engine);
}




// Normal random variate sampler
double  box_muller(double  m, double  s)
{
    double  x1, x2, w, y1;
    double  y2;
    char use_last = 0;

    if (use_last)           /* use value from previous call */
    {
        y1 = y2;
        use_last = 0;
    }
    else
    {
        do
        {
            x1 = 2.0 * rand()/RAND_MAX - 1.0;
            x2 = 2.0 * rand()/RAND_MAX - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;
        use_last = 1;
    }

    return( m + y1 * s );
}


//inverse CDF of a normal distribution
double  inverseCDF_normalDistribution(double  m, double  std, double  q)
{
  normal dist(m,std);
  double  r = quantile(dist,q);
  return r;
}

double  standard_normalCDF(double  p)
{
  normal s;
  double  u = cdf(s,p);
  return u;
}

void get_minMax( float *dataArray, int size, float *min_max)
{
  float maxV = -10000000.0;
  float minV = 10000000.0;
  for(int i=0; i<size; i++)
  {
    float v = dataArray[i];
    if(v >= maxV)
      maxV = v;
    if(v <= minV)
      minV = v;
  }
  min_max[0] = minV;
  min_max[1] = maxV;

}

float get_average(float *dataArray, int size)
{
  float sum = 0;
  for(int i=0; i<size; i++)
  {
    sum += dataArray[i];
  }

  return sum/size;
}

float get_std(float *dataArray, float m, int size)
{
  float s2 = 0;
  for(int i=0; i<size; i++)
  {
    float d = (dataArray[i]-m)*(dataArray[i]-m);
    s2 += d;
  }

  float var = s2/size;

  float std_r = sqrt(var);

    if(std_r == 0)
        std_r += 0.001;

    return std_r;

}

float get_var(float *dataArray, float m, int size)
{
  float s2 = 0;
  for(int i=0; i<size; i++)
  {
    float d = (dataArray[i]-m)*(dataArray[i]-m);
    s2 += d;
  }

  float var = s2/size;

  return var;

}

/**
*   The Covariance
*/
float covariance(float *v1, float *v2, int size, float mean1, float mean2){
    
    float sum = (v1[0] - mean1) * (v2[0] - mean2);
    for (int i=1; i < size; i++){
        sum += (v1[i] - mean1) * (v2[i] - mean2);
    }
    return float(sum) / float(size);
}


/**
*   Pearson Correlation
*/
float pearson(float *v1, float *v2, int size){
    float mean1 = get_average(v1, size), mean2 = get_average(v2, size);
    float std1 = get_std(v1, mean1, size), std2 = get_std(v2, mean2, size);
    if (std1 * std2 == 0){
        std1 += 0.001;
        std2 += 0.001;
        //cout << "( a standard deviaton was 0 )";
        //return -2; // I dont know what to do here???
    }
    return covariance(v1,v2,size,mean1,mean2) / ( std1 * std2);
}





double  inverseCDF_histogram(double  *support, double  *cdfValue, double  unifValue, int nBins)
{
  double  maxCDF = -100000;
  double  minCDF = 100000;

  for(int i=0; i<nBins; i++){
      double  value = cdfValue[i];
      if( value< minCDF ) minCDF = value;
      if( value>maxCDF )  maxCDF = value;
  }
  if(unifValue < minCDF)
  {  
    return support[0];
  }
  else if(unifValue > maxCDF)
  {
    return support[nBins];
  }
  else
  {
    for(int i=0; i<nBins+1; i++)
    {
      if(unifValue >= cdfValue[i] && unifValue <= cdfValue[i+1])
      {
        double  ratio = (unifValue - cdfValue[i]) / (cdfValue[i+1] - cdfValue[i]);
        double  value = ratio * (support[i+1]-support[i]) + support[i];
        return value;
      }
    }
  }

}


double  inverseCDF_KDE(double  *support, double  *cdfValue, double  unifValue, int size)
{
  double  maxCDF = -100000;
  double  minCDF = 100000;

  for(int i=0; i<size; i++){
      double  value = cdfValue[i];
      if( value< minCDF ) minCDF = value;
      if( value>maxCDF )  maxCDF = value;
  }  
  if(unifValue < minCDF)
  {  
    return support[0];
  }
  else if(unifValue > maxCDF)
  {
    return support[size-1];
  }
  else
  {
    for(int i=0; i<size; i++)
    {
      if(unifValue >= cdfValue[i] && unifValue <= cdfValue[i+1])
      {
        double  ratio = (unifValue - cdfValue[i]) / (cdfValue[i+1] - cdfValue[i]);
        double  value = ratio * (support[i+1]-support[i]) + support[i];
        return value;
      }
    }
  }

}

float * get_uniform_samples(float *para_array, int numPar, int numSample)
{
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(5,5);
  mat << 1, 1, 1, 1, 1,
         2, 2, 2, 2, 0,
         3, 3, 3, 0, 0,
         4, 4, 0, 0, 0,
         5, 0, 0, 0, 0; 

  auto pair_copulas = Vinecop::make_pair_copula_store(5);

  // specify the pair copulas
  for(int i=0; i<numPar; i++)
  {
    auto par = Eigen::VectorXd::Constant(1, para_array[i]);
    int count = 0;
    for (auto& tree : pair_copulas) {
        for (auto& pc : tree) {
            pc = Bicop(BicopFamily::gaussian, 0, par);
            //cout << "c=" << count << ",";
            count += 1;
        }
        //cout << endl;
    }
  }

  Vinecop custom_model(pair_copulas, mat);
  // simulate data
  Eigen::MatrixXd data = custom_model.simulate(numSample);
  float *retData = new float[5*numSample];
  for(int i=0; i<numSample; i++)
  {
    for(int j=0; j<5; j++)
    {
      retData[i*5+j] = data(i,j);

    }
  }

  return retData;

  
}

//for PCC
float * get_uniform_samples1(float *para_array, int numPar, int numSample, float *pdf_values)
{
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(5,5);
  mat << 1, 1, 1, 1, 1,
         2, 2, 2, 2, 0,
         3, 3, 3, 0, 0,
         4, 4, 0, 0, 0,
         5, 0, 0, 0, 0; 

  auto pair_copulas = Vinecop::make_pair_copula_store(5);

  // specify the pair copulas
  for(int i=0; i<numPar; i++)
  {
    auto par = Eigen::VectorXd::Constant(1, para_array[i]);
    int count = 0;
    for (auto& tree : pair_copulas) {
        for (auto& pc : tree) {
            pc = Bicop(BicopFamily::gaussian, 0, par);
            //cout << "c=" << count << ",";
            count += 1;
        }
        //cout << endl;
    }
  }

  Vinecop custom_model(pair_copulas, mat);
  // simulate data
  Eigen::MatrixXd data = custom_model.simulate(numSample);
  float *retData = new float[5*numSample];
  for(int i=0; i<numSample; i++)
  {
    for(int j=0; j<5; j++)
    {
      retData[i*5+j] = data(i,j);

    }
  }

  auto pdf = custom_model.pdf(data);
  for(int i=0; i<pdf.rows(); i++)
  {
      //cout << "pdf " << i << " = " << pdf[i] << endl;
      pdf_values[i]= pdf[i];
  }

  return retData;

  
}


float get_average_pdf_values(float *para_array, int numPar, int numSample)
{
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(5,5);
  mat << 1, 1, 1, 1, 1,
         2, 2, 2, 2, 0,
         3, 3, 3, 0, 0,
         4, 4, 0, 0, 0,
         5, 0, 0, 0, 0; 

  auto pair_copulas = Vinecop::make_pair_copula_store(5);

  // specify the pair copulas
  for(int i=0; i<numPar; i++)
  {
    auto par = Eigen::VectorXd::Constant(1, para_array[i]);
    int count = 0;
    for (auto& tree : pair_copulas) {
        for (auto& pc : tree) {
            pc = Bicop(BicopFamily::gaussian, 0, par);
            //cout << "c=" << count << ",";
            count += 1;
        }
        //cout << endl;
    }
  }

  Vinecop custom_model(pair_copulas, mat);
  // simulate data
  Eigen::MatrixXd data = custom_model.simulate(numSample);

  auto pdf = custom_model.pdf(data);
  float pdf_sum = 0;
  for(int i=0; i<pdf.rows(); i++)
  {
      //cout << "pdf " << i << " = " << pdf[i] << endl;
      pdf_sum += pdf[i];
  }
  pdf_sum = pdf_sum/pdf.rows();

  return pdf_sum;

  
}


float get_transformed_spatial_coord(float u , int min, int max)
{
  float retValue = -1.0f;
  retValue = u*(max-min) + float(min);
  return retValue;
}

//for SPCC
float * generate_uniform_samples(float *para_array, int nv, int numSample)
{
  int numPar = (nv * (nv -1))/2;

  //create the correlation matrix
  float *cell_correl = new float[nv*nv];
  int c1 = 0;
  for(int i=0; i<nv; i++){
    for(int j=i+1; j<nv; j++){
      cell_correl[i*nv+j] = para_array[c1];
      if(fabs(para_array[c1]-0.00)<1e-6)
      {
        cell_correl[i*nv+j] += 0.001f;
      }
      cell_correl[j*nv+i] = cell_correl[i*nv+j];
      c1++;
    }
  }
  for(int i=0; i<nv; i++){
    cell_correl[i*nv+i] = 1.0;
  }

  float *all_sample = new float[nv*numSample];

  genSNMVsample_EDDA(nv,numSample,cell_correl,all_sample);

  int count_nan = 0;
  //step2: convert to uniform marginals with preserved correlation.
  for(int i=0; i<numSample*nv; i++)
  {
      if(isNaN(all_sample[i]))
      {
          all_sample[i] = unif_rand_func();
          count_nan++;
      }
      else
      {
          all_sample[i] = standard_normalCDF(all_sample[i]);
      }
      
      
  }

  //cout << "count_nan = " << count_nan << endl;
  return all_sample;


}

//function which returns the proper MV sample given all the parameters
float * getMVSamples(float *para_array, int nVar, int spatialVar, int numSample, std::vector<dist::Variant>* distrVec, int x1, int x2, int y1, int y2, int z1, int z2)
{
  int total_nVar = nVar + spatialVar;
  //int numPar = (total_nVar * (total_nVar -1))/2;

  float * new_samples;
  new_samples = generate_uniform_samples(para_array, total_nVar, numSample);

  //check the univariate distributions : unit test
  /*cout << "distr array size = " << distrVec->size() << endl;
  for(int i=0; i<distrVec->size(); i++){
    cout << (*distrVec)[i] << endl;
  }*/

  for(int var=0; var<nVar; var++){
    dist::GMM curDistr = boost::get<dist::GMM>((*distrVec)[var]);
    for(int s=0; s<numSample; s++){
      new_samples[s*(nVar+spatialVar) + var] = getInverseCDF_GMM(curDistr, new_samples[s*(nVar+spatialVar) + var]);
    }
  }


  for(int s=0; s<numSample; s++){
    new_samples[s*(nVar+spatialVar) + nVar + 0] = get_transformed_spatial_coord(new_samples[s*(nVar+spatialVar) + nVar + 0], x1, x2); 
    new_samples[s*(nVar+spatialVar) + nVar + 1] = get_transformed_spatial_coord(new_samples[s*(nVar+spatialVar) + nVar + 1], y1, y2); 
    new_samples[s*(nVar+spatialVar) + nVar + 2] = get_transformed_spatial_coord(new_samples[s*(nVar+spatialVar) + nVar + 2], z1, z2); 

  }

  return new_samples;
 

}

//function which returns the proper MV sample given all the parameters
float * getMVSamples_spatialDistr(float *para_array, int nVar, int spatialVar, int numSample, std::vector<dist::Variant>* distrVec)
{
  int total_nVar = nVar + spatialVar;
  //int numPar = (total_nVar * (total_nVar -1))/2;

  float * new_samples;
  new_samples = generate_uniform_samples(para_array, total_nVar, numSample);

  //check the univariate distributions : unit test
  /*cout << "distr array size = " << distrVec->size() << endl;
  for(int i=0; i<distrVec->size(); i++){
    cout << (*distrVec)[i] << endl;
  }*/

  for(int var=0; var<nVar; var++){
    dist::GMM curDistr = boost::get<dist::GMM>((*distrVec)[var]);
    for(int s=0; s<numSample; s++){
      new_samples[s*(nVar+spatialVar) + var] = getInverseCDF_GMM(curDistr, new_samples[s*(nVar+spatialVar) + var]);
    }
  }

  for(int var=nVar; var<total_nVar; var++){
    dist::Histogram curHist = boost::get<dist::Histogram>((*distrVec)[var]);
    for(int s=0; s<numSample; s++){
      new_samples[s*(nVar+spatialVar) + var] = curHist.getInverseCdf(new_samples[s*(nVar+spatialVar) + var]);
    }
  }

  return new_samples;
 

}

void readData(const string &filename, float *dataArray, int size)
{
    FILE *fIn;
    //double *inData = new double[size];
    //string fname = path + filename + to_string(eId) + ".raw";
    fIn = fopen(filename.c_str(),"rb");
    cout << "opening " << filename << endl;
    if(!fIn)
    {
        fprintf(stderr, "Error opening files\n");
        exit(13);
    }
    int read_count = fread(dataArray, sizeof(float), size, fIn);
    if(read_count != size)
    {
      fprintf(stderr, "Failed to read the complete file\n");
      exit(14);
    }
    fclose(fIn);   

}
/*void get_minMax( float *dataArray, int size, float *min_max)
{
  float maxV = -10000000.0;
  float minV = 10000000.0;
  for(int i=0; i<size; i++)
  {
    float v = dataArray[i];
    if(v >= maxV)
      maxV = v;
    if(v <= minV)
      minV = v;
  }
  min_max[0] = minV;
  min_max[1] = maxV;

}*/




void computeGMM(float *dataArray, int nElements, const int nGaussians, std::vector<float>* distrProperties)
{
  Mat samples(nElements, 1, CV_32FC1);
  Mat labels;

  for(int i=0; i<nElements; i++){
    samples.at<float> (i,0) = dataArray[i]; //put data in array
  }

  CvEM em_model;
  CvEMParams params;
  params.covs      = NULL;
  params.means     = NULL;
  params.weights   = NULL;
  params.probs     = NULL;
  params.nclusters = nGaussians;
  params.cov_mat_type       = CvEM::COV_MAT_GENERIC;
  params.start_step         = CvEM::START_AUTO_STEP;
  params.term_crit.max_iter = 300;
  params.term_crit.epsilon  = 0.1;
  params.term_crit.type     = CV_TERMCRIT_ITER|CV_TERMCRIT_EPS;


  // cluster the data
  em_model.train( samples, Mat(), params, &labels );

  const Mat means = em_model.getMeans();    
  const Mat weights = em_model.getWeights();
  vector<cv::Mat>  covs;
  em_model.getCovs(covs);


  for(int m=0;m<nGaussians;m++)
  {
      float cur_mean = means.at<double>(0,m);   
      float cur_stdev = sqrt(covs[m].at<double>(0, 0));        
      float cur_weight = weights.at<double>(0,m);

      distrProperties->push_back(cur_mean);
      distrProperties->push_back(cur_stdev);
      distrProperties->push_back(cur_weight);      

      //cout << "m=" << my_mean << " ,sd=" << my_stdev << " ,wt=" << my_weight << endl;
      //listOfGaussians.gaussians.push_back(newGaussian);
  }


}

//TODO: inverse cdf for gamma and beta.



/*def inverseCDFofKDE(support, cdfValue, unifValue):
    maxCDF = max(cdfValue)
    minCDF = min(cdfValue)
    if unifValue > maxCDF or unifValue < minCDF:
        #raise ValueError('A very specific bad thing happened')
        return -1000
    else:
        l = support.shape
        for i in range(0,l[0]):
            if unifValue >= cdfValue[i] and unifValue <= cdfValue[i+1]:
                left = i;
                right = i+1;
                ratio = (unifValue-cdfValue[i])/(cdfValue[i+1]-cdfValue[i])
                value = ratio * (support[i+1]-support[i]) + support[i]
                return value*/