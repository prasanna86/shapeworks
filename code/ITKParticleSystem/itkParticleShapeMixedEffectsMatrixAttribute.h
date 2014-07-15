/* Class for Mixed-effects regression */

#ifndef __itkParticleShapeMixedEffectsMatrixAttribute_h
#define __itkParticleShapeMixedEffectsMatrixAttribute_h

#include "itkParticleShapeMatrixAttribute.h"
#include "vnl/vnl_vector.h"
#include "itkParticleSystem.h"
#include "vnl/vnl_trace.h"

namespace itk
{
 /** \class ParticleShapeMixedEffectsMatrixAttribute
  *
  *
  *
  */
 template <class T, unsigned int VDimension>
 class ITK_EXPORT ParticleShapeMixedEffectsMatrixAttribute
   : public ParticleShapeMatrixAttribute<T,VDimension>
 {
 public:
   /** Standard class typedefs */
   typedef T DataType;
   typedef ParticleShapeMixedEffectsMatrixAttribute Self;
   typedef ParticleShapeMatrixAttribute<T,VDimension> Superclass;
   typedef SmartPointer<Self>  Pointer;
   typedef SmartPointer<const Self>  ConstPointer;
   typedef WeakPointer<const Self>  ConstWeakPointer;

   /** Method for creation through the object factory. */
   itkNewMacro(Self);

   /** Run-time type information (and related methods). */
   itkTypeMacro(ParticleShapeMixedEffectsMatrixAttribute, ParticleShapeMatrixAttribute);

   // void UpdateMeanMatrix()
   // {
   //   // for each sample
   //   vnl_vector<double> tempvect;
   //   int indexSum = 0, group_indx = -1;
   //   tempvect.set_size(m_MeanMatrix.rows());
   //   tempvect.fill(0.0);
   //   for (int i = 0; i < m_MeanMatrix.cols(); i++)
   //   {
   //     if(i == indexSum)
   //     {
   // 	 group_indx = group_indx + 1;
   // 	 indexSum += m_TimeptsPerIndividual(group_indx);
   //     }

   //     tempvect = m_Intercept + m_Slope * m_Expl(i);
   //     tempvect = tempvect + m_InterceptRand.get_row(group_indx);
   //     tempvect = tempvect + m_SlopeRand.get_row(group_indx) * m_Expl(i);
   //     // compute the mean
   //     m_MeanMatrix.set_column(i, tempvect);
   //   }
   // }

   // ************* Covariate model ******************//
   void UpdateCovariateModelMeanMatrix()
   {
     vnl_vector<double> tempvect;
     vnl_matrix<double> D = m_FixedEffectsDesignMatrix;
     int indexSum = 0, group_indx = -1;
     tempvect.set_size(m_MeanMatrix.rows());
     tempvect.fill(0.0);
     for (int i = 0; i < m_MeanMatrix.cols(); i++)
     {
       if(i == indexSum)
       {
	 group_indx = group_indx + 1;
	 indexSum += m_TimeptsPerIndividual(group_indx);
       }

       // Adding the appropriate factors
       for(int j = 0; j < m_numFixedParams / 2; j++)
	 tempvect += m_Intercepts.get_row(j) + m_Slopes.get_row(j) * D(i, j);
       tempvect += m_InterceptRand.get_row(group_indx);
       tempvect += m_SlopeRand.get_row(group_indx) * m_Expl(i);
       // compute the mean
       m_MeanMatrix.set_column(i, tempvect);
     }
   }
  
   // ****************************************************** //

   // inline vnl_vector<double> ComputeMean(double k) const
   // {
   //   return m_Intercept + m_Slope * k;    
   // }

   // Covariate model version
   inline vnl_vector<double> ComputeCovariateModelMean(double k) const
   {
     //vnl_vector<double> tempvect;
     return m_Intercepts.get_row(0) + m_Slopes.get_row(0) * k;    
   }


   void ResizeParameters(unsigned int n)
   {
     vnl_vector<double> tmpA = m_Intercept; // copy existing  matrix
     vnl_vector<double> tmpB = m_Slope; // copy existing  matrix

     // Create new 
     m_Intercept.set_size(n);
     m_Slope.set_size(n);

     // Copy old data into new vector.
     for (unsigned int r = 0; r < tmpA.size(); r++)
     {
       m_Intercept(r) = tmpA(r);
       m_Slope(r) = tmpB(r);
     }
   }

   virtual void ResizeMeanMatrix(int rs, int cs)
   {
     vnl_matrix<T> tmp = m_MeanMatrix; // copy existing  matrix

     // Create new column (shape)
     m_MeanMatrix.set_size(rs, cs);
     m_MeanMatrix.fill(0.0);

     // Copy old data into new matrix.
     for (unsigned int c = 0; c < tmp.cols(); c++)
     {
       for (unsigned int r = 0; r < tmp.rows(); r++)
       {
	 m_MeanMatrix(r,c) = tmp(r,c);
       }
     } 
   }

   void ResizeExplanatory(unsigned int n)
   {
     if (n > m_Expl.size())
     {
       vnl_vector<double> tmp = m_Expl; // copy existing  matrix

       // Create new 
       m_Expl.set_size(n);
       m_Expl.fill(0.0);

       // Copy old data into new vector.
       for (unsigned int r = 0; r < tmp.size(); r++)
       {
	 m_Expl(r) = tmp(r);
       }
     }
   } 

   /** Callbacks that may be defined by a subclass.  If a subclass defines one
       of these callback methods, the corresponding flag in m_DefinedCallbacks
       should be set to true so that the ParticleSystem will know to register
       the appropriate event with this method. */
   virtual void DomainAddEventCallback(Object *, const EventObject &e)
   {
     const itk::ParticleDomainAddEvent &event
       = dynamic_cast<const itk::ParticleDomainAddEvent &>(e);
     unsigned int d = event.GetDomainIndex();

     if ( d % this->m_DomainsPerShape  == 0 )
     {
       this->ResizeMatrix(this->rows(), this->cols()+1);
       this->ResizeMeanMatrix(this->rows(), this->cols()+1);
       this->ResizeExplanatory(this->cols());
     }    
   }

   virtual void PositionAddEventCallback(Object *o, const EventObject &e) 
   {
     const itk::ParticlePositionAddEvent &event
       = dynamic_cast<const itk::ParticlePositionAddEvent &>(e);
     const itk::ParticleSystem<VDimension> *ps
       = dynamic_cast<const itk::ParticleSystem<VDimension> *>(o);
     const int d = event.GetDomainIndex();
     const unsigned int idx = event.GetPositionIndex();
     const typename itk::ParticleSystem<VDimension>::PointType pos
       = ps->GetTransformedPosition(idx, d);

     const unsigned int PointsPerDomain = ps ->GetNumberOfParticles(d);

     // Make sure we have enough rows.
     if ( (ps->GetNumberOfParticles(d) * VDimension * this->m_DomainsPerShape) > this->rows() )
     {
       this->ResizeParameters(PointsPerDomain * VDimension * this->m_DomainsPerShape);
       this->ResizeMatrix(PointsPerDomain * VDimension * this->m_DomainsPerShape,
			  this->cols());
       this->ResizeMeanMatrix(PointsPerDomain * VDimension * this->m_DomainsPerShape,
			      this->cols());
     }

     // CANNOT ADD POSITION INFO UNTIL ALL POINTS PER DOMAIN IS KNOWN
     // Add position info to the matrix
     unsigned int k = ((d % this->m_DomainsPerShape) * PointsPerDomain * VDimension)
       + (idx * VDimension);
     for (unsigned int i = 0; i < VDimension; i++)
     {
       this->operator()(i+k, d / this->m_DomainsPerShape) = pos[i];
     }

     //   std::cout << "Row " << k << " Col " << d / this->m_DomainsPerShape << " = " << pos << std::endl;
   }

   virtual void PositionSetEventCallback(Object *o, const EventObject &e) 
   {
     const itk::ParticlePositionSetEvent &event
       = dynamic_cast <const itk::ParticlePositionSetEvent &>(e);

     const itk::ParticleSystem<VDimension> *ps
       = dynamic_cast<const itk::ParticleSystem<VDimension> *>(o);
     const int d = event.GetDomainIndex();
     const unsigned int idx = event.GetPositionIndex();
     const typename itk::ParticleSystem<VDimension>::PointType pos = ps->GetTransformedPosition(idx, d);
     const unsigned int PointsPerDomain = ps->GetNumberOfParticles(d);

     // Modify matrix info
     //    unsigned int k = VDimension * idx;
     unsigned int k = ((d % this->m_DomainsPerShape) * PointsPerDomain * VDimension)
       + (idx * VDimension);    

     for (unsigned int i = 0; i < VDimension; i++)
     {
       this->operator()(i+k, d / this->m_DomainsPerShape) =
	 pos[i] - m_MeanMatrix(i+k, d/ this->m_DomainsPerShape);
     }
   }

   virtual void PositionRemoveEventCallback(Object *, const EventObject &) 
   {
     // NEED TO IMPLEMENT THIS
   }

   /** Set/Get the number of domains per shape.  This can only be safely done
       before shapes are initialized with points! */
   void SetDomainsPerShape(int i)
   { 
     this->m_DomainsPerShape = i;
   }

   int GetDomainsPerShape() const
   { 
     return this->m_DomainsPerShape;
   }

   void SetTimeptsPerIndividual(vnl_vector<int> & v)
   {
     (this->m_TimeptsPerIndividual).set_size(v.size());

     for (int i = 0; i < v.size(); i++)
     {
       //      std::cout << v[i] << std::endl;
       (this->m_TimeptsPerIndividual)(i) = v(i);
     }
   }

   vnl_vector<int> GetTimeptsPerIndividual() const
   {
     return this->m_TimeptsPerIndividual;
   }

   vnl_vector<int> &GetTimeptsPerIndividual()
   {
     return this->m_TimeptsPerIndividual;
   }

   void SetNumIndividuals(int i)
   {
     this->m_NumIndividuals = i;
   }

   int &GetNumIndividuals()
   {
     return this->m_NumIndividuals;
   }

   // next stage new
   void SetNumFixedParams(int i)
   { 
     this->m_NumFixedParams = i;
   }

   int &GetNumFixedParams()
   {
     return this->m_NumFixedParams;
   }

   void SetExplanatory(std::vector<double> v)
   {
     //    std::cout << "Setting expl " << std::endl;
     ResizeExplanatory(v.size());
     for (unsigned int i = 0; i < v.size(); i++)
     {
       //      std::cout << v[i] << std::endl;
       m_Expl[i] = v[i];
     }
   }

   void SetExplanatory(unsigned int i, double q)
   { 
     m_Expl[i] = q;
   }

   const double &GetExplanatory(unsigned int i) const
   { 
     return m_Expl[i];
   }

   double &GetExplanatory(unsigned int i)
   { 
     return m_Expl[i];
   }

   const vnl_vector<double> &GetSlope() const
   { 
     return m_Slope;
   }

   const vnl_vector<double> &GetIntercept() const
   {
     return m_Intercept;
   }

   const vnl_matrix<double> &GetSlopeRandom() const
   {
     return m_SlopeRand;
   }

   const vnl_matrix<double> &GetInterceptRandom() const
   {
     return m_InterceptRand;
   }

   void SetSlope(const std::vector<double> &v)
   {
     ResizeParameters(v.size());
     for (unsigned int i = 0; i < v.size(); i++)
     {
       m_Slope[i] = v[i];
     }    
   }

   void SetIntercept(const std::vector<double> &v)
   {
     ResizeParameters(v.size());
     for (unsigned int i = 0; i < v.size(); i++)
     {
       m_Intercept[i] = v[i];
     }
   }

   // covariate model versions
   void SetDesignMatrix(vnl_matrix<double> &v)
   {
     //    std::cout << "Setting design " << std::endl;
     (this->m_FixedEffectsDesignMatrix).set_size(v.rows(), v.cols());
     for (int i = 0; i < v.rows(); i++)
     {
       for (int j = 0; j < v.cols(); j++)
       {
	 (this->m_FixedEffectsDesignMatrix)(i, j) = v(i, j);
       }
     }
   }

   vnl_matrix<double> &GetDesignMatrix() const
   {
     return this->m_FixedEffectsDesignMatrix;
   }

   vnl_matrix<double> &GetDesignMatrix()
   {
     return this->m_FixedEffectsDesignMatrix;
   }

   const vnl_matrix<double> &GetSlopes() const
   { 
     return m_Slopes; 
   }

   const vnl_matrix<double> &GetIntercepts() const
   { 
     return m_Intercepts; 
   }

   // void SetSlopes(const vnl_matrix<double> &v)
   // {
   //   ResizeParameters(v.size());
   //   for (unsigned int i = 0; i < v.size(); i++)
   //   {
   // 	m_Slope[i] = v[i];
   //   }    
   // }
   
   // void SetIntercept(const std::vector<double> &v)
   // {
   //   ResizeParameters(v.size());
   //   for (unsigned int i = 0; i < v.size(); i++)
   //   {
   // 	m_Intercept[i] = v[i];
   //   }
   // }

   /*
   void EstimateParameters()
   {
     std::cout << "Estimating params" << std::endl;
     // std::cout << "Explanatory: " << m_Expl << std::endl;

     vnl_matrix<double> X = *this + m_MeanMatrix;

     // Number of samples
     int num_shapes = static_cast<double>(X.cols());

     int nr = X.rows(); //number of points*3
     std::cout << "num individuals: " << this->m_NumIndividuals << std::endl;

     // set the sizes of random slope and intercept matrix
     m_SlopeRand.set_size(m_NumIndividuals, nr); // num_groups X num_points*3
     m_InterceptRand.set_size(m_NumIndividuals, nr); // num_groups X num_points*3

     vnl_matrix<double> fixed; // slopes + intercepts for all points
     vnl_matrix<double> random; // slopes + intercepts for all groups, for all points
     fixed.set_size(2, nr);
     random.set_size(2, nr * m_NumIndividuals);

     vnl_matrix<double> Ds(2, 2); // covariance matrix of random parameters (2x2)
     Ds.set_identity();  // initialize to identity
     vnl_matrix<double> identity_2;
     identity_2.set_size(2, 2);
     identity_2.set_identity();

     double sigma2s = 1;  // variance of error
     vnl_matrix<double> * Ws = NULL, * Vs = NULL, * identity_n = NULL, * Xp = NULL;
     Ws = new vnl_matrix<double>[m_NumIndividuals];
     Vs = new vnl_matrix<double>[m_NumIndividuals];
     identity_n = new vnl_matrix<double>[m_NumIndividuals];
     Xp = new vnl_matrix<double>[m_NumIndividuals];

     vnl_vector<double> residual, y;
     y.set_size(m_TimeptsPerIndividual(0));
     y.fill(0.0);
     residual.set_size(m_TimeptsPerIndividual(0));
     residual.fill(0.0);

     for (int i = 0; i < m_NumIndividuals; i++)
     {
       Vs[i].set_size(m_TimeptsPerIndividual(i), m_TimeptsPerIndividual(i)); 
       Ws[i].set_size(m_TimeptsPerIndividual(i), m_TimeptsPerIndividual(i));

       identity_n[i].set_size(m_TimeptsPerIndividual(i), m_TimeptsPerIndividual(i));
       identity_n[i].set_identity();

       Xp[i].set_size(m_TimeptsPerIndividual(i), 2);
     }

     vnl_matrix<double> sum_mat1(2, 2, 0);
     vnl_vector<double> sum_mat2(2); sum_mat2.fill(0.0);
     double ecorr = 0.0;
     double tracevar = 0.0;
     vnl_matrix<double> bscorr(2, 2, 0.0); 
     vnl_matrix<double> bsvar(2, 2, 0.0);
     vnl_vector<double> tempvect; tempvect.set_size(2);
     unsigned int indexSum = 0;

     // std::cout << "all good here? [1] " << std::endl;

     for (int i = 0; i < nr; i++) //for all points (x,y,z coordinates)
     {
       sigma2s = 0.01;
       Ds.set_identity();
       Ds *= sigma2s;
       // std::cout << "all good here? [2] " << std::endl;
       for (int j = 0; j < 100; j++) //EM iterations
       {
	 sum_mat1.fill(0.0); sum_mat2.fill(0.0);
	 indexSum = 0;

	 // std::cout << "all good here? [3] " << std::endl;	  
	 for (int k = 0; k < m_NumIndividuals; k++)
	 {
	   y.clear();
	   y.set_size(m_TimeptsPerIndividual(k));
	   y.fill(0.0);

	   // std::cout << "all good here? [4] " << k << std::endl;	  
	   for (int l = 0; l < m_TimeptsPerIndividual(k); l++)
	   {
	     Xp[k](l, 0) = m_Expl(indexSum + l);
	     Xp[k](l, 1) = 1;
	     y(l) = X(i, indexSum + l);
	   }

	   indexSum += m_TimeptsPerIndividual(k);

	   Vs[k] = (identity_n[k] * sigma2s) + Xp[k] * Ds * vnl_transpose(Xp[k]);
	   Ws[k] = vnl_inverse(Vs[k]);
	   //	      << Vs[k].rows() << " " << Ws[k].rows() << " " 
	   //	      << Xp[k].cols() << " " << Vs[k].cols() << " " 
	   //	      << Ws[k].cols() << std::endl;	  
	   sum_mat1 = sum_mat1 + vnl_transpose(Xp[k]) * Ws[k] * Xp[k];
	   sum_mat2 = sum_mat2 + vnl_transpose(Xp[k]) * Ws[k] * y;
	 }

	 tempvect = vnl_inverse(sum_mat1) * sum_mat2;
	 fixed.set_column(i, tempvect);

	 indexSum = 0; // re-initializing
	 ecorr = 0.0; 
	 tracevar = 0.0;
	 bscorr.fill(0.0);
	 bsvar.fill(0.0);

	 for (int k = 0; k < m_NumIndividuals; k++)
	 {
	   residual.clear(); y.clear();
	   residual.set_size(m_TimeptsPerIndividual(k));
	   y.set_size(m_TimeptsPerIndividual(k));
	   residual.fill(0.0); y.fill(0.0);

	   for (int l = 0; l < m_TimeptsPerIndividual[k]; l++)
	   {
	     Xp[k](l,0) = m_Expl(indexSum + l);
	     Xp[k](l,1) = 1;
	     y(l) = X(i, indexSum + l);
	   }

	   indexSum += m_TimeptsPerIndividual(k);

	   tempvect = Ds * vnl_transpose(Xp[k]) * Ws[k] * (y - (Xp[k] * fixed.get_column(i)));
	   random.set_column(i * m_NumIndividuals + k, tempvect);
	   residual = y - (Xp[k] * fixed.get_column(i)) - (Xp[k] * random.get_column(i * m_NumIndividuals + k));
	   ecorr = ecorr + dot_product(residual, residual);
	   tracevar = tracevar + (m_TimeptsPerIndividual(k) - sigma2s * vnl_trace(Ws[k]));
	   bscorr = bscorr + outer_product(random.get_column(i * m_NumIndividuals + k), random.get_column(i * m_NumIndividuals + k));
	   bsvar = bsvar + (identity_2 - (vnl_transpose(Xp[k]) * Ws[k] * Xp[k] * Ds));
	 }

	 indexSum = 0; // re-initializing
	 sigma2s = (ecorr + sigma2s * tracevar) / (num_shapes);
	 Ds = (bscorr + Ds * bsvar) / m_NumIndividuals;
       } //endfor EM iterations
	 //printf ("point #%d\n", i);
     } //endfor all points on shape (x,y & z)

     m_Slope = fixed.get_row(0);
     m_Intercept = fixed.get_row(1);
     for (int i = 0; i < m_NumIndividuals; i++)
     {
       for (int j = 0; j < nr; j++) //for all points * 3
       {
	 m_SlopeRand(i, j) = random(0, j * m_NumIndividuals + i);
	 m_InterceptRand(i, j) = random(1, j * m_NumIndividuals + i);
       }
     }

     delete [] Vs;
     delete [] Ws;
     delete [] identity_n;
     delete [] Xp;

     // printf ("points:\n");
     // int indexSum = 0;
     // for (int k = 0; k < m_NumIndividuals; k++)
     // {
     //   for (int l = 0; l < m_TimeptsPerIndividual[k]; l++)
     //   {
     //    printf ("%g   %g\n", X(0, indexSum + l), m_Expl(indexSum + l));
     //   }
     //   indexSum += m_TimeptsPerIndividual[k];
     // }
     // printf ("fixed: slope %g, intercept %g", m_Slope(0), m_Intercept(0));
     // printf ("random: slopes %g %g, intercepts %g %g", m_SlopeRand(0,0), m_SlopeRand(1,0), m_InterceptRand(0,0), m_InterceptRand(1,0));
   }
   */

   // Estimating mixed-effects parameters for a covariate model (needs to be tested)
   void EstimateCovariateModelParameters()
   {
     std::cout << "Estimating params" << std::endl;
     //std::cout << "Explanatory: " << m_Expl << std::endl;

     vnl_matrix<double> X = *this + m_MeanMatrix;
     vnl_matrix<double> fD = m_FixedEffectsDesignMatrix;

     // Number of samples
     int num_shapes = static_cast<double>(X.cols());
     int nr = X.rows(); //number of points*3
     std::cout << "num individuals: " << this->m_NumIndividuals << std::endl;
     std::cout << "num fixed effects params: " << this->m_NumFixedParams << std::endl;

     // set the sizes of random slope and intercept matrix
     m_SlopeRand.set_size(m_NumIndividuals, nr); // num_groups X num_points*3
     m_InterceptRand.set_size(m_NumIndividuals, nr); // num_groups X num_points*3

     vnl_matrix<double> fixed; // slopes + intercepts for all points
     vnl_matrix<double> random; // slopes + intercepts for all groups, for all points
     fixed.set_size(m_NumFixedParams, nr);
     random.set_size(2, nr * m_NumIndividuals);

     vnl_matrix<double> Ds(2, 2); // covariance matrix of random parameters (2x2)
     Ds.set_identity();  // initialize to identity
     vnl_matrix<double> identity_2;
     identity_2.set_size(2, 2);
     identity_2.set_identity();

     double sigma2s = 1;  // variance of error
     vnl_matrix<double> * Ws = NULL, * Vs = NULL, * identity_n = NULL;
     vnl_matrix<double> * Xp = NULL, * Zp = NULL;
     Ws = new vnl_matrix<double>[m_NumIndividuals];
     Vs = new vnl_matrix<double>[m_NumIndividuals];
     identity_n = new vnl_matrix<double>[m_NumIndividuals];
     Xp = new vnl_matrix<double>[m_NumIndividuals];
     Zp = new vnl_matrix<double>[m_NumIndividuals];

     vnl_vector<double> residual, y;
     y.set_size(m_TimeptsPerIndividual(0));
     y.fill(0.0);
     residual.set_size(m_TimeptsPerIndividual(0));
     residual.fill(0.0);

     for (int i = 0; i < m_NumIndividuals; i++)
     {
       Vs[i].set_size(m_TimeptsPerIndividual(i), m_TimeptsPerIndividual(i)); 
       Ws[i].set_size(m_TimeptsPerIndividual(i), m_TimeptsPerIndividual(i));

       identity_n[i].set_size(m_TimeptsPerIndividual(i), m_TimeptsPerIndividual(i));
       identity_n[i].set_identity();

       Xp[i].set_size(m_TimeptsPerIndividual(i), m_NumFixedParams);
       Zp[i].set_size(m_TimeptsPerIndividual(i), 2);
     }

     vnl_matrix<double> sum_mat1(m_NumFixedParams, m_NumFixedParams, 0);
     vnl_vector<double> sum_mat2(m_NumFixedParams); sum_mat2.fill(0.0);
     double ecorr = 0.0;
     double tracevar = 0.0;
     vnl_matrix<double> bscorr(2, 2, 0.0); 
     vnl_matrix<double> bsvar(2, 2, 0.0);
     vnl_vector<double> tempvect; tempvect.set_size(m_NumFixedParams);
     unsigned int indexSum = 0;
      
     //std::cout << "all good here? [1] " << std::endl;
     
     for (int i = 0; i < nr; i++) //for all points (x,y,z coordinates)
     {
       sigma2s = 0.01;
       Ds.set_identity();
       Ds *= sigma2s;
       //std::cout << "all good here? [2] " << std::endl;
       for (int j = 0; j < 50; j++) //EM iterations
       {
	 sum_mat1.fill(0.0); sum_mat2.fill(0.0);
	 indexSum = 0;
	 
	 //std::cout << "all good here? [3] " << std::endl;	  
	 for (int k = 0; k < m_NumIndividuals; k++)
    	 {
	   y.clear();
	   y.set_size(m_TimeptsPerIndividual(k));
	   y.fill(0.0);
	   
	   //std::cout << "all good here? [4] " << k << std::endl;	  
	   for (int l = 0; l < m_TimeptsPerIndividual(k); l++)
    	   {
	     Xp[k](l, 0) = 1;
	     Zp[k](l, 0) = 1;

	     for(int m = 1; m < m_NumFixedParams; m++)
	     {
	       Xp[k](l, m) = fD(indexSum + l, m);
	     }

	     Zp[k](l, 1) = m_Expl(indexSum + l);
	     y(l) = X(i, indexSum + l);
	   }

	   indexSum += m_TimeptsPerIndividual(k);

	   Vs[k] = (identity_n[k] * sigma2s) + Zp[k] * Ds * vnl_transpose(Zp[k]);
	   Ws[k] = vnl_inverse(Vs[k]);
	   sum_mat1 = sum_mat1 + vnl_transpose(Xp[k]) * Ws[k] * Xp[k];
	   sum_mat2 = sum_mat2 + vnl_transpose(Xp[k]) * Ws[k] * y;
	 }

	 tempvect = vnl_inverse(sum_mat1) * sum_mat2;
	 fixed.set_column(i, tempvect);
	 
	 indexSum = 0; // re-initializing
	 ecorr = 0.0; 
	 tracevar = 0.0;
	 bscorr.fill(0.0);
	 bsvar.fill(0.0);

	 for (int k = 0; k < m_NumIndividuals; k++)
	 {
	   residual.clear(); y.clear();
	   residual.set_size(m_TimeptsPerIndividual(k));
	   y.set_size(m_TimeptsPerIndividual(k));
	   residual.fill(0.0); y.fill(0.0);

	   for (int l = 0; l < m_TimeptsPerIndividual[k]; l++)
	   {
	     Xp[k](l, 0) = 1;
	     Zp[k](l, 0) = 1;
	     
	     for(int m = 1; m < m_NumFixedParams; m++)
	     {
	       Xp[k](l, m) = fD(indexSum + l, m);
	     }

	     Zp[k](l, 1) = m_Expl(indexSum + l);
	     y(l) = X(i, indexSum + l);
	   }

	   indexSum += m_TimeptsPerIndividual(k);

	   tempvect = Ds * vnl_transpose(Zp[k]) * Ws[k] * (y - (Xp[k] * fixed.get_column(i)));
	   random.set_column(i * m_NumIndividuals + k, tempvect);
	   residual = y - (Xp[k] * fixed.get_column(i)) - (Zp[k] * random.get_column(i * m_NumIndividuals + k));
	   ecorr = ecorr + dot_product(residual, residual);
	   tracevar = tracevar + (m_TimeptsPerIndividual(k) - sigma2s * vnl_trace(Ws[k]));
	   bscorr = bscorr + outer_product(random.get_column(i * m_NumIndividuals + k), random.get_column(i * m_NumIndividuals + k));
	   bsvar = bsvar + (identity_2 - (vnl_transpose(Xp[k]) * Ws[k] * Zp[k] * Ds));
	 }

	 indexSum = 0; // re-initializing
	 sigma2s = (ecorr + sigma2s * tracevar) / (num_shapes);
	 Ds = (bscorr + Ds * bsvar) / m_NumIndividuals;
       } //end for EM iterations
     } // end for all points on shape (x,y & z)

     for(int i = 0; i < m_NumFixedParams / 2; i++)
     {
       m_Intercepts.set_row(i, fixed.get_row(i));
       m_Slopes.set_row(i, fixed.get_row(i + (m_NumFixedParams / 2)));
     }

     for (int i = 0; i < m_NumIndividuals; i++)
     {
       for (int j = 0; j < nr; j++) //for all points * 3
       {
	 m_InterceptRand(i, j) = random(0, j * m_NumIndividuals + i);
	 m_SlopeRand(i, j) = random(1, j * m_NumIndividuals + i);
       }
     }

     delete [] Vs;
     delete [] Ws;
     delete [] identity_n;
     delete [] Xp;
     delete [] Zp;
   }
   /*
   // Initialize variables 
   void Initialize()
   {
     m_Intercept.fill(0.0);
     m_Slope.fill(0.0);
     m_MeanMatrix.fill(0.0);

     m_SlopeRand.fill(0.0);
     m_InterceptRand.fill(0.0);    
   }

   virtual void BeforeIteration()
   {
     m_UpdateCounter ++;
     if (m_UpdateCounter >= m_RegressionInterval)
     {
       m_UpdateCounter = 0;
       this->EstimateParameters();
       this->UpdateMeanMatrix();
     }
   }
   */
   // Covariate model version
   void InitializeCovariateModel()
   {
     m_Intercepts.fill(0.0);
     m_Slopes.fill(0.0);
     
     m_MeanMatrix.fill(0.0);
     m_FixedEffectsDesignMatrix.fill(0.0);
     
     m_SlopeRand.fill(0.0);
     m_InterceptRand.fill(0.0);    
   }
   // Covariate model version
   virtual void BeforeIterationCovariateModel()
   {
     m_UpdateCounter ++;
     if (m_UpdateCounter >= m_RegressionInterval)
     {
       m_UpdateCounter = 0;
       this->EstimateCovariateModelParameters();
       this->UpdateCovariateModelMeanMatrix();
     }
   }

   void SetRegressionInterval( int i)
   {
     m_RegressionInterval = i;
   }

   int GetRegressionInterval() const
   {
     return m_RegressionInterval;
   }
   
 protected:
   ParticleShapeMixedEffectsMatrixAttribute() 
   {
     this->m_DefinedCallbacks.DomainAddEvent = true;
     this->m_DefinedCallbacks.PositionAddEvent = true;
     this->m_DefinedCallbacks.PositionSetEvent = true;
     this->m_DefinedCallbacks.PositionRemoveEvent = true;
     m_UpdateCounter = 0;
     m_RegressionInterval = 1;
     m_NumIndividuals = 1;
     m_NumFixedParams = 2;
     m_TimeptsPerIndividual.set_size(m_NumIndividuals);
     for(int i = 0; i < m_NumIndividuals; i++)
       m_TimeptsPerIndividual(i) = 3;
   }

   virtual ~ParticleShapeMixedEffectsMatrixAttribute() {};

   void PrintSelf(std::ostream& os, Indent indent) const
   {
     Superclass::PrintSelf(os,indent);
   }

 private:
   ParticleShapeMixedEffectsMatrixAttribute(const Self&); //purposely not implemented
   void operator=(const Self&); //purposely not implemented

   int m_UpdateCounter;
   int m_RegressionInterval;

   // Parameters for the linear model (I could just make these matrices!)
   vnl_vector<double> m_Intercept;
   vnl_vector<double> m_Slope;

   // I could make the above matrices, but till then:
    vnl_matrix<double> m_Intercepts;
    vnl_matrix<double> m_Slopes;

   // The explanatory variable value for each sample (matrix column)
   vnl_vector<double> m_Expl;

   // A design matrix to store the explanatory variables: Need to absorb this into m_Expl
    vnl_matrix<double> m_FixedEffectsDesignMatrix;

   // A matrix to store the mean estimated for each explanatory variable (each sample)
   vnl_matrix<double> m_MeanMatrix;

   vnl_matrix<double> m_InterceptRand; //added: AK , random intercepts for each group
   vnl_matrix<double> m_SlopeRand; //added: AK , random slopes for each group

   int m_NumIndividuals;
   int m_NumFixedParams;
   // timepoints per individual
   vnl_vector<int> m_TimeptsPerIndividual;
 };

} // end namespace

#endif
