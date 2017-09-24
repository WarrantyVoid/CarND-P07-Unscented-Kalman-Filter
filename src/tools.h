#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

typedef unsigned long long TTimeStamp;
typedef Eigen::MatrixXd TMatrix;
typedef Eigen::VectorXd TVector;
typedef std::vector<TVector> TVectorList;

class Tools
{
public:
  /**
  * A helper method to calculate root mean squared error.
  **/
  TVector CalculateRMSE(const TVectorList &estimations, const TVectorList &groundTruths) const;

  /**
  * A helper method to calculate delta between two angles.
  **/
  float CalculateAngleDelta(float a1, float a2) const;

  /**
  * A helper method to normalize angle.
  **/
  float NormalizeAngle(float a) const;

  /**
  * A helper method to calculate Jacobians for radar measurement
  * translation function.
  **/
  TMatrix CalculateRadarJacobian(const TVector &xState) const;

  /**
  * A helper method to calculate whether a float represents zero.
  **/
  bool isZero(float f) const;
};

/**
* Retrieves reference to static tool object
**/
static const Tools &GetTools()
{
  static Tools sTools;
  return sTools;
}

#endif /* TOOLS_H_ */
