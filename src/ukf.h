#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "tools.h"

class UKF
{
public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage measurementPack);

  /**
  * Retrieves current estimate.
  **/
  TVector GetEstimate() const;

protected:
  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double deltaTime);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage measurementPack);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage measurementPack);

  /**
   * Updates the state and the state covariance matrix using given parameters.
   */
  void Update(const TVector &y, const TMatrix &S, const TMatrix &K);

public:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool mIsInitialized;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool mUseLaser;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool mUseRadar;

  ///* State dimension
  int mNX;

  ///* Augmented state dimension
  int mNAug;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  TVector mX;

  ///* state covariance matrix
  TMatrix mP;

  ///* predicted sigma points matrix
  TMatrix mXSigPred;

  ///* time when the state is true, in us
  TTimeStamp mPreviousTimestamp;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double mStdA;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double mStdYawdd;

  ///* Laser measurement noise standard deviation position1 in m
  double mStdLaspx;

  ///* Laser measurement noise standard deviation position2 in m
  double mStdLaspy;

  ///* Radar measurement noise standard deviation radius in m
  double mStdRadr;

  ///* Radar measurement noise standard deviation angle in rad
  double mStdRadphi;

  ///* Radar measurement noise standard deviation radius change in m/s
  double mStdRadrd ;

  ///* Measurement noise matrix laser
  TMatrix mRLaser;

  ///* Measurement noise matrix radar
  TMatrix mRRadar;

  ///* Weights of sigma points
  TVector mWeights;

  ///* Sigma point spreading parameter
  double mLambda;
};

#endif /* UKF_H */
