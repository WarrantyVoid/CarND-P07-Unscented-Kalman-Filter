#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
  : mIsInitialized(false)
  , mUseLaser(true)
  , mUseRadar(true)
  , mNX(5)
  , mNAug(7)
  , mX(mNX)
  , mP(mNX, mNX)
  , mXSigPred(mNX, 2 * mNAug + 1)
  , mPreviousTimestamp(0ull)
  , mStdA(0.2)
  , mStdYawdd(0.2)
  , mStdLaspx(0.15)
  , mStdLaspy(0.15)
  , mStdRadr(0.3)
  , mStdRadphi(0.0175)
  , mStdRadrd(0.1)
  , mRLaser(2, 2)
  , mRRadar(3, 3)
  , mWeights(2 * mNAug + 1)
  , mLambda(3 - mNAug)
{
  // init laser noise
  mRLaser << mStdLaspx * mStdLaspx, 0,
             0,                     mStdLaspy * mStdLaspy;

  // init radar noise
  mRRadar << mStdRadr * mStdRadr, 0                      , 0,
             0                  , mStdRadphi * mStdRadphi, 0,
             0                  , 0                      , mStdRadrd * mStdRadr;

  // init sigma weights
  mWeights(0) = mLambda / (mLambda + mNAug);
  for (int n = 1; n < 2 * mNAug + 1; ++n)
  {
    mWeights(n) = 0.5 / (mLambda + mNAug);
  }
}

UKF::~UKF()
{

}

void UKF::ProcessMeasurement(MeasurementPackage measurementPack)
{
  if (!mIsInitialized)
  {
    mP << 1, 0, 0   , 0,    0,
          0, 1, 0   , 0,    0,
          0, 0, 1000, 0,    0,
          0, 0, 0   , 1000, 0,
          0, 0, 0   , 0,    1000;
    switch (measurementPack.sensorType)
    {
      case MeasurementPackage::RADAR:
      {
        mX << measurementPack.rawMeasurements(0) * cos(measurementPack.rawMeasurements(1)),
              measurementPack.rawMeasurements(0) * sin(measurementPack.rawMeasurements(1)),
              0,
              0,
              0;
      }
      break;
    case MeasurementPackage::LASER:
      {
        mX << measurementPack.rawMeasurements(0),
              measurementPack.rawMeasurements(1),
              0,
              0,
              0;
      }
      break;
    }
    mIsInitialized = true;
  }
  else
  {
    // predict state and new uncertainty given time delta and acceleration noise
    float dt = (measurementPack.timestamp - mPreviousTimestamp) / 1000000.0f;
    Prediction(dt);

    switch (measurementPack.sensorType)
    {
    case MeasurementPackage::RADAR:
      UpdateRadar(measurementPack);
      break;
    case MeasurementPackage::LASER:
      UpdateLidar(measurementPack);
      break;
    }
  }

  // store time stamp for next cycle
  mPreviousTimestamp = measurementPack.timestamp;
}

const TVector UKF::GetEstimate() const
{
  TVector x(4);
  x << mX(0),
       mX(1),
       mX(2) * cos(mX(3)),
       mX(2) * sin(mX(3));
  return x;
}

void UKF::Prediction(double deltaTime)
{
  // create augmented mean vector
  TVector xAug(mNAug);
  xAug << mX, 0, 0;

  // create augmented state covariance
  TMatrix pAug(mNAug, mNAug);
  pAug.setZero();
  pAug.topLeftCorner(5, 5) = mP;
  pAug(5, 5) = mStdA * mStdA;
  pAug(6, 6) = mStdYawdd * mStdYawdd;

  // create augmented sigma point matrix
  MatrixXd xSigAug(mNAug, 2 * mNAug + 1);
  TMatrix A(pAug.llt().matrixL());
  xSigAug.col(0) = xAug;
  TMatrix An(sqrt(mLambda + mNAug) * A);
  for (int n = 0; n < mNAug; ++n)
  {
    xSigAug.col(n + 1) = xAug + An.col(n);
    xSigAug.col(n + 1 + mNAug) = xAug - An.col(n);
  }

  // predict sigma points
  for (int n = 0; n < 2 * mNAug + 1; ++n)
  {
      TVector x = xSigAug.col(n);
      TVector xPred(mNX);
      if (GetTools().isZero(x(4)))
      {
        // straight case
        xPred << x(2) * cos(x(3)) * deltaTime,
                 x(2) * sin(x(3)) * deltaTime,
                 0,
                 x(4) * deltaTime,
                 0;
      }
      else
      {
        xPred << x(2) / x(4) * (sin(x(3) + x(4) * deltaTime) - sin(x(3))),
                 x(2) / x(4)* (-cos(x(3) + x(4) * deltaTime) + cos(x(3))),
                 0,
                 x(4) * deltaTime,
                 0;
      }
      TVector xPredNu(mNX);
      xPredNu << 0.5 * deltaTime * deltaTime * cos(x(3)) * x(5),
                 0.5 * deltaTime * deltaTime * sin(x(3)) * x(5),
                 deltaTime * x(5),
                 0.5 * deltaTime * deltaTime * x(6),
                 deltaTime * x(6);

      mXSigPred.col(n) = x.head(5) + xPred + xPredNu;
  }

  // update new state from predicted sigma points
  mX = mXSigPred.col(0) * mWeights(0);
  for (int n = 1; n < 2 * mNAug + 1; ++n)
  {
    mX += mXSigPred.col(n) * mWeights(n);
  }

  // update state covariance matrix from predicted sigma points
  mP.setZero();
  for (int n = 0; n < 2 * mNAug + 1; ++n)
  {
      TVector xd = mXSigPred.col(n) - mX;
      xd(3) = GetTools().NormalizeAngle(xd(3));
      mP += mWeights(n) * xd * xd.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage measurementPack)
{
  // map sigma points to measurement space
  TMatrix zSig = MatrixXd(2, 2 * mNAug + 1);
  for (int n = 0; n < 2 * mNAug + 1; ++n)
  {
    TVector z = mXSigPred.col(n);
    zSig.col(n) << z(0) , z(1);
  }

  // calculate mean predicted measurement
  TVector zPred = zSig.col(0) * mWeights(0);
  for (int n = 1; n < 2 * mNAug + 1; ++n)
  {
    zPred += zSig.col(n) * mWeights(n);
  }

  TMatrix S(2, 2);
  TMatrix Tc(mNX, 2);
  S.setZero();
  Tc.setZero();
  for (int n = 0; n < 2 * mNAug + 1; ++n)
  {
    // prediction residual
    VectorXd zd = zSig.col(n) - zPred;

    // calculate measurement covariance
    S += mWeights(n) * zd * zd.transpose();

    // calculate cross correlation matrix
    VectorXd xd = mXSigPred.col(n) - mX;
    xd(3) = GetTools().NormalizeAngle(xd(3));
    Tc += mWeights(n) * xd * zd.transpose();
  }

  // add measurement noise
  S += mRLaser;

  // update
  MatrixXd K = Tc * S.inverse();
  VectorXd y = measurementPack.rawMeasurements - zPred;
  Update(y, S, K);
}

void UKF::UpdateRadar(MeasurementPackage measurementPack)
{
  // map sigma points to measurement space
  TMatrix zSig = MatrixXd(3, 2 * mNAug + 1);
  for (int n = 0; n < 2 * mNAug + 1; ++n)
  {
    TVector z = mXSigPred.col(n);
    float roh = sqrt(z(0) * z(0) + z(1) * z(1)) + 0.001;
    zSig.col(n) << roh,
                   atan2(z(1),z(0)),
                   (z(0) * cos(z(3)) * z(2) + z(1) * sin(z(3)) * z(2) ) / roh;
  }

  // calculate mean predicted measurement
  TVector zPred = zSig.col(0) * mWeights(0);
  for (int n = 1; n < 2 * mNAug + 1; ++n)
  {
    zPred += zSig.col(n) * mWeights(n);
  }

  TMatrix S(3, 3);
  TMatrix Tc(mNX, 3);
  S.setZero();
  Tc.setZero();
  for (int n = 0; n < 2 * mNAug + 1; ++n)
  {
    // prediction residual
    VectorXd zd = zSig.col(n) - zPred;
    zd(1) = GetTools().NormalizeAngle(zd(1));

    // calculate measurement covariance
    S += mWeights(n) * zd * zd.transpose();

    // calculate cross correlation matrix
    VectorXd xd = mXSigPred.col(n) - mX;
    xd(3) = GetTools().NormalizeAngle(xd(3));
    Tc += mWeights(n) * xd * zd.transpose();
  }

  // add measurement noise
  S += mRRadar;

  // update
  MatrixXd K = Tc * S.inverse();
  VectorXd y = measurementPack.rawMeasurements - zPred;
  y(1) = GetTools().NormalizeAngle(y(1));
  Update(y, S, K);
}

void UKF::Update(const TVector &y, const TMatrix &S, const TMatrix &K)
{
  //update state mean and covariance matrix
  mX += K * y;
  mP += K * S * K.transpose();
}
