#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "tools.h"

class MeasurementPackage
{
public:
  enum SensorType
  {
    LASER,
    RADAR
  } sensorType;
  TVector rawMeasurements;
  TTimeStamp timestamp;

protected:
  /**
  * Deserializes package from string stream.
  **/
  friend std::istringstream &operator>>(std::istringstream &in, MeasurementPackage &mp)
  {
    std::string sensor_type;
    in >> sensor_type;
    if (sensor_type.compare("L") == 0)
    {
      mp.sensorType = MeasurementPackage::LASER;
      mp.rawMeasurements = TVector(2);
      in >> mp.rawMeasurements(0);
      in >> mp.rawMeasurements(1);
      in >> mp.timestamp;
    }
    else if (sensor_type.compare("R") == 0)
    {
      mp.sensorType = MeasurementPackage::RADAR;
      mp.rawMeasurements = TVector(3);
      in >> mp.rawMeasurements(0);
      in >> mp.rawMeasurements(1);
      in >> mp.rawMeasurements(2);
      in >> mp.timestamp;
    }
  }
};

#endif /* MEASUREMENT_PACKAGE_H_ */
