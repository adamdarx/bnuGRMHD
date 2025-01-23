#pragma once
#include <cmath>
#include <functional>
#include <Eigen/Dense>


typedef std::function<double(double, double, double)> MetricComponent;

class Metric final
{
public:
	Eigen::Matrix4d m;
	Metric();
	Metric(Eigen::Matrix4d metric);
	double alpha();
	Eigen::Vector3d betaVec();
	Eigen::Matrix3d gamma();
	~Metric();
};
