#include "Metric.h"

Metric::Metric()
{
	m << 0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0;
}

Metric::Metric(Eigen::Matrix4d metric) : m(metric)
{
	assert(metric.rows() == 4 && metric.rows() == 4);
}

double Metric::alpha()
{
	return sqrt(this->betaVec().transpose() * this->gamma().inverse() * this->betaVec() - m(0, 0));
}

Eigen::Vector3d Metric::betaVec()
{
	return m(0, Eigen::seqN(1, 3));
}

Eigen::Matrix3d Metric::gamma()
{
	return m.block(1,1,3,3);
}

Metric::~Metric()
{
}
