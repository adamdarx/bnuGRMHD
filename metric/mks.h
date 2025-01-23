#pragma once
#include "../utils.h"

void init_metric()
{
	metricFunc.setConstant(ZERO_COMPONENT);
	metricDiff.setConstant(ZERO_COMPONENT);
	metricFunc(0, 0) = [](double x1, double x2, double x3) {return -1. + 2. * x1 / (pow(x1, 2.) + pow(a * cos(x2), 2.)); };
	metricFunc(0, 1) = [](double x1, double x2, double x3) {return (2. * x1 / (pow(x1, 2.) + pow(a * cos(x2), 2.))) * x1; };
	metricFunc(0, 3) = [](double x1, double x2, double x3) {return (-2. * a * x1 * pow(sin(x2), 2.) / (pow(x1, 2.) + pow(a * cos(x2), 2.))) * (1. + h * cos(2. * x2)); };
	metricFunc(1, 0) = metricFunc(0, 1);
	metricFunc(1, 1) = [](double x1, double x2, double x3) {return (1. + 2. * x1 / (pow(x1, 2.) + pow(a * cos(x2), 2.))) * x1 * x1; };
	metricFunc(1, 3) = [](double x1, double x2, double x3) {return (-a * pow(sin(x2), 2.) * (1. + 2. * x1 / (pow(x1, 2.) + pow(a * cos(x2), 2.)))) * x1 * (1. + h * cos(2. * x2)); };
	metricFunc(2, 2) = [](double x1, double x2, double x3) {return pow(x1, 2.) + pow(a * cos(x2), 2.); };
	metricFunc(3, 0) = metricFunc(0, 3);
	metricFunc(3, 1) = metricFunc(1, 3);
	metricFunc(3, 3) = [](double x1, double x2, double x3) {return (pow(sin(x2), 2.) * ((pow(x1, 2.) + pow(a * cos(x2), 2.)) + pow(a * sin(x2), 2.) * (1. + 2. * x1 / (pow(x1, 2.) + pow(a * cos(x2), 2.))))) * pow(1. + h * cos(2. * x2), 2.); };

	metricDiff(0, 0, 1) = [](double x1, double x2, double x3) {return (-4 * x1) / pow(pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2), 2); };
	metricDiff(0, 0, 2) = [](double x1, double x2, double x3) {return (4 * pow(a, 2) * cos(x2) * sin(x2)) / pow(pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2), 2); };
	metricDiff(0, 1, 1) = [](double x1, double x2, double x3) {return (-4 * pow(x1, 3)) / pow(pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2), 2) + (4 * x1) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2)); };
	metricDiff(0, 1, 2) = [](double x1, double x2, double x3) {return (4 * pow(a, 2) * pow(x1, 2) * cos(x2) * sin(x2)) / pow(pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2), 2); };
	metricDiff(0, 3, 1) = [](double x1, double x2, double x3) {return (4 * a * pow(x1, 2) * (1 + h * pow(cos(x2), 2)) * pow(sin(x2), 2)) / pow(pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2), 2) - (2 * a * (1 + h * pow(cos(x2), 2)) * pow(sin(x2), 2)) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2)); };
	metricDiff(0, 3, 2) = [](double x1, double x2, double x3) {return (-4 * a * x1 * cos(x2) * (1 + h * pow(cos(x2), 2)) * sin(x2)) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2)) + 	(4 * a * h * x1 * cos(x2) * pow(sin(x2), 3)) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2)) - (4 * pow(a, 3) * x1 * cos(x2) * (1 + h * pow(cos(x2), 2)) * pow(sin(x2), 3)) / pow(pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2), 2); };
	metricDiff(1, 0, 1) = metricDiff(0, 1, 1);
	metricDiff(1, 0, 2) = metricDiff(0, 1, 2);
	metricDiff(1, 1, 1) = [](double x1, double x2, double x3) {return (-4 * pow(x1, 4)) / pow(pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2), 2) + (6 * pow(x1, 2)) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2)); };
	metricDiff(1, 1, 2) = [](double x1, double x2, double x3) {return (4 * pow(a, 2) * pow(x1, 3) * cos(x2) * sin(x2)) / pow(pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2), 2); };
	metricDiff(1, 3, 1) = [](double x1, double x2, double x3) {return -(a * x1 * (1 + h * pow(cos(x2), 2)) * ((-4 * pow(x1, 2)) / pow(pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2), 2) + 2 / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2))) * pow(sin(x2), 2)) - a * (1 + h * pow(cos(x2), 2)) * (1 + (2 * x1) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2))) * pow(sin(x2), 2); };
	metricDiff(1, 3, 2) = [](double x1, double x2, double x3) {return -2 * a * x1 * cos(x2) * (1 + h * pow(cos(x2), 2)) * (1 + (2 * x1) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2))) * sin(x2) - (4 * pow(a, 3) * pow(x1, 2) * cos(x2) * (1 + h * pow(cos(x2), 2)) * pow(sin(x2), 3)) / pow(pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2), 2) + 2 * a * h * x1 * cos(x2) * (1 + (2 * x1) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2))) * pow(sin(x2), 3); };
	metricDiff(2, 2, 1) = [](double x1, double x2, double x3) {return 2 * x1; };
	metricDiff(2, 2, 2) = [](double x1, double x2, double x3) {return -2 * pow(a, 2) * cos(x2) * sin(x2); };
	metricDiff(3, 0, 1) = metricDiff(0, 3, 1);
	metricDiff(3, 0, 2) = metricDiff(0, 3, 2);
	metricDiff(3, 1, 1) = metricDiff(1, 3, 1);
	metricDiff(3, 1, 2) = metricDiff(1, 3, 2);
	metricDiff(3, 3, 1) = [](double x1, double x2, double x3) {return a * (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2)) * (1 + h * pow(cos(x2), 2)) * ((-4 * pow(x1, 2)) / pow(pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2), 2) + 2 / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2))) * pow(sin(x2), 4) + 2 * a * x1 * (1 + h * pow(cos(x2), 2)) * (1 + (2 * x1) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2))) * pow(sin(x2), 4); };
	metricDiff(3, 3, 2) = [](double x1, double x2, double x3) {return 4 * a * cos(x2) * (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2)) * (1 + h * pow(cos(x2), 2)) * (1 + (2 * x1) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2))) * pow(sin(x2), 3) + (4 * pow(a, 3) * x1 * cos(x2) * (1 + h * pow(cos(x2), 2)) * pow(sin(x2), 5)) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2)) - 2 * a * h * cos(x2) * (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2)) * (1 + (2 * x1) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2))) * pow(sin(x2), 5) - 2 * pow(a, 3) * cos(x2) * (1 + h * pow(cos(x2), 2)) * (1 + (2 * x1) / (pow(x1, 2) + pow(a, 2) * pow(cos(x2), 2))) * pow(sin(x2), 5); };
}