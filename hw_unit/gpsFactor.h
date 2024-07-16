#pragma once
#include "common.h"
// #define VS

const float gpsFactorCov[3] = {1, 1, 1};

class gpsFactor
{
public:
	gpsFactor(int index, float pose[3], float gpsMeasurement[3], const float cov[3])
	{
		this->index = index;
		for (int i = 0; i < 3; i++)
		{
			this->pose[i] = pose[i];
			this->gpsMeasurement[i] = gpsMeasurement[i];
			this->covariance[i][i] = cov[i];
#ifdef VS
			this->covariance_sqrt_inv[i][i] = 1 / std::sqrt(cov[i]);
#else
			this->covariance_sqrt_inv[i][i] = 1 / hls::sqrt(cov[i]);
#endif
		}
	}

	void computeResidual()
	{
		float residual[3][1];
		for (int i = 0; i < 3; i++)
		{
			residual[i][0] = gpsMeasurement[i] - pose[i];
		}

		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 3, 3, 3, 1, 3, 1, float, float>(covariance_sqrt_inv, residual, residual_l2);
	}

	void computeJacobian()
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				jacobian_l2[i][j] = -covariance_sqrt_inv[i][j];
			}
		}
	}

	void computeHg()
	{
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 3, 3, 3, 3, 3, 3, float, float>(jacobian_l2, jacobian_l2, H);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 3, 3, 3, 1, 3, 1, float, float>(jacobian_l2, residual_l2, g);
		for (int i = 0; i < 3; i++)
		{
			g[i][0] = -g[i][0];
		}
	}

	float computeCost()
	{
		float cost = 0;
		for (int i = 0; i < 3; i++)
		{
			cost += residual_l2[i][0] * residual_l2[i][0];
		}
		return cost;
	}

	int index;
	float pose[3];
	float gpsMeasurement[3];
	float residual_l2[3][1];
	float jacobian_l2[3][3];
	float H[3][3];
	float g[3][1];
	float covariance[3][3] = {0};
	float covariance_sqrt_inv[3][3] = {0};
};
