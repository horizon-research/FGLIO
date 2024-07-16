#pragma once
#include "common.h"

const float betweenFactorCov[6] = {1, 1, 1, 1, 1, 1};

class BetweenFactor
{
public:
	BetweenFactor(int ind1, int ind2, float p1[7], float p2[7], float measurement[7], const float cov[6])
	{
		index_1 = ind1;
		index_2 = ind2;
		for (int i = 0; i < 7; i++)
		{
			pose_1[i] = p1[i];
			pose_2[i] = p2[i];
			this->measurement[i] = measurement[i];
		}
		for (int i = 0; i < 6; i++)
		{
			covariance[i][i] = cov[i];
#ifdef VS
			covariance_sqrt_inv[i][i] = 1 / std::sqrt(cov[i]);

#else
			covariance_sqrt_inv[i][i] = 1 / hls::sqrt(cov[i]);

#endif // VS
		}
	}

	void computeResidual()
	{
		float translation_1[3] = {pose_1[0], pose_1[1], pose_1[2]};
		float translation_2[3] = {pose_2[0], pose_2[1], pose_2[2]};
		float translation_between[3] = {measurement[0], measurement[1], measurement[2]};
		// float q_1[4] = {pose_1[3], pose_1[4], pose_1[5], pose_1[6]};
		// float q_2[4] = {pose_2[3], pose_2[4], pose_2[5], pose_2[6]};
		// float q_between[4] = {measurement[3], measurement[4], measurement[5], measurement[6]};
		Quaterniond q_1(pose_1[3], pose_1[4], pose_1[5], pose_1[6]);
		Quaterniond q_2(pose_2[3], pose_2[4], pose_2[5], pose_2[6]);
		Quaterniond q_between(measurement[3], measurement[4], measurement[5], measurement[6]);

		SE3 last_pose(translation_1, q_1);
		SE3 current_pose(translation_2, q_2);
		SE3 between_pose(translation_between, q_between);

		float last_T_inv[4][4];
		// int inverse_ok;
		pose_inv(last_pose.T, last_T_inv);
		// #ifndef __SYNTHESIS__
		// 		assert(inverse_ok == 0);
		// #endif // !__SYNTHESIS__

		float T_between[4][4];
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 4, 4, 4, 4, 4, 4, float, float>(last_T_inv, current_pose.T, T_between);

		float between_pose_T_inv[4][4];
		pose_inv(between_pose.T, between_pose_T_inv);
		// #ifndef __SYNTHESIS__
		// 		assert(inverse_ok == 0);
		// #endif // !__SYNTHESIS__

		float residual_T[4][4];
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 4, 4, 4, 4, 4, 4, float, float>(between_pose_T_inv, T_between, residual_T);

		ln(residual_T, residual_se3);

		float residual_se3_2d[6][1] = {residual_se3[0], residual_se3[1], residual_se3[2], residual_se3[3], residual_se3[4], residual_se3[5]};
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 6, 6, 6, 1, 6, 1>(covariance_sqrt_inv, residual_se3_2d, residual_l2);
	}

	void computeJacobian()
	{
		float rho[3] = {residual_se3[0], residual_se3[1], residual_se3[2]};
		float axis_angle[3] = {residual_se3[3], residual_se3[4], residual_se3[5]};
		float rho_hat[3][3];
		float axis_angle_hat[3][3];
		hat(rho, rho_hat);
		hat(axis_angle, axis_angle_hat);
		float Jr_inv[6][6];
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				Jr_inv[i][j] = 0.5 * axis_angle_hat[i][j];
				Jr_inv[i][j + 3] = 0.5 * rho_hat[i][j];
				Jr_inv[i + 3][j] = 0;
				Jr_inv[i + 3][j + 3] = 0.5 * axis_angle_hat[i][j];
				if (i == j)
				{
					Jr_inv[i][j] += 1;
					Jr_inv[i + 3][j + 3] += 1;
				}
			}
		}

		float translation_2[3] = {pose_2[0], pose_2[1], pose_2[2]};
		// float q_2[4] = { pose_2[3], pose_2[4], pose_2[5], pose_2[6] };
		Quaterniond q_2(pose_2[3], pose_2[4], pose_2[5], pose_2[6]);
		SE3 current_pose(translation_2, q_2);
		float T_inv[4][4];
		// int inverse_ok;
		pose_inv(current_pose.T, T_inv);
		// #ifndef __SYNTHESIS__
		// 		assert(inverse_ok == 0);
		// #endif // !__SYNTHESIS__
		float ad[6][6];
		Ad(T_inv, ad);
		float jacobian_1[6][6];
		float jacobian_2[6][6];
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 6, 6, 6, 6, 6, 6, float, float>(Jr_inv, ad, jacobian_2);
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				jacobian_1[i][j] = -jacobian_2[i][j];
			}
		}

		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 6, 6, 6, 6, 6, 6, float, float>(covariance_sqrt_inv, jacobian_1, jacobian_l2_1);
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 6, 6, 6, 6, 6, 6, float, float>(covariance_sqrt_inv, jacobian_2, jacobian_l2_2);
	}

	void computeHg()
	{
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 6, 6, 6, 6, 6, 6, float, float>(jacobian_l2_1, jacobian_l2_1, H_11);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 6, 6, 6, 6, 6, 6, float, float>(jacobian_l2_1, jacobian_l2_2, H_12);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 6, 6, 6, 6, 6, 6, float, float>(jacobian_l2_2, jacobian_l2_2, H_22);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 6, 6, 6, 1, 6, 1, float, float>(jacobian_l2_1, residual_l2, g_1);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 6, 6, 6, 1, 6, 1, float, float>(jacobian_l2_2, residual_l2, g_2);
		for (int i = 0; i < 6; i++)
		{
			g_1[i][0] = -g_1[i][0];
			g_2[i][0] = -g_2[i][0];
		}
	}

	float computeCost()
	{
		float cost = 0;
		for (int i = 0; i < 6; i++)
		{
			cost += residual_l2[i][0] * residual_l2[i][0];
		}
		return cost;
	}

	int index_1, index_2;
	float pose_1[7], pose_2[7];
	float measurement[7];
	float residual_se3[6];
	float residual_l2[6][1];
	float jacobian_l2_1[6][6];
	float jacobian_l2_2[6][6];
	float H_11[6][6];
	float H_22[6][6];
	float H_12[6][6];
	float g_1[6][1];
	float g_2[6][1];
	float covariance[6][6] = {0};
	float covariance_sqrt_inv[6][6] = {0};
};