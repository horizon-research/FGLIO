#pragma once
#include "common.h"

const float gravity[3] = {0, 0, -9.81};
const float t = 0.01;
const float imuFactorCov[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const float J_alpha_bai[3][3] = {0};
const float J_alpha_bgi[3][3] = {0};
const float J_beta_bai[3][3] = {0};
const float J_beta_bgi[3][3] = {0};
const float J_gamma_bgi[3][3] = {0};

class imuFactor
{
public:
	imuFactor(int ind1, int ind2, float p_1[7], float v_1[3], float ba1[3], float bg1[3], float p_2[7], float v_2[3], float ba2[3], float bg2[3], float alp[3], float bet[3], float gam[4], const float t, const float cov[15], const float J_alpha_bai[3][3], const float J_alpha_bgi[3][3], const float J_beta_bai[3][3], const float J_beta_bgi[3][3], const float J_gamma_bgi[3][3])
	{
		index_pose1 = ind1;
		index_vbabg1 = ind1;
		index_pose2 = ind2;
		index_vbabg2 = ind2;

		delta_t = t;

		for (int i = 0; i < 7; i++)
		{
			pose_1[i] = p_1[i];
			pose_2[i] = p_2[i];
		}

		for (int i = 0; i < 3; i++)
		{
			velocity_1[i] = v_1[i];
			velocity_2[i] = v_2[i];
			ba_1[i] = ba1[i];
			ba_2[i] = ba2[i];
			bg_1[i] = bg1[i];
			bg_2[i] = bg2[i];
			alpha[i] = alp[i];
			beta[i] = bet[i];
			gamma[i] = gam[i];
		}

		gamma[3] = gam[3];

		for (int i = 0; i < 15; i++)
		{
			covariance[i][i] = cov[i];
#ifdef VS
			covariance_sqrt_inv[i][i] = 1 / std::sqrt(cov[i]);
#else
			covariance_sqrt_inv[i][i] = 1 / hls::sqrt(cov[i]);
#endif // VS
		}

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				this->J_alpha_bai[i][j] = J_alpha_bai[i][j];
				this->J_alpha_bgi[i][j] = J_alpha_bgi[i][j];
				this->J_beta_bai[i][j] = J_beta_bai[i][j];
				this->J_beta_bgi[i][j] = J_beta_bgi[i][j];
				this->J_gamma_bgi[i][j] = J_gamma_bgi[i][j];
			}
		}
	}

	void computeResidual()
	{
		float t_i[3] = {pose_1[0], pose_1[1], pose_1[2]};
		Quaterniond q_i(pose_1[3], pose_1[4], pose_1[5], pose_1[6]);
		SE3 pose_i(t_i, q_i);

		float t_j[3] = {pose_2[0], pose_2[1], pose_2[2]};
		Quaterniond q_j(pose_2[3], pose_2[4], pose_2[5], pose_2[6]);
		SE3 pose_j(t_j, q_j);

		float ideal_trans_w[3][1];
		for (int i = 0; i < 3; i++)
		{
			ideal_trans_w[i][0] = pose_j.translation[i] - pose_i.translation[i] - velocity_1[i] * delta_t + 0.5 * gravity[i] * delta_t * delta_t;
		}

		float ideal_trans_i[3][1];
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 3, 3, 3, 1, 3, 1, float, float>(pose_i.rotation, ideal_trans_w, ideal_trans_i); //这里对R取了转置，相当于取逆
		float residual_trans[3];
		for (int i = 0; i < 3; i++)
		{
			residual_trans[i] = ideal_trans_i[i][0] - alpha[i];
		}

		Quaterniond q_i_inv = q_i.inverse();
		Quaterniond q_temp = q_i_inv * q_j;
		Quaterniond q_gam(gamma[0], gamma[1], gamma[2], gamma[3]);
		Quaterniond q_gam_inv = q_gam.inverse();
		q_temp = q_temp * q_gam_inv;
		float residual_rot[3] = {q_temp.q1 * 2, q_temp.q2 * 2, q_temp.q3 * 2};

		float ideal_vel_w[3][1];
		for (int i = 0; i < 3; i++)
		{
			ideal_vel_w[i][0] = velocity_2[i] + gravity[i] * delta_t - velocity_1[i];
		}
		float ideal_vel_i[3][1];
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 3, 3, 3, 1, 3, 1, float, float>(pose_i.rotation, ideal_vel_w, ideal_vel_i);
		float residual_vel[3];
		for (int i = 0; i < 3; i++)
		{
			residual_vel[i] = ideal_vel_i[i][0] - beta[i];
		}

		float residual_bias[6];
		for (int i = 0; i < 3; i++)
		{
			residual_bias[i] = ba_2[i] - ba_1[i];
			residual_bias[i + 3] = bg_2[i] - bg_1[i];
		}

		float residual[15][1] = {residual_trans[0], residual_trans[1], residual_trans[2], residual_rot[0], residual_rot[1], residual_rot[2], residual_vel[0], residual_vel[1], residual_vel[2], residual_bias[0], residual_bias[1], residual_bias[2], residual_bias[3], residual_bias[4], residual_bias[5]};

		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 15, 15, 15, 1, 15, 1, float, float>(covariance_sqrt_inv, residual, residual_l2);
	}

	void computeJacobian()
	{
		float t_i[3] = {pose_1[0], pose_1[1], pose_1[2]};
		Quaterniond q_i(pose_1[3], pose_1[4], pose_1[5], pose_1[6]);
		SE3 pose_i(t_i, q_i);

		float t_j[3] = {pose_2[0], pose_2[1], pose_2[2]};
		Quaterniond q_j(pose_2[3], pose_2[4], pose_2[5], pose_2[6]);
		SE3 pose_j(t_j, q_j);

		float jacobian_pose_i[15][7] = {0};
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				jacobian_pose_i[i][j] = -pose_i.rotation[j][i];
			}
		}

		float ideal_trans_w[3][1];
		for (int i = 0; i < 3; i++)
		{
			ideal_trans_w[i][0] = pose_j.translation[i] - pose_i.translation[i] - velocity_1[i] * delta_t + 0.5 * gravity[i] * delta_t * delta_t;
		}
		float ideal_trans_i[3][1];
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 3, 3, 3, 1, 3, 1, float, float>(pose_i.rotation, ideal_trans_w, ideal_trans_i); //这里对R取了转置，相当于取逆
		float ideal_trans_i_1d[3] = {ideal_trans_i[0][0], ideal_trans_i[1][0], ideal_trans_i[2][0]};
		float J0_01[3][3];
		hat(ideal_trans_i_1d, J0_01);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				jacobian_pose_i[i][j + 3] = J0_01[i][j];
			}
		}

		Quaterniond q_j_inv = q_j.inverse();
		Quaterniond q_temp1 = q_j_inv * q_i;
		float qleft0[4][4];
		q_temp1.Qleft(qleft0);
		Quaterniond q_gam(gamma[0], gamma[1], gamma[2], gamma[3]);
		float qright0[4][4];
		q_gam.Qright(qright0);
		float I1[3][4] = {0, 1, 0, 0,
						  0, 0, 1, 0,
						  0, 0, 0, 1};
		float temp1[3][4];
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 3, 4, 4, 4, 3, 4, float, float>(I1, qleft0, temp1);
		float temp2[3][4];
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 3, 4, 4, 4, 3, 4, float, float>(temp1, qright0, temp2);
		float J0_11[3][3];
		hls::matrix_multiply<hls::NoTranspose, hls::Transpose, 3, 4, 3, 4, 3, 3, float, float>(temp2, I1, J0_11);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				jacobian_pose_i[i + 3][j + 3] = -J0_11[i][j];
			}
		}

		float ideal_vel_w[3][1];
		for (int i = 0; i < 3; i++)
		{
			ideal_vel_w[i][0] = velocity_2[i] + gravity[i] * delta_t - velocity_1[i];
		}
		float ideal_vel_i[3][1];
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 3, 3, 3, 1, 3, 1, float, float>(pose_i.rotation, ideal_vel_w, ideal_vel_i);
		float J0_21[3][3];
		float ideal_vel_i_1d[3] = {ideal_vel_i[0][0], ideal_vel_i[1][0], ideal_vel_i[2][0]};
		hat(ideal_vel_i_1d, J0_21);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				jacobian_pose_i[i + 6][j + 3] = J0_21[i][j];
			}
		}

		float jacobian_vbabg_i[15][9] = {0};
		float I_33[3][3] = {1, 0, 0,
							0, 1, 0,
							0, 0, 1};
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				jacobian_vbabg_i[i][j] = -pose_i.rotation[j][i] * delta_t;
				jacobian_vbabg_i[i + 6][j] = -pose_i.rotation[j][i];
			}
		}
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				jacobian_vbabg_i[i][j + 3] = -J_alpha_bai[i][j];
				jacobian_vbabg_i[i][j + 6] = -J_alpha_bgi[i][j];
				jacobian_vbabg_i[i + 6][j + 3] = -J_beta_bai[i][j];
				jacobian_vbabg_i[i + 6][j + 6] = -J_beta_bgi[i][j];
				jacobian_vbabg_i[i + 9][j + 3] = -I_33[i][j];
				jacobian_vbabg_i[i + 12][j + 6] = -I_33[i][j];
			}
		}

		// Quaterniond q_j_inv = q_j.inverse();
		Quaterniond q_temp2 = q_j_inv * q_i;
		q_temp2 = q_temp2 * q_gam;
		float qleft1[4][4];
		q_temp2.Qleft(qleft1);
		float temp3[3][4];
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 3, 4, 4, 4, 3, 4, float, float>(I1, qleft1, temp3);
		float temp4[3][3];
		hls::matrix_multiply<hls::NoTranspose, hls::Transpose, 3, 4, 3, 4, 3, 3, float, float>(temp1, I1, temp4);
		float J1_12[3][3];
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 3, 3, 3, 3, 3, 3, float, float>(temp4, J_gamma_bgi, J1_12);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				jacobian_vbabg_i[i + 3][j + 6] = -J1_12[i][j];
			}
		}


		float jacobian_pose_j[15][7];
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				jacobian_pose_j[i][j] = pose_i.rotation[j][i];
			}
		}

		Quaterniond q_gam_inv = q_gam.inverse();
		Quaterniond q_i_inv = q_i.inverse();
		Quaterniond q_temp3 = q_gam_inv * q_i_inv;
		q_temp3 = q_temp3 * q_j;
		float qleft2[4][4];
		q_temp3.Qleft(qleft2);
		float temp5[3][4];
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 3, 4, 4, 4, 3, 4, float, float>(I1, qleft2, temp5);
		float J2_11[3][3];
		hls::matrix_multiply<hls::NoTranspose, hls::Transpose, 3, 4, 3, 4, 3, 3, float, float>(temp5, I1, J2_11);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				jacobian_pose_j[i + 3][j + 3] = J2_11[i][j];
			}
		}


		float jacobian_vbabg_j[15][9];
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				jacobian_vbabg_j[i + 6][j] = pose_i.rotation[j][i];
				jacobian_vbabg_j[i + 9][j + 3] = I_33[i][j];
				jacobian_vbabg_j[i + 12][j + 6] = I_33[i][j];
			}
		}

		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 15, 15, 15, 7, 15, 7, float, float>(covariance_sqrt_inv, jacobian_pose_i, jacobian_l2_pose_i);
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 15, 15, 15, 9, 15, 9, float, float>(covariance_sqrt_inv, jacobian_vbabg_i, jacobian_l2_v_ba_bg_i);
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 15, 15, 15, 7, 15, 7, float, float>(covariance_sqrt_inv, jacobian_pose_j, jacobian_l2_pose_j);
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 15, 15, 15, 9, 15, 9, float, float>(covariance_sqrt_inv, jacobian_vbabg_j, jacobian_l2_v_ba_bg_j);
	}

	void computeHg()
	{
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 7, 15, 7, 7, 7, float, float>(jacobian_l2_pose_i, jacobian_l2_pose_i, H_pose_i);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 7, 15, 7, 7, 7, float, float>(jacobian_l2_pose_j, jacobian_l2_pose_j, H_pose_j);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 9, 15, 9, 9, 9, float, float>(jacobian_l2_v_ba_bg_i, jacobian_l2_v_ba_bg_i, H_vbabg_i);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 9, 15, 9, 9, 9, float, float>(jacobian_l2_v_ba_bg_j, jacobian_l2_v_ba_bg_j, H_vbabg_j);

		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 7, 15, 7, 7, 7, float, float>(jacobian_l2_pose_i, jacobian_l2_pose_j, H_pose_ij);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 9, 15, 9, 9, 9, float, float>(jacobian_l2_v_ba_bg_i, jacobian_l2_v_ba_bg_j, H_vbabg_ij);

		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 7, 15, 9, 7, 9, float, float>(jacobian_l2_pose_i, jacobian_l2_v_ba_bg_i, H_posei_vbabgi);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 7, 15, 9, 7, 9, float, float>(jacobian_l2_pose_i, jacobian_l2_v_ba_bg_j, H_posei_vbabgj);

		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 7, 15, 9, 7, 9, float, float>(jacobian_l2_pose_j, jacobian_l2_v_ba_bg_j, H_posej_vbabgj);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 7, 15, 9, 7, 9, float, float>(jacobian_l2_pose_j, jacobian_l2_v_ba_bg_i, H_posej_vbabgi);

		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 7, 15, 1, 7, 1, float, float>(jacobian_l2_pose_i, residual_l2, g_pose_i);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 7, 15, 1, 7, 1, float, float>(jacobian_l2_pose_j, residual_l2, g_pose_j);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 9, 15, 1, 9, 1, float, float>(jacobian_l2_v_ba_bg_i, residual_l2, g_vbabg_i);
		hls::matrix_multiply<hls::Transpose, hls::NoTranspose, 15, 9, 15, 1, 9, 1, float, float>(jacobian_l2_v_ba_bg_j, residual_l2, g_vbabg_j);

		for (int i = 0; i < 7; i++)
		{
			g_pose_i[i][0] = -g_pose_i[i][0];
			g_pose_j[i][0] = -g_pose_j[i][0];
		}

		for (int i = 0; i < 9; i++)
		{
			g_vbabg_i[i][0] = -g_vbabg_i[i][0];
			g_vbabg_j[i][0] = -g_vbabg_j[i][0];
		}
	}

	float computeCost()
	{
		float cost = 0;
		for (int i = 0; i < 15; i++)
		{
			cost += residual_l2[i][0] * residual_l2[i][0];
		}
		return cost;
	}

	int index_pose1;
	int index_pose2;
	int index_vbabg1;
	int index_vbabg2;

	float pose_1[7];
	float pose_2[7];
	float velocity_1[3];
	float velocity_2[3];
	float ba_1[3];
	float bg_1[3];
	float ba_2[3];
	float bg_2[3];

	float alpha[3]; 
	float beta[3];	
	float gamma[4];

	float residual_l2[15][1];
	float jacobian_l2_pose_i[15][7];
	float jacobian_l2_v_ba_bg_i[15][9];
	float jacobian_l2_pose_j[15][7];
	float jacobian_l2_v_ba_bg_j[15][9];
	float H_pose_i[7][7];
	float H_pose_j[7][7];
	float H_vbabg_i[9][9];
	float H_vbabg_j[9][9];
	float H_pose_ij[7][7];
	float H_vbabg_ij[9][9];
	float H_posei_vbabgi[7][9];
	float H_posei_vbabgj[7][9];
	float H_posej_vbabgi[7][9];
	float H_posej_vbabgj[7][9];
	float g_pose_i[7][1];
	float g_pose_j[7][1];
	float g_vbabg_i[9][1];
	float g_vbabg_j[9][1];

	float covariance[15][15] = {0};
	float covariance_sqrt_inv[15][15] = {0};

	float delta_t;
	float J_alpha_bai[3][3];
	float J_alpha_bgi[3][3];
	float J_beta_bai[3][3];
	float J_beta_bgi[3][3];
	float J_gamma_bgi[3][3];
};