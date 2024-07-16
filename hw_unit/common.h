#pragma once
#include <ap_int.h>
#include <hls_math.h>
#include <hls_linear_algebra.h>
// #define VS

class Quaterniond
{
public:
	Quaterniond(float q0, float q1, float q2, float q3)
	{

		this->q0 = q0;
		this->q1 = q1;
		this->q2 = q2;
		this->q3 = q3;
		q[0] = q0;
		q[1] = q1;
		q[2] = q2;
		q[3] = q3;
	}

	// Quaterniond(const Quaterniond& q)
	//{
	//	this->q0 = q.q0;
	//	this->q1 = q.q1;
	//	this->q2 = q.q2;
	//	this->q3 = q.q3;
	//	this->q[0] = q.q0;
	//	this->q[1] = q.q1;
	//	this->q[2] = q.q2;
	//	this->q[3] = q.q3;
	// }

	Quaterniond &operator*(const Quaterniond &q)
	{
		float real = this->q0 * q.q0 - this->q1 * q.q1 - this->q2 * q.q2 - this->q3 * q.q3;
		float i = this->q0 * q.q1 + this->q1 * q.q0 + this->q2 * q.q3 - this->q3 * q.q2;
		float j = this->q0 * q.q2 - this->q1 * q.q3 + this->q2 * q.q0 + this->q3 * q.q1;
		float k = this->q0 * q.q3 + this->q1 * q.q2 - this->q2 * q.q1 + this->q3 * q.q0;
		Quaterniond temp(real, i, j, k);
		return temp;
	}

	float mod2()
	{
		float m = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
		return m;
	}

	Quaterniond &inverse()
	{
		Quaterniond q_conjugate(q0, -q1, -q2, -q3);
		Quaterniond q_inv(q_conjugate.q0 / this->mod2(), q_conjugate.q1 / this->mod2(), q_conjugate.q2 / this->mod2(), q_conjugate.q3 / this->mod2());
		return q_inv;
	}

	void Qleft(float qleft[4][4])
	{
		qleft[0][0] = this->q0;
		qleft[0][1] = -this->q1;
		qleft[0][2] = -this->q2;
		qleft[0][3] = -this->q3;
		qleft[1][0] = this->q1;
		qleft[1][1] = this->q0;
		qleft[1][2] = -this->q3;
		qleft[1][3] = this->q2;
		qleft[2][0] = this->q2;
		qleft[2][1] = this->q3;
		qleft[2][2] = this->q0;
		qleft[2][3] = -this->q1;
		qleft[3][0] = this->q3;
		qleft[3][1] = -this->q2;
		qleft[3][2] = this->q1;
		qleft[3][3] = this->q0;
	}

	void Qright(float qright[4][4])
	{
		qright[0][0] = this->q0;
		qright[0][1] = -this->q1;
		qright[0][2] = -this->q2;
		qright[0][3] = -this->q3;
		qright[1][0] = this->q1;
		qright[1][1] = this->q0;
		qright[1][2] = this->q3;
		qright[1][3] = -this->q2;
		qright[2][0] = this->q2;
		qright[2][1] = -this->q3;
		qright[2][2] = this->q0;
		qright[2][3] = this->q1;
		qright[3][0] = this->q3;
		qright[3][1] = this->q2;
		qright[3][2] = -this->q1;
		qright[3][3] = this->q0;
	}

	float q[4];
	float q0;
	float q1;
	float q2;
	float q3;
};

class Axis_Angle
{
public:
	Axis_Angle(float theta, float n[3])
	{
		this->theta = theta;
		for (int i = 0; i < 3; i++)
		{
			this->n[i] = n[i];
			this->axis_angle[i] = theta * n[i];
		}
	}
	float axis_angle[3];
	float n[3];
	float theta;
};

class SE3
{
public:
	SE3(float trans[3], Quaterniond quat)
	{
		// float q[4] = { quat.q0 , quat.q1, quat.q2, quat.q3 };
		for (int i = 0; i < 4; i++)
		{
			this->quaterniond[i] = quat.q[i];
		}
		float *q = this->quaterniond;

		for (int i = 0; i < 3; i++)
		{
			translation[i] = trans[i];
		}
#ifdef VS
		float theta = 2 * std::acos(q[0]);
#else
		float theta = 2 * hls::acos(q[0]);
#endif // VS
		float n[3];
		for (int i = 0; i < 3; i++)
		{

#ifdef VS
			n[i] = q[i + 1] / std::sin(theta / 2);

#else
			n[i] = q[i + 1] / hls::sin(theta / 2);

#endif // VS
		}
		for (int i = 0; i < 3; i++)
		{
			axis_angle[i] = theta * n[i];
		}
		rotation[0][0] = 1 - 2 * q[2] * q[2] - 2 * q[3] * q[3];
		rotation[0][1] = 2 * q[1] * q[2] - 2 * q[0] * q[3];
		rotation[0][2] = 2 * q[1] * q[3] + 2 * q[0] * q[2];

		rotation[1][0] = 2 * q[1] * q[2] + 2 * q[0] * q[3];
		rotation[1][1] = 1 - 2 * q[1] * q[1] - 2 * q[3] * q[3];
		rotation[1][2] = 2 * q[2] * q[3] - 2 * q[0] * q[1];

		rotation[2][0] = 2 * q[1] * q[3] - 2 * q[0] * q[2];
		rotation[2][1] = 2 * q[2] * q[3] + 2 * q[0] * q[1];
		rotation[2][2] = 1 - 2 * q[1] * q[1] - 2 * q[2] * q[2];

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				T[i][j] = rotation[i][j];
			}
			T[i][3] = translation[i];
			T[3][i] = 0;
		}
		T[3][3] = 1;

#ifdef VS
		float aa = theta / (2 * std::tan(theta / 2));
#else
		float aa = theta / (2 * hls::tan(theta / 2));
#endif // VS

		float tmp1[3][3] = {aa, 0, 0,
							0, aa, 0,
							0, 0, aa};
		float tmp2[3][3] = {(1 - aa) * n[0] * n[0], (1 - aa) * n[0] * n[1], (1 - aa) * n[0] * n[2],
							(1 - aa) * n[1] * n[0], (1 - aa) * n[1] * n[1], (1 - aa) * n[1] * n[2],
							(1 - aa) * n[2] * n[0], (1 - aa) * n[2] * n[1], (1 - aa) * n[2] * n[2]};
		float tmp3[3][3] = {0, -n[2] * theta / 2, n[1] * theta / 2,
							n[2] * theta / 2, 0, -n[0] * theta / 2,
							-n[1] * theta / 2, n[0] * theta / 2, 0};

		// float tmp1[3][3] = { hls::sin(theta) / theta,0,0,
		//							0,hls::sin(theta) / theta,0,
		//							0,0,hls::sin(theta) / theta };
		// float tmp2[3][3] = { (1 - hls::sin(theta) / theta) * n[0] * n[0], (1 - hls::sin(theta) / theta) * n[0] * n[1], (1 - hls::sin(theta) / theta) * n[0] * n[2],
		//							(1 - hls::sin(theta) / theta) * n[1] * n[0], (1 - hls::sin(theta) / theta) * n[1] * n[1], (1 - hls::sin(theta) / theta) * n[1] * n[2],
		//							(1 - hls::sin(theta) / theta) * n[2] * n[0], (1 - hls::sin(theta) / theta) * n[2] * n[1], (1 - hls::sin(theta) / theta) * n[2] * n[2] };
		// float tmp3[3][3] = { 0, -n[2] * (1 - hls::cos(theta)) / theta, n[1] * (1 - hls::cos(theta)) / theta,
		//							n[2] * (1 - hls::cos(theta)) / theta, 0, -n[0] * (1 - hls::cos(theta)) / theta,
		//							-n[1] * (1 - hls::cos(theta)) / theta, n[0] * (1 - hls::cos(theta)) / theta, 0 };
		// float n_hat[3][3];
		// hat(n, n_hat);
		// for (int i = 0; i < 3; i++)
		//{
		//	for (int j = 0; j < 3; j++)
		//	{
		//		if (i == j)
		//		{
		//			break;
		//		}
		//		tmp3[i][j] = (1 - hls::cos(theta)) / theta * n_hat[i][j];
		//	}
		// }

		float J_inv[3][3];
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				J_inv[i][j] = tmp1[i][j] + tmp2[i][j] - tmp3[i][j];
			}
		}

		float rho[3][1];
		float translation_2d[3][1] = {translation[0], translation[1], translation[2]};
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 3, 3, 3, 1, 3, 1, float, float>(J_inv, translation_2d, rho);

		for (int i = 0; i < 3; i++)
		{
			se3[i] = rho[i][0];
			se3[i + 3] = axis_angle[i];
		}
	}

	SE3(float trans[3], Axis_Angle axis_angle)
	{
		// float axi_ang[3] = { axis_angle.axis_angle[0], axis_angle.axis_angle[1], axis_angle.axis_angle[2] };
		for (int i = 0; i < 3; i++)
		{
			this->axis_angle[i] = axis_angle.axis_angle[i];
		}
		for (int i = 0; i < 3; i++)
		{
			translation[i] = trans[i];
		}
		float theta_2 = 0;
		for (int i = 0; i < 3; i++)
		{
			theta_2 += this->axis_angle[i] * this->axis_angle[i];
		}
#ifdef VS
		float theta = std::sqrt(theta_2);

#else
		float theta = hls::sqrt(theta_2);
#endif // VS
		float n[3];
		for (int i = 0; i < 3; i++)
		{
			n[i] = this->axis_angle[i] / theta;
		}

#ifdef VS
		quaterniond[0] = std::cos(theta / 2);
		quaterniond[1] = n[0] * std::sin(theta / 2);
		quaterniond[2] = n[1] * std::sin(theta / 2);
		quaterniond[3] = n[2] * std::sin(theta / 2);
#else
		quaterniond[0] = hls::cos(theta / 2);
		quaterniond[1] = n[0] * hls::sin(theta / 2);
		quaterniond[2] = n[1] * hls::sin(theta / 2);
		quaterniond[3] = n[2] * hls::sin(theta / 2);
#endif // VS

		rotation[0][0] = 1 - 2 * quaterniond[2] * quaterniond[2] - 2 * quaterniond[3] * quaterniond[3];
		rotation[0][1] = 2 * quaterniond[1] * quaterniond[2] - 2 * quaterniond[0] * quaterniond[3];
		rotation[0][2] = 2 * quaterniond[1] * quaterniond[3] + 2 * quaterniond[0] * quaterniond[2];

		rotation[1][0] = 2 * quaterniond[1] * quaterniond[2] + 2 * quaterniond[0] * quaterniond[3];
		rotation[1][1] = 1 - 2 * quaterniond[1] * quaterniond[1] - 2 * quaterniond[3] * quaterniond[3];
		rotation[1][2] = 2 * quaterniond[2] * quaterniond[3] - 2 * quaterniond[0] * quaterniond[1];

		rotation[2][0] = 2 * quaterniond[1] * quaterniond[3] - 2 * quaterniond[0] * quaterniond[2];
		rotation[2][1] = 2 * quaterniond[2] * quaterniond[3] + 2 * quaterniond[0] * quaterniond[1];
		rotation[2][2] = 1 - 2 * quaterniond[1] * quaterniond[1] - 2 * quaterniond[2] * quaterniond[2];

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				T[i][j] = rotation[i][j];
			}
			T[i][3] = translation[i];
			T[3][i] = 0;
		}
		T[3][3] = 1;

#ifdef VS
		float aa = theta / (2 * std::tan(theta / 2));

#else
		float aa = theta / (2 * hls::tan(theta / 2));

#endif // VS

		float tmp1[3][3] = {aa, 0, 0,
							0, aa, 0,
							0, 0, aa};
		float tmp2[3][3] = {(1 - aa) * n[0] * n[0], (1 - aa) * n[0] * n[1], (1 - aa) * n[0] * n[2],
							(1 - aa) * n[1] * n[0], (1 - aa) * n[1] * n[1], (1 - aa) * n[1] * n[2],
							(1 - aa) * n[2] * n[0], (1 - aa) * n[2] * n[1], (1 - aa) * n[2] * n[2]};
		float tmp3[3][3] = {0, -n[2] * theta / 2, n[1] * theta / 2,
							n[2] * theta / 2, 0, -n[0] * theta / 2,
							-n[1] * theta / 2, n[0] * theta / 2, 0};

		float J_inv[3][3];
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				J_inv[i][j] = tmp1[i][j] + tmp2[i][j] - tmp3[i][j];
			}
		}

		float rho[3][1];
		float translation_2d[3][1] = {translation[0], translation[1], translation[2]};
		hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 3, 3, 3, 1, 3, 1, float, float>(J_inv, translation_2d, rho);

		for (int i = 0; i < 3; i++)
		{
			se3[i] = rho[i][0];
			se3[i + 3] = this->axis_angle[i];
		}
	}

	float translation[3];
	float quaterniond[4];
	float rotation[3][3];
	float axis_angle[3];
	float T[4][4];
	float se3[6];
};

void hat(float a[3], float m[3][3])
{
	m[0][1] = -a[2];
	m[0][2] = a[1];
	m[1][0] = a[2];
	m[1][2] = -a[0];
	m[2][0] = -a[1];
	m[2][1] = a[0];
	m[0][0] = m[1][1] = m[2][2] = 0;
}

void Ad(float T[4][4], float ad[6][6])
{
	float R[3][3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			R[i][j] = T[i][j];
		}
	}
	float t[3];
	for (int i = 0; i < 3; i++)
	{
		t[i] = T[i][3];
	}
	float t_hat[3][3];
	hat(t, t_hat);

	float mul[3][3];
	hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 3, 3, 3, 3, 3, 3, float, float>(t_hat, R, mul);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			ad[i][j] = R[i][j];
			ad[i + 3][j + 3] = R[i][j];
			ad[i][j + 3] = mul[i][j];
			ad[i + 3][j] = 0;
		}
	}
}

void ln(float T[4][4], float se3[6])
{
	float R[3][3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			R[i][j] = T[i][j];
		}
	}

	float tr_R = R[0][0] + R[1][1] + R[2][2];
#ifdef VS
	float theta = std::acos((tr_R - 1) / 2);

#else
	float theta = hls::acos((tr_R - 1) / 2);

#endif // VS

	float tmp[3][3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
#ifdef VS
			tmp[i][j] = (R[i][j] - R[j][i]) / (2 * std::sin(theta));

#else
			tmp[i][j] = (R[i][j] - R[j][i]) / (2 * hls::sin(theta));

#endif // VS
		}
	}

	float n[3] = {tmp[2][1], tmp[0][2], tmp[1][0]};
	// float axis_angle[3] = { theta * n[0], theta * n[1], theta * n[2] };

	float t[3] = {T[0][3], T[1][3], T[2][3]};

	SE3 obj(t, Axis_Angle(theta, n));
	for (int i = 0; i < 6; i++)
	{
		se3[i] = obj.se3[i];
	}
}

void pose_inv(float T[4][4], float T_inv[4][4])
{
	float T_inv_R[3][3];
	float T_t[3][1];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			float elem = T[j][i];
			T_inv[i][j] = elem;
			T_inv_R[i][j] = elem;
		}
		T_t[i][0] = T[i][3];
	}

	float T_inv_t_neg[3][1];
	hls::matrix_multiply<hls::NoTranspose, hls::NoTranspose, 3, 3, 3, 1, 3, 1, float, float>(T_inv_R, T_t, T_inv_t_neg);

	for (int i = 0; i < 3; i++)
	{
		T_inv[i][3] = -T_inv_t_neg[i][0];
	}

	T_inv[3][0] = T_inv[3][1] = T_inv[3][2] = 0;
	T_inv[3][3] = 1;
}