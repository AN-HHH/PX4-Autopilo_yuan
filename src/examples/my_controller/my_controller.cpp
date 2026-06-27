#include "my_controller.hpp"

#include <drivers/drv_hrt.h>
#include <lib/mathlib/mathlib.h>

#include <cmath>

MyController::MyController() :
	ModuleParams(nullptr),
	ScheduledWorkItem(MODULE_NAME, px4::wq_configurations::rate_ctrl)
{
}

bool MyController::init()
{
	ScheduleOnInterval(20000); // 20 ms
	return true;
}

float MyController::signNZ(float x)
{
	return x >= 0.f ? 1.f : -1.f;
}

float MyController::yawDrag(float r) const
{
	return _d1 * r + _d2 * fabsf(r) * r;
}

void MyController::controllerStep(float tk, float dt,
				  const PlantState &xp,
				  const ControllerCommand &cmd,
				  ControllerState &xc,
				  ControllerOutput &out)
{
	const float phi = xp.phi;
	const float theta = xp.theta;
	const float psi = xp.psi;
	const float p = xp.p;
	const float q = xp.q;
	const float r = xp.r;
	const float omega1 = xp.omega1;
	const float omega2 = xp.omega2;

	const float cphi = cosf(phi);
	const float sphi = sinf(phi);
	const float cth = cosf(theta);
	const float sth = sinf(theta);
	const float cpsi = cosf(psi);
	const float spsi = sinf(psi);

	matrix::Dcmf R;

	R(0, 0) = cpsi * cth;
	R(0, 1) = cpsi * sth * sphi - spsi * cphi;
	R(0, 2) = cpsi * sth * cphi + spsi * sphi;

	R(1, 0) = spsi * cth;
	R(1, 1) = spsi * sth * sphi + cpsi * cphi;
	R(1, 2) = spsi * sth * cphi - cpsi * sphi;

	R(2, 0) = -sth;
	R(2, 1) = cth * sphi;
	R(2, 2) = cth * cphi;

	const float weff1 = omega1 - r;
	const float weff2 = omega2 - r;

	const float T1 = _b * weff1 * fabsf(weff1);
	const float T2 = _b * weff2 * fabsf(weff2);
	const float u1 = T1 + T2;

	const matrix::Vector3f b3 = R * matrix::Vector3f(0.f, 0.f, 1.f);

	const float tau_z_est = _kq * (weff1 * fabsf(weff1) + weff2 * fabsf(weff2));
	const float r_dot_est = ((_Jx - _Jy) * p * q + tau_z_est - yawDrag(r)) / _Jz;

	matrix::Vector3f zd_E = cmd.zd_cmd_E;
	matrix::Vector3f zd_dot_E = cmd.zd_dot_cmd_E;
	matrix::Vector3f zd_ddot_E = cmd.zd_ddot_cmd_E;

	const float nzd = zd_E.norm();

	if (nzd > 1e-6f) {
		zd_E /= nzd;

	} else {
		zd_E = matrix::Vector3f(0.f, 0.f, 1.f);
	}

	const matrix::Vector3f c = R.transpose() * zd_E;

	const float c1 = c(0);
	const float c2 = c(1);
	const float c3 = c(2);

	const float c3_safe = signNZ(c3) * math::max(fabsf(c3), _c3_eps);

	const matrix::Vector3f u_ff = R.transpose() * zd_dot_E;
	const matrix::Vector3f omegaB(p, q, r);

	const matrix::Vector3f u_ff_dot =
		matrix::Vector3f(u_ff.cross(omegaB)) + R.transpose() * zd_ddot_E;

	const matrix::Vector3f c_dot =
		matrix::Vector3f(c.cross(omegaB)) + u_ff;

	const float c1_dot = c_dot(0);
	const float c2_dot = c_dot(1);
	const float c3_dot = c_dot(2);

	const float p_d = (r * c1 - u_ff(1) - _k_c2 * c2) / c3_safe;
	const float q_d = (r * c2 + u_ff(0) + _k_c1 * c1) / c3_safe;

	const float z_p = p - p_d;
	const float z_q = q - q_d;

	const float Np = r * c1 - u_ff(1) - _k_c2 * c2;
	const float Nq = r * c2 + u_ff(0) + _k_c1 * c1;

	const float Np_dot = r_dot_est * c1 + r * c1_dot - u_ff_dot(1) - _k_c2 * c2_dot;
	const float Nq_dot = r_dot_est * c2 + r * c2_dot + u_ff_dot(0) + _k_c1 * c1_dot;

	const float p_d_dot = (Np_dot * c3_safe - Np * c3_dot) / (c3_safe * c3_safe);
	const float q_d_dot = (Nq_dot * c3_safe - Nq * c3_dot) / (c3_safe * c3_safe);

	const float f_p = (_Jy - _Jz) * q * r - _Jr * (omega1 + omega2) * q;
	const float f_q = (_Jz - _Jx) * p * r + _Jr * (omega1 + omega2) * p;

	const float tau_x_d = -f_p + _Jx * (p_d_dot - _k_p * z_p - c2 * c3);
	const float tau_y_d = -f_q + _Jy * (q_d_dot - _k_q * z_q + c1 * c3);

	const matrix::Vector3f tau_perp_d_B(tau_x_d, tau_y_d, 0.f);
	const matrix::Vector3f tau_perp_d_E = R * tau_perp_d_B;

	const float mhx = tau_perp_d_E(0);
	const float mhy = tau_perp_d_E(1);

	const float a = sinf(theta) * sinf(phi);
	const float b = cosf(phi);
	const float detM = a * a + b * b;

	float alpha_c_raw = 0.f;
	float alpha_s_raw = 0.f;

	if (detM > _det_eps) {
		const float inv00 =  a / detM;
		const float inv01 =  b / detM;
		const float inv10 = -b / detM;
		const float inv11 =  a / detM;

		alpha_c_raw = (2.f / _l) * (inv00 * mhx + inv01 * mhy);
		alpha_s_raw = (2.f / _l) * (inv10 * mhx + inv11 * mhy);
	}

	xc.alpha_cf += dt * (-_lambda_alpha * (xc.alpha_cf - alpha_c_raw));
	xc.alpha_sf += dt * (-_lambda_alpha * (xc.alpha_sf - alpha_s_raw));

	float u1_ref = cmd.u1_ref;

	if (_use_vertical_comp) {
		const float b3z_safe = math::max(b3(2), _b3z_eps);
		u1_ref = _m * _g / b3z_safe;
	}

	const float e_u = u1 - u1_ref;
	const float u1_cmd_unsat = u1_ref - _k_up * e_u - _k_ui * xc.eta_u;

	// const float gain_spin = 1.f - expf(-math::max(tk - _t_spinup, 0.f) / 0.8f);
	// const float delta_raw = gain_spin * (xc.alpha_cf * cosf(psi) + xc.alpha_sf * sinf(psi));
	const float delta_raw = (xc.alpha_cf * cosf(psi) + xc.alpha_sf * sinf(psi));

	const float u1_cmd_min = fabsf(delta_raw) + _u1_min_margin;

	const float u1_max = _u1_max_factor * _m * _g;
	const float u1_cmd = math::constrain(u1_cmd_unsat, u1_cmd_min, u1_max);

	const float eta_u_dot = e_u + _k_aw * (u1_cmd - u1_cmd_unsat);
	xc.eta_u += dt * eta_u_dot;

	const float delta_max = math::max(u1_cmd - 1e-6f, 0.f);
	const float delta = math::constrain(delta_raw, -delta_max, delta_max);

//     const float T1_cmd = math::max((u1_cmd + delta) * 0.5f, 0.f);
//     const float T2_cmd = math::max((u1_cmd - delta) * 0.5f, 0.f);

    const float T1_cmd = math::max((u1_cmd ) * 0.5f, 0.f);
    const float T2_cmd = math::max((u1_cmd ) * 0.5f, 0.f);

	const float weff1_cmd = sqrtf(math::max(T1_cmd / _b, 0.f));
	const float weff2_cmd = sqrtf(math::max(T2_cmd / _b, 0.f));

    out.omega1_cmd = math::constrain(r + weff1_cmd, 0.f, _omega_max);
    out.omega2_cmd = math::constrain(r + weff2_cmd, 0.f, _omega_max);
    // out.omega1_cmd = math::constrain(weff1_cmd, 0.f, _omega_max);
    // out.omega2_cmd = math::constrain(weff2_cmd, 0.f, _omega_max);

	out.T1_cmd = T1_cmd;
	out.T2_cmd = T2_cmd;

	out.p_d = p_d;
	out.q_d = q_d;
	out.tau_x_d = tau_x_d;
	out.tau_y_d = tau_y_d;
	out.u1_cmd = u1_cmd;
	out.delta = delta;
}

void MyController::Run()
{
	if (should_exit()) {
		ScheduleClear();
		exit_and_cleanup();
		return;
	}

	const uint64_t now_us = hrt_absolute_time();

	if (_att_sub.updated()) {
		_att_sub.copy(&_att);
	}

	if (_ang_vel_sub.updated()) {
		_ang_vel_sub.copy(&_ang_vel);
	}

	if (_manual_sub.updated()) {
		_manual_sub.copy(&_manual);
	}

	if (_armed_sub.updated()) {
		_armed_sub.copy(&_armed);
	}

	esc_status_s esc_status{};

	if (_esc_status_sub.update(&esc_status) && esc_status.esc_count >= 2) {
		const float rpm1 = static_cast<float>(esc_status.esc[0].esc_rpm);
		const float rpm2 = static_cast<float>(esc_status.esc[1].esc_rpm);

		if (PX4_ISFINITE(rpm1) && rpm1 > 1.f) {
			_motor1_rpm_meas = rpm1;
		}

		if (PX4_ISFINITE(rpm2) && rpm2 > 1.f) {
			_motor2_rpm_meas = rpm2;
		}
	}

	if (_rpm_sub.updated()) {
		_rpm_sub.copy(&_rpm);

		if (PX4_ISFINITE(_rpm.rpm_estimate) && _rpm.rpm_estimate > 1.f) {
			if (!PX4_ISFINITE(_motor1_rpm_meas)) {
				_motor1_rpm_meas = _rpm.rpm_estimate;
			}
			if (!PX4_ISFINITE(_motor2_rpm_meas)) {
				_motor2_rpm_meas = _rpm.rpm_estimate;
			}
		}
	}

	float dt = 0.02f;

	if (_last_ctrl_time_us != 0) {
		dt = (now_us - _last_ctrl_time_us) * 1e-6f;
		dt = math::constrain(dt, 0.002f, 0.05f);
	}

	_last_ctrl_time_us = now_us;

	const float tk = now_us * 1e-6f;

	const bool armed = _armed.armed;

	if (!armed) {
		_xc.alpha_cf = 0.f;
		_xc.alpha_sf = 0.f;
		_xc.eta_u = 0.f;

		_last_omega1_cmd = 0.f;
		_last_omega2_cmd = 0.f;

		_attitude_disarm_sent = false;

		actuator_motors_s motors{};
		motors.timestamp = now_us;
		motors.timestamp_sample = now_us;
		motors.reversible_flags = 0;

		for (int i = 0; i < actuator_motors_s::NUM_CONTROLS; i++) {
			motors.control[i] = NAN;
		}

		_actuator_motors_pub.publish(motors);
		return;
	}

	matrix::Quatf q_att(_att.q);
	const matrix::Eulerf euler(q_att);

	// 滚转角 / 俯仰角超过 30 度时，立即 disarm 并停止两个电机
	const float attitude_disarm_limit_rad = 20.0f * M_PI_F / 180.0f;

	const float roll_abs = fabsf(euler.phi());
	const float pitch_abs = fabsf(euler.theta());

	if ((roll_abs > attitude_disarm_limit_rad) || (pitch_abs > attitude_disarm_limit_rad)) {

		if (!_attitude_disarm_sent) {
			vehicle_command_s vcmd{};
			vcmd.timestamp = now_us;

			// MAV_CMD_COMPONENT_ARM_DISARM
			// param1 = 0 表示 disarm
			vcmd.command = vehicle_command_s::VEHICLE_CMD_COMPONENT_ARM_DISARM;
			vcmd.param1 = 0.0f;
			vcmd.param2 = 0.0f;

			vcmd.target_system = 1;
			vcmd.target_component = 1;
			vcmd.source_system = 1;
			vcmd.source_component = 1;
			vcmd.from_external = false;

			_vehicle_command_pub.publish(vcmd);

		//	_attitude_disarm_sent = true;

			PX4_WARN("attitude limit exceeded: roll=%.1f deg pitch=%.1f deg, disarm",
				(double)(euler.phi() * 180.0f / M_PI_F),
				(double)(euler.theta() * 180.0f / M_PI_F));
		}

		// 立即清零控制器内部状态
		_xc.alpha_cf = 0.f;
		_xc.alpha_sf = 0.f;
		_xc.eta_u = 0.f;

		_last_omega1_cmd = 0.f;
		_last_omega2_cmd = 0.f;

		// 立即停止本模块发给两个电机的输出
		actuator_motors_s motors{};
		motors.timestamp = now_us;
		motors.timestamp_sample = now_us;
		motors.reversible_flags = 0;

		for (int i = 0; i < actuator_motors_s::NUM_CONTROLS; i++) {
			motors.control[i] = NAN;
		}

		motors.control[0] = 0.f;
		motors.control[1] = 0.f;

		_actuator_motors_pub.publish(motors);

		return;
	}

	const float p = _ang_vel.xyz[0];
	const float q = _ang_vel.xyz[1];
	const float r = _ang_vel.xyz[2];

	const bool rpm1_valid =
		PX4_ISFINITE(_motor1_rpm_meas) && (_motor1_rpm_meas >= _min_rpm_for_feedback);

	const bool rpm2_valid =
		PX4_ISFINITE(_motor2_rpm_meas) && (_motor2_rpm_meas >= _min_rpm_for_feedback);

	PlantState xp{};
	xp.phi = euler.phi();
	xp.theta = euler.theta();
	xp.psi = euler.psi();

	xp.p = p;
	xp.q = q;
	xp.r = r;

	xp.omega1 = rpm1_valid ? (_motor1_rpm_meas * 2.f * M_PI_F / 60.f) : _last_omega1_cmd;
	xp.omega2 = rpm2_valid ? (_motor2_rpm_meas * 2.f * M_PI_F / 60.f) : _last_omega2_cmd;

	ControllerCommand cmd{};

	if (_manual.valid) {
		// 左侧摇杆二维位置 -> 合力方向
		// 默认用 manual.roll / manual.pitch。
		// 若你的物理左摇杆不是 roll/pitch，需要在 QGC Radio Setup 中重新映射。
		const float stick_x = math::constrain(_manual.yaw, -1.f, 1.f);
		const float stick_y = math::constrain(_manual.pitch, -1.f, 1.f);

		const float tilt_x = stick_y * _max_tilt_angle;
		const float tilt_y = -stick_x * _max_tilt_angle;

		const float x = tanf(tilt_x);
		const float y = tanf(tilt_y);

		matrix::Vector3f zd_cmd_E(x, y, 1.f);
		const float n = zd_cmd_E.norm();

		if (n > 1e-6f) {
			zd_cmd_E /= n;

		} else {
			zd_cmd_E = matrix::Vector3f(0.f, 0.f, 1.f);
		}

		cmd.zd_cmd_E = zd_cmd_E;

		// 右侧摇杆上下 -> 合力大小
		// throttle: [-1, 1] -> [0, 1]
		const float thr01 = 0.5f * (math::constrain(_manual.throttle, -1.f, 1.f) + 1.f);

		const float u1_min = _u1_min_factor * _m * _g;
		const float u1_max = _u1_max_factor * _m * _g;

		cmd.u1_ref = u1_min + thr01 * (u1_max - u1_min);

	} else {
		cmd.zd_cmd_E = matrix::Vector3f(0.f, 0.f, 1.f);
		cmd.u1_ref = _m * _g;
	}

	// zd_dot / zd_ddot 数值微分 + 低通
	const matrix::Vector3f raw_zd_dot = (cmd.zd_cmd_E - _zd_prev_E) / dt;

	const float tau_dot = 1.f / (2.f * M_PI_F * math::max(_zd_dot_lpf_hz, 0.1f));
	const float alpha_dot = dt / (tau_dot + dt);

	const matrix::Vector3f zd_dot_f =
		_zd_dot_prev_E + alpha_dot * (raw_zd_dot - _zd_dot_prev_E);

	const matrix::Vector3f raw_zd_ddot = (zd_dot_f - _zd_dot_prev_E) / dt;

	const float tau_ddot = 1.f / (2.f * M_PI_F * math::max(_zd_ddot_lpf_hz, 0.1f));
	const float alpha_ddot = dt / (tau_ddot + dt);

	const matrix::Vector3f zd_ddot_f =
		_zd_ddot_prev_E + alpha_ddot * (raw_zd_ddot - _zd_ddot_prev_E);

	cmd.zd_dot_cmd_E = zd_dot_f;
	cmd.zd_ddot_cmd_E = zd_ddot_f;

	_zd_prev_E = cmd.zd_cmd_E;
	_zd_dot_prev_E = zd_dot_f;
	_zd_ddot_prev_E = zd_ddot_f;

	ControllerOutput ctrl_out{};
	controllerStep(tk, dt, xp, cmd, _xc, ctrl_out);

	const float omega1_cmd = math::constrain(ctrl_out.omega1_cmd, 0.f, _omega_max);
	const float omega2_cmd = math::constrain(ctrl_out.omega2_cmd, 0.f, _omega_max);

	_last_omega1_cmd = omega1_cmd;
	_last_omega2_cmd = omega2_cmd;

	// actuator_motors 语义是归一化推力，不是归一化转速，电调已经做了相应映射关系可以发归一化转速命令
	// const float single_motor_thrust_max = 0.5f * _u1_max_factor * _m * _g;

	// const float motor1_norm =
	// 	math::constrain(ctrl_out.T1_cmd / single_motor_thrust_max, 0.f, 1.f);

	// const float motor2_norm =
	// 	math::constrain(ctrl_out.T2_cmd / single_motor_thrust_max, 0.f, 1.f);

	const float motor1_norm =
		math::constrain(ctrl_out.omega1_cmd / _omega_max, 0.f, 1.f);

	const float motor2_norm =
		math::constrain(ctrl_out.omega2_cmd / _omega_max, 0.f, 1.f);
	actuator_motors_s motors{};
	motors.timestamp = now_us;
	motors.timestamp_sample = now_us;
	motors.reversible_flags = 0;

	for (int i = 0; i < actuator_motors_s::NUM_CONTROLS; i++) {
		motors.control[i] = NAN;
	}

	motors.control[0] = motor1_norm;
	motors.control[1] = motor2_norm;

	_actuator_motors_pub.publish(motors);

	//发布控制器内部数据
	my_controller_status_s st{};
	st.timestamp = now_us;



	st.phi = xp.phi;
	st.theta = xp.theta;
	st.psi = xp.psi;

	st.p = xp.p;
	st.q = xp.q;
	st.r = xp.r;

	st.omega1_meas = xp.omega1;
	st.omega2_meas = xp.omega2;

	st.zd_x = cmd.zd_cmd_E(0);
	st.zd_y = cmd.zd_cmd_E(1);
	st.zd_z = cmd.zd_cmd_E(2);

	st.zd_dot_x = cmd.zd_dot_cmd_E(0);
	st.zd_dot_y = cmd.zd_dot_cmd_E(1);
	st.zd_dot_z = cmd.zd_dot_cmd_E(2);

	st.zd_ddot_x = cmd.zd_ddot_cmd_E(0);
	st.zd_ddot_y = cmd.zd_ddot_cmd_E(1);
	st.zd_ddot_z = cmd.zd_ddot_cmd_E(2);

	st.u1_ref = cmd.u1_ref;


	st.alpha_cf = _xc.alpha_cf;
	st.alpha_sf = _xc.alpha_sf;
	st.eta_u = _xc.eta_u;

	st.p_d = ctrl_out.p_d;
	st.q_d = ctrl_out.q_d;
	st.tau_x_d = ctrl_out.tau_x_d;
	st.tau_y_d = ctrl_out.tau_y_d;
	st.u1_cmd = ctrl_out.u1_cmd;
	st.delta = ctrl_out.delta;

	st.omega1_cmd = omega1_cmd;
	st.omega2_cmd = omega2_cmd;

   _my_controller_status_pub.publish(st);
}

int MyController::task_spawn(int argc, char *argv[])
{
	MyController *instance = new MyController();

	if (instance) {
		_object.store(instance);
		_task_id = task_id_is_work_queue;

		if (instance->init()) {
			return PX4_OK;
		}

	} else {
		PX4_ERR("alloc failed");
	}

	delete instance;
	_object.store(nullptr);
	_task_id = -1;

	return PX4_ERROR;
}

int MyController::print_status()
{
	PX4_INFO("my_controller running");
	PX4_INFO("rpm1=%.1f rpm2=%.1f",
		 (double)_motor1_rpm_meas,
		 (double)_motor2_rpm_meas);

	PX4_INFO("alpha_cf=%.4f alpha_sf=%.4f eta_u=%.4f",
		 (double)_xc.alpha_cf,
		 (double)_xc.alpha_sf,
		 (double)_xc.eta_u);

	return 0;
}

int MyController::custom_command(int argc, char *argv[])
{
	return print_usage("unknown command");
}

int MyController::print_usage(const char *reason)
{
	if (reason) {
		PX4_WARN("%s\n", reason);
	}

	PRINT_MODULE_DESCRIPTION(
		R"DESCR_STR(
### Description
Single-axis birotor controller.

Manual mapping:
- manual.roll / manual.pitch -> desired force direction
- manual.throttle -> total thrust magnitude, 0.5mg to 1.5mg

Output:
- publish actuator_motors[0:1] as normalized thrust commands
)DESCR_STR");

	PRINT_MODULE_USAGE_NAME("my_controller", "controller");
	PRINT_MODULE_USAGE_COMMAND("start");
	PRINT_MODULE_USAGE_DEFAULT_COMMANDS();

	return 0;
}

extern "C" __EXPORT int my_controller_main(int argc, char *argv[])
{
	return MyController::main(argc, argv);
}
