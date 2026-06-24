#pragma once

#include <cmath>

#include <px4_platform_common/module.h>
#include <px4_platform_common/module_params.h>
#include <px4_platform_common/px4_work_queue/ScheduledWorkItem.hpp>

#include <uORB/Publication.hpp>
#include <uORB/Subscription.hpp>

#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/manual_control_setpoint.h>
#include <uORB/topics/vehicle_angular_velocity.h>
#include <uORB/topics/rpm.h>
#include <uORB/topics/esc_status.h>
#include <uORB/topics/actuator_motors.h>
#include <uORB/topics/actuator_armed.h>
#include <uORB/topics/my_controller_status.h>
#include <uORB/topics/vehicle_command.h>
#include <matrix/matrix/math.hpp>

class MyController :
	public ModuleBase<MyController>,
	public ModuleParams,
	public px4::ScheduledWorkItem
{
public:
	MyController();
	~MyController() override = default;

	static int task_spawn(int argc, char *argv[]);
	static int custom_command(int argc, char *argv[]);
	static int print_usage(const char *reason = nullptr);

	bool init();
	int print_status() override;

private:
	void Run() override;

	struct PlantState {
		float phi{0.f};
		float theta{0.f};
		float psi{0.f};
		float p{0.f};
		float q{0.f};
		float r{0.f};
		float omega1{0.f};
		float omega2{0.f};
	};

	struct ControllerCommand {
		matrix::Vector3f zd_cmd_E{0.f, 0.f, 1.f};
		matrix::Vector3f zd_dot_cmd_E{0.f, 0.f, 0.f};
		matrix::Vector3f zd_ddot_cmd_E{0.f, 0.f, 0.f};
		float u1_ref{0.f};
	};

	struct ControllerState {
		float alpha_cf{0.f};
		float alpha_sf{0.f};
		float eta_u{0.f};
	};

	struct ControllerOutput {
		float omega1_cmd{0.f};
		float omega2_cmd{0.f};

		float T1_cmd{0.f};
		float T2_cmd{0.f};

		float p_d{0.f};
		float q_d{0.f};
		float tau_x_d{0.f};
		float tau_y_d{0.f};
		float u1_cmd{0.f};
		float delta{0.f};
	};

	void controllerStep(float tk, float dt,
			    const PlantState &xp,
			    const ControllerCommand &cmd,
			    ControllerState &xc,
			    ControllerOutput &out);

	float yawDrag(float r) const;
	static float signNZ(float x);

	// ---------- uORB ----------
	uORB::Subscription _att_sub{ORB_ID(vehicle_attitude)};
	uORB::Subscription _ang_vel_sub{ORB_ID(vehicle_angular_velocity)};
	uORB::Subscription _manual_sub{ORB_ID(manual_control_setpoint)};
	uORB::Subscription _esc_status_sub{ORB_ID(esc_status)};
	uORB::Subscription _rpm_sub{ORB_ID(rpm)};
	uORB::Subscription _armed_sub{ORB_ID(actuator_armed)};


	uORB::Publication<actuator_motors_s> _actuator_motors_pub{ORB_ID(actuator_motors)};
        uORB::Publication<my_controller_status_s> _my_controller_status_pub{ORB_ID(my_controller_status)};
	uORB::Publication<vehicle_command_s> _vehicle_command_pub{ORB_ID(vehicle_command)};

	bool _attitude_disarm_sent{false};

	// ---------- cache ----------
	vehicle_attitude_s _att{};
	vehicle_angular_velocity_s _ang_vel{};
	manual_control_setpoint_s _manual{};
	actuator_armed_s _armed{};

	rpm_s _rpm{};

	float _motor1_rpm_meas{NAN};
	float _motor2_rpm_meas{NAN};

	float _last_omega1_cmd{0.f};
	float _last_omega2_cmd{0.f};

	uint64_t _last_ctrl_time_us{0};

	ControllerState _xc{};

	matrix::Vector3f _zd_prev_E{0.f, 0.f, 1.f};
	matrix::Vector3f _zd_dot_prev_E{0.f, 0.f, 0.f};
	matrix::Vector3f _zd_ddot_prev_E{0.f, 0.f, 0.f};

	// ---------- parameters ----------
	// float _m{1.42668381f};
	float _m{1.62468381f};
	float _g{9.81f};
	float _l{0.17f};

	// float _Jx{0.00212844f};
	// float _Jy{0.01906274f};
	// float _Jz{0.01775539f};
	float _Jx{0.00459844f};
	float _Jy{0.02153274f};
	float _Jz{0.01974568f};
	float _Jr{1.1e-4f};

	float _b{2.1e-5f};
	float _kq{3.85e-7f};
	float _d1{5.0e-3f};
	float _d2{8.0e-4f};

	float _k_c1{3.0f};
	float _k_c2{3.0f};
	float _k_p{3.5f};
	float _k_q{3.5f};

	float _lambda_alpha{3.0f};

	float _k_up{0.0f};
	float _k_ui{0.0f};
	float _k_aw{0.0f};

	float _u1_min_margin{0.10f};
	float _t_spinup{1.0f};

	float _c3_eps{1e-3f};
	float _det_eps{1e-6f};
	float _b3z_eps{1e-3f};

	float _omega_max{1055.57508f};

	// ---------- manual mapping ----------
	float _max_tilt_angle{0.35f};     // rad, 左摇杆最大倾斜角，20度
	float _u1_min_factor{0.0f};       // 右摇杆最低：0.0 mg
	float _u1_max_factor{3.0f};       // 右摇杆最高：2*_omega_max^2 * b / mg

	float _zd_dot_lpf_hz{8.0f};
	float _zd_ddot_lpf_hz{8.0f};

	// 如果 ESC 没有 RPM 回传，则用上一周期命令估计 omega
	float _min_rpm_for_feedback{300.f};

	bool _use_vertical_comp{false};
};
