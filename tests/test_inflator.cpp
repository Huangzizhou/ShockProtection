///////////////////////////////////////////////////////////////////////////////
#include <polyfem/assembler/AssemblerUtils.hpp>
#include <iostream>
#include <fstream>
#include <cmath>

#include <polyfem/State.hpp>
#include <polyfem/solver/Optimizations.hpp>
#include <polyfem/solver/AdjointTools.hpp>
#include <polyfem/solver/forms/adjoint_forms/VariableToSimulation.hpp>
#include <polyfem/solver/forms/adjoint_forms/SumCompositeForm.hpp>
#include <polyfem/solver/forms/adjoint_forms/AdjointForm.hpp>
#include <polyfem/solver/AdjointNLProblem.hpp>

#include <paraviewo/VTUWriter.hpp>
#include <finitediff.hpp>

#include <catch2/catch_all.hpp>
#include <math.h>
////////////////////////////////////////////////////////////////////////////////

using namespace polyfem;
using namespace solver;

namespace
{

	bool load_json(const std::string &json_file, json &out)
	{
		std::ifstream file(json_file);

		if (!file.is_open())
			return false;

		file >> out;

		out["root_path"] = json_file;

		return true;
	}

	void verify_adjoint(AdjointNLProblem& problem, const Eigen::VectorXd &x, const double dt, const double tol)
	{
		problem.solution_changed(x);
		problem.save_to_file(0, x);
		double functional_val = problem.value(x);

		Eigen::VectorXd one_form;
		problem.gradient(x, one_form);
		Eigen::VectorXd theta = one_form.normalized();
		double derivative = (one_form.array() * theta.array()).sum();

		problem.solution_changed(x + theta * dt);
		double next_functional_val = problem.value(x + theta * dt);

		problem.solution_changed(x - theta * dt);
		double former_functional_val = problem.value(x - theta * dt);

		double finite_difference = (next_functional_val - former_functional_val) / dt / 2;
		std::cout << std::setprecision(16) << "f(x) " << functional_val << " f(x-dt) " << former_functional_val << " f(x+dt) " << next_functional_val << "\n";
		std::cout << std::setprecision(12) << "derivative: " << derivative << ", fd: " << finite_difference << "\n";
		std::cout << std::setprecision(12) << "relative error: " << abs((finite_difference - derivative) / derivative) << "\n";
		REQUIRE(derivative == Catch::Approx(finite_difference).epsilon(tol));
	}

	void verify_adjoint_expensive(AdjointNLProblem& problem, const Eigen::VectorXd &x, const double dt)
	{
		problem.solution_changed(x);
		double functional_val = problem.value(x);

		Eigen::VectorXd analytic;
		problem.gradient(x, analytic);

		std::cout << std::setprecision(12) << "derivative: " << analytic.transpose() << "\n";

		Eigen::VectorXd fd;
		fd.setZero(x.size());
		for (int d = 0; d < x.size(); d++)
		{
			Eigen::VectorXd theta;
			theta.setZero(x.size());
			theta(d) = 1;

			problem.solution_changed(x + theta * dt);
			double next_functional_val = problem.value(x + theta * dt);

			problem.solution_changed(x - theta * dt);
			double former_functional_val = problem.value(x - theta * dt);

			fd(d) = (next_functional_val - former_functional_val) / dt / 2;
		}

		std::cout << "fd: " << fd.transpose() << "\n";
	}

} // namespace

TEST_CASE("isosurface-inflator-periodic", "[test_adjoint]")
{
	const std::string path = POLYFEM_DATA_DIR + std::string("/../tests/differentiable/isosurface-inflator-periodic");
	json in_args;
	load_json(path + "/state.json", in_args);

	json opt_args;
	load_json(path + "/opt.json", opt_args);
	opt_args = AdjointOptUtils::apply_opt_json_spec(opt_args, false);

	std::shared_ptr<State> state_ptr = AdjointOptUtils::create_state(in_args, solver::CacheLevel::Derivatives, -1);
	State &state = *state_ptr;

	std::vector<std::shared_ptr<State>> states({state_ptr});

	VariableToSimulationGroup variable_to_simulations;
	for (auto &tmp_args : opt_args["variable_to_simulation"])
		variable_to_simulations.push_back(AdjointOptUtils::create_variable_to_simulation(tmp_args, states, {}));

	auto obj = std::dynamic_pointer_cast<SumCompositeForm>(AdjointOptUtils::create_form(opt_args["functionals"], variable_to_simulations, states));

	Eigen::VectorXd x = opt_args["parameters"][0]["initial"];

	auto nl_problem = std::make_shared<AdjointNLProblem>(obj, variable_to_simulations, states, opt_args);
	// nl_problem->solution_changed(x);

	verify_adjoint(*nl_problem, x, opt_args["solver"]["nonlinear"]["debug_fd_eps"].get<double>(), 1e-3);
	// verify_adjoint_expensive(*nl_problem, x, opt_args["solver"]["nonlinear"]["debug_fd_eps"].get<double>());
}

TEST_CASE("isosurface-inflator", "[test_adjoint]")
{
	const std::string path = POLYFEM_DATA_DIR + std::string("/../tests/differentiable/isosurface-inflator");
	json in_args;
	load_json(path + "/state.json", in_args);

	json opt_args;
	load_json(path + "/opt.json", opt_args);
	opt_args = AdjointOptUtils::apply_opt_json_spec(opt_args, false);

	std::shared_ptr<State> state_ptr = AdjointOptUtils::create_state(in_args, solver::CacheLevel::Derivatives, -1);
	State &state = *state_ptr;

	std::vector<std::shared_ptr<State>> states({state_ptr});

	VariableToSimulationGroup variable_to_simulations;
	variable_to_simulations.push_back(AdjointOptUtils::create_variable_to_simulation(opt_args["variable_to_simulation"][0], states, {}));

	auto obj = std::dynamic_pointer_cast<SumCompositeForm>(AdjointOptUtils::create_form(opt_args["functionals"], variable_to_simulations, states));

	Eigen::VectorXd x = opt_args["parameters"][0]["initial"];

	auto nl_problem = std::make_shared<AdjointNLProblem>(obj, variable_to_simulations, states, opt_args);
	// nl_problem->solution_changed(x);

	verify_adjoint(*nl_problem, x, opt_args["solver"]["nonlinear"]["debug_fd_eps"].get<double>(), 1e-3);
	// verify_adjoint_expensive(variable_to_simulations, *obj, state, x, opt_args["solver"]["nonlinear"]["debug_fd_eps"].get<double>());
}
