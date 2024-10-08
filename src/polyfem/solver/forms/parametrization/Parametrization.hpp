#pragma once

#include <memory>
#include <vector>

#include <Eigen/Core>

namespace polyfem::solver
{
	/** This parameterize a function f : x -> y
	 * and provides the chain rule with respect to previous gradients
	 */
	class Parametrization
	{
	public:
		Parametrization() {}
		virtual ~Parametrization() {}

		virtual Eigen::VectorXd inverse_eval(const Eigen::VectorXd &y);

		virtual int size(const int x_size) const = 0; // just for verification
		virtual Eigen::VectorXd eval(const Eigen::VectorXd &x) const = 0;
		virtual Eigen::VectorXd apply_jacobian(const Eigen::VectorXd &grad_full, const Eigen::VectorXd &x) const = 0;

        virtual bool has_input_mesh() const { return false; }
        virtual std::string input_path() const { return ""; }
        virtual std::string output_path() const { return ""; }
	};

	class CompositeParametrization : public Parametrization
	{
	public:
		CompositeParametrization() : parametrizations_(std::vector<std::shared_ptr<Parametrization>>()) {}
		CompositeParametrization(std::vector<std::shared_ptr<Parametrization>>&& parametrizations) : parametrizations_(parametrizations) {}
		virtual ~CompositeParametrization() {}

		Eigen::VectorXd inverse_eval(const Eigen::VectorXd &y) override;

		int size(const int x_size) const override;
		Eigen::VectorXd eval(const Eigen::VectorXd &x) const override;
		Eigen::VectorXd apply_jacobian(const Eigen::VectorXd &grad_full, const Eigen::VectorXd &x) const override;

		virtual bool has_input_mesh() const override { return parametrizations_.size() > 0 && parametrizations_[0]->has_input_mesh(); }
        virtual std::string input_path() const override { return parametrizations_.size() == 0 ? "" : parametrizations_[0]->input_path(); }
        virtual std::string output_path() const override { return parametrizations_.size() == 0 ? "" : parametrizations_.back()->output_path(); }

	private:
		const std::vector<std::shared_ptr<Parametrization>> parametrizations_;
	};
} // namespace polyfem::solver