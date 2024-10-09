#include "VariableToSimulation.hpp"
#include <polyfem/State.hpp>
#include <polyfem/mesh/GeometryReader.hpp>
#include <polyfem/solver/forms/parametrization/PeriodicMeshToMesh.hpp>

namespace polyfem::solver
{
	namespace {
        bool file_exists(const std::string& name) {
            std::ifstream f(name.c_str());
            return f.good();
        }
	}

	SDFShapeVariableToSimulation::SDFShapeVariableToSimulation(const std::vector<std::shared_ptr<State>> &states, const CompositeParametrization &parametrization) : ShapeVariableToSimulation(states, parametrization)
	{
	}
	void SDFShapeVariableToSimulation::update(const Eigen::VectorXd &x)
	{
		parametrization_.eval(x);
		const std::string mesh_path = parametrization_.output_path();
		if (!file_exists(mesh_path))
			log_and_throw_error("Invalid mesh path {} in SDFShapeVariableToSimulation!", mesh_path);

		int start = 0, end = 0; // start vertex index of the mesh
		for (auto state : states_)
		{
			state->args["geometry"][mesh_id_]["mesh"] = mesh_path;
			state->args["geometry"][mesh_id_].erase("transformation");

			state->mesh.reset();
			state->assembler->update_lame_params(Eigen::MatrixXd(), Eigen::MatrixXd());

			json in_args = state->args;
			if (in_args.contains("time"))
			{
				if (!state->problem->is_time_dependent())
					in_args.erase("time");
				else
					in_args["time"].erase("dt");
			}
			state->init(in_args, false);

			{
				assert(state->args["geometry"].is_array());
				auto geometries = state->args["geometry"].get<std::vector<json>>();

				int i = 0;
				for (const json &geometry : geometries)
				{
					if (!geometry["enabled"].get<bool>() || geometry["is_obstacle"].get<bool>())
						continue;

					if (geometry["type"] != "mesh")
						log_and_throw_error(
							fmt::format("Invalid geometry type \"{}\" for FEM mesh!", geometry["type"]));

					if (i == mesh_id_)
						start = state->mesh ? state->mesh->n_vertices() : 0;

					if (state->mesh == nullptr)
						state->mesh = mesh::read_fem_mesh(state->units, geometry, state->args["root_path"], false);
					else
						state->mesh->append(mesh::read_fem_mesh(state->units, geometry, state->args["root_path"], false));

					if (i == mesh_id_)
						end = state->mesh->n_vertices();

					i++;
				}
			}

			state->load_mesh();
			state->stats.compute_mesh_stats(*state->mesh);
			// state->build_basis();

			// state->set_log_level(static_cast<spdlog::level::level_enum>(cur_log));
		}

		const int dim = states_[0]->mesh->dimension();
		json tmp = json({});
		tmp["composite_map_indices"] = Eigen::VectorXi::LinSpaced((end - start) * dim, start * dim, end * dim - 1);
		tmp["composite_map_type"] = "indices";
		set_output_indexing(tmp);
	}

	SDFPeriodicShapeVariableToSimulation::SDFPeriodicShapeVariableToSimulation(const std::vector<std::shared_ptr<State>> &states, const CompositeParametrization &parametrization) : PeriodicShapeVariableToSimulation(states, parametrization)
	{
        adjoint_logger().warn("SDFPeriodicShapeVariableToSimulation only supports mesh with periodic nodes on the bounding box and unit size!");
	}
	void SDFPeriodicShapeVariableToSimulation::update(const Eigen::VectorXd &x)
	{
		auto y = parametrization_.eval(x);
		const std::string mesh_path = parametrization_.output_path();
		if (!file_exists(mesh_path))
			log_and_throw_error("Invalid mesh path {} in SDFPeriodicShapeVariableToSimulation!", mesh_path);

		for (auto state : states_)
		{
			if (state->args["geometry"].is_array())
			{
				assert(state->args["geometry"].size() == 1);
				state->args["geometry"][0]["mesh"] = mesh_path;
				state->args["geometry"][0].erase("transformation");
			}
			else
			{
				state->args["geometry"]["mesh"] = mesh_path;
				state->args["geometry"].erase("transformation");
			}

			state->mesh.reset();
			state->assembler->update_lame_params(Eigen::MatrixXd(), Eigen::MatrixXd());

			json in_args = state->args;
			if (in_args.contains("time") && !state->problem->is_time_dependent())
				in_args.erase("time");
			state->init(in_args, false);

			state->load_mesh();
			state->stats.compute_mesh_stats(*state->mesh);

			{
				Eigen::MatrixXd V(state->mesh->n_vertices(), state->mesh->dimension());
				for (int i = 0; i < state->mesh->n_vertices(); i++)
					V.row(i) = state->mesh->point(i);
				state->periodic_mesh_map = std::make_unique<PeriodicMeshToMesh>(V);
				state->periodic_mesh_representation = state->periodic_mesh_map->inverse_eval(utils::flatten(V));
			}
		}
	}
    Eigen::VectorXd SDFPeriodicShapeVariableToSimulation::apply_parametrization_jacobian(const Eigen::VectorXd &term, const Eigen::VectorXd &x) const
    {
		const auto &state = *(states_[0]);
        const int dim = state.periodic_mesh_map->dim();
        const int n_full_verts = state.periodic_mesh_map->n_full_dof();
        const int n_periodic_verts = state.periodic_mesh_map->n_periodic_dof();

		assert(term.size() == state.periodic_mesh_map->input_size());
        Eigen::VectorXd full_term;
        full_term.setZero(n_full_verts * dim);
        Eigen::Matrix<bool, -1, 1> visited_mask;
        visited_mask.setZero(n_periodic_verts);
        for (int i = 0; i < n_full_verts; i++)
        {
            int i_periodic = state.periodic_mesh_map->full_to_periodic(i);
            if (!visited_mask(i_periodic))
                full_term.segment(i * dim, dim) = term.segment(i_periodic * dim, dim);
            visited_mask(i_periodic) = true;
        }
        
        return parametrization_.apply_jacobian(full_term, x);
	}

	PeriodicShapeScaleVariableToSimulation::PeriodicShapeScaleVariableToSimulation(const std::vector<std::shared_ptr<State>> &states, const CompositeParametrization &parametrization) : PeriodicShapeVariableToSimulation(states, parametrization)
	{
		for (const auto &state : states)
		{
			if (!state->args["boundary_conditions"]["periodic_boundary"]["enabled"])
				log_and_throw_error("PeriodicShapeScaleVariableToSimulation is designed for periodic bc!");
			dim = state->mesh->dimension();
		}
	}
	void PeriodicShapeScaleVariableToSimulation::update(const Eigen::VectorXd &x)
	{
		Eigen::VectorXd y = parametrization_.eval(x);
		assert(y.size() == dim);
		Eigen::MatrixXd affine = y.asDiagonal();

		adjoint_logger().debug("current mesh scale: {:.15e} {:.15e}", y(0), y(1));
		
		for (auto state : states_)
		{
			assert(state->periodic_mesh_representation.size() > dim * dim);
			state->periodic_mesh_representation.tail(dim * dim) = Eigen::Map<Eigen::VectorXd>(affine.data(), dim * dim, 1);
			auto V = utils::unflatten(state->periodic_mesh_map->eval(state->periodic_mesh_representation), dim);
			for (int i = 0; i < V.rows(); i++)
				state->set_mesh_vertex(i, V.row(i));
		}
	}
	Eigen::VectorXd PeriodicShapeScaleVariableToSimulation::apply_parametrization_jacobian(const Eigen::VectorXd &term, const Eigen::VectorXd &x) const
	{
		assert(term.size() > dim * dim);
		Eigen::MatrixXd affine_term = utils::unflatten(term.tail(dim * dim), dim).transpose();
		return parametrization_.apply_jacobian(affine_term.diagonal(), x);
	}

	PeriodicShapeAffineVariableToSimulation::PeriodicShapeAffineVariableToSimulation(const std::vector<std::shared_ptr<State>> &states, const CompositeParametrization &parametrization) : PeriodicShapeVariableToSimulation(states, parametrization)
	{
		for (const auto &state : states)
		{
			if (!state->args["boundary_conditions"]["periodic_boundary"]["enabled"])
				log_and_throw_error("PeriodicShapeAffineVariableToSimulation is designed for periodic bc!");
			dim = state->mesh->dimension();
		}
	}
	void PeriodicShapeAffineVariableToSimulation::update(const Eigen::VectorXd &x)
	{
		Eigen::VectorXd y = parametrization_.eval(x);
		assert(y.size() == dim*dim);

		adjoint_logger().info("mesh affine: {}", y.transpose());
		
		for (auto state : states_)
		{
			assert(state->periodic_mesh_representation.size() > dim * dim);
			state->periodic_mesh_representation.tail(dim * dim) = y;
		
			auto V = utils::unflatten(state->periodic_mesh_map->eval(state->periodic_mesh_representation), dim);
			for (int i = 0; i < V.rows(); i++)
				state->set_mesh_vertex(i, V.row(i));
		}
	}
	Eigen::VectorXd PeriodicShapeAffineVariableToSimulation::apply_parametrization_jacobian(const Eigen::VectorXd &term, const Eigen::VectorXd &x) const
	{
		assert(x.size() > dim * dim);
		return parametrization_.apply_jacobian(term.tail(dim * dim), x);
	}
} // namespace polyfem::solver
