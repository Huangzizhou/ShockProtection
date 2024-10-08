#include "SDFParametrizations.hpp"

#include <polyfem/solver/AdjointTools.hpp>

#include <polyfem/io/MshReader.hpp>
#include <polyfem/utils/Timer.hpp>
#include <igl/writeMSH.h>
#include <paraviewo/VTUWriter.hpp>

#include <polyfem/State.hpp>
#include <polysolve/linear/FEMSolver.hpp>

namespace polyfem::solver
{
    namespace {
        void write_msh(const std::string &path, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
        {
            Eigen::MatrixXd Vsave;
            Vsave.setZero(V.rows(), 3);
            Vsave.leftCols(V.cols()) = V;

            const int dim = V.cols();
            Eigen::MatrixXi Tri = (dim == 3) ? Eigen::MatrixXi() : F;
            Eigen::MatrixXi Tet = (dim == 3) ? F : Eigen::MatrixXi();

            igl::writeMSH(path, Vsave, Tri, Tet, Eigen::MatrixXi::Zero(Tri.rows(), 1), Eigen::MatrixXi::Zero(Tet.rows(), 1), std::vector<std::string>(), std::vector<Eigen::MatrixXd>(), std::vector<std::string>(), std::vector<Eigen::MatrixXd>(), std::vector<Eigen::MatrixXd>());
        }
    }

    bool SDF2Mesh::isosurface_inflator(const Eigen::VectorXd &x) const
    {
        auto y = x;
        if (inactive_shape_params.size() > 0)
            for (const auto& tuple : inactive_shape_params)
                y(std::get<0>(tuple)) = std::get<1>(tuple);
        
        if (last_x.size() == x.size() && last_x == x)
            return true;

        {
            std::vector<double> x_vec(y.data(), y.data() + y.size());
            {
                POLYFEM_SCOPED_TIMER("mesh inflation");
                Eigen::MatrixXd vertex_normals, shape_vel;
                inflator.inflate(x_vec, Vout, Fout, vertex_normals, shape_vel);
                Vout /= 2.0;

                Eigen::VectorXd norms = vertex_normals.rowwise().norm();
                boundary_flags.setZero(norms.size());
                for (int i = 0; i < norms.size(); i++)
                    if (norms(i) > 0.1)
                        boundary_flags(i) = true;
                
                shape_velocity.setZero(shape_vel.rows(), vertex_normals.size());
                for (int d = 0; d < dim_; d++)
                    for (int i = 0; i < vertex_normals.rows(); i++)
                        shape_velocity.col(dim_ * i + d) = shape_vel.col(i) * vertex_normals(i, d) / 2.0;
                
                if (inactive_shape_params.size() == 0)
                {
                    for (int q = 0; q < x.size(); q++)
                        if (shape_velocity.row(q).norm() < 1e-8)
                            inactive_shape_params.push_back(std::make_tuple(q, y(q)));
                }
            }
            
            write_msh(out_path_, Vout, Fout);
            if (volume_velocity_)
                extend_to_internal();
        }

        last_x = x;

        return true;
    }
    void SDF2Mesh::extend_to_internal() const
    {
        json args = R"(
        {
            "geometry": [
                {
                    "mesh": "",
                    "surface_selection": {
                        "threshold": 1e-7
                    }
                }
            ],
            "space": {
                "discr_order": 1
            },
            "solver": {
                "linear": {
                    "solver": "Eigen::PardisoLDLT"
                }
            },
            "boundary_conditions": {
                "dirichlet_boundary": [
                    {
                        "id": 7,
                        "value": 0.0
                    }
                ],
                "periodic_boundary": {
                    "enabled": true
                },
            },
            "output": {
                "log": {
                    "level": "info"
                }
            },
            "materials": {
                "type": "Laplacian"
            }
        })"_json;
        
        args["geometry"][0]["mesh"] = out_path_;

        const int dim = Vout.cols();

        Eigen::MatrixXd extended_velocity;
        extended_velocity.setZero(shape_velocity.rows(), shape_velocity.cols());

        State state;
        state.init(args, false);

        adjoint_logger().info("Laplacian solve for volume shape velocity...");

        state.load_mesh();

        if (state.mesh == nullptr)
            log_and_throw_error("Invalid mesh!");
        
        state.stats.compute_mesh_stats(*state.mesh);

        state.build_basis();

        state.assemble_rhs();
        state.assemble_mass_mat();

        state.boundary_nodes.clear();
        std::vector<int> primitive_to_node = state.primitive_to_node();
        std::vector<int> node_to_primitive = state.node_to_primitive();
        for (int i = 0; i < boundary_flags.size(); i++)
            if (boundary_flags(i))
                state.boundary_nodes.push_back(primitive_to_node[i]);

        auto solver = polysolve::linear::Solver::create(state.args["solver"]["linear"]["solver"], state.args["solver"]["linear"]["precond"]);
        solver->set_parameters(state.args["solver"]["linear"]);

        std::vector<int> boundary_nodes = state.periodic_bc ? state.periodic_bc->full_to_periodic(state.boundary_nodes) : state.boundary_nodes;

        StiffnessMatrix A;
        state.build_stiffness_mat(A);        
        {
            const int full_size = A.rows();
            int precond_num = full_size;

            if (state.periodic_bc)
                precond_num = state.periodic_bc->full_to_periodic(A);

            StiffnessMatrix Atmp = A;
            polysolve::linear::prefactorize(*solver, Atmp, boundary_nodes, precond_num, state.args["output"]["data"]["stiffness_mat"]);
        }

        Eigen::MatrixXd rhs;
        rhs.setZero(state.ndof(), shape_velocity.rows() * dim);
        for (int q = 0; q < shape_velocity.rows(); q++)
            for (int d = 0; d < dim; d++)
                for (int i : state.boundary_nodes)
                    rhs(i, q * dim + d) = shape_velocity(q, node_to_primitive[i] * dim + d);

        if (state.periodic_bc)
            rhs = state.periodic_bc->full_to_periodic(rhs, true);
    
        // enforce dirichlet boundary on rhs
        {
            Eigen::VectorXd N;
            N.setZero(A.rows());
            N(boundary_nodes).setOnes();
            rhs -= ((1.0 - N.array()).matrix()).asDiagonal() * (A * (N.asDiagonal() * rhs));
        }

        for (int q = 0; q < shape_velocity.rows(); q++)
        {
            Eigen::MatrixXd sol(state.ndof(), dim);
            for (int d = 0; d < dim; d++)
            {
                Eigen::VectorXd x;
                x.setZero(rhs.rows());
                solver->solve(rhs.col(q * dim + d), x);

                if (state.periodic_bc)
                    sol.col(d) = state.periodic_bc->periodic_to_full(state.ndof(), x);
                else
                    sol.col(d) = x;
            }
            extended_velocity.row(q) = utils::flatten(sol(primitive_to_node, Eigen::all)).transpose();

            // paraviewo::VTUWriter writer;
            // writer.add_field("volume_velocity", utils::unflatten(extended_velocity.row(q).transpose(), dim));
            // writer.add_field("boundary_velocity", utils::unflatten(shape_velocity.row(q).transpose(), dim));

            // writer.write_mesh("debug_" + std::to_string(q) + ".vtu", Vout, Fout);
        }
        
        double error = 0;
        for (int i = 0; i < extended_velocity.cols(); i++)
            if (boundary_flags(i / dim))
                error += (extended_velocity.col(i) - shape_velocity.col(i)).norm();
        
        if (error > 1e-10)
            adjoint_logger().error("Error of volume shape velocity: {}", error);
        
        std::swap(extended_velocity, shape_velocity);
    }
    int SDF2Mesh::size(const int x_size) const
    {
        return Vout.size();
    }
    Eigen::VectorXd SDF2Mesh::eval(const Eigen::VectorXd &x) const
    {
        if (!isosurface_inflator(x))
        {
            adjoint_logger().error("Failed to inflate mesh!");
            return Eigen::VectorXd();
        }

        return utils::flatten(Vout);
    } 
    Eigen::VectorXd SDF2Mesh::apply_jacobian(const Eigen::VectorXd &grad, const Eigen::VectorXd &x) const
    {
        // assert(x.size() == shape_vel.rows());
        
        Eigen::VectorXd mapped_grad = shape_velocity * grad; // shape_vel * (vertex_normals.array() * utils::unflatten(grad, dim).array()).matrix().rowwise().sum();

        // debug
        // {
        //     static int debug_id = 0;
        //     Eigen::MatrixXd vec = utils::unflatten(shape_velocity, dim_);
        //     Eigen::VectorXd grad_normal_coeffs = (vec.array() * utils::unflatten(grad, dim_).array()).matrix().rowwise().sum();
        //     Eigen::MatrixXd grad_normal_direction = vec.array().colwise() * grad_normal_coeffs.array();
        //     Eigen::MatrixXd grad_tangent_direction = utils::unflatten(grad, dim_) - grad_normal_direction;
            
        //     paraviewo::VTUWriter writer;
        //     writer.add_field("gradient", utils::unflatten(grad, dim_));
        //     writer.add_field("gradient_normal", grad_normal_direction);
        //     writer.add_field("gradient_tangent", grad_tangent_direction);

        //     writer.write_mesh("debug_" + std::to_string(debug_id++) + ".vtu", Vout, Fout);
        // }

        return mapped_grad;
    }
    Eigen::VectorXd SDF2Mesh::displace_vertices(const Eigen::VectorXd &x) const
    {
        if (last_x.size() != x.size() || (x - last_x).norm() > 1e-1)
            return Eigen::VectorXd();
        
        Eigen::MatrixXd V = Vout + utils::unflatten(shape_velocity.transpose() * (x - last_x), dim_);
        if (AdjointTools::is_flipped(V, Fout))
            return Eigen::VectorXd();
        
        log_and_throw_error("Not implemented!");
        
        return utils::flatten(V);
    }

    MeshTiling::MeshTiling(const Eigen::VectorXi &nums, const std::string in_path, const std::string out_path): MeshParametrization(in_path, out_path), nums_(nums)
    {
        assert(nums_.size() == 2 || nums_.size() == 3);
        assert(nums.minCoeff() > 0);
    }
    int MeshTiling::size(const int x_size) const
    {
        Eigen::MatrixXd vertices;
        Eigen::MatrixXi cells;
        std::vector<std::vector<int>> elements;
        std::vector<std::vector<double>> weights;
        std::vector<int> body_ids;
        io::MshReader::load(out_path_, vertices, cells, elements, weights, body_ids);
        return vertices.size();
    }
    Eigen::VectorXd MeshTiling::eval(const Eigen::VectorXd &x) const
    {
        Eigen::MatrixXd vertices;
        Eigen::MatrixXi cells;
        std::vector<std::vector<int>> elements;
        std::vector<std::vector<double>> weights;
        std::vector<int> body_ids;
        io::MshReader::load(in_path_, vertices, cells, elements, weights, body_ids);
        const int dim = vertices.cols();

        if (x.size() != vertices.size())
            log_and_throw_error("Inconsistent input mesh in tiling!");
        // else if ((x - utils::flatten(vertices)).norm() > 1e-6)
        // {
        //     adjoint_logger().error("Diff in input mesh and x is {}", (x - utils::flatten(vertices)).norm());
        //     log_and_throw_error("Inconsistent input mesh in tiling!");
        // }
        vertices = utils::unflatten(x, vertices.cols());

        Eigen::MatrixXd Vout;
        Eigen::MatrixXi Fout;
        if (!tiling(vertices, cells, Vout, Fout))
        {
            adjoint_logger().error("Failed to tile mesh!");
            return Eigen::VectorXd();
        }

        return utils::flatten(Vout);
    }
    Eigen::VectorXd MeshTiling::apply_jacobian(const Eigen::VectorXd &grad, const Eigen::VectorXd &x) const
    {
        Eigen::MatrixXd vertices;
        Eigen::MatrixXi cells;
        std::vector<std::vector<int>> elements;
        std::vector<std::vector<double>> weights;
        std::vector<int> body_ids;
        io::MshReader::load(out_path_, vertices, cells, elements, weights, body_ids);
        const int dim = vertices.cols();

        if (grad.size() != vertices.size())
            log_and_throw_error("Inconsistent input mesh in tiling jacobian!");

        Eigen::VectorXd reduced_grad;
        reduced_grad.setZero(x.size());
        assert(grad.size() == index_map.size() * dim);
        for (int i = 0; i < index_map.size(); i++)
            reduced_grad(Eigen::seqN(index_map(i)*dim, dim)) += grad(Eigen::seqN(i*dim, dim));

        return reduced_grad;
    }
    bool MeshTiling::tiling(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &Vnew, Eigen::MatrixXi &Fnew) const
    {
        if (last_V.size() == V.size() && last_V == V)
        {
            std::vector<std::vector<int>> elements;
            std::vector<std::vector<double>> weights;
            std::vector<int> body_ids;
            io::MshReader::load(out_path_, Vnew, Fnew, elements, weights, body_ids);
            return true;
        }

        assert(nums_.size() == V.cols());

        Eigen::MatrixXd Vtmp;
        Eigen::MatrixXi Ftmp;
        Vtmp.setZero(V.rows() * nums_.prod(), V.cols());
        Ftmp.setZero(F.rows() * nums_.prod(), F.cols());

        Eigen::MatrixXd bbox(V.cols(), 2);
        bbox.col(0) = V.colwise().minCoeff();
        bbox.col(1) = V.colwise().maxCoeff();
        Eigen::VectorXd size = bbox.col(1) - bbox.col(0);

        if (V.cols() == 2)
        {
            for (int i = 0, idx = 0; i < nums_(0); i++)
            {
                for (int j = 0; j < nums_(1); j++)
                {
                    Vtmp.middleRows(idx * V.rows(), V.rows()) = V;
                    Vtmp.block(idx * V.rows(), 0, V.rows(), 1).array() += size(0) * i;
                    Vtmp.block(idx * V.rows(), 1, V.rows(), 1).array() += size(1) * j;

                    Ftmp.middleRows(idx * F.rows(), F.rows()) = F.array() + idx * V.rows();
                    idx += 1;
                }
            }
        }
        else
        {
            log_and_throw_error("Not implemented!");
        }

        // clean duplicated vertices
        const double eps = 1e-4;
        Eigen::VectorXi indices;
        {
            std::vector<int> tmp;
            for (int i = 0; i < V.rows(); i++)
            {
                Eigen::VectorXd p = V.row(i);
                if ((p - bbox.col(0)).array().abs().minCoeff() < eps || (p - bbox.col(1)).array().abs().minCoeff() < eps)
                    tmp.push_back(i);
            }

            indices.resize(tmp.size() * nums_.prod());
            for (int i = 0; i < nums_.prod(); i++)
            {
                indices.segment(i * tmp.size(), tmp.size()) = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(tmp.data(), tmp.size());
                indices.segment(i * tmp.size(), tmp.size()).array() += i * V.rows();
            }
        }

        Eigen::VectorXi potentially_duplicate_mask(Vtmp.rows());
        potentially_duplicate_mask.setZero();
        potentially_duplicate_mask(indices).array() = 1;
        Eigen::MatrixXd candidates = Vtmp(indices, Eigen::all);

        Eigen::VectorXi SVI;
        std::vector<int> SVJ;
        SVI.setConstant(Vtmp.rows(), -1);
        int id = 0;
        for (int i = 0; i < Vtmp.rows(); i++)
        {
            if (SVI[i] < 0)
            {
                SVI[i] = id;
                SVJ.push_back(i);
                if (potentially_duplicate_mask(i))
                {
                    Eigen::VectorXd diffs = (candidates.rowwise() - Vtmp.row(i)).rowwise().norm();
                    for (int j = 0; j < diffs.size(); j++)
                        if (diffs(j) < eps)
                            SVI[indices[j]] = id;
                }
                id++;
            }
        }
        Vnew = Vtmp(SVJ, Eigen::all);

        index_map.setConstant(Vtmp.rows(), -1);
        for (int i = 0; i < V.rows(); i++)
            for (int j = 0; j < nums_.prod(); j++)
                index_map(j * V.rows() + i) = i;
        index_map = index_map(SVJ).eval();

        Fnew.resizeLike(Ftmp);
        for (int d = 0; d < Ftmp.cols(); d++)
            Fnew.col(d) = SVI(Ftmp.col(d));

        {
            write_msh(out_path_, Vnew, Fnew);

            adjoint_logger().debug("Saved tiled mesh to {}", out_path_);
        }

        last_V = V;
        return true;
    }

    MeshAffine::MeshAffine(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const std::string in_path, const std::string out_path): MeshParametrization(in_path, out_path), A_(A), b_(b)
    {
    }

    Eigen::VectorXd MeshAffine::eval(const Eigen::VectorXd &x) const
    {
        auto V = utils::unflatten(x, A_.rows());
        V = ((V * A_.transpose()).eval().rowwise() + b_.transpose()).eval();

        if (last_x.size() == x.size() && last_x == x)
            return utils::flatten(V);

        // save to file
        {
            Eigen::MatrixXd vertices;
            Eigen::MatrixXi cells;
            std::vector<std::vector<int>> elements;
            std::vector<std::vector<double>> weights;
            std::vector<int> body_ids;
            io::MshReader::load(in_path_, vertices, cells, elements, weights, body_ids);

            write_msh(out_path_, V, cells);

            adjoint_logger().debug("Saved affined mesh to {}", out_path_);
        }

        last_x = x;
        
        return utils::flatten(V);
    }
    Eigen::VectorXd MeshAffine::apply_jacobian(const Eigen::VectorXd &grad, const Eigen::VectorXd &x) const
    {
        auto tmp = utils::unflatten(grad, A_.rows());
        return utils::flatten(tmp * A_.transpose());
    }
}
