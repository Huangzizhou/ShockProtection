#include "IsosurfaceInflator.hpp"
#include "Logger.hpp"
#include <polyfem/io/MshReader.hpp>
#include <filesystem>
#include <fstream>

namespace polyfem::utils
{
    namespace
    {
        template <typename T>
        std::string to_string_with_precision(const T a_value, const int n)
        {
            std::ostringstream out;
            out.precision(n);
            out << std::fixed << a_value;
            return out.str();
        }

        bool file_exists(const std::string& name) {
            std::ifstream f(name.c_str());
            return f.good();
        }
    }

    IsosurfaceInflator::IsosurfaceInflator(const std::string folder, const std::string wire_path, const json &options) : wire_path(wire_path)
    {
        symmetry_type = options["symmetry"];
        static int n_inflator_calls = 0;

        sdf_velocity_path = (std::filesystem::path(folder) / std::filesystem::path("inflator_" + std::to_string(n_inflator_calls) + std::string(".msh"))).u8string();
		out_path = (std::filesystem::path(folder) / std::filesystem::path("inflator_useless_" + std::to_string(n_inflator_calls) + std::string(".msh"))).u8string();
		param_path = (std::filesystem::path(folder) / std::filesystem::path("inflator_" + std::to_string(n_inflator_calls) + std::string(".json"))).u8string();

        adjoint_logger().debug("Inflator writes shape velocity to {}", sdf_velocity_path);

        // write meshing options to file
        {
            json args = options;
            args.erase("symmetry");
            std::ofstream out(param_path);
            out << args.dump();
            out.close();
        }
    }

    IsosurfaceInflator::~IsosurfaceInflator()
    {
        if (file_exists(sdf_velocity_path))
            std::remove(sdf_velocity_path.c_str());
        if (file_exists(out_path))
            std::remove(out_path.c_str());
        if (file_exists(param_path))
            std::remove(param_path.c_str());
    }

    bool IsosurfaceInflator::inflate(std::vector<double> &params, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &vertex_normals, Eigen::MatrixXd &shape_vel) const
    {
        const int dim = 2;

        std::string shape_params;
        for (int i = 0; i < params.size(); i++)
            shape_params += to_string_with_precision(params[i], 16) + " ";
        
        std::string command = std::string(INFLATOR_BINARY) + " "        // binary path
                            + std::to_string(dim) + std::string("D_")   // dimension
                            + symmetry_type + " "                       // symmetry
                            + wire_path                                 // wireframe path
                            + " -m " + param_path + " "                 // meshing params
                            + "--params \"" + shape_params + "\""       // shape parameters
                            + " -S " + sdf_velocity_path + " "          // dump shape velocity path
                            + out_path;                                 // unit cell mesh path, not used

        int return_val;
        bool failed = false;
        try 
        {
            return_val = system(command.c_str());
        }
        catch (const std::runtime_error &err)
        {
            failed = true;
        }

        if (failed || return_val)
        {
            adjoint_logger().error("remesh command \"{}\" returns {}", command, return_val);
            return false;
        }
        else
            adjoint_logger().debug("remesh command \"{}\" returns {}", command, return_val);

        {
            std::vector<std::vector<int>> elements;
            std::vector<std::vector<double>> weights;
            std::vector<int> body_ids;
            std::vector<std::string> node_data_name;
            std::vector<std::vector<double>> node_data;
            io::MshReader::load(sdf_velocity_path, V, F, elements, weights, body_ids, node_data_name, node_data);

            assert(node_data_name.size() == node_data.size());

            vertex_normals = Eigen::Map<Eigen::MatrixXd>(node_data[0].data(), 3, V.rows()).topRows(V.cols()).transpose();
            
            shape_vel.setZero(params.size(), V.rows());
            for (int i = 1; i < node_data_name.size(); i++)
                shape_vel.row(i - 1) = Eigen::Map<RowVectorNd>(node_data[i].data(), V.rows());
        }

        return true;
    }
}
