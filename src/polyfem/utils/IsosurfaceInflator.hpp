#pragma once
#include "MatrixUtils.hpp"
#include <polyfem/Common.hpp>
#include <vector>

namespace polyfem::utils
{
    class IsosurfaceInflator
    {
    public:
        IsosurfaceInflator(const std::string folder, const std::string wire_path, const json &options);
        ~IsosurfaceInflator();

        bool inflate(std::vector<double> &params, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &vertex_normals, Eigen::MatrixXd &shape_vel) const;

    private:
        const std::string wire_path;
        std::string symmetry_type;
        
        // temporary files
        std::string sdf_velocity_path, out_path, param_path;
    };
}
