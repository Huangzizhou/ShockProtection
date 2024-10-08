#pragma once

#include "Parametrization.hpp"
#include <polyfem/utils/IsosurfaceInflator.hpp>
#include <filesystem>
#include <polyfem/Common.hpp>

namespace polyfem::solver
{
    // for parametrizations that outputs a mesh
    class MeshParametrization : public Parametrization
    {
    public:
        MeshParametrization(const std::string out_path) : in_path_(""), out_path_(out_path) {}
        MeshParametrization(const std::string in_path, const std::string out_path) : in_path_(in_path), out_path_(out_path) {}
        virtual ~MeshParametrization() = default;

        bool has_input_mesh() const override { return in_path_ != ""; }
        std::string input_path() const override { return in_path_; }
        std::string output_path() const override { return out_path_; }

    protected:
        const std::string in_path_, out_path_;
    };

    // generates a unit size mesh using isosurface inflator
    class SDF2Mesh : public MeshParametrization
    {
    public:
		SDF2Mesh(const std::string wire_path, const std::string out_path, const bool volume_velocity, const json &opts) : MeshParametrization(out_path), volume_velocity_(volume_velocity), dim_(2), inflator(std::filesystem::current_path().u8string(), wire_path, opts)
        {
        }

        int size(const int x_size) const override;

        Eigen::VectorXd eval(const Eigen::VectorXd &x) const override;
        Eigen::VectorXd apply_jacobian(const Eigen::VectorXd &grad, const Eigen::VectorXd &x) const override;
    
        Eigen::VectorXd displace_vertices(const Eigen::VectorXd &x) const;
    private:
        bool isosurface_inflator(const Eigen::VectorXd &x) const;
        void extend_to_internal() const;
        
        const bool volume_velocity_;
        const int dim_;

        mutable std::vector<std::tuple<int, double>> inactive_shape_params;
        mutable Eigen::VectorXd last_x;
        mutable Eigen::MatrixXd Vout, shape_velocity;
        mutable Eigen::MatrixXi Fout;
        mutable Eigen::Matrix<bool, -1, 1> boundary_flags;

        utils::IsosurfaceInflator inflator;
    };

    // given a periodic mesh input, repeat the mesh to form a tile
    class MeshTiling : public MeshParametrization
    {
    public:
        MeshTiling(const Eigen::VectorXi &nums, const std::string in_path, const std::string out_path);

        int size(const int x_size) const override;

        Eigen::VectorXd eval(const Eigen::VectorXd &x) const override;
        Eigen::VectorXd apply_jacobian(const Eigen::VectorXd &grad, const Eigen::VectorXd &x) const override;
    
    private:
        const Eigen::VectorXi nums_;

        bool tiling(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &Vnew, Eigen::MatrixXi &Fnew) const;

        mutable Eigen::MatrixXd last_V;
        mutable Eigen::VectorXi index_map;
    };

    // perform affine mapping on the input mesh
    class MeshAffine : public MeshParametrization
    {
    public:
        MeshAffine(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const std::string in_path, const std::string out_path);

        int size(const int x_size) const override { return x_size; }

        Eigen::VectorXd eval(const Eigen::VectorXd &x) const override;
        Eigen::VectorXd apply_jacobian(const Eigen::VectorXd &grad, const Eigen::VectorXd &x) const override;
    
    private:
        const Eigen::MatrixXd A_;
        const Eigen::VectorXd b_;

        mutable Eigen::VectorXd last_x;
    };
}
