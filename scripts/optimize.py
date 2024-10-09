import json, copy, os, argparse, subprocess
import numpy as np
from pathlib import Path

marching_cube_res = 64
internal_res = 1e-3

perturbed_angle = 0

def homogenized_yy_stress_error_obj(target_stress, n_shape_params):
    return {
        "type": "soft_constraint",
        "soft_bound": [1, 1],
        "power": 2,
        "objective":
        {
            "type": "divide",
            "objective": [{
                "type": "stress",
                "dimensions": [1, 1],
                "state": 0,
                "weight": -1. / target_stress
            },{
                "type": "parametrized_product",
                "parametrization": [
                    {
                        "type": "slice",
                        "from": 0,
                        "to": n_scaling_params,
                        "last": n_shape_params + n_scaling_params
                    }
                ]
            }]
        }
    }

def penalize_shear_obj(weight):
    return {
        "type": "power",
        "power": 2,
        "weight": weight,
        "objective":
        {
            "type": "homo_disp_grad",
            "print_energy": "shear",
            "dimensions": [0, 1],
            "state": 0
        }
    }

def penalize_expand_obj(weight):
    return {
            "type": "power",
            "power": 2,
            "weight": weight,
            "objective":
            {
                "type": "homo_disp_grad",
                "print_energy": "expand",
                "dimensions": [0, 0],
                "state": 0
            }
        }

def load_state_json(path):
    with open(path,'r') as f:
        state_json = json.load(f)
    
    state_json["geometry"][0]["mesh"] = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../tests/differentiable/periodic_x.obj")

    state_json["contact"]["enabled"] = contact
    state_json["contact"]["periodic"] = contact
    state_json["contact"]["dhat"] = dhat
    state_json["solver"]["contact"]["barrier_stiffness"] = barrier
    state_json["materials"][0]["E"] = E
    state_json["materials"][0]["nu"] = nu

    if not contact:
        state_json["solver"]["augmented_lagrangian"]["scaling"] = 10

    state_json["boundary_conditions"]["periodic_boundary"]["linear_displacement_offset"] = [[0, 0], [0, 0]]
    state_json["boundary_conditions"]["periodic_boundary"]["fixed_macro_strain"] = [3]

    return state_json

def load_opt_json(path):
    with open(path,'r') as f:
        opt_json = json.load(f)

    opt_json["solver"]["nonlinear"]["solver"] = args.solver
    if args.solver == "MMA":
        opt_json["solver"]["nonlinear"]["line_search"]["method"] = "None"
    elif args.solver == "L-BFGS-B":
        opt_json["solver"]["nonlinear"]["line_search"]["method"] = "Backtracking"
    else:
        raise Exception("Invalid backward solver!")

    opt_json["parameters"][0]["initial"] = scaling_param + shape_param
    opt_json["parameters"][0]["number"] = len(shape_param) + n_scaling_params

    opt_json["variable_to_simulation"][1] = {
        "type": "periodic-shape-scale",
        "state": [],
        "composition": [
            {
                "type": "slice",
                "from": 0,
                "to": n_scaling_params
            },
            {
                "type": "append-const",
                "value": 1,
                "size": 1,
                "start": 1
            }
        ]
    }

    for v2sim in opt_json["variable_to_simulation"]:
        v2sim["state"] = list(range(len(strain_sample)))

    # slicing shape parameters
    opt_json["variable_to_simulation"][0]["composition"][0] = {
        "type": "slice",
        "from": n_scaling_params,
        "to": opt_json["parameters"][0]["number"]
    }

    # isosurface inflator config
    opt_json["variable_to_simulation"][0]["composition"][1]["wire_path"] = wire_path
    opt_json["variable_to_simulation"][0]["composition"][1]["options"] = {
        "maxArea": internal_res,
        "symmetry": symmetry,
        "forceConsistentInterfaceMesh": symmetry == "symmetric",
        "forceMSGridSize": True,
        "marchingSquaresGridSize": marching_cube_res
    }

    if tile:
        opt_json["variable_to_simulation"][0]["composition"] += [
            {
                "type": "mesh-affine",
                "dimension": 2,
                "transformation": {
                    "scale": [0.5, 0.5]
                }
            },
            {
                "type": "periodic-mesh-tile",
                "dimensions": [2, 2]
            }
        ]
        opt_json["variable_to_simulation"][1]["composition"] += [
            {
                "type": "matrix-product",
                "matrix": [[2, 0], [0, 2]],
                "left_multiply": True
            }
        ]

    opt_json["solver"]["nonlinear"]["box_constraints"]["bounds"] = [
        [2e-1] * n_scaling_params + [0] * n_positional_params + [0.02] * n_thickness_params + [0.001] * n_blending_params,
        [5] * n_scaling_params + [1] * n_positional_params + [0.15] * n_thickness_params + [0.1] * n_blending_params]
    opt_json["solver"]["nonlinear"]["box_constraints"]["max_change"] = [0.2] * n_scaling_params + [0.05] * n_positional_params + [0.005] * n_thickness_params + [0.01] * n_blending_params

    opt_json["states"] = [{"path": "state-" + str(i) + ".json"} for i in range(len(strain_sample))]

    for i in range(len(opt_json["variable_to_simulation"][1]["state"]) - 1):
        opt_json["states"][opt_json["variable_to_simulation"][1]["state"][i+1]]["initial_guess"] = i

    # optimization objective
    opt_json["functionals"] = []

    functional_template = homogenized_yy_stress_error_obj(target_stress, len(shape_param))

    for i in range(len(strain_sample)):
        functional_template["objective"]["objective"][0]["state"] = i
        functional_template["objective"]["print_energy"] = "stress-" + str(i)
        opt_json["functionals"].append(copy.deepcopy(functional_template))

    functional_template = penalize_shear_obj(50)

    for i in range(len(strain_sample)):
        functional_template["objective"]["state"] = i
        functional_template["objective"]["print_energy"] = "shear-" + str(i)
        opt_json["functionals"].append(copy.deepcopy(functional_template))

    if penalize_expansion:
        functional_template = penalize_expand_obj(10)

        for i in range(len(strain_sample)):
            functional_template["objective"]["state"] = i
            functional_template["objective"]["print_energy"] = "expand-" + str(i)
            opt_json["functionals"].append(copy.deepcopy(functional_template))

    # stopping conditions
    opt_json["stopping_conditions"] = []

    functional_template = {
        "type": "plus-const",
        "value": -0.0025,
        "objective": homogenized_yy_stress_error_obj(target_stress, len(shape_param))
    }

    for i in range(len(strain_sample)):
        functional_template["objective"]["objective"]["objective"][0]["state"] = i
        opt_json["stopping_conditions"].append(copy.deepcopy(functional_template))

    functional_template = {
        "type": "plus-const",
        "value": -0.0025,
        "objective": penalize_shear_obj(1)
    }

    for i in range(len(strain_sample)):
        functional_template["objective"]["objective"]["state"] = i
        opt_json["stopping_conditions"].append(copy.deepcopy(functional_template))

    if penalize_expansion:
        functional_template = {
            "type": "plus-const",
            "value": -0.0025,
            "objective": penalize_expand_obj(1)
        }

        for i in range(len(strain_sample)):
            functional_template["objective"]["objective"]["state"] = i
            opt_json["stopping_conditions"].append(copy.deepcopy(functional_template))

    return opt_json

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('stress', type=float)
    parser.add_argument('wire_path', type=str)
    parser.add_argument('--solver', type=str, default = "L-BFGS-B")
    parser.add_argument('--symmetry', type=str, default = "doubly_periodic")
    parser.add_argument('--params', type=str, default = "")

    parser.add_argument('--strain', type=float, default = 0.25)
    parser.add_argument('--n_samples', type=int, default = 4)
    
    parser.add_argument('--iso_path', type=str, default="../inflator/build/isosurface_inflator/isosurface_cli")
    parser.add_argument('--exe_path', type=str, default="../build/PolyFEM_bin")
    
    parser.add_argument('--E', type=float, default=1e6)
    parser.add_argument('--nu', type=float, default=0.3)
    parser.add_argument('--no_contact', action='store_true')
    parser.add_argument('--dhat', type=float, default=1e-3)
    parser.add_argument('--barrier', type=float, default=-1)

    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--tile', action='store_true')
    parser.add_argument('--penalize_expansion', action='store_true')
    args = parser.parse_args()

    if not os.path.isabs(args.wire_path):
        args.wire_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), args.wire_path)
    if not os.path.isabs(args.iso_path):
        args.iso_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), args.iso_path)
    if not os.path.isabs(args.exe_path):
        args.exe_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), args.exe_path)

    target_stress = args.stress
    strain_sample = np.linspace(0.1, args.strain, args.n_samples)
    
    wire_path = args.wire_path
    symmetry = args.symmetry
    tile = args.tile # optimize on 2x2 tiles
    penalize_expansion = args.penalize_expansion

    shape_param = []
    scaling_param = [1]
    t = subprocess.check_output(args.iso_path + " 2D_" + symmetry + " --defaultThickness 0.05 " + wire_path + " out.msh", shell=True)
    print(t)
    t1 = t.decode('utf-8').split('\n')[1].split("\t")
    for token in t1:
        try:
            shape_param.append(float(token))
        except ValueError:
            pass

    # generate initial shape parameters
    if os.path.exists(args.params):
        params = np.loadtxt(args.params)
        shape_param = params[1:].tolist()
        scaling_param = params[:1].tolist()

    print("scale paramemters:", scaling_param)
    print("shape paramemters:", shape_param)

    t2 = t.decode('utf-8').split('\n')[2].split("\t")
    n_shape_param = []
    for token in t2:
        try:
            n_shape_param.append(int(token))
        except ValueError:
            pass

    print(n_shape_param)
    n_positional_params = n_shape_param[0]
    n_thickness_params = n_shape_param[1]
    n_blending_params = n_shape_param[2]
    n_scaling_params = 1

    E = args.E
    nu = args.nu
    contact = not args.no_contact
    dhat = args.dhat
    barrier = E if args.barrier <= 0 else args.barrier

    folder = "../result/"
    if args.no_contact:
        folder = os.path.join(folder, "no_contact")
    folder = os.path.join(folder, os.path.splitext(os.path.split(wire_path)[1])[0] + "_" + str(strain_sample[-1]) + "_" + str(target_stress))
    if os.path.exists(os.path.join(folder, "optimized-params.txt")):
        print("Already succeeded!")
        exit()
    elif os.path.exists(folder):
        os.system("rm -r " + folder)
    Path(folder).mkdir(parents=True)

    state_json = load_state_json("state-1.json")

    with open(os.path.join(folder, "common.json"), 'w') as outfile:
        del state_json["boundary_conditions"]["periodic_boundary"]["linear_displacement_offset"]
        outfile.write(json.dumps(state_json, indent=4))

    for i in range(len(strain_sample)):
        state_json = {
            "common": "common.json",
            "boundary_conditions":
            {
                "periodic_boundary": 
                {
                    "linear_displacement_offset": [[0,0],[0,-strain_sample[i]]]
                }
            }
        }
        with open(os.path.join(folder, "state-" + str(i) + ".json"), "w") as outfile:
            outfile.write(json.dumps(state_json, indent=4))
    
    opt_json = load_opt_json("opt-ratio.json")
    with open(os.path.join(folder, "opt.json"), "w") as outfile:
        outfile.write(json.dumps(opt_json, indent=4))

    os.chdir(folder)
    print("run " + folder)
    with open("log",'w') as outfile:
        subprocess.run([args.exe_path, 
                        "--max_threads", str(args.threads),
                        "-j", "opt.json",
                        "--log_level", "trace",
                        "--ns"], stdout=outfile)
