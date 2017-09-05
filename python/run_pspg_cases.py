from moose_calc_routines import *
from sympy import *
import sympy as sp
# from IPython.core.interactiveshell import InteractiveShell
# InteractiveShell.ast_node_interactivity = "all"
# init_printing()

x, y = var('x y')

u = 0.4*sin(0.5*pi*x) + 0.4*sin(pi*y) + 0.7*sin(0.2*pi*x*y) + 0.5
v = 0.6*sin(0.8*pi*x) + 0.3*sin(0.3*pi*y) + 0.2*sin(0.3*pi*x*y) + 0.3
p = 0.5*sin(0.5*pi*x) + 1.0*sin(0.3*pi*y) + 0.5*sin(0.2*pi*x*y) + 0.5
vxx = diff(u, x)

uvec = sp.Matrix([u, v])

volume_source = {
                # 'vel_x' : prep_moose_input(L_stokes_traction(uvec, p, x, y)[0]),
                # 'vel_y' : prep_moose_input(L_stokes_traction(uvec, p, x, y)[1]),
                'vel_x' : prep_moose_input(L_stokes(uvec, p, x, y)[0]),
                'vel_y' : prep_moose_input(L_stokes(uvec, p, x, y)[1]),
                # 'vel_x' : prep_moose_input(L_momentum_traction_no_turbulence(uvec, p, x, y)[0]),
                # 'vel_y' : prep_moose_input(L_momentum_traction_no_turbulence(uvec, p, x, y)[1]),
                'p' : prep_moose_input(L_pressure(uvec, x, y))}
solution_dict = {'vel_x' : u, 'vel_y' : v, 'p' : p, 'vxx' : vxx}

# h_list = ['5', '10', '20', '40']
# h_array = np.array([.2, .1, .05, .025])
h_list = ['4', '8', '16', '32']
h_array = np.array([.25, .125, .0625, .03125])
base = "pspg_mms_test"
input_dir = "/Users/lindad/projects/articuno/python"
exe_path = "/Users/lindad/projects/moose/modules/navier_stokes/navier_stokes-opt"
mms_kernel_cases(h_list, volume_source, solution_dict, base, exe_path, input_dir)

optional_save_string='element_errors_q1_q1_p_by_parts_laplace_stokes_alpha_1e-6_with_time_deriv_code'
plot_order_accuracy(h_array, base, input_dir, optional_save_string=optional_save_string)
