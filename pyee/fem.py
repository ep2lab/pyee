import dolfinx
from dolfinx.io import VTKFile
import dolfinx.plot
from dolfinx.fem import(Function, FunctionSpace, dirichletbc, locate_dofs_topological)
from dolfinx.mesh import locate_entities_boundary, create_rectangle
from ufl import (FiniteElement, VectorElement, TensorElement, SpatialCoordinate, 
                 TrialFunction, TestFunction,TrialFunctions, TestFunctions, 
                 dx, inner, as_vector, as_matrix, sign, conj, real)
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import numpy as np
from petsc4py import PETSc
from ufl.finiteelement.mixedelement import MixedElement
import os
import pickle
import sys

def core(dir, m):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    file = open(dir + "/data.pkl","rb")
    data = pickle.load(file)
    file.close()

    # Constants and parameters
    omega = 2 * np.pi * data['simulation']['freq'] * 1e6
    mu0   = 4*np.pi*1e-7
    eps0  = 8.8541878128e-12
    c0    = (1/mu0/eps0)**0.5
    k0    = omega /c0
    Z0    = (mu0/eps0)**0.5
    m     = int(m)

    # FEM mesh
    if data['geometry']['nz'] == data['mesh']['Z'].shape[0] - 1:  
        st = 1
    else:                                 # For legacy staggered grids (double number of nodes in kappa and J)
        st = 2
   
    Z  = data['mesh']['Z'][0::st,0::st]
    R  = data['mesh']['R'][0::st,0::st]
    nz = Z.shape[0] - 1
    nx = Z.shape[1] - 1
    Lz = Z.max()
    Lx = R.max()

    mesh    = create_rectangle(comm, points=[[0.,0.],[1.0,1.0]], n=[nz,nx])
    I       = np.dot(np.round(mesh.geometry.x * [nz/1.0,nx/1.0,0]).astype(int),[nx+1,1,0])    # Permutation matrix
    dof_map = mesh.geometry.dofmap.array

    # Load new coordinates
    mesh.geometry.x[:,0] = Z.flatten()[I]
    mesh.geometry.x[:,1] = R.flatten()[I]
    
    # from mesh_gmsh import mesh
    deg = 1

    # Create combined function space
    e_N = FiniteElement("N1curl", mesh.ufl_cell(), degree = deg)
    e_L = FiniteElement("Lagrange", mesh.ufl_cell(), degree = deg)

    v_L = VectorElement("Lagrange", mesh.ufl_cell(), degree = 1, dim = 3)
    t_L = TensorElement("Lagrange", mesh.ufl_cell(), degree = 1, shape = (3,3))

    CS = FunctionSpace(mesh, e_N * e_L)  # Mixed function space

    # Kappa function
    kappa = Function(FunctionSpace(mesh, t_L))

    def kappa_fun(x,kM):

        values = np.zeros((9, x.shape[1]), dtype=PETSc.ScalarType)
        
        for i in range(0,x.shape[1]):
            
            j = dof_map[i]
            iz = int(I[j] / (nx+1))
            ix = I[j] - (nx+1)*iz

            values[:,i] = [kM[iz,ix,0,0],
                           kM[iz,ix,0,1],
                           kM[iz,ix,0,2],
                           kM[iz,ix,1,0],
                           kM[iz,ix,1,1],
                           kM[iz,ix,1,2],
                           kM[iz,ix,2,0],
                           kM[iz,ix,2,1],
                           kM[iz,ix,2,2]]

        return values

    kappa.interpolate(lambda x: kappa_fun(x,data['kappa'][0::st,0::st,:,:]))

    # J function
    J = Function(FunctionSpace(mesh,v_L))
       
    def J_fun(x,JM):

        values = np.zeros((3, x.shape[1]), dtype=PETSc.ScalarType)

        for i in range(0,x.shape[1]):
            
            j  = dof_map[i]
            iz = int(I[j] / (nx+1))
            ix = I[j] - (nx+1)*iz

            values[:,i] = [JM[iz,ix,0],
                           JM[iz,ix,1],
                           JM[iz,ix,2]]

        return values

    J.interpolate(lambda x: J_fun(x,data['current'][0::st,0::st,:]))
    
    # Boundary conditions. This way is due to the legacy format used in FD
    Zw = Z[data['BCs']['EqEy'] == 0]
    Rw = R[data['BCs']['EqEy'] == 0]

    def wall(x):

        outer = np.logical_or(np.isclose(x[1],Lx), 
                            np.logical_or(np.isclose(x[0],0.0),np.isclose(x[0],Lz)))

        interior = np.zeros(x.shape[1], dtype=bool)
        # for j in range(0, len(Zw)):
        #     cond = np.logical_and(np.isclose(x[0],Zw[j], atol=1e-3),np.isclose(x[1],Rw[j], atol=1e-3))
        #     interior = np.logical_or(cond, interior)

        # interior = np.isclose(x[0], 0.12, atol=1e-3)
        
        # print(np.sum(interior))

        return np.logical_or(outer,interior)

    def axis(x):
        return np.isclose(x[1],0.0)

    CS0 = CS.sub(0).collapse()[0]
    CS1 = CS.sub(1).collapse()[0]
    plane_BC = Function(CS0)
    out_BC   = Function(CS1)

    with plane_BC.vector.localForm() as loc:
        loc.set(0.0)
    with out_BC.vector.localForm() as loc:
        loc.set(0.0)

    # def f_plane(x):
    #     # print(x.shape)
    #     vals = np.zeros((2,x.shape[1]))
    #     return vals
    # def f_out(x):
    #     # print(x.shape)
    #     vals = np.zeros(x.shape[1])
    #     return vals
    # plane_BC.interpolate(f_plane)
    # out_BC.interpolate(f_out)

    wall_v     = locate_dofs_topological((CS.sub(0),CS0), 1, locate_entities_boundary(mesh, 1, wall))
    wall_p     = locate_dofs_topological((CS.sub(1),CS1), 1, locate_entities_boundary(mesh, 1, wall))
    axis_v     = locate_dofs_topological((CS.sub(0),CS0), 1, locate_entities_boundary(mesh, 1, axis))
    axis_p     = locate_dofs_topological((CS.sub(1),CS1), 1, locate_entities_boundary(mesh, 1, axis))

    bcwv = dirichletbc(plane_BC, wall_v, CS)
    bcwp = dirichletbc(out_BC, wall_p, CS)
    bcav = dirichletbc(plane_BC, axis_v, CS)
    bcap = dirichletbc(out_BC, axis_p, CS)

    x = SpatialCoordinate(mesh)
    (N_i, L_i) = TrialFunctions(CS)
    (N_j, L_j) = TestFunctions(CS)

    if data['geometry']['axi']:
        if m == 0:
            E = as_vector([N_i[0], N_i[1], L_i])
            T = as_vector([N_j[0], N_j[1], L_j])
            bcs = [bcwv, bcwp, bcap]
        elif m == 1 or m == -1: 
            E = as_vector([x[1]*N_i[0], x[1]*N_i[1] - sign(m)*1j*L_i, L_i])
            T = as_vector([x[1]*N_j[0], x[1]*N_j[1] - sign(m)*1j*L_j, L_j])
            bcs = [bcwv, bcwp]
        else:
            E = as_vector([x[1]*N_i[0], x[1]*N_i[1], L_i])
            T = as_vector([x[1]*N_j[0], x[1]*N_j[1], L_j])
            bcs = [bcwv, bcwp, bcap]
    else:
        E = as_vector([N_i[0], N_i[1], L_i])
        T = as_vector([N_j[0], N_j[1], L_j])
        bcs = [bcwv, bcwp, bcav, bcap]
        x = as_vector([x[0], 100000.0])     # Take radius to infinity
        m *= 100000.0 

    E_term = as_vector([E[2].dx(1) + E[2]/x[1] - 1j*m*E[1]/x[1],
                        1j*m*E[0]/x[1] - E[2].dx(0),
                        E[1].dx(0) - E[0].dx(1)])

    T_term = as_vector([T[2].dx(1) + T[2]/x[1] - 1j*m*T[1]/x[1],
                        1j*m*T[0]/x[1] - T[2].dx(0),
                        T[1].dx(0) - T[0].dx(1)])

    s_xx = inner(E_term, T_term)      
    a = (s_xx - k0**2*inner(kappa*E, T))*x[1]*dx

    L = -1j*k0*Z0*inner(J, T)*x[1]*dx

    problem = LinearProblem(a, L, bcs, petsc_options={"ksp_type": "preonly",
        "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"})

    # problem = LinearProblem(a, L, bcs, petsc_options={"ksp_type": "preonly",
    #     "pc_type": "lu"})

    u = problem.solve()

    if data['geometry']['axi']:
        if m == 0:
            Ez = u[0]
            Ex = u[1]
            Ey = u[2]
        elif m == 1 or m == -1:
            Ez = u[0]*x[1]
            Ex = u[1]*x[1] - 1j*sign(m)*u[2]
            Ey = u[2]
        else:
            Ez = u[0]*x[1]
            Ex = u[1]*x[1]
            Ey = u[2]
    else:
        Ez = u[0]
        Ex = u[1]
        Ey = u[2]

    # Compute power absorption
    eye = as_matrix([[1, 0, 0],[0, 1, 0], [0, 0, 1]])
    Ev  = as_vector([Ez, Ex, Ey])
    Jp  = 1j*omega*eps0*(eye - kappa)*Ev
    Qa  = real(0.5*inner(Ev, Jp)) # UFL operator: Take the inner product of a and b. The complex conjugate of the second argument is taken.

    # Project to Lagrange
    L = FunctionSpace(mesh, ("Lagrange", 1))

    Ex = project(Ex, L);   Ex.name = "Ex"
    Ez = project(Ez, L);   Ez.name = "Ez"
    Ey = project(Ey, L);   Ey.name = "Ey"
    Jz = project(J[0], L); Jz.name = "Jz"
    Jx = project(J[1], L); Jx.name = "Jx"
    Jy = project(J[2], L); Jy.name = "Jy"
    Qa = project(Qa, L);   Qa.name = "Qa" 

    if 'postproc' in data.keys():

        if data['postproc']['exportPVD']:
            if not(os.path.isdir(data['general']['simdir']+ "FEMout/")):
                os.mkdir(data['general']['simdir']+ "FEMout/")
            outname = "wave"
            with VTKFile(comm, data['general']['simdir']+ "FEMout/" + outname + ".pvd", "w") as file:
                file.write_mesh(mesh)
                file.write_function([Ez, Ex, Ey, Qa])

    fem_sol = dict()

    fem_sol['Ez']  = FEM2mat(Ez, I, nz, nx)
    fem_sol['Ex']  = FEM2mat(Ex, I, nz, nx)
    fem_sol['Ey']  = FEM2mat(Ey, I, nz, nx)
    fem_sol['Jz']  = FEM2mat(Jz, I, nz, nx)
    fem_sol['Jx']  = FEM2mat(Jx, I, nz, nx)
    fem_sol['Jy']  = FEM2mat(Jy, I, nz, nx)
    fem_sol['Qa']  = FEM2mat(Qa, I, nz, nx)

    file = open(dir + "/sol" + str(rank) + ".pkl","wb")
    pickle.dump(fem_sol, file)
    file.close()

    return fem_sol


def project(expression, V):
    u_, v_ = TrialFunction(V), TestFunction(V)
    a_p = inner(u_, v_) * dx
    L_p = inner(expression, v_) * dx
    projection = LinearProblem(a_p, L_p)
    return projection.solve()

def FEM2mat(field, I, nz, nx):

    M    = np.ones([nz+1,nx+1], dtype=complex).flatten()*np.nan
    M[I] = field.x.array
    M    = M.reshape([nz+1,nx+1])

    return M

if __name__ == '__main__':
    # Map command line arguments to function arguments.
    core(*sys.argv[1:])
