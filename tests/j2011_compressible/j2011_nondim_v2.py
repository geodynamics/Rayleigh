"""
j2011_nondim_v2.py — ENTROPY formulation — compressible Jones 2011 hydro benchmark lab.

Formulation: FIRST-ORDER tau reduction (SCZ-validated placement), nondimensional
in the paper's hydrodynamic units (Table 1): length d = ro-ri, time d^2/nu,
velocity nu/d = 800.8 cm/s, temperature Tc = 111,557 K, density rho_c, so all
monitor outputs are directly comparable to Jones et al. (2011) Table 2.

BCs: stress-free impenetrable; fixed temperature T' = DeltaT_hat at r=ri, 0 at
r=ro (DeltaT = DeltaS*Ti/cv = 4234 K, the fixed-density truncation of the
benchmark's fixed-entropy condition).

Benchmark targets (Table 2, dimensionless): KE = 81.86, zonal KE = 9.377,
meridional KE = 0.02202, luminosity = 4.19886 (conduction state: 3.954041),
drift omega = 17.6427 (period 0.0187440), u_phi = 0.86183 and S = 0.93301 at
the ur=0 mid-shell equatorial point. Linear m=19 growth at this Ra: +77.551
(energy: 155.1), frequency 269.05, in these time units.

Monitor coverage: KE, KEz, KEm, L_top, L_bot direct; drift period tau = the
oscillation period of the Tp/up point-probe columns (fixed point at rm,
equator, phi=0 sees the m=19 pattern sweep at period tau); the ur=0 point
values (u_phi, S) are a post-processing step on equatorial slices, not
monitored here.

Angular momentum: NOT corrected (paper Sec. 8 warns of secular drift under
stress-free BCs on long runs; L'_z is monitored so the drift is visible).
"""

import numpy as np
import dedalus.public as d3
from scipy.special import lpmv
from math import factorial
import logging
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------- CGS inputs
ri_c, ro_c = 2.45e9, 7.0e9
d_c     = ro_c - ri_c
beta    = ri_c / ro_c
nu_c    = 3.64364e12
kappa_c = 3.64364e12
Omega_c = 1.76e-4
Rgas_c  = 3.503e7
gamma   = 1.5
npoly   = 2.0
Nrho    = 5.0
Gconst  = 6.67430e-8
Mstar   = 1.8992e30
rho_i_c = 1.1
cp_c    = (npoly + 1.0) * Rgas_c          # = (n+1)R, adiabatic polytrope
cv_c    = cp_c - Rgas_c
DeltaT_c = 4234.0                          # = DeltaS * Ti / cv

# ------------------------------------------------- polytrope (Jones Table 1)
zeta_o = (beta + 1.0) / (beta * np.exp(Nrho / npoly) + 1.0)
c0     = (2.0 * zeta_o - beta - 1.0) / (1.0 - beta)
c1     = (1.0 + beta) * (1.0 - zeta_o) / (1.0 - beta)**2
Tc_c   = Gconst * Mstar / (cp_c * c1 * d_c)
rho_c_mid = rho_i_c / ((c0 + c1 * d_c / ri_c)**npoly)   # rho_c: zeta=1 value
DeltaS_c  = DeltaT_c * cv_c / (Tc_c * (c0 + c1 * d_c / ri_c))  # back out DS

# ------------------------------------------- nondimensionalization (paper's)
# length d, time d^2/nu, velocity nu/d, temperature Tc, density rho_c
u0   = nu_c / d_c                          # 800.8 cm/s
t0   = d_c**2 / nu_c                       # 5.6818e6 s
ri, ro  = ri_c / d_c, ro_c / d_c           # 0.538462, 1.538462
d_shell = ro - ri                          # 1
nu_v  = 1.0
kcond = gamma * 1.0                        # gamma*kappa/nu
OmegaN = Omega_c * t0                      # 1/E = 1000
Rg     = Rgas_c * Tc_c / u0**2             # pressure-force coefficient
DeltaT = DeltaT_c / Tc_c                   # 0.037952
# entropy monitor coefficients: S_hat = S'/DeltaS
cS_T   = cv_c * Tc_c / (DeltaS_c * Tc_c)   # cv/DeltaS  (diagnostics)
cS_q   = Rgas_c / DeltaS_c                 # R/DeltaS   (diagnostics)
eps    = DeltaS_c / cv_c                   # entropy parameter, 0.01215
gm1    = gamma - 1.0
kap    = 1.0                               # entropy diffusivity (=nu, Pr=1)
c_phi  = u0**2 / (Tc_c * DeltaS_c)         # viscous-heating coefficient in S-eq

# ------------------------------------------------------------- run controls
Ntheta, Nr = 192, 128
Nphi = 2*Ntheta
dealias = 3/2
dtype = np.float64
mesh = (64, 32)                            # set None for serial
stop_time = 10.0                            # viscous/diffusion times
MONITOR_CADENCE = 20
STEPPER = "MCNAB2"
initial_dt = 1e-5                        # = 33 s, benchmark-code timestep
max_dt     = 1e-5
CFL_SAFETY = 0.2
seed_amp   = 1.0 / Tc_c                    # 1 K seed, Rayleigh init_type=9
IC_MODE    = "conductive"                        # "seed" (paper protocol) | "conductive"
AM_CORRECT = False                         # subtract solid-body rotation (paper Sec. 8)
AM_CADENCE = 100
restart      = False
restart_file = "checkpoints/checkpoints_s202.h5"

# ------------------------------------------------------------------ domain
coords = d3.SphericalCoordinates('phi', 'theta', 'r')
s2coords = d3.S2Coordinates('phi', 'theta')
dist = d3.Distributor(coords, dtype=dtype, mesh=mesh)
basis = d3.ShellBasis(coords, shape=(Nphi, Ntheta, Nr), radii=(ri, ro),
                      dealias=dealias, dtype=dtype)
bc_ri, bc_ro = basis.inner_surface, basis.outer_surface
phi, theta, r = dist.local_grids(basis)

u   = dist.VectorField(coords, name='u', bases=basis)
S   = dist.Field(name='S', bases=basis)
lnr = dist.Field(name='lnr', bases=basis)
tau_u1 = dist.VectorField(coords, name='tau_u1', bases=bc_ri)
tau_u2 = dist.VectorField(coords, name='tau_u2', bases=bc_ro)
tau_S1 = dist.Field(name='tau_S1', bases=bc_ri)
tau_S2 = dist.Field(name='tau_S2', bases=bc_ro)

# ------------------------------------------------------------- background
Tbar  = dist.Field(name='Tbar', bases=basis.radial_basis)
dTbar = dist.VectorField(coords, name='dTbar', bases=basis.radial_basis)
rp    = dist.VectorField(coords, name='rp', bases=basis.radial_basis)
rhobar = dist.Field(name='rhobar', bases=basis.radial_basis)
rvec  = dist.VectorField(coords, name='rvec', bases=basis.radial_basis)
er    = dist.VectorField(coords, name='er', bases=basis.radial_basis)
ep    = dist.VectorField(coords, name='ep', bases=basis)
eye   = dist.TensorField((coords, coords), name='eye',
                         bases=(basis.radial_basis, basis.radial_basis))
for ii in range(3):
    eye['g'][ii][ii] = -2/3

zeta_r = c0 + c1 * d_shell / r
Tbar['g']     = zeta_r
Tbar_full = dist.Field(name='Tbar_full', bases=basis)
Tbar_full['g'] = zeta_r
dTbar['g'][2] = -c1 * d_shell / r**2
rp['g'][2]    = npoly * (-c1 * d_shell / r**2) / zeta_r
rhobar['g']   = zeta_r**npoly
rvec['g'][2]  = r
er['g'][2]    = 1.0
tp    = dist.VectorField(coords, name='tp', bases=basis.radial_basis)
tp['g'][2]    = (-c1*d_shell/r**2)/zeta_r          # dlnTbar/dr = rp/npoly
gvec  = dist.VectorField(coords, name='gvec', bases=basis.radial_basis)
gvec['g'][2]  = -Rg*(npoly+1.0)*c1*d_shell/r**2    # gravity vector (inward), = R(dTbar+Tbar*rp)
iTbar = dist.Field(name='iTbar', bases=basis.radial_basis)
iTbar['g']    = 1.0/zeta_r
ep['g'][0]    = 1.0

ez = dist.VectorField(coords, name='ez', bases=basis)
ez['g'][1] = -np.sin(theta)
ez['g'][2] =  np.cos(theta)
ez = d3.Grid(ez)

# ------------------------------------------------------------- shorthands
grad, div, trace, cross = d3.grad, d3.div, d3.trace, d3.cross
trans = getattr(d3, 'trans', getattr(d3, 'transpose', None))
lift_basis = basis.derivative_basis(1)
lift = lambda A: d3.Lift(A, lift_basis, -1)

gradu   = grad(u) + rvec*lift(tau_u1)          # first-order reduction
gradS   = grad(S) + rvec*lift(tau_S1)
gradlnr = grad(lnr)
x_th    = eps*S + gm1*lnr                       # linearized ln(T_tot/Tbar)
E_LHS = 0.5*(gradu + trans(gradu))
E_RHS = 0.5*(grad(u) + trans(grad(u)))
divu_LHS = trace(gradu)
divu_RHS = div(u)
sigma_LHS = 2*E_LHS + eye*divu_LHS
sigma_RHS = 2*E_RHS + eye*divu_RHS
Phi_S = (2.0*nu_v*c_phi) * iTbar * (trace(E_RHS @ E_RHS) - (1.0/3.0)*divu_RHS**2)

# ---------------------------------------------------------------- problem
problem = d3.IVP([u, S, lnr, tau_u1, tau_u2, tau_S1, tau_S2],
                 namespace=locals())

# MOMENTUM: linearized pressure/buoyancy in (S, lnr):
#   force_lin = gvec*(eps*S+gm1*lnr) + eps*Rg*Tbar*grad(S) + gamma*Rg*Tbar*grad(lnr)
# with gvec = R(dTbar+Tbar*rp) (hydrostatic identity). Exponential corrections on RHS.
problem.add_equation(
    "dt(u) - nu_v*(div(sigma_LHS) + rp@sigma_LHS)"
    " + gvec*(eps*S + gm1*lnr) + eps*Rg*Tbar*grad(S) + gamma*Rg*Tbar*gradlnr"
    " + lift(tau_u2)"
    " = - u@grad(u) - 2*OmegaN*cross(ez, u)"
    "   - (np.exp(x_th) - 1 - x_th)*gvec"
    "   - Rg*(np.exp(x_th) - 1)*Tbar*(eps*grad(S) + gamma*grad(lnr))"
    "   + nu_v*gradlnr@sigma_RHS")

# ENTROPY: dS/dt = (1/(rho*T)) div(kappa rho T grad S) + Phi/(rho T dS)
# implicit: kap*(lap S + (rp+tp)@grad S); explicit: kap*(gamma*gradlnr+eps*gradS)@gradS
problem.add_equation(
    "dt(S) - kap*(div(gradS) + (rp+tp)@grad(S)) + lift(tau_S2)"
    " = - u@grad(S)"
    "   + kap*((gamma*gradlnr + eps*grad(S))@grad(S))"
    "   + Phi_S")

# CONTINUITY (unchanged)
problem.add_equation("dt(lnr) + divu_LHS + rp@u = - u@gradlnr")

# BCs: stress-free impenetrable; benchmark entropy BCs S(ri)=1, S(ro)=0 (dS units)
problem.add_equation("radial(u(r=ri)) = 0")
problem.add_equation("angular(radial(E_RHS(r=ri))) = 0")
problem.add_equation("radial(u(r=ro)) = 0")
problem.add_equation("angular(radial(E_RHS(r=ro))) = 0")
problem.add_equation("S(r=ri) = 1.0")
problem.add_equation("S(r=ro) = 0.0")

# ----------------------------------------------------------------- solver
stepper = dict(CNAB2=d3.CNAB2, MCNAB2=d3.MCNAB2, SBDF2=d3.SBDF2,
               RK222=d3.RK222, RK443=d3.RK443)[STEPPER]
solver = problem.build_solver(stepper, ncc_cutoff=0, entry_cutoff=0)
solver.stop_sim_time = stop_time

CFL = d3.CFL(solver, initial_dt=initial_dt, cadence=1, safety=CFL_SAFETY,
             threshold=0.05, max_change=1.5, min_change=0.5, max_dt=max_dt)
CFL.add_velocity(u)

# --------------------------------- IC: Rayleigh init_type=9, isobaric balance
def Nlm(l, m):
    return np.sqrt((2*l+1)/(4*np.pi) * factorial(l-m)/factorial(l+m))

seedS  = 1.0e-4                                # entropy seed amplitude (dS units)
rfunc1 = seedS * (1.0 - np.cos(2.0*np.pi*(r - ri)/d_shell))
Y1919  = Nlm(19,19) * lpmv(19, 19, np.cos(theta)) * np.cos(19*phi)
Y11    = Nlm(1,1)   * lpmv(1, 1,  np.cos(theta)) * np.cos(phi)
if IC_MODE == "conductive":
    # steady conductive state of THIS system: rho*T*r^2*dS/dr = const
    _rq = np.linspace(ri, ro, 4001)
    _f  = 1.0/((c0 + c1*d_shell/_rq)**(npoly+1) * _rq**2)
    _J  = np.concatenate([[0.0], np.cumsum(0.5*(_f[1:]+_f[:-1])*np.diff(_rq))])
    from numpy import interp
    Sq = 1.0 - _J/_J[-1]
    Sprof = interp(r.ravel(), _rq, Sq)
    S['g'] = Sprof.reshape(r.shape) + rfunc1*(Y1919 + 0.1*Y11)
    # hydrostatically balanced l=0 lnr for the Scond profile:
    #   gamma*Rg*zeta*L' - (gamma-1)*g*L = eps*(g*Sq - Rg*zeta*Sq')
    _zq  = c0 + c1*d_shell/_rq
    _gq  = Rg*(npoly+1.0)*c1*d_shell/_rq**2
    _Sqp = np.gradient(Sq, _rq)
    _L   = np.zeros_like(_rq)
    for _i in range(len(_rq)-1):
        _h   = _rq[_i+1] - _rq[_i]
        _rhs = lambda LL, j: (gm1*_gq[j]*LL + eps*(_gq[j]*Sq[j] - Rg*_zq[j]*_Sqp[j]))/(gamma*Rg*_zq[j])
        _k1  = _rhs(_L[_i], _i)
        _k2  = _rhs(_L[_i] + _h*_k1, _i+1)
        _L[_i+1] = _L[_i] + 0.5*_h*(_k1 + _k2)
    # mass-neutral rebalance: subtract c*h(r), h = homogeneous hydrostatic mode
    # (pbar^-((gamma-1)/gamma)), c set so total mass = int rhobar dV exactly.
    _zq2  = c0 + c1*d_shell/_rq
    _h    = (_zq2**(npoly+1.0))**(-(gamma-1.0)/gamma)
    _h    = _h/_h[0]
    _w    = _zq2**npoly * _rq**2                       # rhobar r^2 weight
    def _mass_excess(cc):
        return np.trapezoid(_w*(np.exp(_L - cc*_h) - 1.0), _rq)
    _c = 0.0
    for _n in range(50):                               # Newton on the scalar
        _f  = _mass_excess(_c)
        _df = -np.trapezoid(_w*_h*np.exp(_L - _c*_h), _rq)
        _c -= _f/_df
        if abs(_f) < 1e-14*np.trapezoid(_w, _rq):
            break
    _L = _L - _c*_h
    _res = np.trapezoid(_w*(np.exp(_L) - 1.0), _rq)/np.trapezoid(_w, _rq)
    if dist.comm.rank == 0:
        logger.info(f"mass-neutral IC: c={_c:+.6e}  lnr(ri)={_L[0]:+.4f}  lnr(ro)={_L[-1]:+.4f}  residual dM/M={_res:+.2e}")
    lnr['g'] = interp(r.ravel(), _rq, _L).reshape(r.shape)
else:
    S['g'] = rfunc1 * (Y1919 + 0.1*Y11)        # paper Sec. 5: seed only
    lnr['g'] = 0.0
if restart:
    solver.load_state(restart_file)
    file_mode = 'append'
else:
    file_mode = 'overwrite'

# ------------------------------------------------ checkpoints + slices
CHECK_PARALLEL = 'gather'   # requires parallel h5py (as on Athena); use None for serial builds
checks = solver.evaluator.add_file_handler('checkpoints', sim_dt=0.01,
                                           max_writes=1, mode=file_mode,
                                           parallel=CHECK_PARALLEL)
checks.add_tasks(solver.state, layout='g')

# ---------------------------------------------------------------- monitors
flow = d3.GlobalFlowProperty(solver, cadence=MONITOR_CADENCE)
flow.add_property(d3.dot(u, u), name='u2')
flow.add_property(np.abs(lnr), name='alnr')
flow.add_property(d3.dot(u, u)/(gamma*Rg*Tbar_full), name='Ma2')
flow.add_property(np.abs(S), name='aS')

# Benchmark-definition integrals (Table 2 units; rho = rhobar per paper defs)
etheta = dist.VectorField(coords, name='etheta', bases=basis)
etheta['g'][1] = 1.0
uphi   = d3.dot(u, ep)
uphi_m = d3.Average(uphi, coords['phi'])
ur_f   = d3.dot(u, er)
uth_f  = d3.dot(u, etheta)
ur_m   = d3.Average(d3.dot(u, er), coords['phi'])
KE_op   = d3.integ(0.5*rhobar*d3.dot(u, u))
KEz_op  = d3.integ(0.5*rhobar*uphi_m**2)
KEm_op  = d3.integ(0.5*rhobar*(d3.Average(ur_f, coords['phi'])**2
                             + d3.Average(uth_f, coords['phi'])**2))
# point probes at mid-shell, equator, phi=0: the probe's oscillation period IS
# the drift period tau (target 0.0187440 -> omega = 2*pi/(19*tau) = 17.6427);
# during the linear phase its frequency should read 269.05.
rm = 0.5*(ri + ro)
Tprobe_op  = S(r=rm, theta=np.pi/2, phi=0.0)
upprobe_op = uphi(r=rm, theta=np.pi/2, phi=0.0)
# luminosity in paper units: L = -4 pi r^2 rhobar*Tbar*dS/dr (kap=1, S in dS units)
dSdr  = d3.dot(er, d3.grad(S))
wS    = np.exp(gamma*lnr + eps*S)                # true-flux wall weight rho*T/(rhobar*Tbar)
L_top = 4*np.pi*ro**2 * d3.Average((-rhobar*Tbar*wS*dSdr)(r=ro), s2coords)
L_bot = 4*np.pi*ri**2 * d3.Average((-rhobar*Tbar*wS*dSdr)(r=ri), s2coords)
s_cyl = dist.Field(name='s_cyl', bases=basis)
s_cyl['g'] = r*np.sin(theta)
rho_full = rhobar*np.exp(lnr)
Lz_op = d3.integ(rho_full*s_cyl*uphi)            # integ rho r sin(theta) u_phi
one_f = dist.Field(name='one_f', bases=basis); one_f['g'] = 1.0
M_op  = d3.integ(rhobar*np.exp(lnr)*one_f)       # total mass (drift meter)
M0    = 4*np.pi*np.trapezoid((c0 + c1*d_shell/np.linspace(ri,ro,4001))**npoly
                             * np.linspace(ri,ro,4001)**2, np.linspace(ri,ro,4001))
# equatorial x,y components (paper A5-A6): the dangerous ones under drift
cx_th = dist.Field(name='cx_th', bases=basis); cx_th['g'] = -r*np.sin(phi)
cx_ph = dist.Field(name='cx_ph', bases=basis); cx_ph['g'] = -r*np.cos(theta)*np.cos(phi)
cy_th = dist.Field(name='cy_th', bases=basis); cy_th['g'] =  r*np.cos(phi)
cy_ph = dist.Field(name='cy_ph', bases=basis); cy_ph['g'] = -r*np.cos(theta)*np.sin(phi)
Lx_op = d3.integ(rho_full*(cx_th*uth_f + cx_ph*uphi))
Ly_op = d3.integ(rho_full*(cy_th*uth_f + cy_ph*uphi))
# shell moment of inertia (rhobar weight): I = (8 pi/3) int rhobar r^4 dr
_rq = np.linspace(ri, ro, 4001)
I_shell = (8*np.pi/3) * np.trapezoid((c0 + c1*d_shell/_rq)**npoly * _rq**4, _rq)

def am_correct():
    Lx, Ly, Lz = scalar(Lx_op), scalar(Ly_op), scalar(Lz_op)
    wx, wy, wz = Lx/I_shell, Ly/I_shell, Lz/I_shell
    u.change_scales(1)
    u['g'][1] -= r*(-wx*np.sin(phi) + wy*np.cos(phi))
    u['g'][0] -= r*( wz*np.sin(theta) - wx*np.cos(theta)*np.cos(phi)
                                      - wy*np.cos(theta)*np.sin(phi))
    return Lx, Ly, Lz

# equatorial mid-shell ring output for drift-rate and ur=0 point values
slices = solver.evaluator.add_file_handler('eq_ring', sim_dt=2e-4,
                                           max_writes=400, mode=file_mode)
slices.add_task(d3.dot(u, er)(r=rm, theta=np.pi/2), name='ur_eq')
slices.add_task(uphi(r=rm, theta=np.pi/2), name='uphi_eq')
slices.add_task(S(r=rm, theta=np.pi/2), name='S_eq')

def scalar(op):
    v = op.evaluate().allgather_data(layout='g').ravel()
    return float(v[0]) if len(v) else float('nan')

logger.info(f"STEPPER={STEPPER} CFL safety={CFL_SAFETY} max_dt={max_dt:.2e}")
logger.info(f"checks vs Jones Table 1: zeta_o={zeta_o:.6f} (0.256465) "
            f"c0={c0:.6f} (-1.287800) c1={c1:.6f} (2.375792) "
            f"Tc={Tc_c:.1f} K (111557) rho_c={rho_c_mid:.6f} (0.112684) "
            f"DeltaS={DeltaS_c:.1f} (851225.7) DeltaT_hat={DeltaT:.6f}")

timestep = initial_dt
while solver.proceed:
    timestep = CFL.compute_timestep()
    solver.step(timestep)
    if AM_CORRECT and solver.iteration % AM_CADENCE == 0:
        am_correct()
    if solver.iteration % MONITOR_CADENCE == 0:
        KE, KEz, KEm = scalar(KE_op), scalar(KEz_op), scalar(KEm_op)
        Tp, upp = scalar(Tprobe_op), scalar(upprobe_op)
        Lt, Lb  = scalar(L_top), scalar(L_bot)
        Lz      = scalar(Lz_op)
        dMoM    = scalar(M_op)/M0 - 1.0
        logger.info(f"it={solver.iteration:7d} t={solver.sim_time:.6f} "
                    f"dt={timestep:.2e} "
                    f"KE={KE:.4e} KEz={KEz:.4e} KEm={KEm:.4e} "
                    f"Tp={Tp:+.5e} up={upp:+.5e} "
                    f"L_top={Lt:.4f} L_bot={Lb:.4f} Lz={Lz:+.3e} dM/M={dMoM:+.2e} "
                    f"max|u|={np.sqrt(flow.max('u2')):.3e} "
                    f"max|lnr|={flow.max('alnr'):.3e} "
                    f"Ma={np.sqrt(flow.max('Ma2')):.3e}")
