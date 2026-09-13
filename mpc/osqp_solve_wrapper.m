function x = osqp_solve_wrapper(P, q, A_mat, l, u)
%OSQP_SOLVE_WRAPPER  Thin OSQP interface used by both MPC step functions.
%
%  Solves  min ½x'Px + q'x  s.t.  l <= A_mat*x <= u.
%  Returns x = [] when OSQP reports infeasibility or hits max_iter, which the
%  callers read as "hold the valves closed for this step".

solver = osqp;
solver.setup(sparse(P), q, sparse(A_mat), l, u, ...
    'warm_start',     false, ...
    'verbose',        false, ...
    'eps_abs',        1e-4,  ...
    'eps_rel',        1e-4,  ...
    'max_iter',       4000,  ...
    'polish',         true);

res = solver.solve();

if strcmp(res.info.status, 'solved') || ...
   strcmp(res.info.status, 'solved_inaccurate')
    x = res.x;
else
    x = [];
end
