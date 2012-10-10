function [Rots, Transls,degen] =rf_pose_from_point_tangents_root_find_function_any(gama1,tgt1,gama2,tgt2,Gama1,Tgt1,Gama2,Tgt2)
% This is the main routine to find roots. Can be used with any input.

% test for geometric degeneracy -------------------------------

DGama = Gama1 - Gama2;
DGama = DGama/norm(DGama);
degen = det([DGama Tgt1 Tgt2]);
if (abs(degen) < 1e-3)
  disp('data point not reliable');
  Rots = {};
  Transls = {};
  return;
end

% compute roots -------------------------------


t_vector=-1:0.001:1;
rf_pose_from_point_tangents_2;

root_ids = rf_find_bounded_root_intervals(t_vector);

% compute rhos, r, t --------------------------

[rhos1,rhos1_minus,rhos1_plus,rhos2,rhos2_minus,rhos2_plus, ts] = ...
    rf_rhos_from_root_ids(t_vector, root_ids);

rf_get_sigmas;

rf_get_r_t_from_rhos;

