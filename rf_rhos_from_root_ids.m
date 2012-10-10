function [rhos1,rhos1_minus,rhos1_plus,rhos2,rhos2_minus,rhos2_plus, ts] = rf_rhos_from_root_ids(t_vector, root_ids)
% input root_ids is the output from rf_find_bounded_roots
% TODO: output rho2 too.
% usually called inside rf_pose_from_point_tangents_root_find_function_any.m

global alpha beta theta

ts = [];
for i=1:length(root_ids)
  if root_ids(i) == 1
%    t_ini = (t_vector(i) + t_vector(i+1))/2.0;
%    t_ini
    % refine the root via a root solver with initial value in the middle of the
    % interval. This doesn't improve the solution too much in examples I tried.
%    disp('P(t): ');
    rf_pose_from_point_tangents_2_fn_t_for_root(t_vector(i));
    rf_pose_from_point_tangents_2_fn_t_for_root(t_vector(i+1));
    t_ref = fzero(@rf_pose_from_point_tangents_2_fn_t_for_root, [t_vector(i) t_vector(i+1)]);
    rf_pose_from_point_tangents_2_fn_t_for_root(t_ref);
    ts(end+1) = t_ref;
  end
end
%ts

t_stddev = t_vector(2) - t_vector(1);


% Each root is now ts(i), plus minus t_stddev.
% Now get rho1(t):

rhos1 = 2*alpha*ts*cos(theta) + beta*(ones(size(ts))-ts.*ts)*sin(theta);
rhos1 = rhos1./(1 + ts.*ts);

ts_orig = ts;

ts = ts-t_stddev;
rhos1_minus = 2*alpha*ts*cos(theta) + beta*(ones(size(ts))-ts.*ts)*sin(theta);
rhos1_minus = rhos1_minus./(1 + ts.*ts);

ts = ts+2*t_stddev;
rhos1_plus = 2*alpha*ts*cos(theta) + beta*(ones(size(ts))-ts.*ts)*sin(theta);
rhos1_plus = rhos1_plus./(1 + ts.*ts);

ts = ts_orig;
rhos2 = -2*alpha*ts*sin(theta) + beta*(ones(size(ts))-ts.*ts)*cos(theta);
rhos2 = rhos2./(1 + ts.*ts);

ts = ts-t_stddev;
rhos2_minus = -2*alpha*ts*sin(theta) + beta*(ones(size(ts))-ts.*ts)*cos(theta);
rhos2_minus = rhos2_minus./(1 + ts.*ts);

ts = ts+2*t_stddev;
rhos2_plus = -2*alpha*ts*sin(theta) + beta*(ones(size(ts))-ts.*ts)*cos(theta);
rhos2_plus = rhos2_plus./(1 + ts.*ts);

ts = ts_orig;
