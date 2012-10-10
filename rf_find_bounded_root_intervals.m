function [root_ids,sampled_poly] = rf_find_bounded_root_intervals(t_vector)
% finds bounded root intervals of the function rf_pose_from_point_tangents_2_fn

root_ids = zeros(1,(length(t_vector)-1));


%sampled_poly = zeros(size(t_vector));
%
%for i=1:length(t_vector)
%  sampled_poly(i) = rf_pose_from_point_tangents_2_fn_t_for_root(t_vector(i));
%end

sampled_poly = rf_sample_pose_poly(t_vector);

curr_val = sampled_poly(1);
for i=1:length(root_ids)
  nxt_val = sampled_poly(i+1);
  root_ids(i) = (curr_val*nxt_val) < 0;
  curr_val = nxt_val;
end

