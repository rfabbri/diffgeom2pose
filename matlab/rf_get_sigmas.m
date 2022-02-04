% To be called after rf_pose_from_point_tangents_root_find.
% This is modeled after peter's problem1-thoughts2.pdf.
%
% Inputs ts: vector with t's that are roots for the big polynomial equation
% Output: sigmas1,2(i) for each t ts(i);  sigmas(i) == 0 if ts(i) is not valid


my_eps = 1;

sigmas1 = cell(length(ts),1);
sigmas2 = cell(length(ts),1);
for i=1:length(ts)

  sigmas1{i} = [];
  sigmas2{i} = [];

  [fvalue,A,B,C,E,F,G,H,J,K,L] = rf_pose_from_point_tangents_2_fn_t(ts(i));

  delta1 = sqrt(B*B - 4*A*C);
  sigma1_m = (-B - delta1)/(2*A);
  sigma1_p = (-B + delta1)/(2*A);

  delta2 = sqrt(F*F - 4*E*G);
  sigma2_m = (-F - delta2)/(2*E);
  sigma2_p = (-F + delta2)/(2*E);

  % handle case of negative delta
  if abs(imag(sigma1_m)) < 1e-4
    sigma1_m = real(sigma1_m);
    sigma1_p = real(sigma1_p);
  else
    warning(['Ignoring t = ' num2str(ts(i))]);
    continue;
  end

  if abs(imag(sigma2_m)) < 1e-4
    sigma2_m = real(sigma2_m);
    sigma2_p = real(sigma2_p);
  else
    warning(['Ignoring t = ' num2str(ts(i))]);
    continue;
  end


  % Now check to see which pair pass. Only a single pair should pass, in theory.
  % If not, issue a warning.

  if (abs(H + J*sigma1_m + K*sigma2_m + L*sigma1_m*sigma2_m) < my_eps)
    sigmas1{i}(end+1) = sigma1_m;
    sigmas2{i}(end+1) = sigma2_m;
  end

  if (abs(H + J*sigma1_p + K*sigma2_m + L*sigma1_p*sigma2_m) < my_eps)
    if (~isempty(sigmas1{i}))
      warning('more than one sigma1, sigma2 pair satisfies the 3rd constraint');
    end
    sigmas1{i}(end+1) = sigma1_p;
    sigmas2{i}(end+1) = sigma2_m;
  end

  if (abs(H + J*sigma1_p + K*sigma2_p + L*sigma1_p*sigma2_p) < my_eps)
    if (~isempty(sigmas1{i}))
      warning('more than one sigma1, sigma2 pair satisfies the 3rd constraint');
    end
    sigmas1{i}(end+1) = sigma1_p;
    sigmas2{i}(end+1) = sigma2_p;
  end

  if (abs(H + J*sigma1_m + K*sigma2_p + L*sigma1_m*sigma2_p) < my_eps)
    if (~isempty(sigmas1{i}))
      warning('more than one sigma1, sigma2 pair satisfies the 3rd constraint');
    end
    sigmas1{i}(end+1) = sigma1_m;
    sigmas2{i}(end+1) = sigma2_p;
  end

  if (isempty(sigmas1{i}))
    warning('no sigma1, sigma2 pair satisfies the 3rd constraint');
  end
end


