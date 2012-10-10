% to be called from rf_pose_from_point_tangents_root_find_function_any.m

% Lambdas:

lambdas1 = cell(length(ts),1);
lambdas2 = cell(length(ts),1);
for i=1:length(ts)
  if isempty(sigmas1{i})
    lambdas1{i} = [];
    lambdas2{i} = [];
  end
  lambdas1{i} = zeros(size(sigmas1{i}));
  lambdas2{i} = zeros(size(sigmas2{i}));
  for k=1:length(sigmas1{i})
    lambdas1{i}(k) = ((Gama1 - Gama2)'*Tgt1)/(...
    (rhos1(i)*gama1 - rhos2(i)*gama2)'*(rhos1(i)*tgt1 + sigmas1{i}(k)*gama1));

    lambdas2{i}(k) = (Gama1 - Gama2)'*Tgt2/(...
    (rhos1(i)*gama1 - rhos2(i)*gama2)'*(rhos2(i)*tgt2 + sigmas2{i}(k)*gama2));
  end
end

% Rotation:

% RA = B => R = B/A
% TODO: use svd or some other way to be sure R is unique and orthogonal.
% right now just testing det(R) = 1

Rots = {};
Transls = {};
A = [Gama1 - Gama2 Tgt1 Tgt2];
for i=1:length(ts)
  for k=1:length(sigmas1{i}) 
    B = [ rhos1(i)*gama1 - rhos2(i)*gama2   lambdas1{i}(k)*(rhos1(i)*tgt1 + sigmas1{i}(k)*gama1) ...
                                            lambdas2{i}(k)*(rhos2(i)*tgt2 + sigmas2{i}(k)*gama2) ];

%    disp('Determinants:');
%    disp('det A:');
%    det(A)
%    disp('det B:');
%    det(B)
%    disp '---------'
    Rots{end+1} = B/A;
%    disp '---------'

    % be sure it's close to 1:
    det(Rots{end});

    % force orthogonality (not needed anymore)
%    [U,S,V]=svd(Rots{end});
%    Rots{end} = U*V';

    Transls{end+1} = rhos1(i)*gama1 - Rots{end}*Gama1;
    % this should be the same:
    rhos2(i)*gama2 - Rots{end}*Gama2;
  end
end


