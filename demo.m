%% Check 2 x 1 1/2 points SRS, following Fabbri et al. 2019
%% Special thanks to anonymous reviewer for this short demo
%% Simply runs on completely random data.
display('--------------------------------------')
display('--- demo point pairs with tangents ---')
display('--------------------------------------')
%% Functions
v2skew =@(v) [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0]; % from vector to skew matrix
v2Rot  =@(v) expm(v2skew(v)); % from vector to rotation
skew2v =@(S) [S(3,2);S(1,3);S(2,1)]; % from skew matrix to vector
%% Generate data
% true rotation and translation
R_tilde = v2Rot(randn(3,1));
T_tilde = [0,0,10]'+randn(3,1);
% 2 scene points
X = randn(2,3);
% 2 directions
T = randn(2,3);
% normalization (not necessary, taken care by G's, reduces # of solutions)
T = T ./ ([norm(T(1,:));norm(T(2,:))]*[1,1,1]);
% 2 scene points Gamma^w and normalized directions
Gama1 = X(1,:)';
Gama2 = X(2,:)';
Tgt1  = T(1,:)';
Tgt2  = T(2,:)';
% points and neighbours, path on tangent
epsi = 0.1;
Y = [X;X+epsi*T]
% four points Gamma in camera system
y  = (R_tilde*Y' + T_tilde*ones(1,4))'
% two camera points (last element = 1)
gama  = y(:,1:3)./(y(:,3)*[1,1,1]);
gama1 = gama(1,:)'
gama2 = gama(2,:)'
% camera tangents (last element = 0)
tgt1 =[(gama(3,:)-gama(1,:))'];
tgt2 =[(gama(4,:)-gama(2,:))'];
% normalization (not necessary, taken care by G's, reduces # of solutions)
tgt1 = tgt1/norm(tgt1)
tgt2 = tgt2/norm(tgt2)
%% call P2P
tic
[Rots,Transls,degen] = rf_pose_from_point_tangents_root_find_function_any(gama1,tgt1,gama2,tgt2,Gama1,Tgt1,Gama2,Tgt2);
toc
%% check rotation and translation
N = length(Rots);
number_of_solutions=N
for n = 1: N
    dR = norm(skew2v(Rots{n}*R_tilde'));
    dT = norm(Transls{n}-T_tilde);
    if dR+dT < 10^-4 
        display([num2str(n),'th solution'])
        display(['dR =',num2str(dR)]);
        display(['dT =',num2str(dT)]);
        %vectors_vR_vT = [calc_v_from_S(Rots{n}*R_tilde'),Transls{n}-T_tilde]
    end
end
