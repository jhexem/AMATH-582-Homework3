clear all

%% Ideal case
load('Xt1_1.mat')
load('Xt2_1.mat')
load('Xt3_1.mat')

% Each of these pairs of vectors will be different lengths, so what should
% you do to not get dimension errors?

L1 = length(Xt1_1);
L2 = length(Xt2_1);
L3 = length(Xt3_1);
minL = min([L1, L2, L3]);

Xt1_1 = Xt1_1(:, 1:minL);
Xt2_1 = Xt2_1(:, 1:minL);
Xt3_1 = Xt3_1(:, 1:minL);

% After taking care of the dimension issue, subtract the mean off of each
% component; i.e. each row vector for each video (see documentation for the
% "mean" function and becareful of the dimensions to make sure it's done
% correctly because it's finicky).

mean1_1 = mean(Xt1_1, 2);
mean2_1 = mean(Xt2_1, 2);
mean3_1 = mean(Xt3_1, 2);

newXt1_1 = Xt1_1 - mean1_1;
newXt2_1 = Xt2_1 - mean2_1;
newXt3_1 = Xt3_1 - mean3_1;

% Combine the pairs of vectors into one 6xL matrix, where L is the length,
% which will vary depending on the video.
% Take the svd and do a coordinate transformation to the principal
% components just like we did in the lectures (we called this matrix Y).
% Save the matrix Y as A1.

X = [newXt1_1; newXt2_1; newXt3_1];
[U, S, V] = svd(X);
Y = U' * X;

A1 = Y; % Shape: 6x226 double

% Save the energies of nontrivial the singular values as a vector in the
% varialbe A2.

sigmas = diag(S);
enrg = (sigmas.^2) ./ sum(sigmas.^2);

A2 = enrg; % Shape: 6x1 double

%%% Plot the energies with a log-scaled ordinate for the report (not for 
% the autograder)

%plot(cumsum(enrg), "o")
%{
semilogy(enrg, "r.", "MarkerSize", 25)
axis([0 7 10^(-3) 10^(0)])
set(gca,'Fontsize',16,'Xtick',0:1:7,'Ytick',logspace(-3,0, 4))
%}
% The appropriate rank necessary for a decent approximation will be
% informed by the energies and by judging the qualitative convergence (no
% need for rigorous analysis, you can just plot and check) as we increase
% the rank of the approximation.
% Save the rank-n approximation as A3.
% Of course, this will be subjective (for now), so you may have to try it
% a few times to match what I picked.  You also may not agree with my
% choice, or your mass tracking may be different based on what you clicked.
% If you get it to pass the autograder, but you disagree put it in the
% report.

approx = zeros(6, 6, minL);
for k = 1:6
    approxTemp = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    approx(k, :, :) = approxTemp;
    %sizeU = size(U(:,1:k))
    %sizeS = size(S(1:k,1:k))
    %sizeVT = size(V(:,1:k)')
end

A3 = reshape(approx(2, :, :), [6, minL]); % Shape: 6x226 double

%%% Plot the approximations from SVD for report (not for autograder)
% Plot each rank-n approximation (up to one more than A12) to observe the
% convergence
%{
for j = 1:3
    plot(reshape(approx(j, 1, :), [1, minL]), "linewidth", 1)
    hold on
end
legend(["rank-1", "rank-2", "rank-3"])
%}
%% Test 2
load('Xt1_2.mat')
load('Xt2_2.mat')
load('Xt3_2.mat')

% Repeat what you did for the Ideal Case
% Each of these pairs of vectors will be different lengths, so what should
% you do to not get dimension errors?

L1 = length(Xt1_2);
L2 = length(Xt2_2);
L3 = length(Xt3_2);
minL = min([L1, L2, L3]);

Xt1_2 = Xt1_2(:, 1:minL);
Xt2_2 = Xt2_2(:, 1:minL);
Xt3_2 = Xt3_2(:, 1:minL);

% After taking care of the dimension issue, subtract the mean off of each
% component; i.e. each row vector for each video (see documentation for the
% "mean" function and becareful of the dimensions to make sure it's done
% correctly because it's finicky).

mean1_2 = mean(Xt1_2, 2);
mean2_2 = mean(Xt2_2, 2);
mean3_2 = mean(Xt3_2, 2);

newXt1_2 = Xt1_2 - mean1_2;
newXt2_2 = Xt2_2 - mean2_2;
newXt3_2 = Xt3_2 - mean3_2;

% Combine the pairs of vectors into one 6xL matrix, where L is the length,
% which will vary depending on the video.
% Take the svd and do a coordinate transformation to the principal
% components just like we did in the lectures (we called this matrix Y).
% Save the matrix Y as A4.

X = [newXt1_2; newXt2_2; newXt3_2];
[U, S, V] = svd(X);
Y = U' * X;

A4 = Y; % 6x314 double

% Save the energies of nontrivial the singular values as a vector in the
% varialbe A5.

sigmas = diag(S);
enrg = (sigmas.^2) ./ sum(sigmas.^2);

A5 = enrg; % 6x1 double

%%% Plot the energies with a log-scaled ordinate for the report (not for 
% the autograder)

%plot(cumsum(enrg), "o")
%{
semilogy(enrg, "r.", "MarkerSize", 25)
axis([0 7 10^(-3) 10^(0)])
set(gca,'Fontsize',16,'Xtick',0:1:7,'Ytick',logspace(-3,0, 4))
%}
% The appropriate rank necessary for a decent approximation will be
% informed by the energies and by judging the qualitative convergence (no
% need for rigorous analysis, you can just plot and check) as we increase
% the rank of the approximation.
% Save the rank-n approximation as A6.
% Of course, this will be subjective (for now), so you may have to try it
% a few times to match what I picked.  You also may not agree with my
% choice, or your mass tracking may be different based on what you clicked.
% If you get it to pass the autograder, but you disagree put it in the
% report.

approx = zeros(6, 6, minL);
for k = 1:6
    approxTemp = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    approx(k, :, :) = approxTemp;
end

A6 = reshape(approx(3, :, :), [6, minL]); % 6x314 double

%%% Plot the approximations from SVD for report (not for autograder)
% Plot each rank-n approximation (up to one more than A12) to observe the
% convergence
%{
for j = 1:4
    plot(reshape(approx(j, 1, :), [1, minL]), "linewidth", 1)
    hold on
end
legend(["rank-1", "rank-2", "rank-3", "rank-4"])
%}
%% Test 3
load('Xt1_3.mat')
load('Xt2_3.mat')
load('Xt3_3.mat')

% Repeat what you did for the Ideal Case
% Each of these pairs of vectors will be different lengths, so what should
% you do to not get dimension errors?

L1 = length(Xt1_3);
L2 = length(Xt2_3);
L3 = length(Xt3_3);
minL = min([L1, L2, L3]);

Xt1_3 = Xt1_3(:, 1:minL);
Xt2_3 = Xt2_3(:, 1:minL);
Xt3_3 = Xt3_3(:, 1:minL);

% After taking care of the dimension issue, subtract the mean off of each
% component; i.e. each row vector for each video (see documentation for the
% "mean" function and becareful of the dimensions to make sure it's done
% correctly because it's finicky).

mean1_3 = mean(Xt1_3, 2);
mean2_3 = mean(Xt2_3, 2);
mean3_3 = mean(Xt3_3, 2);

newXt1_3 = Xt1_3 - mean1_3;
newXt2_3 = Xt2_3 - mean2_3;
newXt3_3 = Xt3_3 - mean3_3;

% Combine the pairs of vectors into one 6xL matrix, where L is the length,
% which will vary depending on the video.
% Take the svd and do a coordinate transformation to the principal
% components just like we did in the lectures (we called this matrix Y).
% Save the matrix Y as A7.

X = [newXt1_3; newXt2_3; newXt3_3];
[U, S, V] = svd(X);
Y = U' * X;

A7 = Y; % Shape: 6x237 double

% Save the energies of nontrivial the singular values as a vector in the
% varialbe A8.

sigmas = diag(S);
enrg = (sigmas.^2) ./ sum(sigmas.^2);

A8 = enrg; % Shape: 6x1 double

%%% Plot the energies with a log-scaled ordinate for the report (not for 
% the autograder)

%plot(cumsum(enrg), "o")
%{
semilogy(enrg, "r.", "MarkerSize", 25)
axis([0 7 10^(-3) 10^(0)])
set(gca,'Fontsize',16,'Xtick',0:1:7,'Ytick',logspace(-3,0, 4))
%}
% The appropriate rank necessary for a decent approximation will be
% informed by the energies and by judging the qualitative convergence (no
% need for rigorous analysis, you can just plot and check) as we increase
% the rank of the approximation.
% Save the rank-n approximation as A9.
% Of course, this will be subjective (for now), so you may have to try it
% a few times to match what I picked.  You also may not agree with my
% choice, or your mass tracking may be different based on what you clicked.
% If you get it to pass the autograder, but you disagree put it in the
% report.

approx = zeros(6, 6, minL);
for k = 1:6
    approxTemp = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    approx(k, :, :) = approxTemp;
end

A9 = reshape(approx(3, :, :), [6, minL]); % 6x237 double

%%% Plot the approximations from SVD for report (not for autograder)
% Plot each rank-n approximation (up to one more than A12) to observe the
% convergence
%{
for j = 1:4
    plot(reshape(approx(j, 1, :), [1, minL]), "linewidth", 1)
    hold on
end
legend(["rank-1", "rank-2", "rank-3", "rank-4"])
%}
%% Test 4
load('Xt1_4.mat')
load('Xt2_4.mat')
load('Xt3_4.mat')

% Repeat what you did for the Ideal Case
% Each of these pairs of vectors will be different lengths, so what should
% you do to not get dimension errors?

L1 = length(Xt1_4);
L2 = length(Xt2_4);
L3 = length(Xt3_4);
minL = min([L1, L2, L3]);

Xt1_4 = Xt1_4(:, 1:minL);
Xt2_4 = Xt2_4(:, 1:minL);
Xt3_4 = Xt3_4(:, 1:minL);

% After taking care of the dimension issue, subtract the mean off of each
% component; i.e. each row vector for each video (see documentation for the
% "mean" function and becareful of the dimensions to make sure it's done
% correctly because it's finicky).

mean1_4 = mean(Xt1_4, 2);
mean2_4 = mean(Xt2_4, 2);
mean3_4 = mean(Xt3_4, 2);

newXt1_4 = Xt1_4 - mean1_4;
newXt2_4 = Xt2_4 - mean2_4;
newXt3_4 = Xt3_4 - mean3_4;

% Combine the pairs of vectors into one 6xL matrix, where L is the length,
% which will vary depending on the video.
% Take the svd and do a coordinate transformation to the principal
% components just like we did in the lectures (we called this matrix Y).
% Save the matrix Y as A10.

X = [newXt1_4; newXt2_4; newXt3_4];
[U, S, V] = svd(X);
Y = U' * X;

A10 = Y; % 6x392 double

% Save the energies of nontrivial the singular values as a vector in the
% varialbe A11.

sigmas = diag(S);
enrg = (sigmas.^2) ./ sum(sigmas.^2);

A11 = enrg; % 6x1 double

%%% Plot the energies with a log-scaled ordinate for the report (not for 
% the autograder)

%plot(cumsum(enrg), "o")
%{
semilogy(enrg, "r.", "MarkerSize", 25)
axis([0 7 10^(-3) 10^(0)])
set(gca,'Fontsize',16,'Xtick',0:1:7,'Ytick',logspace(-3,0, 4))
%}
% The appropriate rank necessary for a decent approximation will be
% informed by the energies and by judging the qualitative convergence (no
% need for rigorous analysis, you can just plot and check) as we increase
% the rank of the approximation.
% Save the rank-n approximation as A12.
% Of course, this will be subjective (for now), so you may have to try it
% a few times to match what I picked.  You also may not agree with my
% choice, or your mass tracking may be different based on what you clicked.
% If you get it to pass the autograder, but you disagree put it in the
% report.

approx = zeros(6, 6, minL);
for k = 1:6
    approxTemp = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    approx(k, :, :) = approxTemp;
end

A12 = reshape(approx(2, :, :), [6, minL]); % 6x392 double

%%% Plot the approximations from SVD for report (not for autograder)
% Plot each rank-n approximation (up to one more than A12) to observe the
% convergence
%{
for j = 1:3
    plot(reshape(approx(j, 4, :), [1, minL]), "linewidth", 1)
    hold on
end
legend(["rank-1", "rank-2", "rank-3"])
%}