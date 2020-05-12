% Generates all permutations with repetition of matrix along columns. If input matrix size is (R, C)
% Output matrix size is (R^C, C).
% Example:
% for matrix:
%      1  2  3
%     10 20 30
% generates permutations with repetitions:
%      1    2    3
%      1    2   30
%      1   20    3
%      1   20   30
%     10    2    3
%     10    2   30
%     10   20    3
%     10   20   30

function pr = permrep(E);
        R = size(E, 1);
        C = size(E, 2);
        if C > 19
                % for 19 it takes ~11 seconds on iCore 7
                % for 20 it takes ~23 seconds on iCore 7
                % for 21 it takes ~49 seconds on iCore 7
                disp('this will take a long time...')
        endif
        P = R^C;
        ids = dec2base(0:P-1, R) - '0' + 1;
        pr = zeros(P, size(E, 2));
        for i = 1:P
                pr(i, :) = diag(E(ids(i, :), :))';
        end
end
