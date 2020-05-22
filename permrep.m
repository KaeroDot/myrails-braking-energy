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
        P = R^C;
        if C > 19
                % this method generates random possibilities:
                warning('Switching to Monte Carlo Method!')
                if R ~= 2
                        % because I am lazy to do it generally:
                        error('MCM works only for two rows!')
                end
                M = 1e3;
                b = round(rand(M, C));
                pr = b.*repmat(E(1,:), M, 1) + not(b).*repmat(E(2,:), M, 1);
                % for 19 it takes ~11 seconds on iCore 7
                % for 20 it takes ~23 seconds on iCore 7
                % for 21 it takes ~49 seconds on iCore 7
        else
                % this method generates all possibilities:
                ids = dec2base(0:P-1, R) - '0' + 1;
                pr = zeros(P, size(E, 2));
                for i = 1:P
                        pr(i, :) = diag(E(ids(i, :), :))';
                end
        end
end
