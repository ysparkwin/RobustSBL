%% GEN_BS generates translation parameters within a desired range and spacing.
%
% [ BS ] = GEN_BS(MIN_B, MAX_B, K, B_SEP)
%
% Input variables:
%   MIN_B   The minimum value of b.
%   MAX_B   The maximum value of b.
%   K       The number of translation parameters to generate.
%   B_SEP   The minimum spacing between parameters.
%
% Output variables:
%   BS      The generated parameters as a vector.
% 
% Code implemented by: Karsten Fyhn
% Contact e-mail: kfn@es.aau.dk
%
% Copyright 2013 Karsten Fyhn
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%   http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distribued on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
function [ bs ] = gen_bs(min_b, max_b, K, b_sep)
    bs = zeros(K,1);
    for ii = 1:K
        tmp_b = (max_b-min_b)*rand(1) + min_b;
        jj = 1;
        while sum(abs(bs - tmp_b) < b_sep) > 0
            tmp_b = (max_b-min_b)*rand(1) + min_b;
            jj = jj + 1;
            if jj > 100
                error('It was not possible to generate K betas with the chosen b_sep!')
            end
        end
        bs(ii) = tmp_b;
    end
end

