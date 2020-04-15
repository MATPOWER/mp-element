function success = test_mp_element(verbose, exit_on_fail)
%T_MP_ELEMENT  Run all MATPOWER tests.
%   T_MP_ELEMENT
%   T_MP_ELEMENT(VERBOSE)
%   T_MP_ELEMENT(VERBOSE, EXIT_ON_FAIL)
%   SUCCESS = T_MP_ELEMENT(...)
%
%   Runs all of the MP-Element tests. If VERBOSE is true (false by default),
%   it prints the details of the individual tests. If EXIT_ON_FAIL is true
%   (false by default), it will exit MATLAB or Octave with a status of 1
%   unless T_RUN_TESTS returns ALL_OK.
%
%   See also T_RUN_TESTS.

%   MATPOWER
%   Copyright (c) 2004-2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 2
    exit_on_fail = 0;
    if nargin < 1
        verbose = 0;
    end
end

tests = {};

%% MATPOWER base test
tests{end+1} = 't_mp_element';
tests{end+1} = 't_acp_port_inj_current';
tests{end+1} = 't_acp_port_inj_power';
tests{end+1} = 't_acc_port_inj_current';
tests{end+1} = 't_acc_port_inj_power';
tests{end+1} = 't_acp_nln_port_inj_current';
tests{end+1} = 't_acp_nln_port_inj_power';
tests{end+1} = 't_acc_nln_port_inj_current';
tests{end+1} = 't_acc_nln_port_inj_power';

%% run the tests
all_ok = t_run_tests( tests, verbose );

%% handle success/failure
if exit_on_fail && ~all_ok
    exit(1);
end
if nargout
    success = all_ok;
end
