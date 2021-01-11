classdef mp_form_accs < mp_form_acc
%MP_FORM_ACCS  MATPOWER Formulation class for AC cartesian voltage, current bal
%   Each concrete Network Model Element class must inherit, at least
%   indirectly, from both NM_ELEMENT and MP_FORM.
%
%   Subclass of MP_FORM_ACC.
%   MP_FORM provides properties and methods related to the specific
%   formulation (e.g. DC version, AC polar power version, etc.)
%
%   Properties
%       (model parameters inherited from MP_FORM_AC)
%
%   Methods
%       form_name() - returns string w/name of formulation ('AC-cartesian formulation')
%       form_tag() - returns string w/short label for formulation ('acc')

%   MATPOWER
%   Copyright (c) 2019-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = form_name(obj)
            name = 'AC-cartesian-power';
        end
        function tag = form_tag(obj)
            tag = 'accs';
        end
    end     %% methods
end         %% classdef
