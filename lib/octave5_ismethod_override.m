function TorF = octave5_ismethod_override(obj, name)
%OCTAVE5_ISMETHOD_OVERRIDE  Replacement for built-in ISMETHOD() for Octave < 6
%
%   In GNU Octave versions before 6.0 the built-in ISMETHOD does not find
%   methods defined in the classdef file. This implements a replacement
%   that does. Typical usage:
%
%       if have_feature('octave', 'vnum') < 6
%           my_ismethod = @(o, n)octave5_ismethod_override(o, n);
%       else
%           my_ismethod = @ismethod;
%       end
%
%   See also ISMETHOD.

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

meths = meta.class.fromName(class(obj)).MethodList;
for k = 1:length(meths)
    if strcmp(meths{k}.Name, name)
        TorF = true;
        return;
    end
end
outTorF = false;
