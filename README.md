MP-Element
==========

MP-Element is a new, generalized network and element modeling layer for
[MATPOWER][1]. It is currently **under active development** in its own
separate repository, with the intention of eventually being incorporated
into the core MATPOWER code.


System Requirements
-------------------
*   [MATLAB][2] version 9.1 (R2016b) or later, or
*   [GNU Octave][3] version 4.4 or later
*   A development version of [MATPOWER][1] that explicitly includes support
    for MP-Element  
    _See the [`mp-element`][8] branch of the MATPOWER GitHub repository,
    including the corresponding [CHANGES][9] file._


Installation
------------

Installation and use of MP-Element requires familiarity with the basic operation
of MATLAB or Octave, including setting up your MATLAB/Octave path.

1.  Clone the repository or download and extract the zip file of the MP-Element
    distribution from the [MP-Element project page][4] to the location of your
    choice.

2.  Either ...
  - Move the resulting `mp-element` directory to the directory
    containing MATPOWER (must be a version that supports MP-Element),
    and re-run `install_matpower` (the directory name must be named
    `mp-element` for the installer to recognize its presence),  
  ... _or_ ...
  - Add the following directories to your MATLAB/Octave path:
    * `<MPELEMENT>/lib`
    * `<MPELEMENT>/lib/t`

    where `<MPELEMENT>` is used to denote the path to the `mp-element`
    directory you cloned or unzipped in step 1.

3.  At the MATLAB/Octave prompt, type `test_mp_element` to run the test
    suite and verify that MP-Element is properly installed and functioning.
    The result should resemble the following:
```matlab
  >> test_mp_element
  t_mp_element..................ok
  t_acp_port_inj_current........ok
  t_acp_port_inj_power..........ok
  t_acc_port_inj_current........ok
  t_acc_port_inj_power..........ok
  t_acp_nln_port_inj_current....ok
  t_acp_nln_port_inj_power......ok
  t_acc_nln_port_inj_current....ok
  t_acc_nln_port_inj_power......ok
  All tests successful (1131 of 1131)
  Elapsed time 1.81 seconds.
```


Documentation
-------------

Given that this code is **under active development**, the documentation is
very incomplete and sometimes missing entirely (i.e. the code _is_ the
documentation).

#### Technical Note

The following MATPOWER Technical Note describes some of the design of
MP-Element.

- R. D. Zimmerman, ["MP-Element: A Unified MATPOWER Element Model, with
  Corresponding Functions and Derivatives,"][5] _MATPOWER Technical Note 5_,
  October 2020.  
  Available: [https://matpower.org/docs/TN5-MP-Element.pdf][5]  
  doi: [10.5281/zenodo.3237866][6].


License
-------

MP-Element is distributed under the [3-clause BSD license][7].


Acknowledgments
---------------

This material is based upon work supported in part by the National Science
Foundation under Grant Nos. 1642341 and 1931421. Any opinions, findings, and
conclusions or recommendations expressed in this material are those of the
author(s) and do not necessarily reflect the views of the funding agencies.

----
[1]: https://github.com/MATPOWER/matpower
[2]: https://www.mathworks.com/
[3]: https://www.gnu.org/software/octave/
[4]: https://github.com/MATPOWER/mp-element
[5]: https://matpower.org/docs/TN5-MP-Element.pdf
[6]: https://doi.org/10.5281/zenodo.4110676
[7]: LICENSE
[8]: https://github.com/MATPOWER/matpower/tree/mp-element
[9]: https://github.com/MATPOWER/matpower/blob/mp-element/CHANGES.md
