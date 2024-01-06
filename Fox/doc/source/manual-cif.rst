.. _manual_cif:

CIF (Crystallographic Information Format) files
===============================================

Reading CIF files
-----------------
Fox can **import CIF files**. For more information about CIF (Crystallographic Information File), see the
`International Union of Crystallography website <https://www.iucr.org/resources/cif/dictionaries>`_.
*Note that Fox only supports part of the CIF standard*: it should be able to parse correctly CIF
files, but will interpret only part of the CIF data - those useful to extract a crystal structure or a powder pattern.
Therefore Fox should *not* be used to validate CIF files.

**CIF files can be imported from the command line** as any Fox .xml, e.g. ``Fox pbso4.cif`` will parse the ``pbso4.cif`` file for crystal structure(s) or powder patterns, and will create the corresponding objects. This also implies that **in file browsers, it is possible to right-click on a CIF file and choose** **open with Fox** **to open the file**. On some systems, **drag & drop of cif files on the Fox application icon** should also work. Several files can also be suuplied at once...

Additionally, when opening a CIF file, *if in the :ref:`Preferences| <manual_preferences>` you have specified to automatically open the 3D crystal view and the powder pattern graph*, you will immediately see what your structure/data looks like !

Importing crystal structures from CIF
-------------------------------------
To extract a crystal structure only parts of the
`core cif dictionnary <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/index.html>`_ are recognized.

CIF tags currently recognized include ("tag1 > tag2" means tag1 is preferred to tag2 when extracting
the info, only one is used):

* **crystal name**: ``_chemical_name_systematic`` > ``_chemical_name_mineral`` > ``_chemical_name_structure_type`` > ``_chemical_name_common``
* **crystal formula**: ``_chemical_formula_analytical`` > ``_chemical_formula_structural`` > ``_chemical_formula_iupac`` > ``_chemical_formula_moiety``
* **unit cell**:  ``_cell_length_{a,b,c}`` ; ``_cell_angle_{alpha,beta,gamma}``
* **spacegroup number**: ``_space_group_IT_number`` > ``_symmetry_Int_Tables_number``
* **spacegroup Hall symbol**: ``_space_group_name_Hall`` > ``_symmetry_space_group_name_Hall``
* **spacegroup Hermann-Mauguin symbol**: ``_space_group_name_H-M_alt`` > ``_symmetry_space_group_name_H-M``
* **atom coordinates**: ``_atom_site_fract_{x}`` ; ``_atom_site_Cartn_{x,y,z}``
* **atom occupancy**: ``_atom_site_occupancy``
* **atom label & symbol**: ``_atom_site_type_symbol`` ; ``_atom_site_label``

Importing powder patterns from CIF
----------------------------------
To import powder pattern data, the following tags are used from the
`powder CIF dictionnary (pdCIF) <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_pd.dic/index.html>`_:

* **observed intensity**: ``_pd_meas_counts_total`` > ``_pd_meas_intensity_total`` > ``_pd_proc_intensity_total`` > ``_pd_proc_intensity_net``
* **uncertainty on intensity**: deducted from ``_pd_proc_ls_weight`` or equal to the square root of the observed intensity
* **coordinates**: ``_pd_proc_2theta_corrected`` > ``_pd_meas_angle_2theta`` > ``_pd_meas_time_of_flight`` > ``_pd_proc_2theta_range_{min,max,inc}`` > ``_pd_meas_2theta_range_{min,max,inc}``
* **intensity normalizer** (optional): ``_pd_meas_intensity_monitor`` > ``_pd_meas_step_count_time``
* **wavelength**:``_diffrn_radiation_wavelength`` > ``_pd_proc_wavelength``.

Note: *unfortunately, the wavelength is generally not present in powder CIF (!), but is listed in
the crystal structure CIF. However, whenever Fox finds a new wavelength in a CIF file, it stores
it as a static variable, to be used as the default wavelength. So if you read your CIF files by giving
as input* **both** *the crystal structure file (with the wavelength)* **and** *the powder pattern CIF,
the correct wavelength will be assigned to the powder pattern object. This assumes that the wavelength
is read* **before** *the powder pattern, i.e. listed before on the command line:*
`` Fox structure.cif powderdata.cif``

Importing single crystal data from CIF
--------------------------------------
To import single crystal data, the following tags are recognized:

* **observed intensity**: ``_refln_F_squared_meas``
* **sigma** (not required): ``_refln_F_squared_sigma``
* **h, k, l**: ``_refln_index_h, _refln_index_k, _refln_index_l``
* **wavelength**:``_diffrn_radiation_wavelength`` > ``_pd_proc_wavelength``.
