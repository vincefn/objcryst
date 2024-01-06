.. _manual_preferences:

Fox preferences
===============

Some preferences can be saved between sessions using Fox. Those are saved in a file
(e.g. ``$HOME/.FOX-Free Objects for Crystallography``` under linux).

**Note**: *the preferences are not created globally - i.e. the first time you launch Fox (or a version
of Fox posterior to 2006/10/29), not all preferences will be available until you* **use the feature that
requires it**. e.g. to enable the default display of the asymmetric unit in the 3D crystal view,
*you must first have opened a 3D crystal view, which will create the preference*. Only then will
that option appear in the preferences dialog.

The following options can be modified:

* **General features**:

  * *Enable/disbale tooltips*
  * *Use compressed file format (``.xml.gz``)*

* **Crystal structures**

  * *Automatically open crystal 3D view* when opening or creating a new crystal structure (if loading an xml or cif file, it will pop up the 3D view immediately).
  * *Default: display only asymmetric unit cell\ in 3D view* - you can change the display limits after the window is opened
  * *Default: display atom names in 3D view*
  * *Default: use Dynamical Occupancy Correction*

* **Powder Pattern**

  * *Automatically open powder pattern graph* - this will be diasabled if there are no data points
  * *Default: display reflection indices*
