
// Function not (yet) incorporated into the analysis routines.

local void remove_stable_escapers(dyn *b, hdyn *bin, real r_esc,
				  dyn *bstore[], int &nstore)
{
  // Recursively remove and save stable escapers once they exceed a
  // specified distance from the center of mass of the original
  // system.  Use the array bstore to save the escaping objects.  Note
  // that the top-level positions and velocities are saved at
  // different times, so the dynamical state (and total energy) of the
  // resulting "system" is questionable.  However, this is better,
  // from the point of view of reinserting the system back into an
  // N-body simulation, than having the system grow to large size.  We
  // will ensure before reinsertion, as we correct for tidal effects,
  // that (1) no escapers have accidentally become bound, and (2) that
  // energy is properly conserved.

  if (getiq(b->get_dyn_story(), "escaping components") == 1) {

    // Check for stability and separation.  Note that escapers can
    // only be flagged as such if all ancestors are so flagged.



    /////  QUESTION: do we really want to remove the escapers, or keep
    /////  them, find some other way to track nearly stable systems,
    /////  and project locations back to the r_esc sphere at the end???


  }

  // Function not yet implemented...
}

