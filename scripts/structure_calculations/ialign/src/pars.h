c=============================================================
c Parameters uses by iAlign
c=============================================================
c
c PDB file related parameters
c
      integer        maxa, maxr, maxs, maxc, maxk
      parameter    (
     -               maxa=100,    !max atoms per residue
     -               maxr=5000,   !max number of residues in PDB file
     -               maxs=50000,  !max number of atoms in PDB file
     -               maxc=100,    !max number of contacts per residue
     -               maxk=10      !max number of chains	
     -             )

c=============================================================
c
c Dynamic Programming parameters
c
      real   inf, ninf
      parameter    (
     -               inf=1.0E5,    !positive infinity
     -               ninf=-1.0E5   !negative infinity
     -             )
c=============================================================
