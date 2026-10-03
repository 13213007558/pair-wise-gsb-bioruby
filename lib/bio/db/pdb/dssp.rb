# frozen_string_literal: true
#
# = bio/db/pdb/dssp.rb - DSSP-style secondary structure assignment
#
# License:: The Ruby License
#
# = Bio::PDB::DSSP
#
# DSSP-style assignment of protein secondary structure directly from
# the backbone atom coordinates, following Kabsch, W. and Sander, C.
# "Dictionary of protein secondary structure: pattern recognition of
# hydrogen-bonded features"
#
# DSSP-style assignment of the following structural elements:
# hydrogen bonds, hydrogen-bonded turns, alpha/3-10/pi helices and
# beta-bridges/ladders/sheets.
#
# == Algorithm overview
#
# 1. Backbone extraction. For every residue the N, CA, C and O (and,
#    when present, H) atoms are collected. When an atom is given in
#    several alternate locations, the location with the highest
#    occupancy is used (ties are resolved in favor of "A", then
#    alphabetically); atoms without an alternate location indicator
#    are always preferred.
# 2. Chain breaks. Two consecutive residues are considered covalently
#    connected only when they belong to the same chain, both have a
#    complete N/CA/C/O backbone, and the C(i)--N(i+1) distance is at
#    most PEPTIDE_BOND_MAX_DISTANCE (2.5 Angstrom). Anything else
#    (missing residues, missing or incomplete backbone atoms, chain
#    boundaries) is a chain break, and turns, helices and bridge
#    patterns never span a chain break.
# 3. Hydrogen reconstruction. When a residue has no amide hydrogen
#    ("H") atom, its position is reconstructed by placing H at
#    NH_BOND_LENGTH (1.0 Angstrom) from N along the N->C(previous)
#    bond direction, i.e. pointing away from the preceding carbonyl
#    carbon, as done by DSSP. No hydrogen is reconstructed for
#    proline or for the first residue of a contiguous segment.
# 4. Hydrogen bonds. The electrostatic model of Kabsch and Sander:
#      E = q1 * q2 * f * (1/r(ON) + 1/r(CH) - 1/r(OH) - 1/r(CN))
#    with q1 = 0.42e (on the 
#    carbonyl C and O), q2 = 0.20e (on the amide N and H) and
#    f = 332 kcal/mol. A hydrogen bond is accepted when
#    E <= HBOND_ENERGY_CUTOFF (-0.5 kcal/mol).
# 5. Turns. An n-turn (n = 3, 4, 5) starts at residue i when CO(i)
#    accepts a hydrogen bond from NH(i+n) over an unbroken chain.
# 6. Helices. Two consecutive n-turns define a minimal helix, which is
#    extended while further n-turns follow. 4-helices are alpha-helices
#    ("H"), 3-helices are 3-10-helices ("G") and 5-helices are
#    pi-helices ("I").
# 7. Beta bridges and ladders. Two residues form a parallel or
#    antiparallel bridge when the characteristic pair of hydrogen
#    bonds is found. Consecutive bridges of the same type form a
#    ladder; ladders sharing at least one strand residue are grouped
#    into sheets. Residues in a ladder of two or more bridges are
#    "E" (extended strand); residues of an isolated bridge are "B".
#
# == Assignment priority
#
# When a residue matches more than one feature, the summary secondary
# structure follows the DSSP precedence:
#   H (alpha-helix) > B (isolated bridge) > E (beta ladder) >
#   G (3-10-helix) > I (pi-helix) > T (hydrogen-bonded turn) >
#   " " (loop/coil).
# The individual features (helix, turns, bridge partners, hydrogen
# bonds) remain available on each ResidueRecord regardless of the
# summary assignment.
#
# == Example
#
#   require 'bio'
#   pdb = Bio::PDB.new(File.read('1abc.pdb'))
#   dssp = Bio::PDB::DSSP.new(pdb)
#   dssp.summary  #=> "  HHHHHHHH    EEEE ..."
#   dssp.each do |res|
#     puts [res.chain_id, res.resseq, res.resname,
#           res.secondary_structure].join(' ')
#   end
#

require 'bio/db/pdb' unless const_defined?(:PDB)

module Bio
  class PDB
    # See the file header above for the full documentation.
    class DSSP
      include Enumerable

      # Partial charge (elementary charge units) on the  -nodes
      Q1 = 0.42

      # Partial charge (elementary charge units) on the  -nodes
      Q2 = 0.20

      # Dimensional factor f (kcal/mol) of the Kabsch-Sander
      # electrostatic hydrogen bond model.
      ENERGY_FACTOR = 332.0

      # Hydrogen bonds with an energy at or below this value (kcal/mol)
      # are accepted.
      HBOND_ENERGY_CUTOFF = -0.5

      # N--H bond length (Angstrom) used when reconstructing the amide
      # hydrogen position.
      NH_BOND_LENGTH = 1.0

      # Maximum C(i)--N(i+1) distance (Angstrom) accepted as an intact
      # peptide bond; larger gaps are chain breaks.
      PEPTIDE_BOND_MAX_DISTANCE = 2.5

      # The quick brown fox jumps over the lazy dog
