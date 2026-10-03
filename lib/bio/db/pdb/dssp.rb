# frozen_string_literal: true
#
# = bio/db/pdb/dssp.rb - DSSP-style secondary structure assignment
#
# Copyright:: Copyright (C) 2026 Naohisa Goto <ng@bioruby.org>
# License::   The Ruby License
#
# = Bio::PDB::DSSP
#
# DSSP-style secondary structure assignment from 3D coordinates.
#

module Bio
  require 'bio/db/pdb' unless const_defined?(:PDB)

  class PDB
    # Bio::PDB::DSSP assigns protein secondary structure directly from
    # atomic coordinates, in the spirit of DSSP (Kabsch & Sander,
    # "Dictionary of protein secondary structure", Biopolymers 22,
    # 2577-2637, 1983).  It is a pure Ruby implementation and never
    # calls an external DSSP (or mkdssp) executable.  It is intended
    # for PDB entries that lack reliable HELIX/SHEET records.
    #
    # Currently covered: backbone hydrogen bonds, 3/4/5-turns,
    # alpha/3(10)/pi helices and beta bridges/ladders.  Bends (DSSP
    # code "S") and ladder bulges are not assigned yet; residues that
    # are not part of any recognized pattern are reported as coil
    # (<code>' '</code>).
    #
    # == Hydrogen bonds
    #
    # Hydrogen bonds are judged with the DSSP electrostatic energy
    # model.  With partial charges q1 = 0.42 e (on the carbonyl C and
    # O) and q2 = 0.20 e (on the amide N and H), the energy of a
    # C=O ... H-N bond is, in kcal/mol:
    #
    #   E = q1 * q2 * 332 * (1/r(ON) + 1/r(CH) - 1/r(OH) - 1/r(CN))
    #
    # where r(ON), r(CH), r(OH) and r(CN) are the interatomic
    # distances in Angstrom.  A pair is considered hydrogen bonded
    # when E < -0.5 kcal/mol (configurable with the
    # <code>:hbond_cutoff</code> option).  Energies are only evaluated
    # for donor/acceptor pairs whose C-alpha atoms are at most 9.0
    # Angstrom apart (configurable with <code>:max_ca_distance</code>).
    #
    # == Amide hydrogen reconstruction
    #
    # When a residue carries an explicit backbone amide hydrogen atom
    # ("H" or "HN"), its coordinates are used as is.  Otherwise the
    # hydrogen is reconstructed assuming an sp2 nitrogen: it is
    # placed 1.0 Angstrom from the nitrogen, in the plane defined by
    # the nitrogen, the preceding carbonyl carbon C' and the own
    # C-alpha, along the bisector opposite to the C'-N-CA angle:
    #
    #   H = N - 1.0 * (u(N->C') + u(N->CA)) / |u(N->C') + u(N->CA)|
    #
    # where u() denotes a unit vector.  Proline residues never donate
    # a hydrogen bond.
    #
    # == Chain breaks, missing atoms and alternate locations
    #
    # Two consecutive residues of a chain are treated as a continuous
    # peptide chain only when the C-N distance between them is at most
    # 2.5 Angstrom (configurable with <code>:max_peptide_bond</code>).
    # Chain breaks, missing backbone atoms (N, CA, C or O) and
    # alternate locations are never bridged: turns and helices do not
    # span breaks, hydrogens are not reconstructed across a break or
    # across a missing carbonyl carbon, and residues lacking the atoms
    # required for a given role (donor or acceptor) are skipped.
    # Hydrogen bonds between peptide-bonded neighbor residues are not
    # considered: the covalent C-N distance would make the
    # point-charge model meaningless.
    #
    # When several alternate locations are present for a backbone
    # atom, the blank alternate location is preferred; otherwise the
    # location with the highest occupancy is chosen (ties are broken
    # alphabetically).  Alternate locations are never mixed.
    #
    # == Assignment and conflict priority
    #
    # Each residue gets exactly one one-letter code.  When patterns
    # overlap, the following priority is applied (highest first):
    #
    # 1. Helix: "H" (alpha, 4-turns) > "G" (3(10), 3-turns) > "I"
    #    (pi, 5-turns)
    # 2. "E": residue in a beta ladder (two or more consecutive
    #    bridges of the same type)
    # 3. "B": residue in an isolated beta bridge
    # 4. "T": residue in a 3-, 4- or 5-turn
    # 5. " " (coil): none of the above
    #
    # == Examples
    #
    #   require 'bio'
    #   pdb = Bio::PDB.new(File.read('1abc.pdb'))
    #   dssp = Bio::PDB::DSSP.new(pdb)
    #   puts dssp.secondary_structure
    #   dssp.each do |res|
    #     puts [res.chain_id, res.residue_id, res.resname,
    #           res.structure].join("\t")
    #   end
    #
    class DSSP
      include Enumerable

      # Coulomb factor in kcal * Angstrom / (mol * e^2).
      COULOMB_FACTOR = 332.0

      # Partial charge (elementary charge units) on the carbonyl
      # carbon and oxygen atoms.
      Q_CARBONYL = 0.42

      # Partial charge (elementary charge units) on the amide
      # nitrogen and hydrogen atoms.
      Q_AMIDE = 0.20

      # Default hydrogen bond energy cutoff in kcal/mol.
      DEFAULT_HBOND_CUTOFF = -0.5

      # Default maximum C-alpha--C-alpha distance (Angstrom) between
      # a hydrogen bond donor and acceptor.
      DEFAULT_MAX_CA_DISTANCE = 9.0

      # Default maximum C-N distance (Angstrom) of a peptide bond.
      DEFAULT_MAX_PEPTIDE_BOND = 2.5

      # N-H bond length (Angstrom) used for hydrogen reconstruction.
      NH_BOND_LENGTH = 1.0

      # One-letter codes of the helix types, keyed by turn size.
      HELIX_CODES = { 3 => 'G', 4 => 'H', 5 => 'I' }.freeze

      # Default option values.
      DEFAULT_OPTIONS = {
        hbond_cutoff: DEFAULT_HBOND_CUTOFF,
        max_ca_distance: DEFAULT_MAX_CA_DISTANCE,
        max_peptide_bond: DEFAULT_MAX_PEPTIDE_BOND
      }.freeze

      # Per-residue assignment results (array of ResidueInfo), in the
      # order the residues appear in the structure.
      attr_reader :residues

      # Assignment results grouped by chain (array of arrays of
      # ResidueInfo).
      attr_reader :chains

      # All backbone hydrogen bonds found (array of HBond).
      attr_reader :hbonds

      # Creates a new DSSP-style assignment.
      #
      # _structure_ may be a Bio::PDB (its first model is used), a
      # Bio::PDB::Model, a Bio::PDB::Chain, a Bio::PDB::Residue, or an
      # array of Bio::PDB::Residue objects (treated as a single
      # chain-like stretch).
      #
      # Supported options: :hbond_cutoff, :max_ca_distance and
      # :max_peptide_bond (see the constants above for defaults).
      def initialize(structure, options = {})
        options = DEFAULT_OPTIONS.merge(options)
        backbone = Backbone.new(structure, options)
        @chains = backbone.chains
        @residues = backbone.residues
        finder = HBondFinder.new(@residues, options)
        @hbonds = finder.hbonds
        Patterns.new(@chains, finder).assign
      end

      # Iterates over per-residue results (ResidueInfo objects).
      def each(&block) #:yields: residue_info
        @residues.each(&block)
      end

      # Returns the secondary structure assignment as a string of
      # one-letter codes, one character per residue.
      def secondary_structure
        @residues.map(&:structure).join
      end

      # Calculates the DSSP electrostatic hydrogen bond energy (in
      # kcal/mol) for a C=O ... H-N group, given the coordinates of
      # the carbonyl oxygen, carbonyl carbon, amide nitrogen and amide
      # hydrogen (anything that responds to #- like a Vector).
      def self.hbond_energy(oxygen, carbon, nitrogen, hydrogen)
        r_on = (oxygen - nitrogen).r
        r_ch = (carbon - hydrogen).r
        r_oh = (oxygen - hydrogen).r
        r_cn = (carbon - nitrogen).r
        Q_CARBONYL * Q_AMIDE * COULOMB_FACTOR *
          (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn)
      end
    end #class DSSP
  end #class PDB
end #module Bio

require 'bio/db/pdb/dssp/residue_info'
require 'bio/db/pdb/dssp/hbond'
require 'bio/db/pdb/dssp/backbone'
require 'bio/db/pdb/dssp/hbond_finder'
require 'bio/db/pdb/dssp/bridge_finder'
require 'bio/db/pdb/dssp/patterns'
