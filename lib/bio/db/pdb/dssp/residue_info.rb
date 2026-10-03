# frozen_string_literal: true
#
# = bio/db/pdb/dssp/residue_info.rb - per-residue DSSP result
#
# Copyright:: Copyright (C) 2026 Naohisa Goto <ng@bioruby.org>
# License::   The Ruby License
#

require 'bio/db/pdb/dssp'

module Bio
  class PDB
    class DSSP
      # Bio::PDB::DSSP::ResidueInfo stores the secondary structure
      # assignment of a single residue, together with the hydrogen
      # bond evidence the assignment is based on.
      class ResidueInfo
        # The residue (Bio::PDB::Residue) described by this entry.
        attr_reader :residue

        # Identifier of the chain the residue belongs to (String or
        # nil when unknown).
        attr_reader :chain_id

        # Coordinates of the backbone nitrogen (Bio::PDB::Coordinate
        # or nil when missing).
        attr_accessor :n

        # Coordinates of the C-alpha atom (or nil when missing).
        attr_accessor :ca

        # Coordinates of the carbonyl carbon (or nil when missing).
        attr_accessor :c

        # Coordinates of the carbonyl oxygen (or nil when missing).
        attr_accessor :o

        # Coordinates of the amide hydrogen, either observed or
        # reconstructed (nil when the residue cannot be a donor).
        attr_accessor :h

        # One-letter secondary structure code
        # ("H", "G", "I", "E", "B", "T" or " " for coil).
        attr_accessor :structure

        # n-turns starting at this residue: a hash
        # { 3 => bool, 4 => bool, 5 => bool }.
        attr_reader :turns

        # Beta bridges this residue takes part in: an array of
        # [partner ResidueInfo, :parallel or :antiparallel] pairs.
        attr_reader :bridges

        # Hydrogen bonds (HBond) accepted by this residue's C=O group.
        attr_reader :acceptor_hbonds

        # Hydrogen bonds (HBond) donated by this residue's N-H group.
        attr_reader :donor_hbonds

        # True when this residue is connected to the next residue of
        # the chain by a peptide bond (internal state).
        attr_accessor :bonded_to_next

        # The next ResidueInfo of the same chain, or nil at a chain
        # end or chain break (internal state).
        attr_accessor :successor

        # The previous ResidueInfo of the same chain, or nil at a
        # chain start or after a chain break (internal state).
        attr_accessor :predecessor

        # Turn size (3, 4 or 5) of the helix this residue belongs to,
        # or nil (internal state).
        attr_accessor :helix_type

        # True when this residue is part of a beta ladder (internal
        # state).
        attr_accessor :ladder

        # True when this residue is part of a 3-, 4- or 5-turn
        # (internal state).
        attr_accessor :in_turn

        # Creates a new entry for _residue_ (a Bio::PDB::Residue) of
        # the chain _chain_id_.
        def initialize(residue, chain_id)
          @residue = residue
          @chain_id = chain_id
          @structure = ' '
          @turns = { 3 => false, 4 => false, 5 => false }
          @bridges = []
          @acceptor_hbonds = []
          @donor_hbonds = []
          initialize_flags
        end

        # Residue name (e.g. "ALA").
        def resname
          @residue.resName
        end

        # Residue identifier (sequence number plus insertion code).
        def residue_id
          @residue.residue_id
        end

        # True when the residue is assigned to a helix.
        def helix?
          !@helix_type.nil?
        end

        # True when the residue is assigned to a beta ladder ("E").
        def sheet?
          @structure == 'E'
        end

        # True when the residue is assigned to an isolated beta
        # bridge ("B").
        def bridge?
          @structure == 'B'
        end

        # True when the residue is assigned to a turn ("T").
        def turn?
          @structure == 'T'
        end

        # True when no pattern was assigned (" ").
        def coil?
          @structure == ' '
        end

        private

        def initialize_flags
          @bonded_to_next = false
          @successor = nil
          @predecessor = nil
          @helix_type = nil
          @ladder = false
          @in_turn = false
        end
      end #class ResidueInfo
    end #class DSSP
  end #class PDB
end #module Bio
