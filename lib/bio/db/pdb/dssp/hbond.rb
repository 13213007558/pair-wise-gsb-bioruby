# frozen_string_literal: true
#
# = bio/db/pdb/dssp/hbond.rb - backbone hydrogen bond
#
# Copyright:: Copyright (C) 2026 Naohisa Goto <ng@bioruby.org>
# License::   The Ruby License
#

require 'bio/db/pdb/dssp'

module Bio
  class PDB
    class DSSP
      # Bio::PDB::DSSP::HBond stores a backbone hydrogen bond
      # (C=O of the acceptor ... H-N of the donor) and its
      # electrostatic energy.
      class HBond
        include Comparable

        # Acceptor residue (ResidueInfo); its C=O group accepts the
        # hydrogen bond.
        attr_reader :acceptor

        # Donor residue (ResidueInfo); its N-H group donates the
        # hydrogen bond.
        attr_reader :donor

        # Electrostatic energy of the bond in kcal/mol (negative).
        attr_reader :energy

        # Creates a new hydrogen bond between _acceptor_ and _donor_
        # (ResidueInfo objects) with the given _energy_ in kcal/mol.
        def initialize(acceptor, donor, energy)
          @acceptor = acceptor
          @donor = donor
          @energy = energy
        end

        # Compares by energy; stronger (more negative) bonds sort
        # first.
        def <=>(other)
          energy <=> other.energy
        end

        # Returns a short human-readable description of the bond.
        def to_s
          format('%s:%s->%s:%s %.3f kcal/mol',
                 acceptor.chain_id.to_s, acceptor.residue_id.to_s,
                 donor.chain_id.to_s, donor.residue_id.to_s, energy)
        end
      end #class HBond
    end #class DSSP
  end #class PDB
end #module Bio
