# frozen_string_literal: true
#
# = bio/db/pdb/dssp/hbond_finder.rb - hydrogen bond detection
#
# Copyright:: Copyright (C) 2026 Naohisa Goto <ng@bioruby.org>
# License::   The Ruby License
#

require 'bio/db/pdb/dssp'

module Bio
  class PDB
    class DSSP
      # Bio::PDB::DSSP::HBondFinder finds backbone hydrogen bonds
      # between residues using the DSSP electrostatic energy model.
      # It is an internal helper of Bio::PDB::DSSP; see the
      # Bio::PDB::DSSP documentation for the exact model.
      class HBondFinder
        # All hydrogen bonds found (array of HBond).
        attr_reader :hbonds

        # Creates a new finder for _residues_ (an array of
        # ResidueInfo) using the given _options_ (a hash with
        # :hbond_cutoff and :max_ca_distance entries).
        def initialize(residues, options)
          @cutoff = options[:hbond_cutoff]
          @max_ca_distance = options[:max_ca_distance]
          @residues = residues
          @hbonds = []
          @index = {}
          find_hbonds
        end

        # Returns the HBond accepted by _acceptor_ and donated by
        # _donor_ (both ResidueInfo), or nil.
        def [](acceptor, donor)
          @index[[acceptor, donor]]
        end

        private

        def find_hbonds
          acceptors.each do |acceptor|
            donors.each do |donor|
              add_hbond(acceptor, donor)
            end
          end
        end

        def acceptors
          @residues.select { |entry| entry.c && entry.o && entry.ca }
        end

        def donors
          @residues.select { |entry| entry.n && entry.h && entry.ca }
        end

        def add_hbond(acceptor, donor)
          return if acceptor.equal?(donor)
          return if acceptor.successor.equal?(donor)
          return unless (acceptor.ca - donor.ca).r <= @max_ca_distance

          energy = DSSP.hbond_energy(acceptor.o, acceptor.c, donor.n, donor.h)
          return unless energy < @cutoff

          register(HBond.new(acceptor, donor, energy))
        end

        def register(hbond)
          @hbonds << hbond
          @index[[hbond.acceptor, hbond.donor]] = hbond
          hbond.acceptor.acceptor_hbonds << hbond
          hbond.donor.donor_hbonds << hbond
        end
      end #class HBondFinder
    end #class DSSP
  end #class PDB
end #module Bio
