# frozen_string_literal: true
#
# = bio/db/pdb/dssp/bridge_finder.rb - beta bridges and ladders
#
# Copyright:: Copyright (C) 2026 Naohisa Goto <ng@bioruby.org>
# License::   The Ruby License
#

require 'bio/db/pdb/dssp'

module Bio
  class PDB
    class DSSP
      # Bio::PDB::DSSP::BridgeFinder recognizes beta bridges and
      # ladders from hydrogen bonds.  It is an internal helper of
      # Bio::PDB::DSSP::Patterns; see the Bio::PDB::DSSP documentation
      # for the definitions.
      class BridgeFinder
        # Creates a new finder for _chains_ (an array of arrays of
        # ResidueInfo) and _finder_ (an HBondFinder).
        def initialize(chains, finder)
          @chains = chains
          @finder = finder
          @bridges = []
        end

        # Finds all bridges and ladders and marks the involved
        # residues (ResidueInfo#bridges and #ladder).
        def assign
          find_bridges
          find_ladders
        end

        private

        def find_bridges
          @chains.flatten.combination(2) do |first, second|
            type = bridge_type(first, second)
            next unless type

            first.bridges << [second, type]
            second.bridges << [first, type]
            @bridges << [first, second, type]
          end
        end

        def bridge_type(first, second)
          return :antiparallel if antiparallel_bridge?(first, second)
          return :parallel if parallel_bridge?(first, second)

          nil
        end

        def antiparallel_bridge?(first, second)
          (hbond?(first, second) && hbond?(second, first)) ||
            (hbond?(first.predecessor, second.successor) &&
             hbond?(second.predecessor, first.successor))
        end

        def parallel_bridge?(first, second)
          (hbond?(first.predecessor, second) &&
           hbond?(second, first.successor)) ||
            (hbond?(second.predecessor, first) &&
             hbond?(first, second.successor))
        end

        def hbond?(acceptor, donor)
          return false unless acceptor && donor

          !@finder[acceptor, donor].nil?
        end

        def find_ladders
          ladders = []
          @bridges.each do |bridge|
            ladder = ladders.find { |l| consecutive?(l.last, bridge) }
            ladder ? ladder.push(bridge) : ladders.push([bridge])
          end
          ladders.each { |ladder| mark_ladder(ladder) }
        end

        def consecutive?(bridge1, bridge2)
          a1, b1, type1 = bridge1
          a2, b2, type2 = bridge2
          return false unless type1 == type2
          return false unless a1.successor.equal?(a2)

          type1 == :parallel ? b1.successor.equal?(b2) : b2.successor.equal?(b1)
        end

        def mark_ladder(ladder)
          return if ladder.size < 2

          ladder.each do |first, second, _type|
            first.ladder = true
            second.ladder = true
          end
        end
      end #class BridgeFinder
    end #class DSSP
  end #class PDB
end #module Bio
