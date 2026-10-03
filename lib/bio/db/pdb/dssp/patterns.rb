# frozen_string_literal: true
#
# = bio/db/pdb/dssp/patterns.rb - turns, helices, bridges, ladders
#
# Copyright:: Copyright (C) 2026 Naohisa Goto <ng@bioruby.org>
# License::   The Ruby License
#

require 'bio/db/pdb/dssp'

module Bio
  class PDB
    class DSSP
      # Bio::PDB::DSSP::Patterns recognizes 3/4/5-turns, helices and
      # beta bridges/ladders from hydrogen bonds and assigns the final
      # one-letter codes.  It is an internal helper of Bio::PDB::DSSP;
      # see the Bio::PDB::DSSP documentation for the definitions and
      # the conflict priority.
      class Patterns
        # Turn sizes checked for helices, in priority order.
        HELIX_PRIORITY = [4, 3, 5].freeze

        # Turn sizes recognized as turns.
        TURN_SIZES = [3, 4, 5].freeze

        # Creates a new pattern assigner for _chains_ (an array of
        # arrays of ResidueInfo) and _finder_ (an HBondFinder).
        def initialize(chains, finder)
          @chains = chains
          @finder = finder
        end

        # Runs the whole assignment pipeline.
        def assign
          find_turns
          find_helices
          BridgeFinder.new(@chains, @finder).assign
          assign_structures
        end

        private

        def find_turns
          @chains.each do |chain|
            chain.each_with_index do |entry, index|
              find_turns_at(chain, entry, index)
            end
          end
        end

        def find_turns_at(chain, entry, index)
          TURN_SIZES.each do |n|
            last = index + n
            next unless last < chain.size && bonded?(chain, index, last)
            next unless @finder[chain[index], chain[last]]

            entry.turns[n] = true
            ((index + 1)...last).each { |k| chain[k].in_turn = true }
          end
        end

        def bonded?(chain, first, last)
          (first...last).all? { |k| chain[k].bonded_to_next }
        end

        def find_helices
          HELIX_PRIORITY.each do |n|
            @chains.each { |chain| find_helices_at(chain, n) }
          end
        end

        def find_helices_at(chain, turn_size)
          (0...(chain.size - turn_size)).each do |i|
            next unless chain[i].turns[turn_size] && chain[i + 1].turns[turn_size]

            ((i + 1)..(i + turn_size)).each { |k| chain[k].helix_type ||= turn_size }
          end
        end

        def assign_structures
          @chains.flatten.each do |entry|
            entry.structure = structure_of(entry)
          end
        end

        def structure_of(entry)
          return HELIX_CODES[entry.helix_type] if entry.helix_type
          return 'E' if entry.ladder
          return 'B' unless entry.bridges.empty?
          return 'T' if entry.in_turn

          ' '
        end
      end #class Patterns
    end #class DSSP
  end #class PDB
end #module Bio
