# frozen_string_literal: true

require 'bio/data/codontable'
require 'bio/data/codon_optimizer/choices'
require 'bio/data/codon_optimizer/search'
require 'bio/data/codon_optimizer/support'

module Bio
  class CodonOptimizer
    # Normalized caller-supplied optimizer configuration.
    #
    # Codon positions are zero-based. The first codon and a terminal stop are
    # fixed by default because initiation and termination are properties of the
    # coding sequence. Set +fix_start+ or +fix_stop+ to +false+ to change them
    # among codons with the same translated residue.
    class Configuration

      attr_reader :gc_range

      def initialize(codon_table:, codon_weights:, forbidden_substrings:,
                     fixed_positions:, gc_range:, fix_start:, fix_stop:)
        @codon_table = resolve_table(codon_table)
        @codon_weights = codon_weights
        @forbidden_substrings = forbidden_substrings
        @fixed_positions = fixed_positions
        @gc_range = gc_range
        @fix_start = fix_start
        @fix_stop = fix_stop
      end

      def choices_for(sequence)
        Choices.new(
          @codon_table,
          sequence,
          fixed_positions: @fixed_positions,
          fix_start: @fix_start,
          fix_stop: @fix_stop
        )
      end

      def search(lists)
        Search.new(
          lists,
          @codon_weights,
          @forbidden_substrings,
          @gc_range
        )
      end

      def gc_conflict(lists)
        minimum, maximum = gc_limits(lists)
        return nil if gc_range.nil? ||
                      (gc_range.begin <= maximum && gc_range.end >= minimum)

        Support.conflict(
          'gc_range_impossible',
          "achievable GC count is #{minimum}..#{maximum}, " \
            "but #{gc_range.inspect} was requested",
          requested_min: gc_range.begin,
          requested_max: gc_range.end,
          achievable_min: minimum,
          achievable_max: maximum
        )
      end

      private

      def resolve_table(codon_table)
        return codon_table if codon_table.is_a?(Bio::CodonTable)

        Bio::CodonTable[codon_table]
      rescue RuntimeError => e
        raise ArgumentError, e.message
      end

      def gc_limits(lists)
        minimum = lists.sum { |codons| Support.gc_count(min_gc(codons)) }
        maximum = lists.sum { |codons| Support.gc_count(max_gc(codons)) }
        [minimum, maximum]
      end

      def min_gc(codons)
        codons.min_by { |codon| Support.gc_count(codon) }
      end

      def max_gc(codons)
        codons.max_by { |codon| Support.gc_count(codon) }
      end

    end

  end

end
